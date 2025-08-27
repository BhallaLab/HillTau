# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
# 

'''
*******************************************************************
 * File:            hillTau.py
 * Description:
 * Author:          Upinder S. Bhalla
 * E-mail:          bhalla@ncbs.res.in
 ********************************************************************/
 '''
from __future__ import print_function
import sys
import json
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import ht

lookupQuantityScale = { "M": 1000.0, "mM": 1.0, "uM": 1e-3, "nM": 1e-6, "pM": 1e-9 }

SIGSTR = "{:.4g}" # Used to format floats to keep to 4 sig fig. Helps when dumping JSON files.

def loadHillTau( fname ):
    with open( fname ) as json_file:
        model = json.load( json_file )
        ft = model.get( "FileType" )
        if not ft:
            print( 'Warning: Is "{}" a HillTau file? It lacks "FileType": "HillTau" specification.'.format(fname) )
        else:
            if ft != "HillTau":
                raise( ValueError( "Error: FileType should be 'HillTau', was {}.".format( ft ) ) )
    return model

def subsetModel( model, subsetInfo ):
    return

class Stim():
    def __init__( self, stim, model, off = False ):
        self.objname = stim[0]
        self.mol = model.molInfo.get( stim[0] )
        if not self.mol:
            print( "Stimulus Molecule '{}' not found".format( stim[0] ) )
            quit()
        self.value = stim[1]
        self.isOff = off
        if off:
            self.time = stim[3]
            self.value = self.mol.concInit
        else:
            self.time = stim[2]

    @staticmethod
    def stimOrder( stim ):
        return stim.time

def getQuantityScale( jsonDict ): 
    qu = jsonDict.get( "QuantityUnits" )
    qs = 1.0
    if not qu:
        qu = jsonDict.get( "quantityUnits" )
    if qu:
        qs = lookupQuantityScale[qu]
    return qs

def scaleConst( holder, name, qs, consts, constDone ):
    # This scales values with concentration units. It checks if the value
    # is a string, in which case it tries to scale the reference constant.
    # Otherwise it scales the entry in situ.
    val = holder.get( name )
    if isinstance( val, str ):
        constName = val
        ci = consts.get( constName )
        if not ci:
            raise( ValueError( "Error: Constant {} not found.".format( constName ) ) )
        if not constDone[ constName ]:
            consts[ constName] = float( SIGSTR.format( ci * qs ) )
            constDone[ constName ] = 1
    else:
        holder[name] = float( SIGSTR.format( val * qs ) )


def scaleDict( jsonDict, qs ):
    # This acts before parsing the model, so it should leave the model
    # definition layout intact. That means it should scale the
    # constant definitions rather than fill in values where the the constants are
    # cited in the model. This means we have to take care not to scale
    # a given constant more than once, as it may be used many times.
    consts = jsonDict.get("Constants")
    if not consts:
        consts = {}
    constDone = { i:0 for i in consts.keys()} # dict of all const names.
    for grpname, grp in jsonDict["Groups"].items():
        sp = grp.get( "Species" )
        if sp:
            for m in sp:
                scaleConst( sp, m, qs, consts, constDone )
                #sp[m] = float( SIGSTR.format( sp[m] * qs ) )
        if "Reacs" in grp:
            for reacname, reac in grp['Reacs'].items():
                # Check if it is a single substrate reac
                if len( reac["subs"] ) != 1:
                    scaleConst( reac, "KA", qs, consts, constDone )
                '''
                reac["tau"] = float( SIGSTR.format( reac["tau"] ) )
                tau2 = reac.get( "tau2" )
                if tau2:
                    reac["tau2"] = float( SIGSTR.format( tau2 ) )
                '''
                bl = reac.get( "baseline" )
                if bl:
                    scaleConst( reac, "baseline", qs, consts, constDone )
                kmod = reac.get( "Kmod" )
                if kmod:
                    scaleConst( reac, "Kmod", qs, consts, constDone )

def extractSubs( expr, consts ):
    # This function extracts the molecule names from a math expression.
    isInMol = 0
    mathFns = ["exp", "log", "ln", "log10", "abs", "sin", "cos", "tan", "sinh", "cosh", "tanh", "sqrt", "pow"]
    molname = ""
    lastch = ''
    mols = []
    for ch in expr:
        if isInMol:
            if ch.isalnum() or ch == '_':
                molname += ch
                continue
            else:
                isInMol = 0
                if not molname in mathFns:
                    mols.append( molname )
        else:
            if ch in "eE" and (lastch.isdigit() or lastch == '.'):
                # This is the e in sci notation
                lastch = ch
                continue
            if ch.isalpha():
                molname = ch
                isInMol = 1
            else:
                lastch = ch
    if isInMol:
        mols.append( molname )
    s = []
    c = []
    for key in mols:
        if key in consts:
            c.append( key )
        else:
            s.append( key )
    return s, c

def convConst( consts, value ):
    if isinstance( value, str ):
        ret  = consts.get( value )
        if ret:
            return ret
        else:
            raise( ValueError( "Error: Const '{}' not found.".format( value ) ) )

    return value

def parseModel( jsonDict ):
    model = ht.Model()
    consts = jsonDict.get( "Constants" )
    if not consts:
        consts = {}
    model.namedConsts = consts
    eqnSubs = {}
    # First, pull together all the species names. They crop up in
    # the Species, the Reacs, and the Eqns. They should be used as
    # an index to the conc and concInit vector.
    # Note that we have an ordering to decide which mol goes in which group:
    # Species; names of reacs, First term of Eqns, substrates.
    # This assumes that every quantity term has already been scaled to mM.
    for grpname, grp in jsonDict['Groups'].items():
        model.addGrp( grpname )
        # We may have repeats in the species names as they are used 
        # in multiple places.
        if "Reacs" in grp:
            for reacname, reac in grp['Reacs'].items():
                for subname in reac["subs"]:
                    model.makeMol( subname, grpname, isSub = 1 )
                    #mi[subname] = ht.MolInfo( subname, grpname, order=0)
        if "Eqns" in grp:
            for lhs, expr in grp["Eqns"].items():
                subs, cs = extractSubs( expr, consts )
                eqnSubs[ lhs ] = [ subs, cs ]
                for subname in subs:
                    model.makeMol( subname, grpname, isSub = 1 )

    for grpname, grp in jsonDict['Groups'].items():
        if "Eqns" in grp:
            for lhs, expr in grp["Eqns"].items():
                model.makeMol( lhs, grpname )
        if "Reacs" in grp:
            for reacname, reac in grp['Reacs'].items():
                model.makeMol( reacname, grpname )

    for grpname, grp in jsonDict['Groups'].items():
        if "Species" in grp:
            for molname, conc in grp['Species'].items():
                conc = convConst( consts, conc )
                model.makeMol( molname, grpname, concInit = conc )
                #mi[molname] = ht.MolInfo( molname, grpname, order=0, concInit = conc )
                grp['Species'][molname] = conc

    # Then assign indices to these unique molnames, and build up the
    # numpy arrays for concInit and conc.
    model.allocConc();

    # Now set up the reactions. we need the mols all defined first.
    for grpname, grp in jsonDict['Groups'].items():
        if "Reacs" in grp:
            for reacname, reac in grp['Reacs'].items():
                subs = reac['subs']
                convReac = { key: convConst( consts, val ) for key,val in reac.items() }
                # hack to interface with model::makeReac, which
                # expects all args in reac to be floats.
                convReac['subs'] = 0.0
                model.makeReac( reacname, grpname, subs, convReac )

    # Now set up the equation, again, we need the mols defined.
    for grpname, grp in jsonDict['Groups'].items():
        if "Eqns" in grp:
            for lhs, expr in grp["Eqns"].items():
                e = expr
                subs = eqnSubs[ lhs ][0]
                cs = eqnSubs[lhs][1]
                # Replace the constant names with their values
                for name in cs:
                    val = consts.get( name )
                    if val:
                        e = e.replace( name, str( val ) )
                    else:
                        raise( ValueError( "Error: unknown const '{}' in equation '{    }'".format( name, expr ) ) )
                model.makeEqn( lhs, grpname, e, subs )

    model.allocConc()
    sortReacs( model )
    model.reinit()
    return model

def breakReacLoop( sri, model, maxOrder ):
    for reacname, reac in sri:
        if model.updateMolOrder( maxOrder, reacname ):
            return
        '''
        if model.molInfo[reacname].order < 0:
            model.molInfo[reacname].order = maxOrder
            #print( " BREAK LOOP on ", reacname, " ", maxOrder)
            #print("Warning; Reaction order loop. Breaking {} loop for {}, assigning order: {}".format( numLoopsBroken, reacname, maxOrder ) )
            return
        '''

'''
def breakEqnLoop( model, maxOrder, numLoopsBroken  ):
    for eqname, eqn in model.eqnInfo.items():
        if model.molInfo[eqname].order < 0:
            model.molInfo[eqname].order = maxOrder
            #print( "    FIX_Eqn ORDER = ", eqname, " ", maxOrder)
            return
'''

def sortReacs( model ):
    # Go through and assign levels to the mols and reacs within a group.
    # This will be used later for deciding evaluation order.
    maxOrder = 0
    numLoopsBroken = 0
    numOrdered = 0
    numReac = len( model.reacInfo )
    sri = sorted( model.reacInfo.items() )
    while numOrdered < numReac: 
        numOrdered = 0
        stuck = True
        #for reacname, reac in sorted( model.reacInfo.items() ):
        for reacname, reac in sri:
            #prevOrder = model.molInfo[reacname].order
            prevOrder = model.getMolOrder(reacname)
            maxOrder = max( maxOrder, prevOrder )
            if prevOrder >= 0:
                numOrdered += 1
            else:
                order = reac.getReacOrder()
                # As a side effect it assigns model.molInfo[reacname].order
                if order >= 0:
                    maxOrder = max( maxOrder, order )
                    numOrdered += 1
                    stuck = False
                '''
                order = [ model.molInfo[i].order for i in reac.subs ]
                if min( order ) >= 0:
                    mo = max( order ) + 1
                    model.molInfo[reacname].order = mo
                    maxOrder = max( maxOrder, mo )
                    numOrdered += 1
                    stuck = False
                '''
        #print ( "numOrdered = ", numOrdered, " / ", numReac, " max = ", maxOrder )
        if stuck:
            breakReacLoop( sri, model, maxOrder+1 )
            numLoopsBroken += 1

    # We don't need to sort equations, because they do not cascade.
    # They are all executed in a bunch after the reacs, at which point
    # there should be no unknowns

    maxOrder += 1
    model.setReacSeqDepth( maxOrder )
    for name, reac in model.reacInfo.items():
        order = model.molInfo[name].order
        model.assignReacSeq( name, order )

def writeOutput( fname, model, plotvec, x ):
    with open( fname, "w" ) as fd:
        olist = sorted([ i for i in model.molInfo])
        header = "Time\t"
        outvec = [[str(v) for v in x]]
        rx = range( len( x ) )
        for name in olist:
            header += name + "\t"
            idx = model.molInfo[name].index
            outvec.append( [str(v) for v in plotvec[idx] ] )
        ry = range( len( outvec ) )
        fd.write( header + "\n" )
        for i in rx:
            for j in ry:
                fd.write( outvec[j][i] + "\t" )
            fd.write( "\n" )


def main():
    parser = argparse.ArgumentParser( description = 'This is the hillTau simulator.\n'
    'This program simulates abstract kinetic/neural models defined in the\n'
    'HillTau formalism. HillTau is an event-driven JSON form to represent\n'
    'dynamics of mass-action chemistry and neuronal activity in a fast, \n'
    'reduced form. The hillTau program loads and checks HillTau models,\n'
    'and optionally does simple stimulus specification and plotting\n')
    parser.add_argument( 'model', type = str, help='Required: filename of model, in JSON format.')
    parser.add_argument( '-r', '--runtime', type = float, help='Optional: Run time for model, in seconds. If flag is not set the model is not run and there is no display', default = 0.0 )
    parser.add_argument( '-dt', '--dt', type = float, help='Optional: Time step for model calculations, in seconds. If this argument is not set the code calculates dt to be a round number about 1/100 of runtime.', default = -1.0 )
    parser.add_argument( '-s', '--stimulus', type = str, nargs = '+', action='append', help='Optional: Deliver stimulus as follows: --stimulus molecule conc [start [stop]]. Any number of stimuli may be given, each indicated by --stimulus. By default: start = 0, stop = runtime', default = [] )
    parser.add_argument( '-p', '--plots', type = str, help='Optional: plot just the specified molecule(s). The names are specified by a comma-separated list.', default = "" )
    parser.add_argument( '-o', '--output', type = str, metavar = "fname", help='Optional: Generate an output tab-separated text file with columns of time conc1 conc2 and so on.' )
    args = parser.parse_args()
    jsonDict = loadHillTau( args.model )
    qs = getQuantityScale( jsonDict )
    scaleDict( jsonDict, qs )
    model = parseModel( jsonDict )

    runtime = args.runtime
    if runtime <= 0.0:
        return

    if args.dt < 0:
        model.dt = 10 ** (np.floor( np.log10( runtime )) - 2.0)
        if runtime / model.dt > 500:
            model.dt *= 2
    else:
        model.dt = args.dt

    stimvec = []
    
    stimMolNames = []
    for i in args.stimulus:
        if len( i ) < 2:
            print( "Warning: need at least 2 args for stimulus, got {i}".format( i ) )
            continue
        stimMolNames.append( i[0] )
        i[1] = float( i[1] ) * qs # Assume stim units same as model units.
        if len(i) == 2:
            i.extend( [0.0, runtime] )
        if len(i) == 3:
            i.extend( [runtime] )
        i[2] = float( i[2] )
        i[3] = float( i[3] )
        runtime = max( runtime, i[3] )
        stimvec.append( Stim( i, model ) )
        stimvec.append( Stim( i, model, off = True ) )

    stimvec.sort( key = Stim.stimOrder )
    model.modifySched( saveList = [], deleteList = list( set( stimMolNames ) ) )
    clPlots = args.plots.split(',')
    if len( args.plots ) > 0 :
        #clPlots = [ i.strip() for i in clPlots if i in model.molInfo]
        for pp in clPlots:
            model.addPlot( pp )
    
    model.reinit()
    currTime = 0.0
    for s in stimvec:
        model.advance( s.time - currTime )
        model.conc[s.mol.index] = s.value
        currTime = s.time
    if runtime > currTime:
        model.advance( runtime - currTime )

    x = model.plotTimePoints

    if args.output:
        writeOutput( args.output, model, plotvec, x )

    qu = jsonDict.get( "QuantityUnits" )
    if not qu:
        qu = jsonDict.get( "quantityUnits" )
    if qu:
        ylabel = 'Conc ({})'.format( qu )
        qs = lookupQuantityScale[qu]
    else:
        ylabel = 'Conc (mM)'
        qs = 1

    if len( args.plots ) > 0:
        for idx, name in enumerate( clPlots ):
            plt.plot( x, np.array(model.plotvec[idx])/qs, label = name )
    else: # do all plots
        for name, val in model.molInfo.items():
            plt.plot( x, np.array(model.getConcVec(val.index))/qs, label = name )

    plt.xlabel('Time (s)')
    plt.ylabel(ylabel)
    plt.title( args.model )
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
