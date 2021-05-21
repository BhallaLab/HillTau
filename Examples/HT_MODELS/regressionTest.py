import os
import sys
#sys.path.insert(1, '/home/bhalla/homework/HILLTAU/REPO/HillTau/PythonCode')
import hillTau

ERR_LIMIT = 1e-6

stimVec = [
    ["exc", "input", [1e-3, 10, 0, 20], "output", [0.2638e-3, 10.75, 0.5e-3, 20, 0.0, 30]],
    ["inh", "input", [1e-3, 10, 0, 20], "output", [0.7362e-3, 10.75, 0.5e-3, 20, 1e-3, 30]],
    ["osc", "mol", [], "output", [0.3968e-3, 500, 0.2948e-3, 2000, 0.1008e-3, 3000]],
    ["bcm", "Ca", [0.01e-3, 0, 0.5e-3, 20, 2e-3, 40, 10e-3, 60], "synAMPAR", [0.384e-3, 15, 0.189e-3, 35, 0.391e-3, 55, 0.503e-3, 75]],
    ["eqn", "input", [1e-3, 10, 0, 20], "eq", [3.46e-3, 10.75, 3.7e-3, 15, 1.2e-3, 30]],
    ["eqn_with_constants", "input", [1e-3, 10, 0, 20], "eq", [3.46e-3, 10.75, 3.7e-3, 15, 1.2e-3, 30]],
    ["exc2ndOrder", "input", [1e-3, 10, 0, 20], "output", [0.2638e-3, 10.75, 0.5e-3, 20, 0.0, 30]],
    ["conv", "input", [1e-3, 10, 0, 20], "output", [0.5276e-3, 10.75, 1e-3, 20, 0.0, 30]],
    ["conv2ndOrder", "input", [1e-3, 10, 0, 20], "output", [0.5276e-3, 10.75, 1e-3, 20, 0.0, 30]],
    ["exc_tau_baseline", "input", [1e-3, 10, 0, 20], "output", [0.7638e-3, 10.75, 1.0e-3, 20, 0.684e-3, 25]],
    ["ff_inhib", "input", [1e-3, 10, 0, 20], "output", [0.855e-3, 12, 0.2085e-3, 20, 0, 30]],
    ["fb_inhib", "input", [1e-3, 20, 0, 60], "output", [0.776e-3, 24, 0.188e-3, 60, 0.5e-6, 100]],
    ["gain", "input", [1e-3, 10, 0, 20], "output", [0.5276e-3, 10.75, 1e-3, 20, 0.0, 30]],
    ["modifier", "input", [1e-3, 10, 0, 20], "output", [0.474e-3, 11, 0.75e-3, 19, 0, 30]],
    #["bistable", "mod", [100, 20, 1, 25, 0.01, 40, 1, 45, 0.01, 60, 1, 70], "output", [0.03, 20, 1.024, 40, 1.024, 60, 0.03, 80]],
    ["bistable", "input", [10, 20, 1, 25, 0.1, 40, 1, 45, 0.1, 60, 1, 70], "output", [0.031, 20, 0.978, 40, 0.978, 60, 0.031, 80]],
    ["bcm_bistable", "Ca", [2e-3, 15, 0.08e-3, 16, 0.4e-3, 35, 0.08e-3, 45], "synAMPAR", [0.207e-6, 15, 0.8205e-3, 20, 0.8331e-3, 30, 0.21e-6, 60 ]],
]

class Event():
    def __init__( self, mol, conc, t, isStim ):
        self.mol = mol
        self.isStim = isStim
        self.t = float(t)
        self.conc = float( conc )

def parseEvents( mol, f, isStim ):
    ret = []
    idx = range( 0, len(f), 2 )
    for i in idx:
        ret.append( Event( mol, f[i], f[i+1], isStim ) )
    return ret


def runit( f ):
    jsonDict = hillTau.loadHillTau( f[0] + ".json" )
    qs = hillTau.getQuantityScale( jsonDict )
    hillTau.scaleDict( jsonDict, qs )
    model = hillTau.parseModel( jsonDict )
    model.dt = 1.0

    ev = []
    ev.extend( parseEvents( f[1], f[2], 1 ) )
    ev.extend( parseEvents( f[3], f[4], 0 ) )
    sev = sorted( ev, key = lambda x: x.t )

    model.reinit()
    lastt = 0.0 
    ans = []
    maxconc = 0.0
    for s in sev:
        delta = s.t - lastt
        if delta > 0.0:
            model.advance( delta )
        if s.isStim:
            model.conc[ model.molInfo[ s.mol ].index ] = s.conc
        else:
            c = model.conc[ model.molInfo[ s.mol ].index ]
            #print( "Conc @ {} = {}".format( s.t, c ) )
            maxconc = max( maxconc, c )
            ans.append( s.conc - model.conc[ model.molInfo[ s.mol ].index ] )
        lastt = s.t
    err = 0.0
    #print( "ANS = ", ans )
    for d in ans:
        err += d*d / ( maxconc * maxconc )
    return err / len( ans ), model

def checkGroups( model ):
    # Assumes that the last loaded model was bcm_bistable
    grps = {"internal": "ampar_g", "synAMPAR": "ampar_g", "on_CaMKII":"CaMKII_g", "CaN":"CaN_g", "fb": "CaMKII_g"}
    print( "Checking groups{:20s}".format( "" ), end = "....     " )
    OK = True
    for mol, grp in grps.items():
        if model.molInfo[mol].grp != grp:
            print( "failed, group of mol '{}' should be '{}', but is '{}'".format( mol, grp, model.molInfo[mol].grp ) )
            OK = False
    if OK:
        print( "OK, all objects are in correct group" )


def main():
    for s in stimVec:
        print( "Checking model {:20s}".format( s[0] ), end = "....     " )
        err, model = runit( s )
        if err > ERR_LIMIT:
            print( "failed, err = {:.5g}".format( err ) )
        else:
            print( "OK, err = {:.5g}".format( err ) )
    checkGroups( model )

if __name__ == '__main__':
    main()
