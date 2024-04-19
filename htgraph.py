#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.

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
 * File:            htgraph.py
 * Description:
 * Author:          G.V. Harsha Rani, Upinder S. Bhalla
 * E-mail:          hrani@ncbs.res.in, bhalla@ncbs.res.in
 ********************************************************************/

/**********************************************************************
** This program converts HILLTAU models defined in JSON format to 
** reaction diagrams. It draws on the 'dot' program for graphical layout
** of networks.
**           copyright (C) 2021 Harsha Rani, Upinder S. Bhalla. and NCBS
**********************************************************************/
'''
'''
2024
Apr
After subset if equation has direct molecule (which name become as concInit) and summation, then that is removed in the layout
2022
Nov
no-group -ng (showGroup default True)
specific_group -sg (Pass Group name in quotes separated by comma -sg "Ras_g","CaMKIII_g". If -sg without arugment then all groups are taken)
block-diagram -bd (block_diag )
 -- (-no-group -specific_group) Allowed 
 -- block-diagram is given preference if
    (-no-group, -block-diagram) or 
    (-specific_group, block-diagram) or
    (-specific_group, no-group,block-diagram) is passed 

Sep
box type just for sub is changed
fix to arrow color to co-ordinate the parent
2021
Mar 31 
sub->prd is connected with double arrow,
ligand with single arrow, 
inhibit with tee and 
modifier with diamond
legends and constant are added

Apr 7: group is added

Apr 8: eqns added with pluse and sigma 

May 4: HillTau API is called for reading json file

May 10:
added matplotlib for getting colors
validation of input file type is done
output image can be saved as png or svg

May 15: set function is remove to get Unique items
May 24: Margin for the cluster is increased from 8 to 22. 

May 31: line flags for eliminating the legend and for changing colors to bw.

June 1: added features for adjusting fontsize and height on command line

June 3: Group colors and node colors added

June 15: more option which are Optional 
	-ranksep'   : set rank separation (vertical spacing) in output.
	-group'     : Display group pass multiple group name with comma seperator
	-fontsize'  : set font size for node labels.
	-no_legend' : Turns off generation of legend'
	-bw'		: Goes into black-and-white plotting mode

Jun 19: Order of molecules
	#mol Kmod inhibit First-element second-element third-element
	 2	  0	     0 	     Input       Activator         --
	 2    0      1       Input       Inhibitor         --
	 3    1      0       Input       Modifier 		Activator
	 3    1      1       Input       Modifier       Inhibitor

Jun 30: with option -sg or --specific group, one can display specific group from the big model
python htgraph.py model.json -sg "group1_g","group2_g"
- If group name doesn't exist then it just ignores that specific group and display rest 
- If no group, specified in the command line exist then entire model is display like wise if no group is specified then
also entire model is displayed. 

'''

import sys,os
#sys.path.insert(1, 'PythonCode/')
from subprocess import call
import matplotlib
from collections import OrderedDict

from hillTau import *


use_bw = False

matplotcolors = []
for name,hexno in matplotlib.colors.cnames.items():
	matplotcolors.append(name)

# matplotcolors = ['#377eb8', '#ff7f00', '#4daf4a',
#                   '#f781bf', '#a65628', '#984ea3',
#                   '#999999', '#e41a1c', '#dede00']
def countX(lst, x):
	return lst.count(x)

def unique(list1):
	output = []
	for x in list1:
		if x not in output:
			output.append(x)
	return output
	#return list(set(list1))

def checkdigit_mol(sp):
	if sp.startswith(tuple('0123456789')):
		sp = "s"+sp
	
	return sp

def checkdigit(startstringdigit,grp,sp):
	if sp.startswith(tuple('0123456789')):
		if grp in startstringdigit:
			startstringdigit[grp][sp] = "s"+sp
		else:
			startstringdigit[grp] ={sp:"s"+sp}
			

def checkdigitEqu(startstringdigit,grp,sp):
	if grp in startstringdigit:
		grpitems = startstringdigit[grp]
		for k,v in grpitems.items():
			if k == sp:
				sp = v
	return sp

def getColor(gIndex):
	ignorecolors= ["chartreuse","peachpuff","paleturquoise","palegreen","olive","slategray","slategrey","lavenderblush","lemonchiffon","lightblue","lightcyan","lightgoldenrodyellow","lavender","khaki","seashell","gainsboro","burlywood","darkgrey","darkgray","palegoldenrod","linen","silver","darkkhaki","lightpink","mediumpurple","lightgreen","thistle","papayawhip","preachpuff","pink","tan","powderBlue","navajowhite","moccasin","mistyrose","lightgrey","lightgray","grey","gray","aquamarine","cadetblue","white","wheat","aqua","whitesmoke","mintcream","oldlace","black","snow","aliceblue","azure","cornsilk","beige","bisque","blanchedalmond","antiquewhite","lightyellow","lightsteelblue","ghostwhite","floralwhite","ivory","honeydew"];
	
	if use_bw:
		return( "black", gIndex )
	if gIndex < len(matplotcolors):
		grpcolor = matplotcolors[gIndex]
		gIndex = gIndex+1
		if grpcolor not in ignorecolors:
			return(grpcolor,gIndex)
		else:
			return getColor(gIndex)
	else:
		gIndex = 0
		grpcolor = matplotcolors[gIndex]
		gIndex = gIndex+1
		if grpcolor not in ignorecolors:

			return(grpcolor,gIndex)
		else:
			return getColor(gIndex)

def writeFuncedge(nIndex,ps,reac_edgelist,specieslist_rec,edge_weight,edge_arrowsize,fontsize,lig_exist,func_edgekey):
	
	color,nIndex = getColor(nIndex)
	pscolor = checkdigit_mol(ps[0])
	reaction_color = node_color[pscolor]
	arrowtype = "vee"
	
	ps0 = checkdigit_mol(ps[0])
	ps1 = checkdigit_mol(ps[1])

	#plusesize = "pluse"+str(equ_p_s)
	if ps[2].startswith('plus'):
		reac_edgelist = reac_edgelist+"\n"+ps[2]+"[label=\"+\",shape=circle,width=0, height=0, margin=0]"
	elif ps[2].startswith("sigma"):
		reac_edgelist = reac_edgelist+"\n"+ps[2]+"[label=<&Sigma;>,shape=circle,width=0, height=0, margin=0]"
	reac_edgelist =reac_edgelist+"\n"+ps0+"->"+ps[2]+"[arrowhead ="+str(arrowtype)+" weight = "+str(edge_weight)+ " minlen = 1 arrowsize = "+str(edge_arrowsize)+" color=\""+reaction_color+"\" fontsize="+str(fontsize)
	if ps[0] not in specieslist_rec:
		specieslist_rec.append(ps0)
	if ps[1] not in specieslist_rec:
		specieslist_rec.append(ps1)
	if ps[3] != 1:
		reac_edgelist =reac_edgelist+" label=\""+str(ps[3])+"\"]"
	else:
		reac_edgelist =reac_edgelist+"]"
	if ps[2] not in func_edgekey.keys():
		func_edgekey[ps[2]]=[ps1]
		reac_edgelist = reac_edgelist+"\n"+ps[2]+"->"+ps1+"[arrowhead ="+str(arrowtype)+" weight = "+str(edge_weight)+ " arrowsize = "+str(edge_arrowsize)+"]"
	
	return reac_edgelist,nIndex,lig_exist,func_edgekey	
def writeReacedge(nIndex,ps,reac_edgelist,specieslist_rec,edge_weight,edge_arrowsize,fontsize,input_exist,lig_exist,inhibit_exist,kmod_exist):

	con_arrow={"inhibit":"tee","Modifier":"odiamond","ligant":"vee","input":"normal"}
	
	pscolor = checkdigit_mol(ps[0])
	reaction_color = node_color[pscolor]
	if ps[2] in con_arrow:
	 	arrowtype = con_arrow[ps[2]]

	if ps[2] == "ligant":
		lig_exist = True
	elif ps[2] == "input":
		input_exist = True
	elif ps[2] == "inhibit":
		inhibit_exist = True
	elif ps[2] == "Modifier":
		kmod_exist = True
	ps0 = checkdigit_mol(ps[0])
	ps1 = checkdigit_mol(ps[1])

	reac_edgelist =reac_edgelist+"\n"+ps0+"->"+ps1+"[arrowhead ="+str(arrowtype)+" weight = "+str(edge_weight)+ " minlen = 1 arrowsize = "+str(edge_arrowsize)+" color=\""+reaction_color+"\" fontsize="+str(fontsize)
	if ps[0] not in specieslist_rec:
		specieslist_rec.append(ps0)
	if ps[1] not in specieslist_rec:
		specieslist_rec.append(ps1)
	if ps[3] != 1:
		reac_edgelist =reac_edgelist+" label=\""+str(ps[3])+"\"]"
	else:
		reac_edgelist =reac_edgelist+"]"
	return reac_edgelist,nIndex,input_exist,lig_exist,inhibit_exist,kmod_exist

def jsontoPng(modelpath, outputfile, ranksep = 0.1, hasLegend = True, fontsize = 18, showGroups = True,specific_group = [],show_blockdiagram = False):
	group_no = 0
	groupmap = OrderedDict()
	global startstringdigit
	startstringdigit= OrderedDict()
	global groupTogroup
	groupTogroup = OrderedDict()
	global node_color,input_exist,lig_exist,inhibit_exist,kmod_exist
	node_color = {}
	input_exist = True
	lig_exist = False
	kmod_exist = False
	inhibit_exist = False
	complex_exist = False
	edge_arrowsize = 1.5
	edge_weight = 1
	con_arrow={"inhibit":"tee","Modifier":"odiamond","ligant":"vee","input":"normal"}
	reac_edgelist = ""
	specieslist_rec =[]
	if show_blockdiagram:
		specific_group = []
		
	s = ""
	st = os.path.splitext(outputfile)
	outputfilename = st[0]
	if len( st ) > 1: 
		outputfiletype = st[1][1:]
	else:
		outputfiletype = "png"
	f_graph = open(outputfilename+".dot", "w")
	f_graph.write("digraph mygraph {\n\trank=TB; compound=true;\n")
	if ranksep > 0.0:
		f_graph.write("\trank = same, ranksep={};\n".format( ranksep ))
	f_graph.write("node [shape=box, penwidth=2,width=0, height=0, margin=0,fontsize={}];".format( fontsize ) )
	displayGroups = []
	if specific_group == None:
		displayGroups = modelpath.grpInfo
	else:
		if any(i in specific_group for i in modelpath.grpInfo):
			displayGroups = specific_group
		else:
			displayGroups = modelpath.grpInfo
		
	
	specieslist,node_color,spelist = writeSpecies(groupmap)
	funclist = writeFunc(groupmap)
	
	reaclist,node_color = writeReac(groupmap)
	nIndex = len(matplotcolors)-1
	species = (list(set(spelist) - set(reaclist+funclist)))
	grp_cluster_info = {}
	
	for grp,items in groupmap.items():
		if grp in displayGroups:
			color,nIndex = getColor(nIndex)
			if show_blockdiagram or showGroups:
				s = s + "\nsubgraph cluster_"+str(group_no)+"i\n{"
				s = s+"\nsubgraph cluster_"+str(group_no)+"\n{"+"\n"+"label=\""+grp+"\";\npenwidth=4; margin=10.0\ncolor=\""+color+"\";\nfontsize="+str(fontsize + 2)+";\n"			

			sps = ""
			items = list(unique(items))
			if not show_blockdiagram:
				for sp in items:
					if items.index(sp) != 0:
						sps = sps+','+sp
					else:
						sps = sps+sp
				s = s+sps+"\n"
					
			else:
				s = s+"grp"+str(group_no)+" [shape=point style=invis] "
				grp_cluster_info[grp] ={"grpObj":"grp"+str(group_no),"cluster":"cluster_"+str(group_no),"color":color}


			if show_blockdiagram or showGroups:
					s = s +"\n} style=invisible\n}style=rounded"

			group_no += 1;
	
	f_graph.write(s)

	if not show_blockdiagram:
		'''Blockdiagram is false'''
		''' group and output from group'''
		layedoutGroup = dict()
		func_edgekey = dict()
		for k,v in groupTogroup.items():
			if k in displayGroups:
				if k not in layedoutGroup:
					layedoutGroup[k]=[]
				for key,val in v.items():
					layedoutGroup[k].append(key)
					for n in val:
						if not n[2].startswith('pluse') and not n[2].startswith('sigma'):
							reac_edgelist,nIndex,input_exist,lig_exist,inhibit_exist,kmod_exist = writeReacedge(nIndex,n,reac_edgelist,specieslist_rec,edge_weight,edge_arrowsize,fontsize,input_exist,lig_exist,inhibit_exist,kmod_exist)
						else:
							reac_edgelist,nIndex,lig_exist,func_edgekey = writeFuncedge(nIndex,n,reac_edgelist,specieslist_rec,edge_weight,edge_arrowsize,fontsize,lig_exist,func_edgekey)
							
		if specific_group != None:
			''' when specific group is asked to display, info of input's coming from differnt
			    Group need to identified'''
			for spg in specific_group:
				for k,v in groupTogroup.items():
					if spg in list(v.keys()):
						if k not in layedoutGroup:
							for ps in v[spg]:
								if not ps[2].startswith('pluse') and not ps[2].startswith('sigma'):
									reac_edgelist,nIndex,input_exist,lig_exist,inhibit_exist,kmod_exist = writeReacedge(nIndex,ps,reac_edgelist,specieslist_rec,edge_weight,edge_arrowsize,fontsize,input_exist,lig_exist,inhibit_exist,kmod_exist)
								else:
									reac_edgelist,nIndex,lig_exist,func_edgekey = writeFuncedge(nIndex,ps,reac_edgelist,specieslist_rec,edge_weight,edge_arrowsize,fontsize,lig_exist,func_edgekey)
						
		f_graph.write(reac_edgelist)
	else:
		''' Display Block diagram '''
		test1 = dict()
		for k,v in groupTogroup.items():
			test = dict()		
			for kv,l in v.items():
				for y in l:
					if k != kv:
						if kv in test:
							if (y[2].startswith("pluse") or y[2].startswith("sigma")):
								test[kv] = ["ligant"]
							else:
								test[kv].append(y[2])
						else:
							if (y[2].startswith("pluse") or y[2].startswith("sigma")):
								test[kv] = ["ligant"]
							else:
								test[kv] = [y[2]]
						if len(test) != 0:
							test1[k] = [test]

		for gk,gv in test1.items():

			for vl in gv:
				for vlk,vll in vl.items():
					reac_edgelist =reac_edgelist+"\n"+grp_cluster_info[gk]["grpObj"]+"->"+grp_cluster_info[vlk]["grpObj"]+"[arrowhead ="
					if len(unique(vll)) == 1:
						
						arrowtype = con_arrow[vll[0]]
						if vll[0] == "input":
							input_exist = True
						if vll[0] == "ligant":
							lig_exist = True
						elif vll[0] == "inhibit":
							inhibit_exist = True
						elif vll[0] == "Modifier":
							kmod_exist = True
						kg_color = grp_cluster_info[gk]["color"]
					else:
						arrowtype = "dot"
						kg_color = "black"
						complex_exist = True
						
					reac_edgelist = reac_edgelist+str(arrowtype)+" weight = "+str(edge_weight)+ " minlen = 1 arrowsize = "+str(edge_arrowsize)+" color=\""+kg_color+"\" fontsize="+str(fontsize)+ " ltail="+grp_cluster_info[gk]["cluster"]+",lhead="+grp_cluster_info[vlk]["cluster"]+"]"
		f_graph.write(reac_edgelist)		
		
	if not show_blockdiagram:
		
		splist =specieslist.values()
		for le in specieslist_rec:
			concValue = 0
			le = checkdigit_mol(le)
			if le != 'func':
				v = node_color[le]
				for s_splist in splist:
					for t in s_splist:
						if t['name'] == le:
							concValue = t['value']
				if le in species:
					f_graph.write("\n"+le+"[color=\""+v+"\",shape=Mrecord, width=0, height=0, margin=0.1,tooltip = \"concInit = "+str(float("{:.5f}".format(concValue)))+"\"]")
				else:
					f_graph.write("\n"+le+"[color=\""+v+"\", width=0, height=0, margin=0.1, tooltip = \"concInit = "+str(float("{:.5f}".format(concValue)))+"\"]")	

		for p,q in startstringdigit.items():
			if p in displayGroups:
				for m,n in q.items():
					f_graph.write("\n"+n+"[label=\""+m+"\"]")	
	
	if hasLegend:
		start = ""
		center = ""
		end = ""
		f_graph.write("\nnode [shape=plaintext]\nsubgraph cluster_01 {\n\tlabel = \"Legend\";\n\t{ rank=sink;\n\tkey [label=<<table border=\"0\" cellpadding=\"2\" cellspacing=\"0\" cellborder=\"0\">\n")
		if input_exist:
			start = start+"\t<tr><td align=\"right\" port=\"i1\">Input</td></tr>\n"
			center = center+"<tr><td port=\"i1\">&nbsp;</td></tr>\n"
			end = end+"\tkey:i1 -> key2:i1 [arrowhead=normal color=\"black:black\" style=bold]\n"
		if lig_exist:
			start = start+"\t<tr><td align=\"right\" port=\"i2\">ligand-Act</td></tr>\n"
			center = center+"\t<tr><td port=\"i2\">&nbsp;</td></tr>\n"
			end = end+"\tkey:i2 -> key2:i2 [arrowhead=vee]\n"
		if kmod_exist:
			start = start+"\t<tr><td align=\"right\" port=\"i3\">Modifier</td></tr>\n"
			center = center+"\t<tr><td port=\"i3\">&nbsp;</td></tr>\n"
			end = end+"\tkey:i3 -> key2:i3 [arrowhead=odiamond]\n"
		if inhibit_exist:
			start = start+"\t<tr><td align=\"right\" port=\"i4\">ligand-Inh</td></tr>\n"
			center = center+"\t<tr><td port=\"i4\">&nbsp;</td></tr>\n"
			end = end+"\tkey:i4 -> key2:i4 [arrowhead=tee]\n"
		if complex_exist:
			start = start+"\t<tr><td align=\"right\" port=\"i5\">Complex</td></tr>\n"
			center = center+"\t<tr><td port=\"i5\">&nbsp;</td></tr>\n"
			end = end+"\tkey:i5 -> key2:i5 [arrowhead=dot]\n"
		f_graph.write(start)
		f_graph.write("\t</table>>]\n\tkey2 [label=<<table border=\"0\" cellpadding=\"2\" cellspacing=\"0\" cellborder=\"0\">\n\t")
		f_graph.write(center)
		f_graph.write("\t</table>>]\n")
		f_graph.write(end)
		f_graph.write("\t}\n\t}")

	f_graph.write("\n}")
	f_graph.close()
	
	command = "dot -T"+ outputfiletype + " "+ outputfilename+".dot -o "+outputfile
	call([command], shell=True)
	print("file written ",outputfile)
	return(outputfile)
def writeSpecies(groupmap):
	# getting all the species
	specieslist = {}
	mIndex = 0
	spelist = []
	for molname, mol in ( modelpath.molInfo.items() ):
		checkdigit(startstringdigit,mol.grp,molname)
		molname = checkdigitEqu(startstringdigit,mol.grp,molname)
		if molname not in node_color:
			spe_color,mIndex = getColor(mIndex)
			node_color[molname] = spe_color
		spelist.append(molname)
		if molname != 'concInit':
			if mol.grp in groupmap:
				groupmap[mol.grp].append(molname)
				specieslist[mol.grp].append({"name":molname,"value":mol.concInit})
			else:
				groupmap[mol.grp] = [molname]
				specieslist[mol.grp] = [{"name":molname,"value":mol.concInit}]
	return specieslist,node_color,spelist
		
def writeFunc(groupmap):
	equation_pluse = 0
	equation_sigma = 0
	funclist = []
	for e,t in modelpath.eqnInfo.items():
		checkdigit(startstringdigit,t.grp,t.name)
		t.name = checkdigitEqu(startstringdigit,t.grp,t.name)
		allpluse = True
		mathSym = []
		funclist.append(t.name)
		for i in t.eqnStr:
			if i in ["*","-","/","+"]:
				mathSym.append(i)
		if len(unique(mathSym)) == 1 and mathSym[0] == "+":
			allpluse = True
		else:
			allpluse = False
		if allpluse:
			plusesize = "pluse"+str(equation_pluse)
			equation_pluse+=1
			#groupmap[t.grp].append(plusesize)
			
		else:
			plusesize = "sigma"+str(equation_sigma)
			equation_sigma+=1
			#groupmap[t.grp].append(plusesize)
		for tsubs in unique(t.subs):
			c = countX(t.subs,tsubs)
			if tsubs != 'concInit':
				if allpluse:
					groupmap[t.grp].append(plusesize)
				else:
					groupmap[t.grp].append(plusesize)
				subgroup = [key for key,value in groupmap.items() if tsubs in value][0]
				if subgroup not in groupTogroup:
					groupTogroup[subgroup] = {t.grp:[((tsubs,t.name,plusesize,c))]}		
				else:
					'''  
						A: {A:((list of values)())
							B : ((list of values))
							}
					'''
					if t.grp in groupTogroup[subgroup].keys():
						groupTogroup[subgroup][t.grp].append(((tsubs,t.name,plusesize,c)))
					else:
						''' under subgroup, reac.grp exist then update the connect to the same list
							A :{A:((list of values),( list of values))}
						'''
						groupTogroup[subgroup].update({t.grp:[((tsubs,t.name,plusesize,c))]})
	
	return funclist

def writeReac(groupmap):
	sIndex = 0
	reaclist = []
	#groupTogroup = OrderedDict()
	#groupTogroup1 = OrderedDict()
	for reacname, reac in ( modelpath.reacInfo.items() ):
		checkdigit(startstringdigit,reac.grp,reacname)
		reacname = checkdigitEqu(startstringdigit,reac.grp,reacname)
		reaclist.append(reacname)
		if reac.grp in groupmap:
			groupmap[reac.grp].append(reacname)
		else:
			groupmap[reac.grp] = [reacname]

		sublist = reac.subs
		sublistU = unique(reac.subs)
		prd = reacname
		for sub in sublistU:
			newsub = sub
			subS =checkdigit_mol(sub)
			subgroup = [key for key,value in groupmap.items() if subS in value][0]
			if reac.grp in startstringdigit:
				t = startstringdigit[reac.grp]
				if newsub in t:
					newsub = t[newsub]
			''' if string starting with number, then replace with s+string'''
			if newsub in node_color:
				reaction_color = node_color[newsub]
			else:
				reaction_color,sIndex = getColor(sIndex)
				node_color[newsub] = reaction_color
			checkdigit(startstringdigit,reac.grp,sub)
			c = countX(sublist,sub)
			#sub = checkdigitEqu(startstringdigit,reac.grp,sub)
			if (reac.inhibit == 1.0 and sublistU.index(sub) == len(sublistU)-1 ) :
				''' inhibit  ligant activator tee  '''
				ctype = "inhibit"
			elif len(sublistU) == 3 and sublist.index(sub) == 1:
				''' kmod Modulator odiamond '''
				ctype = "Modifier"
			else:
				if sublist.index(sub) >= 1:
					''' ligand vee '''
					ctype = "ligant"						
				else:
					if  sublist.index(sub) == 0:
						''' input '''
						ctype = "input"
			'''Creating dictionary with
				key : {key:((value)(value))
					   }
				A : {A :((list of values )(list of values))
					 B : ((list of values))
					}
			'''
			if subgroup not in groupTogroup:
				groupTogroup[subgroup] = {reac.grp:[((sub,prd,ctype,c))]}		
			else:
				'''  
					A: {A:((list of values)())
						B : ((list of values))
						}
				'''
				if reac.grp in groupTogroup[subgroup].keys():
					groupTogroup[subgroup][reac.grp].append(((sub,prd,ctype,c)))
				else:
					''' under subgroup, reac.grp exist then update the connect to the same list
						A :{A:((list of values),( list of values))}
					'''
					groupTogroup[subgroup].update({reac.grp:[((sub,prd,ctype,c))]})
				
	return(reaclist,node_color)
		

def file_choices(choices,fname,iotype):
	ext = (os.path.splitext(fname)[1][1:]).lower()
	if iotype == "imagetype":
		if fname not in choices:
			parser.error("Requires output filetype {}".format(choices))
	# elif iotype == "outputfile":
	# 	if ext not in choices:
	# 		parser.error("Requires output filetype {}".format(choices))
	else:
		if ext != "json":
			parser.error("Requires HillTau file in JSON format ")
			
	return fname

if __name__ == "__main__":

	parser = argparse.ArgumentParser( description = 'This program generates a reaction diagram for a HillTau model. It converts the specified HillTau file in JSON format, to the dot format. The dot file is further converted to an image in png/svg format\n')
	parser.add_argument('model',type=lambda s:file_choices(("json"),s,"input"), help='Required: filename of model, in JSON format.')
	#parser.add_argument( 'model', type = str, help='Required: filename of model, in JSON format.')
	parser.add_argument( '-o', '--output', type = str,help='Optional: writes out the png model into named file. default takes json filename')
	
	parser.add_argument( '-T', '--type', type=lambda imagetype:file_choices(("png","svg"),imagetype,"imagetype"), help='Optional: writes out the image in png or svg is format, default takes png format')
	parser.add_argument( '-r', '--ranksep', type=float, default = 0, help='Optional: set rank separation (vertical spacing) in output.')
	parser.add_argument( '-fs', '--fontsize', type=float, default = 18, help='Optional: set font size for node labels.')
	parser.add_argument( '-nl', '--no_legend', action='store_true', help='Optional: Turns off generation of legend')
	parser.add_argument( '-ng', '--no_groups', action='store_true', help='Optional: Removes grouping. All molecules and reactions sit together.')
	parser.add_argument( '-bw', '--bw', action='store_true', help='Optional: Goes into black-and-white plotting mode')
	parser.add_argument('-sg', '--specific_group', help='Optional: Specfiy group names for display,delimited groupname seprated by comma.',type=lambda s:s.split(","))
	parser.add_argument( '-bd', '--block-diagram', action='store_true', help='Optional: Block diagram of the model with connection are retained.')
	
	args = parser.parse_args()
	use_bw = args.bw
	if args.type != None:
		exttype = args.type

	if args.output == None:
		dirpath = os.path.dirname(args.model)
		basename = os.path.basename(args.model)
		if args.type == None:
			exttype = "png"
		if dirpath:
			outputfile = os.path.dirname(args.model)+"/"+os.path.splitext(os.path.basename(args.model))[0]+"."+exttype	
		else:
			outputfile = os.path.splitext(args.model)[0]+"."+exttype
	else:
		basename = os.path.basename(args.model)
		basestr = os.path.splitext(basename)
		outst = args.output
		dirpath = os.path.dirname(args.output)
		basename1 = os.path.basename(args.output)
		if basename1 != "":
			st = os.path.splitext(basename1)
			if len(st[0]) == 1:	
				if not st[0][-1].isalpha():
					outputfilename = st[0][0:len(st[0])-1]
					if len(outputfilename) ==0:
						outputfilename = basestr[0]
				else:
					outputfilename = st[0]
			else:
				outputfilename = st[0]
			#type	
			if len( st ) > 1:
				if (st[1][1:] != "" and st[1][1:] in ('png','svg')):
					if args.type == None:
						exttype = st[1][1:]
				else:
					if st[1][1:] == "":
						exttype = "png"
		else:
			outputfilename = os.path.splitext(os.path.basename(args.model))[0]
			if args.type == None:
				exttype = "png"

		if dirpath != '/':	
			outputfile = dirpath+"/"+outputfilename+"."+exttype
		else:
			outputfile = os.path.splitext(args.model)[0]+"."+exttype
	
	jsonDict = loadHillTau( args.model )
	modelpath = parseModel( jsonDict )
	jsontoPng(modelpath, outputfile, ranksep = args.ranksep, hasLegend = not args.no_legend, fontsize = args.fontsize, showGroups = not args.no_groups,specific_group = args.specific_group,show_blockdiagram = args.block_diagram )
	#return (outputfile)
