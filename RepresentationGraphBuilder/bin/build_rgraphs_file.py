################################################################################
# build_rgraphs_file.py
# Python script for building representation graphs for HSM model state space
# Usage: 
# build_rgraphs_file.py -i <source file> -o <destination file (graph)> [-draw <destination file (graphml)>] [-ini] [-s]
# Requirements:
# Python 3.8.1 or newer
# Dependencies:
# numpy-1.24.4, netwrokx-3.1, igraph-0.11.4, N2G-0.3.3
# Contributors: 
# Institute of Mathematics and Computer Science, University of Latvia
# v_1.0.15, 15.03.2024
# Distributed under GPLv3 license
# Copyright (c) 2024 Juris Viksna
################################################################################

import sys
import argparse
import networkx as nx
import numpy as np
from N2G import yed_diagram

################################################################################
# check and assign the input arguments
################################################################################

parser = argparse.ArgumentParser(description="Computes representation graph for HSM model state space")
parser.add_argument("-i",dest="fname_in",help="input state space file",metavar="<input file>",required=True)
parser.add_argument("-o",dest="fname_outgraph",help="output representation graph file",metavar="<output file>",required=False)
parser.add_argument("-draw",dest="fname_outdraw",help="output graphml file",metavar="<graphml file>",required=False)
parser.add_argument("-ini",dest="ini",action="store_true",help="analyse only part reachable from INI state",required=False)
parser.add_argument("-s",dest="s",action="store_true",help="run in silent mode",required=False)
args = vars(parser.parse_args())

fname_in = args['fname_in']
fname_outgraph = "out.txt" # default
if args['fname_outgraph'] != None:
	fname_outgraph = args['fname_outgraph']
fname_outdraw = "out.graphml" # default
if args['fname_outdraw'] != None:
	fname_outdraw = args['fname_outdraw']

if not args['s']:
	print("Input file: "+fname_in)
	print("Output graph file: "+fname_outgraph)
	if args['fname_outdraw']:
			print("Output graphml file: "+fname_outdraw)
	if args['ini']:
		print("Running with 'ini' flag - analyse only part of space reachable from INI state\n")


# assign filre names from provided arguments

################################################################################
# Read the input file and save it in simple data structures
################################################################################

# grab the file
inf = open(fname_in, "r")
f_contents = ""
while True:
	f_line = inf.readline()
	if not f_line:
		break
	strip_line = f_line.strip()
	# ignore lines marked as comments by '#'
	if strip_line != "" and strip_line[0] != '#':  
		f_contents = f_contents + " " + strip_line
inf.close()

# parse contents
# read header 
l_contents = str.split(f_contents)
r_idx = 0
N_genes = int(l_contents[r_idx]) # number of genes
r_idx += 1
N_bsites = int(l_contents[r_idx]) # number of bsites
r_idx += 1
N_states  = int(l_contents[r_idx]) # number of states
r_idx += 1
S_ini = int(l_contents[r_idx]) # the initial 0 state
r_idx += 1
Gnames = l_contents[r_idx:r_idx+N_genes] # list of gene names
r_idx += N_genes
BSnames = l_contents[r_idx:r_idx+N_bsites] # list of bsite names
r_idx += N_bsites

# read list of states and transitions

class State:
	def __init__(self,id,genes,bs,type,en,eList):
		self.id = id
		self.name = genes
		self.bs = bs
		self.type = type
		self.en = en # number of outgoing edges
		self.edges = eList
		self.ini = -1 # set to 1 for INI state
		self.cini = -1 # set to 1 for INI state and states that get contracted from it
		self.inireach = -1 # set to 1 for vertices reachable from INI state
		self.grey = -1 # set to 1 for vertices that are no attr of dec regions, but will be kept in representation graph as "other"
		self.attr = -1 # set to attractor id vertex (>=0) if state belongs to such
		self.dec = -1 # set to decision state group id for decision states (as defined by dcond)
		self.reach = [] # list of reachability states
		self.dcond = set() # set object of decision conditions
		self.contr = -1 # set to 1 if deleted by contraction
		self.discard = -1 # set to 1 if discarded as unreachable from INI state
		
N_trans = 0 # number of transitions
SGraph = [] # state graph (raw)
IDs = {} # dictionary mapping 0...N_states range to provided state ids  
rIDs = {} # reverse map
	
for i in range(N_states):
	# read state data
	r_id = int(l_contents[r_idx])
	r_gene = l_contents[r_idx+1]
	r_bs = l_contents[r_idx+2]
	r_type = int(l_contents[r_idx+3])
	r_en = int(l_contents[r_idx+4])
	r_idx += 5
	eList = []
	IDs[i] = r_id # assign id to dictionary
	rIDs[r_id] = i
	# read transition data
	for j in range(r_en):
		r_tdid = int(l_contents[r_idx]) # destination state id
		r_tgid = int(l_contents[r_idx+1]) # label gene id
		r_ttype = l_contents[r_idx+2] # label gene type (+/-)		
		r_idx += 3
		eList.append([r_tdid,r_tgid,r_ttype]) 
	N_trans += r_en
	nState= State(i,r_gene,r_bs,r_type,r_en,eList)
	if r_id == S_ini: # mark INI state
		nState.ini = 1
		nState.cini = 1
	SGraph.append(nState)	

# re-traverse graph and remap edge destination id-s
for i in range(N_states):
	for j in range(SGraph[i].en):
		SGraph[i].edges[j][0] = rIDs[SGraph[i].edges[j][0]]

################################################################################
# Transfer to NetworkX digraph objects
################################################################################
		
# built NetworkX diGraph object for SGraph 
GX = nx.DiGraph()
for i in range(N_states):
	GX.add_node(i,id=i) # preserve reference to the original id
for i in range(N_states):
	State = SGraph[i]
	en = State.en
	edges = State.edges
	for j in range(en):
		edge = edges[j]
		GX.add_edge(i,edge[0],l=edge[1],t=edge[2]) # label with triggering gene name and type 

# check the reachability for INI state, if analysis only of part reachable from INI state is asked, delete unreachable nodes 
reach_list = list(nx.descendants(GX,S_ini))
reach_list.append(S_ini)
for i in range(N_states):
	if i in reach_list:
		SGraph[i].inireach = 1 # mnark as reachable from INI state
	if args['ini']: # discard node if the corresponding attribute flag is set
		if i not in reach_list:
			SGraph[i].discard = 1 # set as discarded
			GX.remove_node(i) # remove from graph

###############################################################################
# Define graph output functions
################################################################################

# write attractor
def draw_attr(numb,vlist,outfg):
	outfg.write("### Attractor region "+str(numb+1)+" \n\n") # to be consistent with yEd number attractors from 1, not 0
	outfg.write("# Attractor id | int\n")
	vr = SGraph[attr[0]].attr # define as representative vertex
	outfg.write(str(vr)+"\n\n")	# identify attractor with a randomly selected vertex from SCC
	outfg.write("# Number of states | int\n")
	outfg.write(str(len(vlist))+"\n\n")	
	outfg.write("# List of states: id genes binding_sites | int str str\n")	
	for vert in vlist:
		outfg.write(str(IDs[vert])+" "+SGraph[vert].name+" "+SGraph[vert].bs+" "+"\n")
	outfg.write("\n")
	outfg.write("# List of transitions: source_id dest_id gene_id gene_name {+,-} | int int int str def\n")	
	for vert in vlist:
		en = SGraph[vert].en
		for j in range(en):
			dest_l = SGraph[vert].edges[j]
			dest_v = dest_l[0]
			if dest_v in vlist and dest_v != vert: # edge is within attractor, discard self-loops
				outfg.write(str(IDs[vert])+" "+str(IDs[dest_v])+" "+str(dest_l[1])+" "+Gnames[dest_l[1]]+" "+dest_l[2]+"\n")
	outfg.write("\n")
	
################################################################################
# Compute the components of interest
################################################################################

# find SCCs		
SCC = list(nx.strongly_connected_components(GX))

# locate non-transitional SCCs and copy them to a list (relying on input labels here) 
SCC_list = []
for i in range(len(SCC)):
	clist = list(SCC[i])
	if SGraph[clist[0]].type != 0:
		SCC_list.append(clist)

#print(SCC_list)

# set attractor attributes in SGraph
for i in range (len(SCC_list)):
	attr = SCC_list[i]
	for j in range (len(attr)):
		SGraph[attr[j]].attr = attr[0]
	
# compute attractor reachability
for i in range (len(SCC_list)):
	attr = SCC_list[i]
	vert = attr[0]
	r_list = list(nx.ancestors(GX,vert))
	r_list.append(vert) # include start vertex just in case (could be a duplicate)
	# assign reachability values in SGraph	
	for j in range(len(r_list)):
		SGraph[r_list[j]].reach.append(vert)

# find decision nodes and mark them with decision conditions	
for i in range(N_states):
	if SGraph[i].attr == -1 and SGraph[i].discard == -1: # skip attractors and discarded states
		nb_list = list(GX.neighbors(i))
		#print(i,nb_list)
		for j in range(len(nb_list)):
			nb = nb_list[j]
			if set(SGraph[i].reach) != set(SGraph[nb].reach):
				SGraph[i].dec = 0
				label = GX[i][nb]['l']
				type = GX[i][nb]['t']
				SGraph[i].dcond.add((label,type))
				
# partition decision nodes according to decision conditions
DConds_set = set()
for i in range(N_states):
	if SGraph[i].dec != -1: # decision node
		DConds_set.add(tuple(SGraph[i].dcond))
DConds_list = list(DConds_set)
for i in range(N_states):
	if SGraph[i].dec != -1: # decision node
		SGraph[i].dec = DConds_list.index(tuple(SGraph[i].dcond))
#print(len(DConds_list))
#print(DConds_list)

################################################################################
# Proceed with contractions in NetworkX object according to component structure
################################################################################


# at start delete all self-loops for simplicity (will remain in SGraph and will be used for attractors from here)
Lloop_nodes = list(GX.nodes) # mark the full node list 
for v in Lloop_nodes:
	if GX.has_edge(v,v):
		GX.remove_edge(v,v)

# contract attractors to single nodes
# do not bother withg attributes here - in attractor region should be the same for all vertices
for i in range(len(SCC_list)): 
	attr = SCC_list[i]
	# find all pairs of nodes within attractor (this could be done more efficiently, but attractors are really small)
	a_pairs = []
	for v1 in attr:
		for v2 in attr:
			a_pairs.append((v1,v2))
			a_pairs.append((v2,v1))
	# go through the edge list and contract these that are still available
	for (s,d) in a_pairs:
		if (s,d) in GX.edges: # if still available, then contract
			if SGraph[s].cini == 1 or SGraph[d].cini == 1: # mark as contraction from INI state
				SGraph[s].cini = 1
				SGraph[d].cini = 1	
			SGraph[d].contr = 1
			nx.contracted_edge(GX,(s,d),self_loops=False, copy=False) # contract d to s 
			#print("contr: ",v,d,IDs[v],IDs[d]," aa")
	# go through the edge list and contract these that are still available (repeat 2x - a quick fix instead fo recomputing and adding all newly created edges)
	for (s,d) in a_pairs:
		if (s,d) in GX.edges: # if still available, then contract
			if SGraph[s].cini == 1 or SGraph[d].cini == 1: # mark as contraction from INI state
				SGraph[s].cini = 1
				SGraph[d].cini = 1	
			SGraph[d].contr = 1
			nx.contracted_edge(GX,(s,d),self_loops=False, copy=False) # contract d to s 
			#print("contr: ",v,d,IDs[v],IDs[d]," aa")
			
# the part below is intended as O(m) iteration though all edges, and adjusting this edge set accordnig to contractions
# until speed becomes the issue the implementation uses simplified approach with re-setting the whole edge set after each iteration and rechecking the edges again
not_completed = True
while not_completed:
	not_completed = False # reset to false, will be set back to true if contractions will be performed
	L_edges = list(GX.edges)
	#print(L_edges)
	for e in L_edges:
		(s,d) = e
		# the first case - edge between two boring vertices, always contract
		if e in GX.edges and SGraph[s].attr == -1 and SGraph[s].dec == -1 and SGraph[d].attr == -1 and SGraph[d].dec == -1: #check that boring and that edge is still available
			if SGraph[s].cini == 1 or SGraph[d].cini == 1: # mark as contraction from INI state
				SGraph[s].cini = 1
				SGraph[d].cini = 1	
			SGraph[d].contr = 1
			nx.contracted_edge(GX,(s,d),self_loops=False, copy=False) # contract d to s 
			#print("contr: ",IDs[s],IDs[d]," bb")
			not_completed = True
		# the second case - edge between two identical switching states, always contract	
		if e in GX.edges and SGraph[s].dec != -1 and SGraph[d].dec != -1: # check  that between decision states
			if SGraph[s].reach == SGraph[d].reach: # reachability regions match
				if SGraph[s].dcond <= SGraph[d].dcond: # decision state of destination is a superset of decision set of source
					if SGraph[s].cini == 1 or SGraph[d].cini == 1: # mark as contraction from INI state
						SGraph[s].cini = 1
						SGraph[d].cini = 1	
					SGraph[s].contr = 1
					nx.contracted_nodes(GX,d,s,self_loops=False, copy=False) # contract s to d
					# nx.contracted_edge(GX,(s,d),self_loops=False, copy=False) # contract d to s 
					#print("contr: ",IDs[d],IDs[s]," cc")
					not_completed = True
		# the third case - boring state to attractor
		if e in GX.edges and SGraph[s].attr == -1 and SGraph[s].dec == -1 and SGraph[d].attr != -1: # check that boring to attractors and that edge is still available
			pr_list = list(GX.predecessors(s))
			if (len(pr_list) == 0): # possibly a single vertex region, keep it as grey
					SGraph[s].grey = 1 
			else: # otherwise contract
				if SGraph[s].cini == 1 or SGraph[d].cini == 1: # mark as contraction from INI state
					SGraph[s].cini = 1
					SGraph[d].cini = 1	
				SGraph[s].contr = 1
				nx.contracted_nodes(GX,d,s,self_loops=False, copy=False) # contract s to d - keep attractor intact
				#print("contr: ",IDs[d],IDs[s]," dd")
				not_completed = True
		# the fourth case - boring state to switching state
		if e in GX.edges and SGraph[s].attr == -1 and SGraph[s].dec == -1 and SGraph[d].dec != -1: # check that boring to switching and that edge is still available
			pr_list = list(GX.predecessors(s))
			if (len(pr_list) == 0): # possibly a single vertex region, keep it as grey
					SGraph[s].grey = 1 
			else: # otherwise contract
				if SGraph[s].cini == 1 or SGraph[d].cini == 1: # mark as contraction from INI state
					SGraph[s].cini = 1
					SGraph[d].cini = 1	
				SGraph[s].contr = 1
				nx.contracted_nodes(GX,d,s,self_loops=False, copy=False) # contract s to d - keep swiching state intact
				#print("contr: ",IDs[d],IDs[s]," ee")
				not_completed = True

				

# reprocess grey nodes and join into a single node groups with the same destination vertices
Lall_edges = list(GX.edges) # mark the full node list 
for (s,d) in Lall_edges:
	if (s,d) in GX.edges and SGraph[s].grey == 1: # edge is still available and source is grey edge
		pr_list = list(GX.predecessors(d)) # predecessors of destination vertex
		for p in pr_list:
			if p!=s and p in GX.nodes and SGraph[p].grey == 1 and set(GX.neighbors(s)) == set(GX.neighbors(p)): # a distinct predecessor that is also grey and with the same reachability modes
				if GX[s][d]['l'] == GX[p][d]['l'] and GX[s][d]['t'] == GX[p][d]['t']: # and transition labels to d are also the same, then contract
					if SGraph[s].cini == 1 or SGraph[d].cini == 1: # mark as contraction from INI state
						SGraph[s].cini = 1
						SGraph[d].cini = 1	
					SGraph[p].contr = 1	
					nx.contracted_nodes(GX,s,p,self_loops=False, copy=False) # contract p to s 
					#print("join: ",IDs[s],IDs[p]," ff")

				
################################################################################
# Some test printouts - uncomment as needed
################################################################################	
			
"""			
print(N_genes)
print(N_bsites)
print(N_states)
print(S_ini)
print(N_trans)
print(Gnames)
print(BSnames)
"""

def vert_print(i):
	print(i)
	print("id",SGraph[i].id)
	print("name",SGraph[i].name)
	print("bs",SGraph[i].bs)
	print("type",SGraph[i].type)
	print("en",SGraph[i].en)
	print("edges",SGraph[i].edges)
	print("ini",SGraph[i].ini)
	print("cini",SGraph[i].cini)
	print("inireach",SGraph[i].inireach)
	print("grey",SGraph[i].grey)
	print("attr",SGraph[i].attr)
	print("dec",SGraph[i].dec)
	print("reach",SGraph[i].reach)
	print("dcond",SGraph[i].dcond)
	print("contr",SGraph[i].contr)
	print("discard",SGraph[i].discard)
	print("")

"""
for i in range(N_states):
	#if i in GX.nodes:
	print(IDs[i])
	vert_print(i)
"""	
#print(len(list(GX.nodes)))
#print(GX.nodes)
#print(len(list(GX.edges)))
#print(GX.edges)

################################################################################
# Write results in files
################################################################################

# print representation graph if .graphmlx format, if requested by call arguments 
if args['fname_outdraw']:
	colour_list = ['#45CE2A','#A2ADD0','#FB7EFD','#1974D2','#FFA089','#DA70D6','#EE82EE','#B2EC5D','#C364C5','#17806d']
	diagram = yed_diagram()
	# create representation graph part
	# RG vertices
	for v in GX.nodes:
		vid = IDs[v]
		if SGraph[v].attr != -1: # add attractor
			attr_colour = colour_list[SGraph[v].attr % len(colour_list)]  # rotate through the defined attractor colour list
			if SGraph[v].cini != - 1:  # state includes INI position
				attr_color = '#FADA5E'
			diagram.add_node(str(vid),label=str(vid),attributes={'Shape': {'type': 'roundrectangle'}, 'NodeLabel' : {'alignment' : 'center', 'autoSizePolicy' : 'content', 'fontFamily' : 'Dialog', 'fontSize' : '20', 'fontStyle' : 'bold'}, 'Fill' : {'color': attr_colour}, 'BorderStyle' : {'width' : '3.0'}}) 
		if SGraph[v].dec != -1: # add switching condition
			dec_colour = '#D80073' #default
			if SGraph[v].cini != - 1:  # state includes INI position
				dec_colour = '#FADA5E'
			diagram.add_node(str(vid),label=str(vid),attributes={'Shape': {'type': 'diamond'}, 'NodeLabel' : {'alignment' : 'center', 'autoSizePolicy' : 'content', 'fontFamily' : 'Dialog', 'fontSize' : '20', 'fontStyle' : 'bold'}, 'Fill' : {'color': dec_colour}, 'BorderStyle' : {'width' : '3.0'}}) 
		if SGraph[v].grey != -1: # add grey
			grey_colour = '#E0E0E0' #default
			if SGraph[v].cini != - 1:  # state includes INI position
				grey_colour = '#FADA5E'
			diagram.add_node(str(vid),label=str(vid),attributes={'Shape': {'type': 'rectangle'}, 'NodeLabel' : {'alignment' : 'center', 'autoSizePolicy' : 'content', 'fontFamily' : 'Dialog', 'fontSize' : '20', 'fontStyle' : 'bold'}, 'Fill' : {'color': grey_colour}, 'BorderStyle' : {'width' : '1.0'}}) 
	# RG edges
	for (s,d) in GX.edges:
		sid = IDs[s]
		did = IDs[d]
		label_str = Gnames[GX[s][d]['l']]+str(GX[s][d]['t'])
		diagram.add_link(str(sid),str(did),label_str,attributes={'Arrows' : {'source' : 'none', 'target' : 'standard'}, 'Label' : {'alignment' : 'center', 'fontFamily' : 'Dialog', 'fontSize' : '16'}, 'LineStyle' : {'width' : '2.0'}})
	# create attractor parts 
	for i in range(len(SCC_list)): 
		attr = SCC_list[i]
		vr = SGraph[attr[0]].attr # define as representative vertex
		attr_colour = colour_list[SGraph[vr].attr % len(colour_list)]
		# attractor vertices
		for v in attr:
			vid = IDs[v]
			node_id = str(i+1)+":"+str(vid) # need different id from representation graph vertices
			diagram.add_node(node_id,label=node_id,attributes={'Shape': {'type': 'rectangle'}, 'NodeLabel' : {'alignment' : 'center', 'autoSizePolicy' : 'node_width', 'fontFamily' : 'Dialog', 'fontSize' : '16'}, 'Fill' : {'color': attr_colour}, 'BorderStyle' : {'width' : '2.0'}}) 
		# attractor edges
		for v in attr:
			vid = IDs[v]
			en = SGraph[v].en
			for j in range(en):
				dest_l = SGraph[v].edges[j]
				dv = dest_l[0]
				dvid = IDs[dv]
				dl = dest_l[1]
				dt = dest_l[2]
				label_str = Gnames[dl]+dt
				if dv in attr and dv != v: # edge is within attractor, discard self-loops
					dvid = IDs[dv]
					from_id = str(i+1)+":"+str(vid) # need different id from representation graph vertices
					to_id = str(i+1)+":"+str(dvid) # need different id from representation graph vertices
					diagram.add_link(from_id,to_id,label_str,attributes={'Arrows' : {'source' : 'none', 'target' : 'standard'}, 'Label' : {'alignment' : 'center', 'fontFamily' : 'Dialog', 'fontSize' : '16'}, 'LineStyle' : {'width' : '2.0'}})
	# call layout and save in file
	diagram.layout(algo="kk")
	diagram.dump_file(filename=fname_outdraw, folder="./")

# output representation graph and attractors in files 
outfg = open(fname_outgraph,"w")
#write RG
outfg.write("####### Representation graph #######\n\n")
outfg.write("### Input file: "+fname_in+"\n")
outfg.write("### Output file: "+fname_outgraph+"\n")
if args['ini']:
	outfg.write("### Included only part reachable from INI state: True\n\n")
else:
	outfg.write("### Included only part reachable from INI state: False\n\n")
# get attractor, switching state and grey lists
l_nodes = list(GX.nodes)
l_attr = []
l_dec = []
l_grey = []
for v in l_nodes:
	if SGraph[v].attr != -1:
		l_attr.append(v)
	if SGraph[v].dec != -1:
		l_dec.append(v)		
	if SGraph[v].grey != -1:
		l_grey.append(v)
# write nodes
outfg.write("### Number of nodes: "+str(len(l_nodes))+"\n")		
outfg.write("### Number of attractor nodes: "+str(len(l_attr))+"\n")	
outfg.write("### Number of switching nodes: "+str(len(l_dec))+"\n")	
outfg.write("### Number of other nodes: "+str(len(l_grey))+"\n\n")	
outfg.write("# List of attractor nodes: id genes binding_sites | int  str str\n")
for vert in l_attr:
	outfg.write(str(IDs[vert])+" "+SGraph[vert].name+" "+SGraph[vert].bs+" "+"\n")
outfg.write("\n")
outfg.write("# List of switching nodes: id genes binding_sites | int  str str\n")
for vert in l_dec:
	outfg.write(str(IDs[vert])+" "+SGraph[vert].name+" "+SGraph[vert].bs+" "+"\n")
outfg.write("\n")
outfg.write("# List of other nodes: id genes binding_sites | int  str str\n")
for vert in l_grey:
	outfg.write(str(IDs[vert])+" "+SGraph[vert].name+" "+SGraph[vert].bs+" "+"\n")
outfg.write("\n")
# write edges
outfg.write("# List of transitions: source_id dest_id gene_id gene_name {+,-} | int int int str def\n")	
for (s,d) in GX.edges:
	sid = IDs[s]
	did = IDs[d]
	outfg.write(str(sid)+" "+str(did)+" "+str(GX[s][d]['l'])+" "+Gnames[GX[s][d]['l']]+" "+str(GX[s][d]['t'])+"\n")
outfg.write("\n")
# write attractors
outfg.write("####### Attractor graphs #######\n\n")
outfg.write("##### Number of attractors: \n"+str(len(SCC_list))+"\n\n")
for i in range(len(SCC_list)): 
	draw_attr(i,SCC_list[i],outfg)
outfg.write("\n####### End #######\n")
outfg.close()

################################################################################
################################################################################
################################################################################