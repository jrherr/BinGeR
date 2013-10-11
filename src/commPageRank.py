#!/usr/bin/python
# -*- coding: utf-8 -*- 

"""
BinGeR (Binner for Genome Recovery): 
	in-situ Genome recovery tool for series metagenomes

Copyright(c) 2013 Chengwei Luo (luo.chengwei@gatech.edu)

	This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

https://github.com/luo-chengwei/BinGeR

for help, type:
python BinGeR.py --help
"""

import sys
import re
import os
import random
import collections
import networkx as nx
from operator import itemgetter

def commCrunch(initCore, options):
	# make a local copy of the input initCore
	Core = initCore.copy()
	
	# get weighted node degrees for each node
	weight = {}
	for node in Core.nodes(data = True):
		print node
	return {}
	nodeDegrees = Core.degree(weight = 'weight')
	
	# sorted nodes using degree in a descending order
	sortedNodeDegrees = sorted(nodeDegrees.iteritems(), key = lambda x: x[1], reverse = True)
	
	# get percentile indices over the sortedNodeDegrees
	percentileIndices = []
	percentile_60_index = int(len(sortedNodeDegrees) * 0.6)
	percentile_85_index = int(len(sortedNodeDegrees) * 0.85)
	percentile_90_index = int(len(sortedNodeDegrees) * 0.9)
	percentile_95_index = int(len(sortedNodeDegrees) * 0.95)
	
	# remove last 5% nodes
	nodes_to_remove = map( itemgetter(0), sortedNodeDegrees[percentile_60_index:] )
	initCore.remove_nodes_from(nodes_to_remove)
	
	n = 0
	for component in nx.connected_components(initCore):
		if len(component) < 100: continue
		print len(component)
		n+=1
	print '======================'
	print n
	return {}
	
# End of commRank

def commPageRank(core, seedNodes, tightNodes, coreIndex, options):
	if not options.quiet:
		sys.stdout.write('Running community personalized PageRank to evaluate the core set.\n')
		
	alpha = options.cpr_alpha
	tol = options.cpr_tol
	maxiter = options.cpr_maxiter
	
	sets = {}
	for lca in seedNodes:
		print lca
		# for each seed set we pick 10 at random
		seedNum = min(20, len(seedNodes[lca]))
		seeds = set(random.sample(seedNodes[lca], seedNum))
		
		nodes = []
		for seed in seeds:
			print seed
			# reduce the search space by searching only subgraph with a certain depth from seeds
			nodes += nx.ego_graph(initCore, seed, radius = 4).nodes()
			print '#nodes:', len(nodes)
		
		# extract the subgraph
		subgraph = initCore.subgraph(set(nodes))
		print 'subgraph extracted, node number:', len(subgraph.nodes())
		
		# run pprc here
		contigSet = pprc(subgraph, seeds, alpha, tol, maxiter)
		sets[lca] = contigSet
		
	return sets
	
# end of commPageRank

def pprc(G, seeds, alpha, tol, maxiter):
	"""
	This personalized PageRank clustering algorithm was originally designed by
	David F. Gleich at Purdue University. Here I tweak it to suit the networkx 
	module and the weighted edge scenario with multiple seeds option.
	"""
	

	Gvol = 2 * len(G.edges())
	
	# initialize personalization with seeds
	personalizationDict = {}
	for node in G.nodes():
		if node in seeds:
			personalizationDict[node] = 1
		else:
			personalizationDict[node] = 0
			
	pr = nx.pagerank(G, alpha = alpha, max_iter = maxiter, 
			personalization = personalizationDict, tol = tol)
	
	sys.stdout.write('Finished init node values.\n')
	
	# find cluster 
	# normalized by weighted degree
	for v in pr:
		pr[v] = pr[v]/G.degree(v, weight = 'weight')
		
	# sort x's keys by value in decreasing order
	sv = sorted(pr.iteritems(), key = lambda x: x[1], reverse = True)
	S = set()
	volS = 0.
	cutS = 0.
	bestcond = 1.
	bestset = sv[0]
	for p in sv:
		s = p[0]  # get the vertex
		volS += G.degree(s, weight = 'weight')
		for v in G[s]:
			if v in S:
				cutS -= 1
			else:
				cutS += 1
		S.add(s)
		
		if Gvol == volS: continue
		
		if cutS/min(volS, Gvol - volS) < bestcond:
			bestcond = cutS/min(volS, Gvol - volS)
			bestset = set(S)
			
	print "set size: ", len(bestset)
	
	return bestset
	
# End of pprc