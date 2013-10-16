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
import cPickle
import community
from multiprocessing import Pool

def commCrunch(initCore, coreIndex, projInfo, options):
	# pickfile that holds the partition results
	partitionFile = projInfo.out_dir + '/initCores/initCore.' + str(coreIndex+1) + '.partition'
	
	if os.path.exists(partitionFile):
		try:
			if not options.quiet:
				sys.stdout.write('[initCore %i] Loading partitions...\n'%(coreIndex+1))
			pfh = open(partitionFile, 'rb')
			partition = cPickle.load(pfh)
			pfh.close()
		except:
			sys.stderr.write('FATAL: error in unpickling partition file: %s' % partitionFile)
	
	else:
		if not options.quiet:
			sys.stdout.write('[initCore %i] Partitioning initCore...\n' % (coreIndex+1))
		partition = community.best_partition(initCore)
		try:
			pfh = open(partitionFile, 'wb')
			cPickle.dump(partition, pfh)
			pfh.close()
		except:
			sys.stderr.write('FATAL: error in pickling partition file: %s' % partitionFile)
	
	# transform the partition dict into a dict keyed by commID and valued by lists of contigIDs
	if not options.quiet:
		sys.stdout.write('[initCore %i] Extracting subgraphs...\n' % (coreIndex+1))
	
	contigSet = []
	for commID in set(partition.values()):
		contigList = [nodes for nodes in partition.keys() if partition[nodes] == commID]
		contigSet.append(contigList)
	
	# get the subgraphs
	subgraphs = []
	for contigList in contigSet:
		subgraph = initCore.subgraph(contigList)
		subgraphs.append(subgraph.copy())
		
	if not options.quiet:
		sys.stdout.write('[initCore %i] Done.\n' % (coreIndex+1))
		
	return subgraphs
	
# End of commCrunch

def commPageRank(cores, coreIndex, pTree, seedNodes, tightNodes, options):
	if not options.quiet:
		sys.stdout.write('[initCore %i] Evaluating initCore using community personalized PageRank...\n' % (coreIndex+1))
		
	alpha = options.cpr_alpha
	tol = options.cpr_tol
	maxiter = options.cpr_maxiter
	
	# create node->lca mapping for seeds
	nodeMap = {}
	for lca in seedNodes:
		for seed in seedNodes[lca]:
			nodeMap[seed] = lca
	
	# create node->lca mapping for tight seeds
	tightMap = {}
	for lca in tightNodes:
		for seed in tightNodes[lca]:
			tightMap[seed] = lca
		
	# core->lca mapping
	subIndex = 0
	coreMap = {}
	coreSeeds = {}
	lcaType = {}
	sets = {}
	tempSets = {}
	# iterate every core
	number_of_cores = len(cores)
	
	for index, core in enumerate(cores):
		nodes = core.nodes()
		if len(nodes) < 100:
			continue
		seeds = {}
		S = {}
		for node in nodes:
			if node not in nodeMap and node not in tightMap:
				continue
			if node in nodeMap:
				lca = nodeMap[node]
				lcaType[lca] = 1
			else:
				lca = tightMap[node]
				lcaType[lca] = 0
				
			if lca not in seeds:
				seeds[lca] = []
			
			seeds[lca].append(node)
		
		if len(seeds) <= 1:
			subIndex += 1
			if len(seeds) == 0:
				coreID = str(coreIndex) + '.' + str(subIndex) + '.unknown'
			else:
				lca = seeds.keys()[0]
				coreID = str(coreIndex) + '.' + str(subIndex) + '.' + str(lca) + '.regular'
			sets[coreID] = core.nodes()
		
		else:
			pprcCMDs = [[core, seeds[lca], alpha, tol, maxiter, lca] for lca in seeds]
			pool = Pool(options.num_proc)
			results = pool.map_async(pprc, pprcCMDs)
			pool.close()
			pool.join()
			
			for result in results.get():
				print result
				lca = result[0]
				contigs = result[1]
				S[lca] = contigs
				
			# calculate the overlaps between LCAs
			lcas = S.keys()
			lcaGraph = nx.Graph()
			lcaGraph.add_nodes_from(lcas)
			lcaLinks = []
			for i, lca1 in enumerate(lcas):
				for j, lca2 in enumerate(lcas):
					if i <= j:
						continue
					minLength = min(len(S[lca1]), len(S[lca2]))
					overlapLength = len(S[lca1] & S[lca2])
					overlapPercentage = float(overlapLength)/minLength
					if overlapPercentage > 0.5:
						lcaLinks.append((lca1, lca2))
			lcaGraph.add_edges_from(lcaLinks)
			
			for component in nx.connected_components(lcaGraph):
				contigSet = set()
				for lca in component:
					contigSet &= S[lca]
				commLCA = commLCA(S)
				sets[commLCA] = contigSet
				
				print component, len(contigSet)
		
		
			print '=============='
	
		if not options.quiet:
			sys.stdout.write('[initCore %i] %i out of %i subcores finished.\n' \
					% (coreIndex+1, index+1, number_of_cores))
			sys.stdout.flush()
	
	if not options.quiet:
		sys.stdout.write('[initCore %i] Done.                               \n'%(coreIndex+1))
		
	return sets
	
# end of commPageRank

def pprc(args):
	"""
	This personalized PageRank clustering algorithm was originally designed by
	David F. Gleich at Purdue University. Here I tweaked it to suit the networkx 
	module and the weighted edge scenario with multiple seeds option.
	"""
	
	(G, seeds, alpha, tol, maxiter, lca) = args

	Gvol = 2 * len(G.edges())
	
	# initialize personalization with seeds
	personalizationDict = {}
	for node in G.nodes():
		if node in seeds:
			personalizationDict[node] = 1
		else:
			personalizationDict[node] = 0
	
	try:		
		pr = nx.pagerank(G, alpha = alpha, max_iter = maxiter, 
				personalization = personalizationDict, tol = tol)
	
	except nx.exception.NetworkXError:	
		pr = {}
		r = {}
		Q = collections.deque()
		for s in seeds:
			r[s] = 1/len(seeds)
			Q.append(s)
			
		iterNum = 0
		while len(Q) > 0 and iterNum <= maxiter:
			v = Q.popleft()
			if v not in pr:
				pr[v] = 0
			pr[v] += (1-alpha)*r[v]
			mass = alpha*r[v]/(2*len(G[v]))
			for u in G[v]:
				if u not in r:
					r[u] = 0.
				if r[u] < len(G[u])*tol and r[u] + mass >= len(G[u])*tol:
					Q.append(u)
				r[u] = r[u] + mass
			r[v] = mass*len(G[v])
			if r[v] >= len(G[v])*tol:
				Q.append(v)
			iterNum += 1
		
	# sys.stdout.write('Finished init node values.\n')
	
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
	
	return lca, bestset
	
# End of pprc