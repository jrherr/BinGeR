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

def commPageRank(initCore, seedNodes, options):
	alpha = options.cpr_alpha
	tol = options.cpr_tol
	maxiter = options.cpr_maxiter
	sets = {}
	for lca in seedNodes:
		print lca
		# for each seed set we pick 20 at random
		seedNum = max(20, len(seedNodes[lca]))
		seeds = set(random.sample(seedNodes[lca], seedNum))
		
		contigCounts = {}
		for seed in seeds:
			print seed
			# reduce the search space by searching only subgraph with a certain depth from seeds
			nodes = nx.ego_graph(initCore, seed, radius = 6).nodes()
			print '#nodes:', len(nodes)
			subgraph = initCore.subgraph(nodes)
			print 'subgraph extracted'
			contigSet = pprc(subgraph, seed, alpha, tol, maxiter)
			for contig in contigSet:
				if contig not in contigCounts:
					contigCounts[contig] = 0
				contigCounts[contig] += 1
		contigs = []
		for contig in contigCounts:
			if contigCounts[contig] >= 0.5*seedNum:
				contigs.append(contig)
		
		sets[lca] = contigs
		
	return sets
	
# end of commPageRank

def pprc(G, seed, alpha, tol, maxiter):
	"""
	This personalized PageRank clustering algorithm was originally designed by
	David F. Gleich at Purdue University. Here I tweak it to suit the networkx 
	module and the weighted edge scenario with multiple seeds option.
	"""
	
	
	"""
	x = {}
	r = {}
	Q = collections.deque()
	
	
	## initialized the seed weights
	r[seed] = 1
	Q.append(seed)
	iter = 0
	while len(Q) > 0 and iter <= maxiter:
		iter += 1
		v = Q.popleft()
		if v not in x:
			x[v] = 0
		x[v] += (1-alpha) * r[v]
		mass = alpha*r[v]/(2*len(G[v]))
		
		for u in G[v]: # for neighbors of v
			if u not in r:
				r[u] = 0.
			if r[u] < G.degree(u, weight = 'weight') * tol and \
				(r[u] + mass) >= G.degree(u, weight = 'weight') * tol:
				Q.append(u)
			r[u] = r[u] + mass
		r[v] = mass * G.degree(v, weight = 'weight')
		if r[v] >= G.degree(v, weight = 'weight') * tol:
			Q.append(v)
	"""
	
	Gvol = 2 * len(G.edges())
	personalizationDict = {}
	for node in G.nodes():
		if node == seed:
			personalizationDict[node] = 1
		else:
			personalizationDict[node] = 0
	pr = nx.pagerank(G, alpha = alpha, max_iter = maxiter, 
			personalization = personalizationDict, tol = tol, weight = 'weight')
	
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