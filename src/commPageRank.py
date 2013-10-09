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
#import networkx as nx

def commPageRank(initCore, seedNodes, options):
	alpha = options.cpr_alpha
	tol = 1. - alpha
	
	sets = {}
	for lca in seedNodes:
		# for each seed set we pick 20 at random
		seeds = random.sample(seedNodes[lca], 20)
		contigCounts = {}
		for seed in seeds:
			contigSet = pprc(initCore, seed, alpha)
			for contig in contigSet:
				if contig not in contigCounts:
					contigCounts[contig] = 0
				contigCounts[contig] += 1
		contigs = []
		for contig in contigCounts:
			if contigCounts[contig] >= 10:
				contigs.append(contig)
		
		sets[lca] = contigs
	
	return sets
	
# end of commPageRank

def pprc(G, seed, alpha):
	"""
	This personalized PageRank clustering algorithm was originally designed by
	David F. Gleich at Purdue University. Here I tweak it to suit the networkx 
	module and the weighted edge scenario with multiple seeds option.
	"""
	tol = 1 - alpha
	Gvol = 2 * len(G.edges())
	
	x = {}
	r = {}
	Q = collections.deque()
	
	## initialized the seed weights
	r[seed] = 1
	Q.append(seed)
		
	while len(Q) > 0:
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
			
	# find cluster 
	# normalized by weighted degree
	for v in x:
		x[v] = x[v]/G.degree(v, weight = 'weight')
		
	# sort x's keys by value in decreasing order
	sv = sorted(x.iteritems(), key = lambda x: x[1], reverse = True)
	S = set()
	volS = 0.
	cutS = 0.
	bestond = 1.
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
		
		if cutS/min(volS, Gvol - volS) < bestcond:
			bestcond = cutS/min(volS, Gvol - volS)
			bestset = set(S)
	return bestset
	
# End of pprc