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

def commRank(initCore, seedNodes, options):
	# get weighted node degrees for each node
	nodeDegrees = initCore.degree(weight = 'weight')
	# sorted nodes using degree in a descending order
	sortedNodeDegrees = sorted(nodeDegrees.iteritems(), key = lambda x: x[1], reverse = True)
	
	# get 85%, 90%, 95% percentile index
	percentile_85_index = int(len(sortedNodeDegrees) * 0.85)
	percentile_90_index = int(len(sortedNodeDegrees) * 0.9)
	percentile_95_index = int(len(sortedNodeDegrees) * 0.95)
	
	# remove last 5% nodes
	nodes_to_remove = map( itemgetter(0), sortedNodeDegrees[percentile_95_index:] )
	initCore.remove_nodes_from(nodes_to_remove)
	
	n = 0
	for component in nx.connected_components(initCore):
		print len(component)
		n+=1
	print '======================'
	print n
	return {}