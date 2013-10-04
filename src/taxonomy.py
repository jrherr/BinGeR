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
from operator import itemgetter

class TaxonNode:
	def __init__(self):
		self.taxonID = 0
		self.prevNode = None
		self.rank = None
		self.sciName = None
		self.commName = None
	
	def isRoot(self):
		if self.taxonID > 1:
			return False
		return True	

class TaxonEdge:
	def __init__(self):
		self.fromNode = 0
		self.toNode = 0

class TaxonTree:
	def __init__(self):
		self.nodes = {}
		self.edges = {}
	
	def loadTreeFromNodeLib(self, nodeLibFile, sciNameFile):
		nodefh = open(nodeLibFile,'r')
		sciNamefh = open(sciNameFile,'r')
		
		# load taxonID->sciName map
		sciNameMap = {}
		while 1:
			line=sciNamefh.readline().rstrip('\n')
			if not line:
				break
			col = line.split('\t')
			taxonID = int(col[0])
			name = col[1]
			sciNameMap[taxonID] = name
		sciNamefh.close()
		
		# construct the taxonomy tree here
		while 1:
			line=nodefh.readline().rstrip('\n')
			if not line:
				break
			col = line.split('\t')
			currentNodeID = int(col[0])
			prevNodeID = int(col[1])
			rank = col[2]
			
			if currentNodeID not in self.nodes:
				self.nodes[currentNodeID] = TaxonNode()
			if prevNodeID not in self.nodes:
				self.nodes[prevNodeID] = TaxonNode()
			
			self.nodes[currentNodeID].taxonID = currentNodeID
			self.nodes[currentNodeID].rank = rank
			if currentNodeID in sciNameMap:
				self.nodes[currentNodeID].sciName = sciNameMap[currentNodeID]

			if currentNodeID not in self.edges:
				self.edges[currentNodeID] = TaxonEdge()
				self.edges[currentNodeID].fromNode = currentNodeID
				self.edges[currentNodeID].toNode = prevNodeID
		nodefh.close()
		
	def getTaxonomyPath(self, nodeID):
		taxonIDs = []
		Ranks = []
		sciNames = []
		try:
			currentNode = self.nodes[nodeID]
		except:
			return taxonIDs, Ranks, sciNames
			
		while (not currentNode.isRoot()):
			taxonIDs.append(currentNode.taxonID)
			Ranks.append(currentNode.rank)
			sciNames.append(currentNode.sciName)
			
			#print currentNode.taxonID, currentNode.rank, currentNode.name
			prevNodeID = self.edges[currentNode.taxonID].toNode
			currentNode = self.nodes[prevNodeID]
			
		return taxonIDs, Ranks, sciNames

## end of loadTreeFromNodeLib
