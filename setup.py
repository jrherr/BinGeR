#! /usr/bin/env python

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

from distutils.core import setup

setup(
	name = 'BinGeR',
	
	description = 'BinGeR (Binner for Genome Recovery): \
		in-situ Genome recovery tool for series metagenomes',
	
	scripts = ['src/BinGeR.py', 'src/contigSpace.py',
				'src/commPageRank.py', 'src/taxonomy.py',
				'src/utilities.py', 'src/phylo.py'],
	
	data_files = [('db', ['db/HMM.txt', 'db/ncbiNodes.lib', 'db/ncbiSciNames.lib'])],
	
	author="Chengwei Luo",
    author_email="luo.chengwei@gatech.edu",
    version="0.0.1",
    license="GNU GPL v3.0",
)