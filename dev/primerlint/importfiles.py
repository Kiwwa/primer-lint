#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

importfiles.py is used within PrimerLint to import data from files.
"""

class ImportCSV(object):
	def __init__(self, filename):
		self.filename = filename
		self.inputfile = open(self.filename, 'r')
		self.inputfiledata = self.inputfile.read().split('\n')

	def process(self):
		csvlist = []
		for splitline in self.inputfiledata:
			csvlist.append(splitline.split(','))
		return filter(None, self.inputfiledata)

	def close(self):
		self.inputfile.close()