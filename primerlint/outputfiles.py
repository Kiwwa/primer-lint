#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

outputfiles.py is used within PrimerLint to import data from files.
"""

class OutputLog(object):
	def __init__(self, filename):
		self.filename = filename
		self.outputfile = open(self.filename, "w").close()

	def open(self):
		self.outputfile = open(self.filename, "a")

	def write(self, outputtext):
		### require newline?
		self.outputfile.write(outputtext)

	def close(self):
		self.outputfile.close()