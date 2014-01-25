#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

importfiles.py is used within PrimerLint to import data from files.
"""

class ImportedData(object):
	"""
	All imported data should conform to this structure and usage;
		> row 1 = exon name
		> row 2 = sequence
	"""
	def __init__(self, records):
		self.records = records


class ImportHiplex(ImportedData):
	def __init__(self, filename):
		super(ImportHiplex, self).__init__([])
		self.filename = filename

	def open(self):
		self.inputfile = open(self.filename, "r")

	def process(self):
		inputfiledata = self.inputfile.read().split('\n')
		for splitline in inputfiledata:
			delimited = splitline.split(',')
			self.records.append([delimited[0], delimited[1]])
		return self.records

	def close(self):
		self.inputfile.close()


class ImportFASTA(object):
	def __init__(self, filename):
		from Bio import SeqIO
		self.filename = filename
		
		seqrecordslist = []
		for seqrecord in SeqIO.parse(filename, "fasta"):
			seqrecordslist.append(seqrecord)

	def printrecords(self):
		for record in self.seqrecordslist:
			print record