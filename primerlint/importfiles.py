#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

importfiles.py is used within PrimerLint to import data from files.
"""
import sys
from Bio import SeqIO

def validatesequences(listofseqs):
	for sequence in listofseqs:
		for base in sequence[1]:
			if base not in 'ACGT':
				sys.exit("Bad base found in sequence '%s'." % sequence[1])
				

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
		ImportedData.__init__(self, [])
		self.filename = filename

	def process(self):
		try:
			with open(self.filename, "r") as self.inputfile:
				inputfiledata = self.inputfile.read().split('\n')
			
			for splitline in inputfiledata:
				delimited = splitline.split(',')
				self.records.append([delimited[0], delimited[1]])
				validatesequences(self.records)
			return self.records
		except:
			raise


class ImportFASTA(ImportedData):
	def __init__(self, filename):
		ImportedData.__init__(self, [])
		self.filename = filename

	def process(self):
		try:
			for seqrecord in SeqIO.parse(self.filename, "fasta"):
				self.records.append([seqrecord.id, str(seqrecord.seq)])
				validatesequences(self.records)
			return self.records
		except:
			raise