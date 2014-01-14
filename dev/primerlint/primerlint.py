#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

primerlint.py is the point-of-execution for the analysis software.
"""

import importfiles as impfile
import analysismodules as anlymod
import argparse
import sys

def arg_parsing():
	parser = argparse.ArgumentParser(description="A tool for analysis of"\
	"primers (distributed with hiplex-primer.")

	parser.add_argument("--importCSV",
						metavar="ICSV",
						type=str,
						required=True,
						help="Primer import file location.")
	parser.add_argument("--cg",
						metavar="CG",
						action="store_true",
						help="Outputs cg-content of all imported primers.")
	parser.add_argument("--hairpin",
						metavar="HAIR",
						type=int,
						help="Computes potential hairpins of input size or"\
						"larger.")
	parser.add_argument("--hpminsize",
						metavar="HPMIN",
						type=int,
						help="Minimum size of hairpin to detect in imported"\
						"primers.")
	parser.add_argument("--primerdimer",
						metavar="PD",
						type=int,
						help="Detects potential primer-dimer interactions"\
						"of input size between set of all primers.")
	return parser


def runsimplehairpincsv(minhpsize, processedcsv):
	"""
	Run a simple (naive/literal) hairpin search over a group of sequences, 
	imported from a CSV file (foratted in FIXME format.)

	Args:
	minhpsize --> size of minimum hairpin to search for
	csv --> imported and processed csv file (use importfiles.ImportCSV function)
	"""
	print "\nHairpin", "-" * 70
	for line in processedcsv:
		linesplit = line.split(',')
		
		currentseq = anlymod.Hairpin(linesplit[1]).simpledetecthairpin(minhpsize)
		if sum(1 for item in currentseq) > 0:
			print "Gene: \t\t", linesplit[0]
			print "Sequence:\t", linesplit[1]
	
		# generator defined twice as sum above apparently destroys the data
		currentseq = anlymod.Hairpin(linesplit[1]).simpledetecthairpin(minhpsize)
		for item in currentseq:
			print item

### CSV Import 
print "CSVRead", "-" * 69
csvimport = impfile.ImportCSV("testcsv.csv")
csvimport.close()

try:
	testlist = csvimport.process()
except:
	print "Error processing file...", sys.exc_info()[0]
	raise
	sys.exit("Quitting due to error.")

print "CSV imported and processed successfully."


### Primer-Dimer
print "\nPrimerDimer", "-" * 66
primerdimertest = anlymod.PrimerDimer(
				'GGGGGGGGGGTGGGGGGGGG',
				'ACCCCCCCCCCCCCCCCCCC',
				10)

primerdimertest.basiccompare()

# imported and processed CSV file prior to calling this method
runsimplehairpincsv(4, testlist)

### GC Content
print "\nGCPercent", "-" * 68
print anlymod.Sequence('AAATTTTGGGGGAAAAAAAAAAAAAACCCC').gcpercent() * 100