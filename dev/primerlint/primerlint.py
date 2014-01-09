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
						"larger")
	parser.add_argument("--primerdimer",
						metavar="PD",
						type=int,
						help="Detects potential primer-dimer interactions"\
						"of input size between set of all primers.")
	return parser


print "CSVRead", "-" * 69
csvimport = impfile.ImportCSV("testcsv.csv")
csvimport.close()
testlist = csvimport.process()

for item in testlist:
	print item

print "\nPrimerDimer", "-" * 66
primerdimertest = anlymod.PrimerDimer(
				'GGGGGGGGGGTGGGGGGGGG',
				'ACCCCCCCCCCCCCCCCCCC',
				10)

primerdimertest.basiccompare()

print "\nHairpin", "-" * 70

newobject = anlymod.Hairpin("AAGGGGAAAAAAAAAAAACCCC")
for item in newobject.simpledetecthairpin(4):
	print item

print "\nGCPercent", "-" * 68
print anlymod.Sequence('AAGGGGAAAAAAAAAAAACCCC').gcpercent() * 100