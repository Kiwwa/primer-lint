#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

primerlint.py is the point-of-execution for the analysis software.
"""
import importfiles as impfile
import analysismodules as amod
import argparse
import sys

def arg_parsing():
	parser = argparse.ArgumentParser(description="A tool for analysis of"\
	"primers (distributed with hiplex-primer.")

	parser.add_argument("--importcsv",
						metavar="CSV",
						type=str,
						required=True,
						help="Primer import file location.")
	parser.add_argument("--outputfile",
						metavar="O",
						type=str,
						required=True,
						help="Name and location of output file.")	
	parser.add_argument("--cg",
						action="store_true",
						help="Output cg-content of all imported primers.")
	parser.add_argument("--hairpin",
						metavar="H",
						type=int,
						help="Compute potential hairpins of input size or"\
						"larger.")
	parser.add_argument("--primerdimer",
						metavar="PD",
						type=int,
						help="Detect potential primer-dimer interactions"\
						"of input size between set of all primers.")
	parser.add_argument("--htmloutput",
						metavar="HTML",
						type=str,
						help="If used, outputs a HTML version of the"\
						"analysed primer information.")
	return parser

def runcsvimport(arguments):
	csvimport = impfile.ImportCSV(args.importcsv)
	csvimport.close()

	try:
		importlist = csvimport.process()
	except:
		print "Error processing file...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	
	print "CSV imported and processed successfully."
	return importlist

def runsimplehairpincsv(minhpsize, processedcsv):
	"""
	Run a simple (naive/literal) hairpin search over a group of sequences, 
	imported from a CSV file (formatted in FIXME format.)

	Args:
	minhpsize --> size of minimum hairpin to search for
	csv --> imported and processed csv file (use importfiles.ImportCSV 
	function)
	"""
	print "\nHairpin", "-" * 70
	for line in processedcsv:
		linesplit = line.split(',')
		
		currentseq = amod.Hairpin(linesplit[1]).simpledetecthairpin(minhpsize)
		if sum(1 for item in currentseq) > 0:
			print "Gene: \t\t", linesplit[0]
			print "Sequence:\t", linesplit[1]
	
		# generator defined twice as sum above apparently destroys the data
		currentseq = amod.Hairpin(linesplit[1]).simpledetecthairpin(minhpsize)
		for item in currentseq:
			print item

def localprimerdimercsv(minscore, processedcsv, arguments):
	"""
	Find the local highest score (score equates to complementary bases) 
	between two primers for all primers in a csv-import-object.

	No return, outputs a visual representation of the best-scoring local 
	alignments to a file.

	###FIXME: needs to return the object so it can be output by a more general
	function upstream (get rid of the args argument)
	"""
	output = open(args.outputfile, "w")

	for numitemone, itemone in enumerate(processedcsv):
		linesplitone = itemone.split(",")
		print " %02d started compares on %s" % (numitemone, linesplitone[0])
		for itemtwo in processedcsv:
			linesplittwo = itemtwo.split(",")

			primerdimertest = amod.PrimerDimer(linesplitone[1],
				amod.Sequence(linesplittwo[1]).complement(), minscore)
			primerdimeroutput = primerdimertest.pdlocal()

			output.write("Comparing '%s' and '%s':\n" %\
				(linesplitone[0], linesplittwo[0]))
			for item in primerdimeroutput:
				output.write(amod.PrimerDimer.\
					format_alignment_compl(item[0],
					item[1],
					item[2],
					item[3],
					item[4]))

def localprimerdimer(minscore, baseprimer, compareprimer):
	"""
	Find the local highest score (score equates to complementary bases) 
	between baseprimer and compareprimer.

	No return, outputs a visual representation of the best-scoring local 
	alignments to a file.

	FIXME: should return the written data (or something) and not have arg
	as an argument.
	"""
	output = open("workfile.log", "w")
	print "started compares on %s" % (baseprimer)

	primerdimertest = amod.PrimerDimer(baseprimer,
		amod.Sequence(compareprimer).complement(), minscore)
	primerdimeroutput = primerdimertest.pdlocal()

	output.write("Comparing '%s' and '%s':\n" %\
		(baseprimer, compareprimer))
	for item in primerdimeroutput:
		output.write(amod.\
			PrimerDimer.\
				format_alignment_compl(item[0],
					item[1],
					item[2],
					item[3],
					item[4]))


### CSV Import ------------------------------------------
argparsing = arg_parsing()
args = argparsing.parse_args()

print "CSVRead", "-" * 69
importedcsv = runcsvimport(args)


### Primer-Dimer ----------------------------------------
print "\nPrimerDimer", "-" * 66
localprimerdimer(1, "GATCAAGAG", "CATGGATTCC")
localprimerdimercsv(1, importedcsv, args)

### FIXME: basiccompare functionality not complete
amod.PrimerDimer("AAAAAAA", "TTTTTTT", 4).basiccompare()


### Hairpin ---------------------------------------------
print "Hairpin", "-" * 66
for item in importedcsv:
	linesplit = item.split(",")
	
	print linesplit[0]
	for item in amod.Hairpin(linesplit[1]).simpledetecthairpin(5):
		print item


### GC Content ------------------------------------------
print "\nGCPercent", "-" * 68
for item in importedcsv:
	linesplit = item.split(",")

	print "%s = %02d%%" % (linesplit[0],
		amod.Sequence(linesplit[1]).gcpercent() * 100)