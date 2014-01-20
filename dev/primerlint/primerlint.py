#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

primerlint.py is the point-of-execution for the analysis software.
"""
import importfiles as infile
import outputfiles as outfile
import analysismodules as amod
import argparse
import sys

def wait():
    raw_input("Press enter to continue...")

def arg_parsing():
	parser = argparse.ArgumentParser(description="A tool for analysis of"\
	"primers (distributed with hiplex-primer.")

	parser.add_argument("--importhiplex",
						metavar="CSV",
						type=str,
						required=True,
						help="Import primers from a hiplex-output "\
						"CSV-formatted file.")
	parser.add_argument("--importfasta",
						metavar="FAS",
						type=str,
						help="Import primers from a FASTA-formatted file.")
	parser.add_argument("--outputlog",
						metavar="O",
						type=str,
						required=True,
						help="Name and location of output file.")	
	parser.add_argument("--cg",
						action="store_true",
						help="Output cg-content of all imported primers.")
	parser.add_argument("--cgmin",
						metavar="CGMIN",
						type=int,
						help="Minimum percent of CG content per primer "\
						"(default: 40).")
	parser.add_argument("--cgmax",
						metavar="CGMAX",
						type=int,
						help="Maximum percent of CG content per primer "\
						"(default: 60).")
	parser.add_argument("--hairpin",
						metavar="H",
						type=int,
						help="Compute potential hairpins of input size or "\
						"larger.")
	parser.add_argument("--primerdimer",
						metavar="PD",
						type=int,
						help="Detect potential primer-dimer interactions "\
						"of input size between set of all primers.")
	parser.add_argument("--htmloutput",
						metavar="HTML",
						type=str,
						help="If used, outputs a HTML version of the "\
						"analysed primer information.")
	return parser

def runhipleximport(arguments):
	hipleximport = infile.ImportHiplex(args.importhiplex)

	try:
		hipleximport.open()
		importlist = hipleximport.process()
		hipleximport.close()	
	except:
		print "runcsvimport: Error importing CSV-file...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	return importlist

def runfastaimport(arguments):
	fastaimport = infile.ImportFASTA(args.importfasta)

def runsimplehairpincsv(minhpsize, processedcsv):
	"""
	Run a simple (naive/literal) hairpin search over a group of sequences, 
	imported from a CSV file (formatted in FIXME format.)

	Args:
	minhpsize --> size of minimum hairpin to search for
	csv --> imported and processed csv file (use importfiles.ImportHiplex 
	function)
	"""
	print "\nHairpin", "-" * 70
	for num, line in enumerate(processedcsv):
		linesplit = line.split(',')
		
		currentseq = amod.Hairpin(linesplit[1]).simpledetecthairpin(minhpsize)
		if sum(1 for item in currentseq) > 0:
			print "Gene: \t\t", linesplit[0]
			print "Sequence:\t", linesplit[1]
	
		# generator defined twice as sum above apparently destroys the data
		currentseq = amod.Hairpin(linesplit[1]).simpledetecthairpin(minhpsize)
		for item in currentseq:
			print item

def localprimerdimer(minscore, processeddata, arguments):
	"""
	Find the local highest score (score equates to complementary bases) 
	between two primers for all primers in a csv-import-object.

	No return, outputs a visual representation of the best-scoring local 
	alignments to a file.

	###FIXME: needs to return the object so it can be output by a more general
	function upstream (get rid of the args argument)
	"""

	for numitemone, itemone in enumerate(processeddata):
		print " %02d. Started compares on %s" % (numitemone, itemone[0])
		for itemtwo in processeddata[numitemone:]:
			primerdimertest = amod.PrimerDimer(itemone[1],
				amod.Sequence(itemtwo[1]).complement(), minscore)
			primerdimeroutput = primerdimertest.pdlocal()

			logoutput.write("Comparing '%s' and '%s':\n" %\
				(itemone[0], itemtwo[0]))
			for item in primerdimeroutput:
				yield amod.PrimerDimer.\
					format_alignment_compl(item[0],
					item[1],
					item[2],
					item[3],
					item[4])

def singlelocalprimerdimers(minscore, baseprimer, compareprimer):
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

def gccontent(processedcsv, cgmin=40, cgmax=60):
	for item in processedcsv:
		exonname = item[0]
		gccontent = amod.Sequence(item[1]).gcpercent() * 100

		if gccontent < cgmin:
			yield "%s = %02d%% (is lower than min = %d)\n" % (exonname,
				gccontent, cgmin)
		elif gccontent > cgmax:
			yield "%s = %02d%% (is higher than max = %d)\n" % (exonname,
				gccontent, cgmax)
		else:
			yield "%s = %02d%%\n" % (exonname, gccontent)

def importfileblock():
	"""
	Check whether multiple import actions have been requested; fail if so 
	and inform user why only one import can be performed at a time.
	"""

	### Test for Multiple Import -----------------------
	
	### add new import functions to this list
	multipleimporttest = [args.importhiplex, args.importfasta]

	importargcounter = 0
	for item in multipleimporttest:
		if item is not None:
			importargcounter += 1

	if importargcounter != 1:
		print "You have selected multiple import locations; please only "\
		"select one at a time."
		sys.exit("Quitting due to multiple import locations.")


	### Hiplex CSV Import -------------------------------
	print "\nHiplex Import", "-" * 67
	try:
		importedhiplex = runhipleximport(args)
	except:
		print "Error importing hiplex-CSV-file...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		print "Hiplex-CSV file imported successfully"
		return importedhiplex

	### FASTA Import ------------------------------------
	print "\nFASTA Import", "-" * 66

### Argument parsing ------------------------------------
argparsing = arg_parsing()
args = argparsing.parse_args()


### Attempt data import ---------------------------------
importeddata = importfileblock()


### Open a log-file -------------------------------------
print "\nLog File Open", "-" * 65
try:
	logoutput = outfile.OutputLog(args.outputlog)
	logoutput.open()
except:
	print "Error opening log-file...", sys.exc_info()[0]
	raise
	sys.exit("Quitting due to error.")
else:
	print "Log file opened at %s" % args.outputlog


### Primer-Dimer ----------------------------------------
print "\nPrimerDimer", "-" * 66
if args.primerdimer != None:	
	try:
		logoutput.write("")
		for item in localprimerdimer(1, importeddata, args):
			logoutput.write(item)
	except:
		print "Error running primer-dimer...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		pass


### Hairpin ---------------------------------------------
print "\nHairpin", "-" * 66
if args.hairpin != None:
	try:
		for num, row in enumerate(importeddata):
			# print ("%d. Testing %s for hairpin...") % (num, row[0])
			for item in amod.Hairpin(row[1]).simpledetecthairpin(5):
				logoutput.write(row[0] + ": " + row[1] + "\n")
	except:
		print "Error running hairpin...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		print "Potential hairpins of greater than size %d calculated"\
			% args.hairpin


### GC Content ------------------------------------------
print "\nGCPercent", "-" * 68
if args.cg == True:
	try:
		if args.cgmin != None and args.cgmax != None:
			gccontent(importexdata, args.cgmin, args.cgmax)
		else:
			print "You have not given cg-minimum and cg-maximum values; "\
			"defaulting to 40 and 60"
			wait()
			for item in gccontent(importeddata):
				logoutput.write(item)
	except:
		print "Error equating GC Content...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		print "GC Content calculated."


### Close open files --------------------------------------
logoutput.close()