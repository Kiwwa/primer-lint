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
	parser.add_argument("--cgbrief",
						action="store_true",
						help="Less verbose cg-content output of all imported primers.")
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
						metavar="HP",
						type=int,
						help="Compute potential hairpins of input size or "\
						"larger.")
	parser.add_argument("--hairpinmaxsizes",
						action="store_true",
						help="Output maximum hairpin size for each primer input.")
	parser.add_argument("--primerdimer",
						metavar="PD",
						type=int,
						help="Detect potential primer-dimer interactions "\
						"of input size between set of all primers.")
	return parser

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
		gccontent = amod.Sequence(item[1]).gcpercent() * 100

		if gccontent < cgmin:
			yield "%02d%% (is lower than min = %d)\n" % (
				gccontent, cgmin)
		elif gccontent > cgmax:
			yield "%02d%% (is higher than max = %d)\n" % (
				gccontent, cgmax)
		else:
			yield "%02d%%\n" % (gccontent)

def gcbrief(processedcsv):
	for item in processedcsv:
		gccontent = amod.Sequence(item[1]).gcpercent() * 100
		yield "%02d\n" % (gccontent)

def importfiles():
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

	if importargcounter > 1:
		print "You have selected multiple import locations; please only "\
		"select one at a time."
		sys.exit("Quitting due to multiple import locations.")
	elif importargcounter == 0:
		print "You have not selected any import locations; please select "\
		"data to import."
		sys.exit("Quitting due to no import locations.")


	### Hiplex CSV Import -------------------------------

	if args.importhiplex != None:
		try:
			print "\nHiplex Import", "-" * 67
			hipleximport = infile.ImportHiplex(args.importhiplex)
			importedhiplex = hipleximport.process()
		except IOError:
			print "IOError: Cannot locate file at '%s'" % args.importhiplex
			sys.exit("Quitting due to import error.")
		except:
			print "Error importing hiplex-CSV file...", sys.exc_info()[0]
			sys.exit("Quitting due to error.")
		else:
			print "Hiplex-CSV file imported successfully."
			return importedhiplex


	### FASTA Import ------------------------------------
	if args.importfasta != None:
		try:
			print "\nFASTA Import", "-" * 66
			fastaimport = infile.ImportFASTA(args.importfasta)
			importedfasta = fastaimport.process()
		except:
			print "Error importing FASTA-formatted-file...", sys.exc_info()[0]
			raise
			sys.exit("Quitting due to error.")
		else:
			print "FASTA-formatted file imported successfully."
			return importedfasta


### Argument parsing ------------------------------------
argparsing = arg_parsing()
args = argparsing.parse_args()


### Attempt data import ---------------------------------
importeddata = importfiles()


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
	print "Log file opened at %s." % args.outputlog


### Primer-Dimer ----------------------------------------
if args.primerdimer != None:
	print "\nPrimerDimer", "-" * 66

if args.primerdimer != None:	
	try:
		logoutput.write("PRIMER-DIMER DETECTION --------------\n")
		for item in localprimerdimer(args.primerdimer, importeddata, args):
			logoutput.write(item)
	except:
		print "Error running primer-dimer...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		pass


### Hairpin ---------------------------------------------
if args.hairpin != None or args.hairpinmaxsizes == True:
	print "\nHairpin", "-" * 66

if args.hairpin != None:
	try:
		logoutput.write("HAIRPIN DETECTION: --------------------\n")
		for num, row in enumerate(importeddata):
			# print ("%d. Testing %s for hairpin...") % (num, row[0])
			for item in amod.Hairpin(row[1]).simpledetecthairpin(args.hairpin):
				logoutput.write(row[0] + " - " + row[1] + "\n")
				logoutput.write(str(len(item[4])) + "\n")
	except:
		print "Error running hairpin...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		print "Potential hairpins of greater than size %d calculated"\
			% args.hairpin

if args.hairpinmaxsizes == True:
	try:
		logoutput.write("MAX HAIRPIN OUTPUT: --------------------\n")
		hpcounter = 0
		for num, row in enumerate(importeddata):
			for item in amod.Hairpin(row[1]).simpledetecthairpin(1):
				pass
			logoutput.write(str(len(item[4])) + "\n")
			if hpcounter % 50 == 0:
				print "Done... %d" % hpcounter
			hpcounter += 1
	except:
		print "Error running hairpin max sizes...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		pass


### GC Content ------------------------------------------
if args.cg == True or args.cgbrief == True:
	print "\nGCPercent", "-" * 68

if args.cg == True:
	try:
		logoutput.write("GC-CONTENT DETECTION -------------------\n")
		if args.cgmin != None and args.cgmax != None:
			for item in gccontent(importeddata, args.cgmin, args.cgmax):
				logoutput.write(item)
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

if args.cgbrief == True:
	try:
		logoutput.write("BRIEF GC-CONTENT DETECTION -------------------\n")
		if args.cgmin != None and args.cgmax != None:
			for item in gcbrief(importeddata):
				logoutput.write(item)
		else:
			print "You have not given cg-minimum and cg-maximum values; "\
			"defaulting to 40 and 60"
			wait()
			for item in gcbrief(importeddata):
				logoutput.write(item)
	except:
		print "Error equating brief GC Content...", sys.exc_info()[0]
		raise
		sys.exit("Quitting due to error.")
	else:
		print "GC Content calculated."


### Close open files --------------------------------------
logoutput.close()