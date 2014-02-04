#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

analysismodules.py contains the classes that represent different analysis 
tasks computed by the PrimerLint software.
"""
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

class Sequence(object):
	def __init__(self, sequence):
		self.sequence = str(sequence).upper()
		self.name = "seq_name_not_set"

	def setname(self, name):
		self.name = name

	def complement(self):
		return ''.join(map(complementbase, self.sequence))

	def reversecomplement(self):
		return ''.join(map(complementbase, self.sequence[::-1]))

	def reversesequence(self):
		return self.sequence[::-1]

	def gcpercent(self):
		countg = self.sequence.count("G")
		countc = self.sequence.count("C")

		if countc == 0 and countg == 0:
			return float(0)
		else:
			return 	float(self.sequence.count('G') +\
				self.sequence.count('C')) / float(len(self.sequence))


class PrimerDimer(Sequence):
	def __init__(self, sequence, comparesequence, pdlength):
		super(PrimerDimer, self).__init__(sequence)
		self.comparesequence = comparesequence
		self.pdlength = pdlength

	def basiccompare(self):
		"""
		FIXME: this needs both a docstring AND to be fixed, currently only 
		comparing the first n number of characters for complementarity. 

		Currently it needs:
			> windowing for different-sized primers AND 
			> to properly count complementary bases (not just to n)
		"""
		p1len = len(self.sequence)
		p2len = len(self.comparesequence)

		if self.complement()[p1len - self.pdlength:p1len] ==\
		self.comparesequence[0:self.pdlength]:
			print "primer-dimer detected between primers..." 
			print "fwd 3': %s" % (self.sequence[p1len - self.pdlength:p1len])
			print "rev 3': %s" % (self.comparesequence[0:self.pdlength])

	def pdlocal(self): 
		"""
		Compares sequence and comparesequence for the highest number of 
		identical characters within a single alignment of the two strings.

		Output is a list with: 
			[sequence, 
			 complement(comparesequence), 
			 score, 
			 beginning-of-alignment, 
			 end-of-alignment]

		Use of the static method 'format_alignment_compl' to format the 
		output is highly recommended to display the alignments.
		"""
		# The negative penalties for gaps in the sequence being set to the
		# size of the sequence means local sequences with gaps can NEVER
		# have higher scores than those without gaps (exactly the behaviour
		# we are looking for).
		for a in pairwise2.align.localxs(self.sequence,
										self.comparesequence,
										-len(self.sequence),
										-len(self.sequence)):
			if int(a[2]) >= self.pdlength:
				yield(a)

	@staticmethod
	def format_alignment_compl(align1, align2, score, begin, end):
		"""
		A modified version of biopython.pairwise2.format_alignment.

		Instead of returning the align2 argument, the complement of align2
		is returned. 

		This more correctly matches the expectation of the user
		who calls the primer-dimer method as they are looking for 
		complementarity rather than literal character matching (which is what
		the pairwise2 algorithm ACTUALLY bases it's returned scores off).
		"""
		### FIXME: If this type of object setup and method-call are expensive
		###	then perhaps passing the align2 complement in a parameter would
		### be better?
		s = []
		compalign2 = Sequence(align2).complement()

		s.append("%s\n" % align1)
		s.append("%s%s\n" % (" "*begin, "|"*(end-begin)))
		s.append("%s\n" % compalign2)
		s.append("  Score=%g\n" % score)
		return ''.join(s)


class Hairpin(Sequence):
	def __init__(self, sequence):
		super(Hairpin, self).__init__(sequence)

	def simpledetecthairpin(self, minlength):
		"""
		Detects potential hairpins of n-size within a given sequence using
		literal matching of substrings within the sequence.

		Args:
		minlength --> size of minimum hairpin to search for.
		"""
		revcomseq = self.reversecomplement()
		revseq = self.reversesequence()
		"""
		Starts hairpinsize at minlength (user defined minimum length of 
		hairpins) and increases till window is len(sequence) - 1
		"""
		### FIXUP: Window should never have to be more than half the size of
		### the sequence(?) as it needs to stringmatch...
		###		> Perhaps it is (windowsize / 2) + minlength as the max 
		###	  	  for hairpinsize???
		hairpinsize = minlength
		while hairpinsize < len(self.sequence):
			"""
			Creates a list of substrings of self.sequence of size minlength.
			"""
			i = 1
			subseqlist 	= []
			while i <= ((len(self.sequence) - hairpinsize) + 1):
				subseqlist.append(self.sequence[(i - 1):
					(i + hairpinsize) - 1])
				i += 1
			"""
			For each substring in subseqlist, window over sequence to find 
			matches.
			"""
			for numinlist, seqstring in enumerate(subseqlist):
				"""
				j is number of hairpinsize-sized windows over complement of
				sequence.
				"""
				j = 0
				while j < ((len(self.sequence) - hairpinsize) + 1):
					subseq 	= revseq[(j - 1):(j + hairpinsize) - 1]
					subrevseq = revcomseq[(j - 1):(j + hairpinsize) - 1]
					seqminone = len(self.sequence) - 1
					fwdrange = range(numinlist, numinlist + hairpinsize)
					revrange = range(seqminone - ((j + hairpinsize) - 1),
						seqminone - (j - 1))
					
					if seqstring == subrevseq and\
					not list(set(fwdrange) & set(revrange)):
						yield [numinlist,
							numinlist + hairpinsize,
							seqminone - ((j + hairpinsize) - 1),
							seqminone - (j - 1),
							seqstring,
							subseq]
					j += 1
			hairpinsize += 1

	def detectlocalhairpin(self, minlength):
		for a in pairwise2.align.localxs(self.sequence,
										self.comparesequence,
										-len(self.sequence),
										-len(self.sequence)):
			if int(a[2]) >= self.pdlength:
				yield(a)

	@staticmethod
	def complement_positions(top, bottom):
		positions = []
		for i, (t, b) in enumerate(zip(top, bottom)):
			if iscomplement(t, b):
				positions.append(i)
		return positions


# GLOBAL VARS
_complements = {'-': '-', 'N': 'N', 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def iscomplement(x, y):
	return _complements[x] == y

def complementbase(base):
	if base in _complements:
		return _complements[base]
	else:
		exit("found bad base %s" % base)