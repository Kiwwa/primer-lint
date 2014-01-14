#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

analysismodules.py contains the classes that represent different analysis 
tasks computed by the PrimerLint software.
"""

class Sequence(object):
	def __init__(self, sequence):
		self.sequence = str(sequence).upper()

	def complement(self):
		return ''.join(map(complementbase, self.sequence))

	def reversecomplement(self):
		return ''.join(map(complementbase, self.sequence[::-1]))

	def reversesequence(self):
		return self.sequence[::-1]

	def gcpercent(self):
		return 	float(self.sequence.count('G') +\
				self.sequence.count('C')) / float(len(self.sequence))


class PrimerDimer(Sequence):
	def __init__(self, primersequence, compareprimer, pdlength):
		super(PrimerDimer, self).__init__(primersequence)
		self.compareprimer = compareprimer
		self.pdlength = pdlength

	def basiccompare(self):
		p1len = len(self.sequence)
		p2len = len(self.compareprimer)

		if self.complement()[p1len - self.pdlength:p1len] ==\
		self.compareprimer[0:self.pdlength]:
			print "primer-dimer detected between primers..." 
			print "fwd 3': %s" % (self.sequence[p1len - self.pdlength:p1len])
			print "rev 3': %s" % (self.compareprimer[0:self.pdlength])


class Hairpin(Sequence):
	def __init__(self, sequence):
		super(Hairpin, self).__init__(sequence)

	def simpledetecthairpin(self, minlength):
		"""
		Detects potential hairpins of n-size within a given sequence using
		literal matching of substrings within the sequence.

		Args:
		minlength --> size of minimum hairpin to search for 
		"""
		revcomseq 	= self.reversecomplement()
		revseq 		= self.reversesequence()
		"""
		Starts hairpinsize at minlength (user defined minimum length of 
		hairpins) and increases till window is len(sequence) - 1
		
		FIXUP: Window should never have to be more than half the size of
		the sequence(?) as it needs to stringmatch...
			> Perhaps it is (windowsize / 2) + minlength as the max 
			  for hairpinsize???
		"""
		hairpinsize = minlength
		while hairpinsize < len(self.sequence):
			"""
			Creates a list of substrings of self.sequence of size minlength.
			"""
			i = 1
			subseqlist 	= []
			while i <= ((len(self.sequence) - hairpinsize) + 1):
				subseqlist.append(	self.sequence[(i - 1):
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
					subseq 		= revseq[(j - 1):(j + hairpinsize) - 1]
					subrevseq 	= revcomseq[(j - 1):(j + hairpinsize) - 1]
					seqminone 	= len(self.sequence) - 1
					fwdrange 	= range(numinlist,
										numinlist + hairpinsize)
					revrange 	= range(seqminone - ((j + hairpinsize) - 1),
										seqminone - (j - 1))
					
					if seqstring == subrevseq and\
					not list(set(fwdrange) & set(revrange)):
						yield [	numinlist,
								numinlist + hairpinsize,
								seqminone - ((j + hairpinsize) - 1),
								seqminone - (j - 1),
								seqstring,
								subseq]
					j += 1
			hairpinsize += 1

	@staticmethod
	def complement_positions(top, bottom):
		positions = []
		for i, (t, b) in enumerate(zip(top, bottom)):
			if iscomplement(t, b):
				positions.append(i)
		return positions


# GLOBAL VARS
_complements = {'N': 'N', 'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

def iscomplement(x, y):
	return _complements[x] == y

def complementbase(base):
	if base in _complements:
		return _complements[base]
	else:
		exit("found bad base %s" % base)