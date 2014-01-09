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
		self.sequence = sequence

	def complement(self):
		return ''.join(map(complementbase, self.sequence))

	def reversecomplement(self):
		return ''.join(map(complementbase, self.sequence[::-1]))

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
		complmentseq = self.complement()
		subseqlist = []
		
		""" 
		Creates a list of substrings of self.sequence of size minlength.
		FIXUP: Requires changes to make size of hairpin count up/down (???)
		"""
		i = 1
		while i <= ((len(self.sequence) - minlength) + 1):
			subseqlist.append(self.sequence[(i - 1):(i + minlength) - 1])
			i += 1

		"""
		For each substring in subseqlist, window over sequence to find 
		matches.
		"""
		for numinlist, seqstring in enumerate(subseqlist):
			"""
			j is number of minlength-sized windows over complement of
			sequence.
			"""
			j = ((len(self.sequence) - minlength) + 1)
			while j >= 0:
				if seqstring == complmentseq[(j - 1):(j + minlength) - 1]:
					compseq = Sequence(complmentseq[(j - 1):(j + minlength) - 1]).complement()
					yield [	numinlist,
							numinlist + 4,
							(j - 1),
							((j + minlength) - 1)]
				j -= 1

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