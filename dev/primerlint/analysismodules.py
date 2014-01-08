#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

analysismodules.py contains the classes that represent different analysis tasks
computed by the PrimerLint software.
"""

class PrimerDimer(object):
	def __init__(self, primerone, primertwo):
		self.primerone = primerone
		self.primertwo = primertwo

	def basiccompare(self):
