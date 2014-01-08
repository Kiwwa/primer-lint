#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

primerlint.py is the point-of-execution for the analysis software.
"""

import importfiles

csvimport = importfiles.ImportCSV("testcsv.csv")
csvimport.process()

testlist = csvimport.process()

for listitem in testlist:
	print listitem

