#!/usr/bin/env python

"""
Author: Luke Shillabeer

Description: PrimerLint is a primer analysis tool with focus on 
determining potentially problematic primers. 

PrimerLint is distributed as part of the hiplex-primer software bundle.

primerlint.py is the point-of-execution for the analysis software.
"""

import importfiles
import analysismodules

csvimport = importfiles.ImportCSV("testcsv.csv")
csvimport.close()
testlist = csvimport.process()

primerdimertest = analysismodules.PrimerDimer(
					'GGGGGGGGGGTGGGGGGGGG',
					'ACCCCCCCCCCCCCCCCCCC',
					10)

primerdimertest.basiccompare()