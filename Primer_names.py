#!/bin/env python

'''
Primer design program.

Authors: Bernie Pope, Danny Park, Tu Ng.
'''

from Bio import SeqIO
import sys
import math
import os
import csv
import itertools
import getopt

# print a usage message
def usage():
    print("""Usage: %s
    [-h | --help]
    --refdir=<reference directory>
    --genes=<gene list CSV file>
    --blocksize=<block size>
    --maxprimersize=<primer size>
    --primervar=<primer size variance>
    --splicebuffer=<splice site buffer size>
    --melt=<ideal Tm>
    --idtcsv=<CSV output file for IDT>
    --log=<log file name>""") % sys.argv[0]

# check that a command line argument is a positive integer
# given the name of the argument and its value as a string
def getPositiveIntArg(name, str):
    if str.isdigit():
        val = int(str)
        if val > 0:
            return val
    print ("%s must be a positive integer" % name)
    sys.exit(1)

longOptionsFlags = ["help", "refdir=", "genes=", "blocksize=",
                    "maxprimersize=", "primervar=", "splicebuffer=", "melt=", "log=",
                    "idtcsv="]
shortOptionsFlags = "h"

# A place to store command line arguments.
class Options(object):
    def __init__(self):
        self.reference_dirname = None
        self.gene_filename = None
        self.block_size = None
        self.max_primer_size = None
        self.primer_variance = None
        self.splice_buffer = None
        self.ideal_melt_temp = None
        self.log_filename = None
        self.idt_filename = None
    def check(self):
        return (
            self.reference_dirname and
            self.gene_filename and
            self.block_size and
            self.max_primer_size and
            self.primer_variance and
            self.splice_buffer and
            self.ideal_melt_temp and
            self.log_filename and
            self.idt_filename
       )

def main():
    # process command line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptionsFlags, longOptionsFlags)
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    options = Options()
    for o, a in opts:
        if o == "--refdir":
            options.reference_dirname = a
        elif o == "--genes":
            options.gene_filename = a
        elif o == "--blocksize":
            options.block_size = getPositiveIntArg("block size", a)
        elif o == "--maxprimersize":
            options.max_primer_size = getPositiveIntArg("primer size", a)
        elif o == "--primervar":
            options.primer_variance = getPositiveIntArg("primer variance", a)
        elif o == "--splicebuffer":
            options.splice_buffer = getPositiveIntArg("splice buffer", a)
        elif o == "--melt":
            options.ideal_melt_temp = getPositiveIntArg("ideal melting temp", a)
        elif o == "--log":
            options.log_filename = a
        elif o == "--idtcsv":
            options.idt_filename = a
        elif o in ('-h', '--help'):
            usage()
            sys.exit(0)
    if not options.check():
        print('Incorrect arguments')
        usage()
        exit(2)
    with open(options.log_filename, 'w') as options.log_file:
        # process each exon in the gene file
        csv_file = open(options.idt_filename, 'w')
        for gene_name, exon_id, coords in process_gene_file(options.gene_filename):
            try:
                chr, exon_start, exon_end = parse_coords(coords)
            except:
                exit('bad coordinates: %s' % coords)
            else:
                # score all the sliding windows over this exon
                window_scores = score_exon_windows(options, chr, gene_name, exon_id, exon_start, exon_end)
                # find the best scoring window for this exon
                best_window = get_best_window(options, window_scores)
                # print out the primers for the best window
		print_best_primers(csv_file, gene_name, exon_id, chr, exon_start, exon_end, window_scores[best_window])
	csv_file.close()

# parse the coordinates of a chromosome position in the format:
# chr:start-end
def parse_coords(coords):
    chr, positions = coords.split(':')
    start, end = positions.split('-')
    return (chr, int(start), int(end))

# read the contents of the exon file.
# transpose columns to rows and parse the important fields.
# correct the exon numbering so that they are sequential from the first numbered exon in a gene
def process_gene_file(gene_filename):
    with open(gene_filename, 'U') as gene_file:
        reader = csv.reader(gene_file)
        for row in itertools.izip_longest(*reader):
            if len(row) >= 1:
                gene_name = row[0]
                index = 1
                row_len = len(row)
                # set exon_number to None to indicate that we haven't yet seen an exon for this gene
                exon_number = None
                while index < row_len:
                    item = row[index]
                    # find the start of the block of text that defines an exon for this gene
                    if type(item) == str and item.startswith('Exon number:'):
                        exon_number_words = item.split()
                        if len(exon_number_words) == 3:
                            # if this is the first exon in the gene, then use its number from the file
                            if exon_number == None:
                                exon_number = int(exon_number_words[2])
                            # otherwise just increment the number from the previous exon to correct for possible mistakes in the file
                            else:
                                exon_number += 1
                            if index+2 < row_len:
                                coords = row[index+2]
                                yield (gene_name, exon_number, coords)
                    index += 1

# write a string to the log file
def log_msg(options, msg):
    options.log_file.write('%s\n' % msg)

# width of banner message separator
banner_width = 80

class ScoredBlock(object):
    def __init__(self, block_num, start, end, primer_forward, primer_reverse):
        self.block_num = block_num # sequential number of the block in the window
        self.start = start # first coordinate in the block
        self.end = end     # last coordinate in the block
        self.primer_forward = primer_forward
        self.primer_reverse = primer_reverse

class ScoredPrimer(object):
    def __init__(self, score, direction, start, end, bases):
        self.score = score
        self.direction = direction
        self.start = start
        self.end = end
        self.bases = bases

# Compute the best primers for an exon, using a sliding window.
#
# |.................Region...................|
#
#   MP   SL   SB      Exon      SB   SL   MP
# |----|----|----|------------|----|----|----|
#
#      ^========^========^=========^
#        block1   block2   block3
#
#      |..........Window...........|
#
# MP = max primer size
# SL = slack
# SB = splice buffer
#
# Window is num_blocks * block_size.
#
# Region = exon + (2 * slack) + (2 * max_primer)
#
# A region is the total chunk of DNA that we read from the fasta file.
#
# We slide the window along slack+1 steps. For each window position
# we find the best primers for each block in that position.

def score_exon_windows(options, chr, gene_name, exon_id, exon_start, exon_end):
    exon_size = (exon_end - exon_start) + 1
    # add the splice buffer on to the start and end of the exon coordinates
    exon_buffer_start = exon_start - options.splice_buffer
    exon_buffer_end = exon_end + options.splice_buffer
    exon_buffer_size = (exon_buffer_end - exon_buffer_start) + 1
    num_blocks = int(math.ceil(float(exon_buffer_size) / float(options.block_size)))
    window_size = num_blocks * options.block_size
    slack = window_size - exon_buffer_size
    # the region must include the exon, plus the splice buffer, plus slack on either side, plus maximum primers on either side
    region_start = exon_buffer_start - slack - options.max_primer_size
    region_end = exon_buffer_end + slack + options.max_primer_size

    log_msg(options, '*' * banner_width)
    log_msg(options, 'chrom:\t\t%s' % chr)
    log_msg(options, 'exon:\t\t%s' % exon_id)
    log_msg(options, 'exon start:\t%d' % exon_start)
    log_msg(options, 'exon buff start:%d' % exon_buffer_start)
    log_msg(options, 'exon end:\t%d' % exon_end)
    log_msg(options, 'exon buff end:\t%d' % exon_buffer_end)
    log_msg(options, 'exon size:\t%d' % exon_size)
    log_msg(options, 'exon buff size:\t%d' % exon_buffer_size)
    log_msg(options, 'block_size:\t%d' % options.block_size)
    log_msg(options, 'num_blocks:\t%d' % num_blocks)
    log_msg(options, 'window_size:\t%d' % window_size)
    log_msg(options, 'slack:\t\t%d' % slack)
    log_msg(options, 'region_start:\t%d' % region_start)
    log_msg(options, 'region_end:\t%d' % region_end)
    log_msg(options, '')

    region = get_region(options, chr, region_start, region_end)

    # get the first and last 5 bases in the region so we can print it out for diagnostic purposes.
    if len(region.bases) > 10:
       region_start_bases = region.bases[0:5]
       region_end_bases = region.bases[-5:]
       log_msg(options, 'region starts with: %s, region ends with: %s' % (region_start_bases, region_end_bases))

    window_scores = {}
    window_count = 1
    # slide the window over the exon, starting "slack" number of positions before the start of the exon (with splice buffer added).
    # the last position of the window is when it lines up with the start of the exon.
    for window_start in range (exon_buffer_start - slack, exon_buffer_start + 1):
        log_msg(options, '=' * banner_width)
        log_msg(options, "Scoring window starting at %d (%d/%d)" % (window_start, window_count, slack+1))
        window_scores[window_start] = []
        # score each block in the window, trying the different primer sizes at each end
        for block_num in range(0, num_blocks):
            block_start = window_start + (block_num * options.block_size)
            block_end = block_start + options.block_size - 1

            log_msg(options, '-' * banner_width)
            log_msg(options, 'block:\t\t%d' % block_num)
            log_msg(options, 'block start:\t%d' % block_start)
            log_msg(options, 'block end:\t%d' % block_end)

            # find the best primer for either the start or the end of a block
            def find_best_primer(direction, primer_end):

                log_msg(options, '~' * (banner_width / 2))
                log_msg(options, "Primer dir:\t%s" % direction)
                if direction == "forward":
                    primer_start = (primer_end - options.max_primer_size) + 1
                    primer = get_primer(region, primer_start, primer_end)
                elif direction == "reverse":
                    primer_start = (primer_end + options.max_primer_size) - 1
                    # the reverse primer comes from the opposite strand
                    primer = get_primer(region, primer_end, primer_start)
                    primer = reverse_complement(primer)

                log_msg(options, 'Primer start:\t%d' % primer_start)
                log_msg(options, 'Primer end:\t%d' % primer_end)
                log_msg(options, '')

                primer_len = len(primer)

                # check that the primer length is what we believe it to be
                assert primer_len == (abs(primer_end - primer_start) + 1)

                # score all suffixes of the primer (we try the first N suffixes, where N = primer_variance)
                best_primer = None
                for suffix_start in range(0, options.primer_variance):
                    primer_suffix = primer[suffix_start:]
                    primer_score = score_primer(options.ideal_melt_temp, primer_suffix)
                    log_msg(options, "Score: %4d, 5'> %s <3'" % (primer_score, primer_suffix))
                    # calculate the start position of this primer candidate
                    if direction == 'forward':
                        candidate_start = primer_start + suffix_start
                    elif direction == 'reverse':
                        candidate_start = primer_start - suffix_start
                    candidate = ScoredPrimer(primer_score, direction, candidate_start, primer_end, primer_suffix)
                    if best_primer == None:
                        # this is the first one we've seen
                        best_primer = candidate
                    elif primer_score == best_primer.score:
                        # we have a tie for best primer, choose the shortest one
                        if len(candidate.bases) < len(best_primer.bases):
                            best_primer = candidate
                    elif primer_score < best_primer.score:
                        # this is the new best score
                        best_primer = candidate
                log_msg(options, '')
                log_msg(options, 'Best primer is: %s' % best_primer.bases)
                return best_primer

            # find the best forward and reverse primers for this block
            # using the coordinate of the last position in the primer to
            # indicate the primer position
            primer_forward = find_best_primer('forward', block_start-1)
            primer_reverse = find_best_primer('reverse', block_end+1)
            if primer_forward != None and primer_reverse != None:
                #block_score = (block, block_start, block_end, primer_forward, primer_reverse)
                block_score = ScoredBlock(block_num, block_start, block_end, primer_forward, primer_reverse)
                window_scores[window_start].append(block_score)
            else:
               print("** Warning **: did not find best primers for window start: %d, block start: %d" %
                         (window_start, block_start))
        window_count += 1
    return window_scores

# given all the scores for each position of the sliding window, find the best one
def get_best_window(options, window_scores):
    best_window = None
    best_score = None
    for window, scoredBlocks in window_scores.items():
        total_score = get_total_window_score(scoredBlocks)
        log_msg(options, "Window: %d, total score: %d" % (window, total_score))
        if best_window == None:
            # this is the first one we've seen
            best_window = window
            best_score = total_score
        elif total_score == best_score:
           # XXX we have a tie for the best, what to do?
           log_msg(options, "Tie for best window between %d and %d" % (best_window, window))
        elif total_score < best_score:
            best_window = window
            best_score = total_score
    log_msg(options, "Best window: %d, score: %d" % (best_window, best_score))
    return best_window

def get_total_window_score(scoredBlocks):
    score = 0
    for block in scoredBlocks:
        score += block.primer_forward.score + block.primer_reverse.score
    return score

def print_best_primers(csv_file, gene_name, exon_id, chr, exon_start, exon_end, scored_blocks):
    primer_name_prefix = gene_name + '_X' + str(exon_id) + '_'
    print('-' * banner_width)
    print('gene: %s, exon: %s, %s:%d-%d' % (gene_name, exon_id, chr, exon_start, exon_end))
    for block in scored_blocks:
        print('block %d, %d-%d' % (block.block_num, block.start, block.end))
        forward = block.primer_forward
        print('forward: %d-%d, %s' % (forward.start, forward.end, forward.bases))
        reverse = block.primer_reverse
        print('reverse: %d-%d, %s' % (reverse.start, reverse.end, reverse.bases))
        primer_name_forward = primer_name_prefix + 'F' + str(block.block_num + 1) 
        primer_name_reverse = primer_name_prefix + 'R' + str(block.block_num + 1) 
        csv_file.write(','.join([primer_name_forward, str(forward.bases), '250nM', 'HPLC']) + '\n')
	csv_file.write(','.join([primer_name_reverse, str(reverse.bases), '250nM', 'HPLC']) + '\n')
        csv_file.flush()

# a (possibly large) chunk of DNA from the reference
class Region(object):
    def __init__(self, chr, start, end, file, bases):
        self.chr = chr
        self.start = start
        self.end = end
        self.file = file
        self.bases = bases

# SeqIO uses 0 based indexes, but our input data uses 1 based indexes,
# so we need to convert here
def get_region(options, chr, start, end):
    ref_filename = os.path.join(options.reference_dirname, chr + ".fa")
    log_msg(options, "Reading region on %s, start: %d, end: %d from file %s" % (chr, start, end, ref_filename))
    with open(ref_filename) as ref_file:
        # XXX should really cache this
        bases_list = list(SeqIO.parse(ref_file, "fasta"))
        if len(bases_list) == 1:
            bases = bases_list[0]
            if start > 0 and end <= len(bases):
                # here is where we do index conversion from 1 based to 0 based
                segment = bases[start-1:end].seq
                normalised_segment = segment.upper()
                validate_sequence(normalised_segment)
                return Region(chr, start, end, ref_filename, normalised_segment)
            else:
                exit("chromosome %s does not span %d %d" % (chr, start, end))
        exit("got wrong number of sequences in bases_list: %d" % len(bases_list))

def get_primer(region, start, end):
    new_start = start - region.start
    new_end = end - region.start
    return region.bases[new_start:new_end+1]

# check our sequence only has A,T,G,C
def validate_sequence(sequence):
    for base in sequence:
        if base not in "ATGC":
            exit("bad base found: %s" % base)
    return sequence

def reverse_complement(sequence):
    return ''.join(map (complement,sequence[::-1]))

complements = { 'N': 'N', 'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G' }
def complement(base):
    if base in complements:
        return complements[base]
    else:
        exit("found bad base %s" % base)

def score_primer(ideal_melting_temp, sequence):
    t_m = melting_temp(sequence)
    return (ideal_melting_temp - t_m) ** 2

base_temp = { 'A' : 2, 'T' : 2, 'G' : 4, 'C' : 4 }
def melting_temp(sequence):
    t_m = 0
    for base in sequence:
        t_m += base_temp[base]
    return t_m

if __name__ == '__main__':
    main()
