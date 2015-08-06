#! /usr/bin/env python

import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)
args = parser.parse_args()  
print >>sys.stderr, "%d sequences written" % SeqIO.write(SeqIO.parse(args.infile,'fastq'),args.outfile,'fasta')

