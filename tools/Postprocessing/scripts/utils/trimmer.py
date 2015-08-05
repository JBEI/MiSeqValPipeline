#! /usr/bin/env python

import sys
import argparse

def main(parser):
  args = parser.parse_args()

  from Bio import SeqIO
  type = 'fastq' if not args.fasta else 'fasta'
  trim = args.trim
  print >>sys.stderr, "%d sequences written" % SeqIO.write((s[trim:-trim] for s in SeqIO.parse(args.infile,type)),args.outfile,type)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(prog='generate_barcode_fasta', description='Create FASTA file with barcode sequences')
  parser.add_argument('--trim', type=int, default=20)
  parser.add_argument('--fasta', action='store_true')
  parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),default=sys.stdin)
  parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)
  main(parser)
