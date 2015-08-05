#! /usr/bin/env python

from pbcore.io import FastqWriter, FastqRecord
from pbcore.io import FastaWriter, FastaRecord
from pbcore.io import BasH5Reader
from pbcore.io.BarcodeH5Reader import BarcodeH5Reader
from collections import defaultdict

def main(parser):
  args = parser.parse_args()

  # Get outfile name
  if args.outFile is None:
    outfile = 'nobarcode.fasta' if args.fasta else 'nobarcode.fastq'
  else:
    outfile = args.outFile
  
  # Input files
  barcodeFofn = (l.strip('\n') for l in args.barcode_fofn)
  baxFofn = (l.strip('\n') for l in args.bax_fofn)
  
  # Get the read names that are not barcoded
  no_barcode = defaultdict(set)
  for barcodeFile in barcodeFofn:
    bcH5 = BarcodeH5Reader(barcodeFile)
    for row in bcH5.bestDS:
      if row[3] / row[1] < args.minAvgBarcodeScore:
        no_barcode[bcH5.movieName].add(row[0])

  if args.fasta:
    outh = FastaWriter(outfile)
  else:
    outh = FastqWriter(outfile)

  for baxFile in baxFofn:
    baxH5 = BasH5Reader(baxFile)
    for holeNum in baxH5.sequencingZmws:
      if holeNum in no_barcode[baxH5.movieName]:
        zmw = baxH5[holeNum]
        if len(zmw.subreads) and max(len(sr.basecalls()) for sr in zmw.subreads) >= args.minMaxInsertLength:
          for subread in zmw.subreads:
            if len(subread.basecalls()) >= args.minSubreadLength:
              if args.fasta:
                outh.writeRecord(FastaRecord(subread.readName,subread.basecalls()))
              else:
                outh.writeRecord(FastqRecord(subread.readName,subread.basecalls(),subread.QualityValue()))

  outh.close()

if __name__ == '__main__':
  import sys
  import argparse
  parser = argparse.ArgumentParser(prog='fasta', description='Create FASTQ file with barcode sequences')
  parser.add_argument('--outFile', type=str)  
  parser.add_argument('--minMaxInsertLength', type=int, default=0)
  parser.add_argument('--minAvgBarcodeScore', type=int, default=0)
  parser.add_argument('--minSubreadLength', type=int, default=0) 
  parser.add_argument('--fasta', action='store_true')  
  parser.add_argument('bax_fofn', type=argparse.FileType('r'))
  parser.add_argument('barcode_fofn', type=argparse.FileType('r'))
  main(parser)
