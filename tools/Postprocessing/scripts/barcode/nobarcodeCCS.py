#! /usr/bin/env python

from pbcore.io import FastaWriter, FastaRecord
from pbcore.io import FastqWriter, FastqRecord
from pbcore.io import BasH5Reader
from pbcore.io.BarcodeH5Reader import BarcodeH5Reader

def main(parser):
  args = parser.parse_args()

  # Get outfile name
  if args.outFile is None:
    outfile = 'nobarcode.fasta' if args.fasta else 'nobarcode.fastq'
  else:
    outfile = args.outFile
  
  # Input files
  barcodeFofn = (l.strip('\n') for l in args.barcode_fofn)
  ccsFofn = (l.strip('\n') for l in args.ccs_fofn)  
  
  # Get the read names that are not barcoded
  no_barcode = set()
  for barcodeFile in barcodeFofn:
    bcH5 = BarcodeH5Reader(barcodeFile)
    for row in bcH5.bestDS:
      if row[3] / row[1] < args.minAvgBarcodeScore:
        no_barcode.add('%s/%d' % (bcH5.movieName,row[0]))
  
  if args.fasta:
    outh = FastaWriter(outfile)
  else:
    outh = FastqWriter(outfile)
  
  for ccsFile in ccsFofn:
    ccsH5 = BasH5Reader(ccsFile)
    for ccsRead in ccsH5.ccsReads():
      if ccsRead.zmw.zmwName in no_barcode:
        basecalls = ccsRead.basecalls()
        if len(basecalls) >= args.minMaxInsertLength:
          if args.fasta:
            outh.writeRecord(FastaRecord(ccsRead.zmw.zmwName, ccsRead.basecalls()))
          else:
            outh.writeRecord(FastqRecord(ccsRead.zmw.zmwName, ccsRead.basecalls(), ccsRead.QualityValue()))
  outh.close()

if __name__ == '__main__':
  import sys
  import argparse
  parser = argparse.ArgumentParser(prog='fasta', description='Create FASTQ file with barcode sequences')
  parser.add_argument('--outFile', type=str)  
  parser.add_argument('--minMaxInsertLength', type=int, default=0)
  parser.add_argument('--minAvgBarcodeScore', type=int, default=0)
  parser.add_argument('--fasta', action='store_true')
  parser.add_argument('ccs_fofn', type=argparse.FileType('r'))
  parser.add_argument('barcode_fofn', type=argparse.FileType('r'))
  main(parser)
