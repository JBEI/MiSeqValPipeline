#! /usr/bin/env python

def main(parser):
  args = parser.parse_args()

  quals = {}
  lines = (l.strip('\n') for l in args.fastq)
  try:
    while lines:
      name = lines.next()
      seq  = lines.next()
      blank = lines.next()
      qual = lines.next()
      quals[name.strip('@')] = qual
  except StopIteration:
    pass

  lines = (l.strip('\n').split('\t') for l in args.samfile)
  for l in lines:
    if not l[0].startswith('@'):
      readid =  '/'.join(l[0].split('/')[:-1])
      assert readid in quals
      assert len(l[9])==len(quals[readid])
      l[10] = quals[readid]
    print >>args.outfile, '\t'.join(l)    

if __name__ == '__main__':
  import sys
  import argparse
  parser = argparse.ArgumentParser(prog='add_qual_to_sam', description='Add quality scores from FASTQ file to SAM file')
  parser.add_argument('fastq', type=argparse.FileType('rU'), help='fastq with qualities')
  parser.add_argument('samfile', type=argparse.FileType('rU'), help='samfile to be added to')
  parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),default=sys.stdout)
  main(parser)