##2-1) generating snippy input from simbac simulations (we need reference and individual fasta sequences)
#2-1-1)since simbac doesn't generate referene we will arbitrarily choose sequence 0 as reference
#2-1-2)Snippy accepts inputs as contig fasta file but for them to be recognizable the extenstion must be '.fa',
    #otherwise, (e.g. with .fasta extension), they will be treated as reads

def get_options():
    import argparse

    description = 'Convert simBac simulations to reference and fasta sequences'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--infile',
                        help='simbac simulation outfile')
    parser.add_argument('--outdir',
                        help='directory for the output files')
    return parser.parse_args()

options = get_options()

import os
from subprocess import call

def simbac2snippyPlugger(infile,outdir):
  txt=open(os.path.join(outdir,'reference.fa'),'w')
  with open(infile) as file:
    l_=file.readline()
    txt.write('>25\n')
    l_=file.readline()
    txt.write(l_)
  with open(infile)as file:
    line=file.readline()
    line=file.readline()
    line=file.readline()
    while line:
      id=line[1:].strip()
      txt=open(os.path.join(outdir,id+'.fa'),'w')
      txt.write('>25\n') #chromosome name should be same for all samples for varcall compatibility, '25' for mitochondrial
      txt.write(file.readline())
      txt.close()    
      line=file.readline()
    
simbac2snippyPlugger(options.infile,options.outdir)