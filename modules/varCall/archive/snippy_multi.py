#2-2) running snippy-multi on simulations

def get_options():
    import argparse

    description = 'Running snippy-multi on simulations'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--infile',
                        help='directory to tab-separated simulation files')
    parser.add_argument('--ref',
                        help='reference genome')
    parser.add_argument('--cpus',
                        default=8,
                        help='number of cpus for parallelization')
    parser.add_argument('--outdir',
                        help='directory for the output files')
    return parser.parse_args()

options = get_options()

import os
from os import chdir
from subprocess import call
def snippy(infile,ref,outdir,cpus):
    #run snippy-multi
  callString='snippy-multi %s --ref %s --cpus %s > %s' % (infile,ref,cpus,os.path.join(outdir,'runme.sh'))
  call(callString, shell=True)
  #changing directory to outdir
  chdir(outdir)
  #run all
  call('sh ./runme.sh', shell=True)
    
snippy(options.infile,options.ref,options.outdir,options.cpus)