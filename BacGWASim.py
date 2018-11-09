# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 17:55:45 2018

@author: Morteza.M Saber

"""
import subprocess

def BacGWASim(Snakefile,CoreNumber=1,DirectedAcyclicGraph=True,Rulegraph=False,quiet=False):
    command='snakemake --snakefile %s --cores %s --latency-wait 120 --rerun-incomplete'%(Snakefile,CoreNumber)
    if DirectedAcyclicGraph:
        command+=' --dag | dot -Tpng > %s.png'%(Snakefile+'_DAG')
    if Rulegraph:
        command+=' --rulegraph| dot -Tpng > %s.png'%(Snakefile+'_RuleGraph')
    if quiet:
        command+=' --quiet'
    SettingPhylogeny=Snakefile+'Phylogeny'
    CommandPhylogeny='snakemake --snakefile %s --cores %s '%(SettingPhylogeny,CoreNumber)
    unlocker='snakemake --snakefile %s --unlock'%Snakefile
    import timeit
    start = timeit.default_timer()
    subprocess.call(CommandPhylogeny,shell=True)
    subprocess.call(unlocker,shell=True)
    subprocess.call(command,shell=True)
    stop = timeit.default_timer()
    total_time = stop - start
    mins, secs = divmod(total_time, 60)
    hours, mins = divmod(mins, 60)
    print("Simulation Completed in: %d:%d:%d.\n" % (hours, mins, secs))

    
    
from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-S','--Snakefile', dest='Snakefile',help='Absolute path to the setting file [required]',type='str',default='')
parser.add_option('-J','--CoreNumber', dest='CoreNumber',help= 'Use at most N cores in parallel (default: 1)',type='int',default=1)
parser.add_option('-D',action="store_true", dest='DirectedAcyclicGraph',help= 'Print the directed acyclic\
                        graph of jobs. (Default:False)',default=False)
parser.add_option('-R',action="store_true", dest='Rulegraph',help="Print the dependency graph of rules. *Caution*Rulegraph can not be set to True if DirectedAcyclicGraph is set to True (default: False)",default=False)
parser.add_option('-Q',action="store_true", dest='quiet',help= 'Do not output any progress or rule information (default: False)',default=False)

(options,args)=parser.parse_args()
Snakefile=options.Snakefile
CoreNumber=options.CoreNumber
DirectedAcyclicGraph=options.DirectedAcyclicGraph
Rulegraph=options.Rulegraph
quiet=options.quiet

if DirectedAcyclicGraph==True  and Rulegraph==True:
	print ("\nSimulatorPipe For Simulating Bacterail Whole Genomes \nOnly one of 'DirectedAcyclicGraph' and 'Rulegraph' can be used\n")
	parser.print_help()
	sys.exit(1)
if (Snakefile=='') :
	print ('\nSimulatorPipe For Simulating Bacterail Whole Genomes \nRequired filed(s) not supplied\n')
	parser.print_help()
	sys.exit(1)
import ast
input_option_dict=ast.literal_eval(options.__str__())
print ('Entered arguments are...')
for one in input_option_dict:
	print (one,input_option_dict[one])
print ('')

BacGWASim(Snakefile,CoreNumber,DirectedAcyclicGraph,Rulegraph,quiet)

