configfile: "/home/masih/BacterialSimulator/Pipeline/ConfigFile.yaml"
outfile=expand("{outputDIR}/{PhylogeneticTreeDIR}/RealTree.nwk",outputDIR=config["outputDIR"],PhylogeneticTreeDIR=config["PhylogeneticTreeDIR"])

#import sys
#print(sys.executable) #To test which python is executing the job
from Bio import Phylo
tree= Phylo.read(str(outfile[0]), 'newick')
names=sorted([item.split(')')[0].replace('=','').replace("'",'') for item in str(tree).split('name')[1:]])