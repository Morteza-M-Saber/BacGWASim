'''Transfer HGT genes to their destination nodes'''


def get_options():
    import argparse

    description = 'Transfer HGT genes to their destination nodes'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--HGTNode',
                        help='The nodes randomly chosen as the destination of HGT')
    parser.add_argument('--HGTRate',
                        default='10',
                        help="Frequency of Horizontal Gene Transfer events [Default: 10]")
    parser.add_argument('--HGTSource',
                        help="Absolute path to the directory where HGTs have been simulated")
    parser.add_argument('--HGTDestination',
                        help="Absolute path to the Assembly directory where the HGTs will be transfered to the corresponding fasta genome")


    return (parser.parse_args())


options = get_options()


#configfile: "/home/masih/Baterial_simulator/Pipeline/ConfigFile.yaml"
import os
import json
from Bio import SeqIO
import re

def HGTTransfer (HGTNode,HGTRate,HGTSource,HGTDestination):
    #Importing randomly chosen nodes as destination of HGT
    with open (HGTNode) as f:
        NodeDict=json.load(f)
    NodeList=[item for item in NodeDict]
    #Transfer the simulated genes to their destination
    
    for HGTs in range(int(HGTRate)):
        TransferDict={}
        for SelectedSpecs in NodeDict[NodeList[HGTs]]:
            TransferDict[SelectedSpecs]=[]
            with open (os.path.join(HGTSource,str(HGTs),'ProtSim','DB',str(SelectedSpecs)+'_dna.fa'),'r') as file:
                for records in SeqIO.parse(file,'fasta'):
                    TransferDict[SelectedSpecs].append([records.id,str(records.seq)])
        for TargetSpecs in NodeDict[NodeList[HGTs]]:
            with open (os.path.join(HGTDestination,str(TargetSpecs)+'.fsa'),'a') as txt:
                for Genes in TransferDict[TargetSpecs]:
                    txt.write('>%s\n%s\n'%(Genes[0],Genes[1]))
        #Writing HGTS in human-readable format
    txt=open(os.path.join(HGTDestination,'HGTTransfer.txt'),'w')
    txt.write('Horizontal Gene Transfer events summary\n')
    txt.write('TranferedGene(Uniprot ID)\tSpeciesTransfered\n')
    for HGTs in range(int(HGTRate)):
        with open (os.path.join(HGTSource,str(HGTs)+'.db'),'r') as file:
            line=file.readline()
            result = re.search('<ID>(.*)</ID>', line)
            GeneName=result.group(1)
        txt.write('%s\t%s\n'%(GeneName,str(NodeDict[NodeList[HGTs]])))
    txt.close()
            
    
                    
                

HGTTransfer(options.HGTNode,options.HGTRate,options.HGTSource,options.HGTDestination)
        
        
        
    
    
        
        