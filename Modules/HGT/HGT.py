'''Select randomly genes for HGT from a UNIPROT database'''


def get_options():
    import argparse

    description = 'Select randomly genes for HGT from a UNIPROT database'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--DBdir',
                        help='UNIPROT bacterial gene database in fasta format')

    parser.add_argument('--HGTRate',
                        default='10',
                        help="Frequency of Horizontal Gene Transfer events [Default: 10]")
    parser.add_argument('--Output',
                        help="Output file including the Protein/DNA sequences of HGT genes in darwin format")


    return (parser.parse_args())


options = get_options()

#Make a back codon-table to translate protein into DNA
codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

backcodon={}
for codon in codontable:
    if codontable[codon] not in backcodon:
        backcodon[codontable[codon]]=[]
    backcodon[codontable[codon]].append(codon)
backcodon['X']=list(codontable.keys()) #Undetermined codons in the database
    


def HGT (DBdir,HGTRate,Output):
    #Creating a DB of bacterial gene along with DNA and PROTEIN sequences
    from Bio import SeqIO
    import random
    DB={}
    with open (DBdir,'r') as file:
        for records in SeqIO.parse(file,'fasta'):
            DNA=str()
            for AA in str(records.seq):
                DNA+=random.choice(backcodon[AA])
            DB[records.id]=[str(records.seq),DNA]
            
    #Selecting genes equal to the number of HGTRate from database
    import numpy as np
    Genes=np.random.choice(list(DB.keys()),int(HGTRate),replace=False)
    
    GeneSelect={}
    for gene in Genes:
        GeneSelect[gene]=DB[gene]
    
    #writing the genes in darwin format usable by ALFSIM
    import os
    MyIter=0
    for gene in GeneSelect:
        txt=open(os.path.join(Output,str(MyIter)+'.db'),'w')
        txt.write('<E><ID>%s</ID><TP>%s</TP><LC>%s</LC><SEQ>%s</SEQ><DNA>%s</DNA></E>\n'%(gene.split('|')[1],'Protein-Coding',gene.split('|')[1],GeneSelect[gene][0],GeneSelect[gene][1]))
        txt.close()
        MyIter+=1
    
HGT (options.DBdir,options.HGTRate,options.Output)