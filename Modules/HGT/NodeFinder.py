'''Extract a set of random internal nodes and their leaves from a newick tree file'''


def get_options():
    import argparse

    description = 'Extract a set of random internal nodes and their leaves'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('phylogeny',
                        help='Tree file')

    parser.add_argument('--NodeLen',
                        default='5',
                        help="Node leaves size [Default: 5]")
    parser.add_argument('--HGTRate',
                        default='10',
                        help="Frequency of Horizontal Gene Transfer events [Default: 10]")
    parser.add_argument('--Output',
                        help="Output file including selected nodes and the corresponding leaves in binary format")


    return (parser.parse_args())


options = get_options()

def NodeFinder(phylogeny,NodeLen,HGTRate,Output):
    from ete3 import ClusterTree
    t=ClusterTree(phylogeny)
    
    TreeDict={}
    
    for node in t.traverse():
        TreeDict[node]=[]
        for leaf in node:
            TreeDict[node].append(leaf.name)
    
    #Selecting nodes with leaf size of selected size    
    Iter=0
    NodeList=[]
    TreeSelect={}
    for node in TreeDict:
        if len(TreeDict[node])>= int(NodeLen):
            NodeList.append(Iter)
            TreeSelect[Iter]=TreeDict[node]
            Iter+=1
            
    #Randomly selecting nodes without repalcement,
    #if number of eligible nodes are less than HGTRate, Randomly choose with replacement
    import numpy as np
    if int(options.HGTRate)<=len(NodeList):
        NodeSelect=np.random.choice(NodeList,int(HGTRate),replace=False)
    elif int(options.HGTRate)>len(NodeList):
        NodeSelect=np.random.choice(NodeList,int(HGTRate),replace=True)
    
    NodeFinal={}
    for node in NodeSelect:
        NodeFinal[str(node)]=TreeSelect[node]
    #Return Selected Node List
    print(NodeFinal)
    import json
    with open(str(Output), 'w') as f:
            json.dump(NodeFinal, f)


NodeFinder(options.phylogeny,options.NodeLen,options.HGTRate,options.Output)  
        
        
        
        