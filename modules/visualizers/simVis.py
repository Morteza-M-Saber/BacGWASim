
"""
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش
#Visualizing the results of simulation.
#Inputs: Simulated Tree, Simulated Phenotypes, Causal Variant list
#Outputs: Graph of simulated tree along with color-coded phenotype status and heatmap of causal variant presence/absence
"""

def get_options():
    import argparse

    description = 'Visualize simulation'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--vcfIn',
                        help='path to vcf file ')
    parser.add_argument('--phen',
                        help='path to phenotype in gcta format ')
    parser.add_argument('--phenType',
                        help='Type of phenotype cc/quant (cc:case-control&quant:quantitative) ')
    parser.add_argument('--par',
                        help='path to causalmarkers in gcta format ')
    parser.add_argument('--phylo',
                        help='path to phylogenetic tree in newick format ')
    parser.add_argument('--out',
                        help='path to output')

    return parser.parse_args()

options = get_options()

import os
os.environ['QT_QPA_PLATFORM']='offscreen' #ete3 relies on qt that requires Xserver connection!!
from ete3 import Tree, TextFace, NodeStyle, TreeStyle, faces, AttrFace, ClusterTree, ProfileFace, CircleFace
import pandas as pd
import numpy, vcf

def simViz(vcfIn,phen,phenType,par,phylo,out):
    #1) Input the simulated phenotypes and generate a dictinary out of it.
    PhenoDict={}
    if phenType =='cc':
      with open(phen,'r') as file:
          line=file.readline()
          while line:
              line=line.split()
              PhenoDict[line[0]]=line[2]
              line=file.readline()
    elif phenType =='quant':
      data=pd.read_csv(phen,sep=' ',header=None,index_col=0)
      #removing samples with unknwon phenotype
      for sample_phen in zip (data[1],data[2]):
        if float(sample_phen[1]) == -9.0 :
          PhenoDict[sample_phen[0]]='ukwn'
          data.drop(labels=sample_phen[0],axis=0,inplace=True)
      #bin the phen values
      val=data[2]
      bins = numpy.linspace(min(val),max(val),7)
      digitized_val = numpy.digitize(val, bins)
      data['bin']=digitized_val
      data['bin'].replace(7,6,inplace=True) #include last right-side value into final bin 
      for sample_phen in zip (data[1],data['bin']):
        PhenoDict[sample_phen[0]]=str(sample_phen[1])
    print('Visuazlier Step1 completed!')
    #2) Input the causal variant list and generate  an array out of it
    #Getting the list of all samples
    VCFile=vcfIn
    vcf_reader = vcf.Reader(open(VCFile, 'r'))
    record = next(vcf_reader)
    Samples=vcf_reader.samples
    print('Visuazlier Step2 completed!')
    #Getting the list of causal variants
    MyCausalList=pd.read_table(par) 
    MyCausalList.set_index('QTL',inplace=True)
    CausalList=MyCausalList.index.tolist()
    print('Visuazlier Step3 completed!')
    #Genearte an empty array (list for now)
    ProfileList=[]
    #Add CausalMarker names as column headers
    HeaderList=['#Names']
    for CausalMarkers in CausalList:
        HeaderList.append(CausalMarkers)
    ProfileList.append(HeaderList)
    print('Visuazlier Step4 completed!')
    #Iterate over samples present in the phylogenetic tree
    ProfileListDic={}
    for samp in Samples:
        ProfileListDic[samp]=[samp.split('_')[0]]+(['Null']*len(CausalList)) #Modify this with acual results 
    vcf_reader = vcf.Reader(open(VCFile, 'r'))
    FoundHit=0
    for records in vcf_reader:
        if records.ID in CausalList:
            FoundHit+=1
            for SamplesNow in Samples:
                    GT=records.genotype(SamplesNow)['GT']
                    if GT == '0':
                        GTUpdate=-1
                    elif GT =='1':
                        GTUpdate=1
                    else:
                        raise Exception ('The identified status of the marker %s is not appropriate for bacterial polymorphism'%GT)
                    ProfileListDic[SamplesNow][CausalList.index(records.ID)+1]=GTUpdate
        if FoundHit==len(CausalList):
            break
    for item in Samples:
        ProfileList.append(ProfileListDic[item])
    print('Visuazlier Step5 completed!')
    # generate an Textarray out of the CausalLoci PresenceAbsence
    Matrix=str()
    for Profiles in ProfileList:
        for words in Profiles:
            Matrix+=str(words)+'\t'
        Matrix=Matrix[:-1]+'\n'
    #Generate ClusterTree using ETE toolkit
    t=ClusterTree( phylo , text_array=Matrix)   
    #Define a Tree visualizaiton layout        
    def ColorCodedNodecc (node):
        if node.is_leaf():
            ColorCode=PhenoDict[node.name]
            if ColorCode == '1':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Blue")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace('name',fsize=20,fgcolor='blue'), node, column=0,aligned=True)
                #faces.add_face_to_node(TextFace(text='marker1',fsize=10,fgcolor='black'), node, column=1,position='aligned')
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == '2':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Red")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='red'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == '-9':
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='black'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
    def ColorCodedNodequant (node):
        if node.is_leaf():
            ColorCode=PhenoDict[node.name]
            if ColorCode == '1':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Blue")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace('name',fsize=20,fgcolor='#0302FC'), node, column=0,aligned=True)
                #faces.add_face_to_node(TextFace(text='marker1',fsize=10,fgcolor='black'), node, column=1,position='aligned')
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == '2':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Red")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='#2A00D5'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == '3':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Red")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='#63009E'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == '4':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Red")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='#A1015D'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == '5':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Red")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='#D80027'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == '6':
                #Name=faces.AttrFace('name',fsize='20',fgcolor="Red")
                #NameFace=TextFace(Name)
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='#FE0002'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')
            elif ColorCode == 'ukwn':
                faces.add_face_to_node(AttrFace("name",fsize=20,fgcolor='black'), node, column=0,aligned=True)
                faces.add_face_to_node(ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2),node,column=1,position='aligned')     
                    
    ts = TreeStyle()
    #coloring the labels based on phenotype
    if phenType =='cc':
      ts.layout_fn= ColorCodedNodecc
    elif phenType == 'quant':
      ts.layout_fn= ColorCodedNodequant

    ts.mode='c' #Circular tree
    ts.show_scale = True
    ts.show_leaf_name = False
    ts.draw_guiding_lines=True
    #Adding color legend
    if phenType == 'cc':
      ts.legend.add_face(CircleFace(40, "red"),column=0)
      ts.legend.add_face(TextFace("Case",fsize=40), column=1)
      ts.legend.add_face(CircleFace(40, "blue"),column=0)
      ts.legend.add_face(TextFace("Control",fsize=40), column=1)
    elif phenType =='quant':
      ts.legend.add_face(CircleFace(30, "#0302FC"),column=0)
      ts.legend.add_face(TextFace("%s_%s"%(round(bins[0],1),round(bins[1],1)),fsize=30), column=1)
      ts.legend.add_face(CircleFace(30, "#2A00D5"),column=0)
      ts.legend.add_face(TextFace("%s_%s"%(round(bins[1],1),round(bins[2],1)),fsize=30), column=1)
      ts.legend.add_face(CircleFace(30, "#63009E"),column=0)
      ts.legend.add_face(TextFace("%s_%s"%(round(bins[2],1),round(bins[3],1)),fsize=30), column=1)
      ts.legend.add_face(CircleFace(30, "#A1015D"),column=0)
      ts.legend.add_face(TextFace("%s_%s"%(round(bins[3],1),round(bins[4],1)),fsize=30), column=1)
      ts.legend.add_face(CircleFace(30, "#D80027"),column=0)
      ts.legend.add_face(TextFace("%s_%s"%(round(bins[4],1),round(bins[5],1)),fsize=30), column=1)
      ts.legend.add_face(CircleFace(30, "#FE0002"),column=0)
      ts.legend.add_face(TextFace("%s_%s"%(round(bins[5],1),round(bins[6],1)),fsize=30), column=1)

    t.render(out, dpi=300, tree_style=ts)
    print('Visuazlier Successfully completed!')        


simViz(options.vcfIn,options.phen,options.phenType,options.par,options.phylo,options.out)