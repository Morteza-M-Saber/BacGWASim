
#       07/17/2018
#       Purpose
#Visualizing the results of simulation.
#Inputs: Simulated Tree, Simulated Phenotypes, Causal Variant list
#Outputs: Graph of simulated tree along with color-coded phenotype status and heatmap of causal variant presence/absence

from ete3 import Tree, TextFace, NodeStyle, TreeStyle, faces, AttrFace, ClusterTree, ProfileFace

#       Method
#1) Input the simulated phenotypes and generate a dictinary out of it.
PhenoDict={}
with open('/home/masih/Projects/BacterialSimulator/gatc.phen','r') as file:
    line=file.readline()
    while line:
        line=line.split()
        PhenoDict[line[0]]=line[2]
        line=file.readline()

print('step1 complete')
#2) Input the causal variant list and generate  an array out of it
import pandas as pd

#Getting the list of all samples
import vcf
VCFile='/home/masih/Projects/BacterialSimulator/BCFToolsResMAF0.01'
vcf_reader = vcf.Reader(open(VCFile, 'r'))
record = next(vcf_reader)
Samples=vcf_reader.samples



print('step2 complete')
#Getting the list of causal variants
MyCausalList=pd.read_table("/home/masih/Projects/BacterialSimulator/gatc.par") 
MyCausalList.set_index('QTL',inplace=True)
CausalList=MyCausalList.index.tolist()




#Genearte an empty array (list for now)
ProfileList=[]
#Add CausalMarker names as column headers
HeaderList=['#Names']
for CausalMarkers in CausalList:
    HeaderList.append(CausalMarkers)
ProfileList.append(HeaderList)
print(CausalList)
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
                
    
            
         

print('step5 complete')
# generate an Textarray out of the CausalLoci PresenceAbsence
Matrix=str()
for Profiles in ProfileList:
    for words in Profiles:
        Matrix+=str(words)+'\t'
    Matrix=Matrix[:-1]+'\n'

print(Matrix)
#print(Matrix)
#Generate ClusterTree using ETE toolkit
t=ClusterTree( '/home/masih/Projects/BacterialSimulator/RealTree.nwk' , text_array=Matrix)   

#Define a Tree visualizaiton layout        
def ColorCodedNode (node):
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
            
                 
ts = TreeStyle()
ts.layout_fn= ColorCodedNode
ts.mode='c' #Circular tree
ts.show_scale = False
ts.show_leaf_name = False
ts.draw_guiding_lines=True
t.render("mytree50.png", w=3840, units="px",tree_style=ts)        
    
    


