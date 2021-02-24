
"""
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش
#Visualizing the results of simulation.
#Inputs: Simulated Tree, Simulated Phenotypes, Causal Variant list
#Outputs: Graph of simulated tree along with color-coded phenotype status and heatmap of causal variant presence/absence
"""
import vcf
import numpy
import pandas as pd
from ete3 import Tree, TextFace, NodeStyle, TreeStyle, faces, AttrFace, ClusterTree, ProfileFace, CircleFace
import os
# ete3 relies on qt that requires Xserver connection!!
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

# Snakemake variables
input_vcf = snakemake.input.vcf
input_par = snakemake.input.par
input_phen = snakemake.input.phen
input_phylo = snakemake.input.phylo
output_simVis = snakemake.output.simVis
phenType = snakemake.config["phenType"]

# Plot variables
quant_colors = [
    "#0302FC", "#2A00D5", "#63009E",
    "#A1015D", "#D80027", "#FE0002"
]

node_colors_cc = {
    "1": "blue",
    "2": "red",
    "-9": "black",
}

node_colors_quant = {
    "1": quant_colors[0],
    "2": quant_colors[1],
    "3": quant_colors[2],
    "4": quant_colors[3],
    "5": quant_colors[4],
    "6": quant_colors[5],
    "ukwn": "black",
}


def ColorCodedNodecc(node):
    if node.is_leaf():
        ColorCode = PhenoDict[node.name]
        faces.add_face_to_node(
            AttrFace("name", fsize=20, fgcolor=node_colors_cc[ColorCode]),
            node, column=0, aligned=True
        )
        faces.add_face_to_node(
            ProfileFace(1, -1, 0, width=200, height=40,style='heatmap', colorscheme=2),
            node, column=1, position='aligned'
        )


def ColorCodedNodequant(node):
    if node.is_leaf():
        ColorCode = PhenoDict[node.name]
        faces.add_face_to_node(
            AttrFace("name", fsize=20, fgcolor=node_colors_quant[ColorCode]), 
            node, column=0, aligned=True
        )
        faces.add_face_to_node(
            ProfileFace(1, -1, 0, width=200, height=40, style='heatmap', colorscheme=2), 
            node, column=1, position='aligned'
        )


# 1) Input the simulated phenotypes and generate a dictinary out of it.
PhenoDict = {}
if phenType == 'cc':
    with open(input_phen, 'r') as file:
        line = file.readline()
        while line:
            line = line.split()
            PhenoDict[line[0]] = line[2]
            line = file.readline()
elif phenType == 'quant':
    data = pd.read_csv(input_phen, sep=' ', header=None, index_col=0)
    # removing samples with unknwon phenotype
    for sample_phen in zip(data[1], data[2]):
        if float(sample_phen[1]) == -9.0:
            PhenoDict[sample_phen[0]] = 'ukwn'
            data.drop(labels=sample_phen[0], axis=0, inplace=True)
    # bin the phen values
    val = data[2]
    bins = numpy.linspace(min(val), max(val), 7)
    digitized_val = numpy.digitize(val, bins)
    data['bin'] = digitized_val
    # include last right-side value into final bin
    data['bin'].replace(7, 6, inplace=True)
    for sample_phen in zip(data[1], data['bin']):
        PhenoDict[sample_phen[0]] = str(sample_phen[1])
print('Visualizer Step1 completed!')

# 2) Input the causal variant list and generate  an array out of it
# Getting the list of all samples
vcf_reader = vcf.Reader(open(input_vcf, 'r'))
record = next(vcf_reader)
Samples = vcf_reader.samples
print('Visualizer Step2 completed!')

# 3) Getting the list of causal variants
MyCausalList = pd.read_table(input_par)
MyCausalList.set_index('QTL', inplace=True)
CausalList = MyCausalList.index.tolist()
print('Visualizer Step3 completed!')

# 4) Generate an empty array (list for now)
ProfileList = []
# Add CausalMarker names as column headers
HeaderList = ['#Names']
for CausalMarkers in CausalList:
    HeaderList.append(CausalMarkers)
ProfileList.append(HeaderList)
print('Visualizer Step4 completed!')

# 5) Iterate over samples present in the phylogenetic tree
ProfileListDic = {}
for samp in Samples:
    # Modify this with acual results
    ProfileListDic[samp] = [samp.split('_')[0]]+(['Null']*len(CausalList))
vcf_reader = vcf.Reader(open(input_vcf, 'r'))
FoundHit = 0
for records in vcf_reader:
    if records.ID in CausalList:
        FoundHit += 1
        for SamplesNow in Samples:
            GT = records.genotype(SamplesNow)['GT']
            if GT == '0':
                GTUpdate = -1
            elif GT == '1':
                GTUpdate = 1
            else:
                raise Exception(
                    'The identified status of the marker %s is not appropriate for bacterial polymorphism' % GT)
            ProfileListDic[SamplesNow][CausalList.index(
                records.ID)+1] = GTUpdate
    if FoundHit == len(CausalList):
        break
for item in Samples:
    ProfileList.append(ProfileListDic[item])
print('Visualizer Step5 completed!')
# Generate an Textarray out of the CausalLoci PresenceAbsence
Matrix = str()
for Profiles in ProfileList:
    for words in Profiles:
        Matrix += str(words)+'\t'
    Matrix = Matrix[:-1]+'\n'
# Generate ClusterTree using ETE toolkit
t = ClusterTree(input_phylo, text_array=Matrix)
# Define a Tree visualizaiton layout
ts = TreeStyle()
# Coloring the labels based on phenotype
if phenType == 'cc':
    ts.layout_fn = ColorCodedNodecc
elif phenType == 'quant':
    ts.layout_fn = ColorCodedNodequant

ts.mode = 'c'  # Circular tree
ts.show_scale = True
ts.show_leaf_name = False
ts.draw_guiding_lines = True
# Adding color legend
if phenType == 'cc':
    ts.legend.add_face(CircleFace(40, "red"), column=0)
    ts.legend.add_face(TextFace("Case", fsize=40), column=1)
    ts.legend.add_face(CircleFace(40, "blue"), column=0)
    ts.legend.add_face(TextFace("Control", fsize=40), column=1)
elif phenType == 'quant':
    for i in range(6):
        ts.legend.add_face(CircleFace(30, quant_colors[i]), column=0)
        ts.legend.add_face(TextFace("%s_%s" % (
            round(bins[i], 1), round(bins[i+1], 1)), fsize=30), column=1)

t.render(output_simVis, dpi=300, tree_style=ts)
print('Visualizer Successfully completed!')