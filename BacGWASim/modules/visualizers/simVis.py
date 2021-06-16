from numpy.lib.function_base import select
import dendropy
import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.cluster.hierarchy as sch
import vcf

mpl.rc('font', weight='bold')

# Snakemake variables
input_vcf = snakemake.input.vcf
input_par = snakemake.input.par
input_phen = snakemake.input.phen
input_phylo = snakemake.input.phylo
output_simVis = snakemake.output.simVis
phen_type = snakemake.config["phen_type"]
dpi = snakemake.config["simvis_dpi"]
max_genomes = 500

# Visual parameters
if max_genomes <= 200:
    linewidth = 2
else:
    linewidth = 1

""" 
    Functions fix_verts, smoothsegment and part of the polar dendrogram plotting 
    logic is from: chthonicdaemon
    https://stackoverflow.com/questions/51936574/how-to-plot-scipy-hierarchy-dendrogram-using-polar-coordinates
"""

def font_size(num):
    """ Discrete selection of font size for plot scalling. """
    if num < 40:
        return 16
    elif num < 100:
        return 12
    elif num < 150:
        return 10
    elif num < 350:
        return 8
    else:
        return 5

def fix_verts(ax, orient=1):
    for coll in ax.collections:
        for pth in coll.get_paths():
            vert = pth.vertices
            vert[1:3,orient] = scipy.average(vert[1:3,orient]) 

def smoothsegment(seg, Nsmooth=100):
    return np.concatenate([[seg[0]], np.linspace(seg[1], seg[2], Nsmooth), [seg[3]]])


def polar2z(theta, r):
    """ Transform polar coordinates to planar. """
    z = r * np.exp(1j*theta)
    return (np.real(z), np.imag(z))

def plot_mutations(ax, mutations, angle, x_rect, y_rect, height, width):
    angle_rad = angle/180*np.pi
    mut_height = height / len(mutations)
    for i, mutation in enumerate(mutations):
        if mutation == 1:
            x = x_rect - np.sin(angle_rad) * mut_height*i
            y = y_rect + np.cos(angle_rad) * mut_height*i
            rect = patches.Rectangle(
                (x, y), width=width, height=mut_height, angle=angle,
                transform=ax.transData._b, facecolor='r', edgecolor='k'
            )
            ax.add_artist(rect)

def plot_dendrogram(icoord, dcoord, labels, ax):
    """ Main plotting function for the dendrogram. """

    # Coord for polar plotting
    dcoord = -np.log(dcoord+1)
    gap = 1 / len(icoord)
    imax = icoord.max()
    imin = icoord.min()
    icoord = ((icoord - imin)/(imax - imin)*(1-gap) + gap/2)*2*np.pi

    # Plotting the dendrogram
    for xs, ys in zip(icoord, dcoord):
        xs = smoothsegment(xs)
        ys = smoothsegment(ys)
        ax.plot(xs, ys, color="black", linewidth=3)

    # Axe manipulations
    ax.spines['polar'].set_visible(False)
    ax.grid(False)
    ax.set_rlabel_position(0)

    # Axe labels
    ax.set_yticklabels([])
    Nxticks = len(icoord) + 1
    xticks = np.linspace(gap/2, 1-gap/2, Nxticks)
    ax.set_xticks(xticks*np.pi*2)
    ax.set_xticklabels(labels)
    ax.tick_params(axis='x', which='major', labelsize=font_size(len(labels)))

    # Setting limits to include outer polar plots
    ylim = ax.get_ylim()
    y_min, y_max = ylim[0], -ylim[0]/2
    y_curr = - y_min
    ax.set_ylim(y_min, y_max)

    # Drawing class circles
    circum = -y_min * 2 * np.pi
    width = min(circum / (len(labels)*1.25), circum/360*5)
    radius = min(width/2, circum/360*5/2)
    y_curr += radius * 1.5
    for tick, label in zip(xticks, labels):
        circle = patches.Circle(
            polar2z(tick*np.pi*2, y_curr), radius=radius, 
            transform=ax.transData._b,
            linewidth=2, edgecolor='k', 
            facecolor=cmap(norm(float(PhenoDict[label]))),
        )
        ax.add_patch(circle)

    # Drawing genomes
    y_curr += radius * 2
    for tick, label in zip(xticks, labels):
        angle = tick*360-90
        angle_rad = angle/180*np.pi
        x_rect, y_rect = polar2z(tick*np.pi*2, y_curr)
        height = ((y_max - y_min) - y_curr) * (1 - font_size(len(labels))/100/2)
        x_rect -= np.cos(angle_rad) * width/2
        y_rect -= np.sin(angle_rad) * width/2

        rect = patches.Rectangle(
            (x_rect, y_rect), width=width, height=height, angle=angle,
            edgecolor='k', facecolor=cmap(norm(float(PhenoDict[label]))),
            transform=ax.transData._b, linewidth=linewidth,
        )
        ax.add_patch(rect)

        plot_mutations(ax, pheno_dict[label], angle, x_rect, y_rect, height, width)


def plot_legend(cmap, phen_type):
    """ Adding a colormap as a legend. """
    a = np.array([[0, 1]])
    img = plt.imshow(a, cmap=cmap)
    img.set_visible(False)
    cax = plt.axes([0.05, 0.075, 0.2, 0.030])
    cax.set_title('Phenotype', fontsize=16, fontweight='bold')
    cbar = plt.colorbar(orientation="horizontal", cax=cax)
    if phen_type == "quant":
        cbar.set_ticks([0, 0.5, 1])
        cbar.set_ticklabels([0, 0.5, 1])
    elif phen_type == "cc":
        cbar.set_ticks([0.25, 0.75])
        cbar.set_ticklabels([0, 1])
    cax.tick_params(axis='x', labelsize=14)
    for tick in cax.get_xticklabels():
        tick.set_fontweight("normal")

    cbar.outline.set_linewidth(2)

# 1) Input the simulated phenotypes and generate a dictinary out of it.
PhenoDict = {}
if phen_type == 'cc':
    with open(input_phen, 'r') as file:
        line = file.readline()
        while line:
            line = line.split()
            PhenoDict[line[0]] = line[2]
            line = file.readline()
    color_min, color_max = 1, 2
    
elif phen_type == 'quant':
    data = pd.read_csv(input_phen, sep=' ', header=None, index_col=0)
    # Removing samples with unknown phenotype
    data = data[data[2] != -9.0]
    # Normalize the data
    data[2]=(data[2]-data[2].mean())/data[2].std()
    color_max = data[2].max()
    color_min = data[2].min()
    # Data to dict
    PhenoDict = data[2].to_dict()

# 2) Getting the list of all samples
samples = list(PhenoDict.keys())
if len(samples) > max_genomes:
    samples = np.random.choice(samples, max_genomes, replace=False)

# 3) Getting the list of causal variants
CausalList = pd.read_csv(input_par, sep="\t")["QTL"].to_list()

# 4) Iterate over samples present in the phylogenetic tree
pheno_dict = {sample: ["Null"]*len(CausalList) for sample in samples}
vcf_reader = vcf.Reader(open(input_vcf, 'r'))
for records in vcf_reader:
    for SamplesNow in samples:
        GT = records.genotype(SamplesNow)['GT']
        if GT == '0':
            GTUpdate = -1
        elif GT == '1':
            GTUpdate = 1
        else:
            raise Exception('The identified status of the marker %s is not appropriate for bacterial polymorphism' % GT)
        pheno_dict[SamplesNow][CausalList.index(records.ID)] = GTUpdate


# Creating tree and pdm
tree = dendropy.Tree.get(path=input_phylo, schema='newick')
pdm = tree.phylogenetic_distance_matrix()

# Selecting genomes (limit of max_genomes is plotted)
taxons = list(tree.taxon_namespace)
if len(taxons) > max_genomes:
    taxons = [taxon for taxon in taxons if taxon._label in samples]

# Labels
labels = [taxon.label for taxon in taxons]

# Creating color map
norm = mpl.colors.Normalize(vmin=color_min,vmax=color_max)
cmap = plt.cm.cool
if phen_type == "cc":
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'binary_cmap', [cmap(i) for i in range(cmap.N)], 2
    )

# Distance matrix
dist = np.zeros(shape=(len(labels)**2 - len(labels))//2)
indx = 0
for i, t1 in enumerate(taxons):
    for t2 in taxons[i+1:]:
        dist[indx] = pdm(t1, t2)
        indx += 1

# Compute
Y = sch.linkage(dist, method='single')
Z2 = sch.dendrogram(Y, labels=labels)

# Getting coordinates
icoord = np.array(Z2['icoord'])
dcoord = np.array(Z2['dcoord'])

# Creating figure handle
fig = plt.figure(figsize=(16,16), dpi=int(dpi))
ax = fig.add_subplot(111, polar=True)

# Plotting main dendrogram
plot_dendrogram(icoord, dcoord, Z2['ivl'], ax)

# Adding legend
plot_legend(cmap, phen_type)

# Saving figure
plt.savefig(output_simVis)