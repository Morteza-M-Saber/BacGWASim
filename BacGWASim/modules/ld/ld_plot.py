import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ld = snakemake.input.ld
vcf = snakemake.input.vcf_subset
plot = snakemake.output.plot
genome_length = snakemake.config["genome_length"]
ld_maf = snakemake.config["ld_maf"]

# Loading LD data
data = np.genfromtxt(ld)
N = data.shape[0]

# Loading VCF data
vcf = pd.read_csv(vcf, sep="\t", skiprows=13)
snps = vcf["POS"].values.tolist()

# Masking the upper triangle
mask = np.tri(N, k=-1).T
data_masked = np.ma.array(data, mask=mask)

# Plotting the lower triangle
im_width = 20
im_height = 12
circle_area_width = np.pi*(im_width*72/2)**2

fig = plt.figure(figsize=(im_width, im_height))
ax = fig.add_axes([0, 0, 1, 1])
ax.axis('equal')
ax.set(xlim=(0, 1), ylim=(0, 0.6))


pad_x = 0.025
width = (1 - 2*pad_x) / (N)
square_area = (im_width*72*(1-pad_x*2)/N / 2)**2 * np.pi

x_grid = pad_x
y_grid = 0.5
gap = 0.04
y_gen = y_grid + gap

norm = mpl.colors.Normalize(vmin=0, vmax=1)
cmap = plt.cm.Reds
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, 11)
bounds = np.linspace(0, 1, 11)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)


# Creating grid of squares
r2_mat_coord = [
    [li, col] 
    for line in range(N) 
    for li, col in zip(range(line, N), range(0, N-line))
]

r2_pos = [
    (x_grid + (i+j+1)*width/2, y_grid - (i-j+1)*width/2) 
    for (i, j) in r2_mat_coord
]

colors = [cmap(norm(data[i, j])) for (i, j) in r2_mat_coord]

collection = mpl.collections.RegularPolyCollection(
    numsides=4,
    rotation=0, sizes=(square_area,),
    facecolors=colors,
    linewidths=(1,),
    offsets=r2_pos,
    transOffset=ax.transData,
    )

ax.add_collection(collection)

# Creating genome
genome_width = 1-2*pad_x
rect = mpl.patches.Rectangle(
    (pad_x, y_gen), height=0.035, width=genome_width,
    facecolor='w', edgecolor='k', linewidth=2,
)
ax.add_artist(rect)

# Adding mutations to genome
x_pos = pad_x - width / 2
mutation_step = genome_width / genome_length
mutations = [
    [
        (pad_x + (snp+1.5)*mutation_step, y_gen), 
        (pad_x + (snp+1.5)*mutation_step, y_gen+0.035)
    ] 
    for snp in snps
]
mutation_collection = mpl.collections.LineCollection(
    mutations, linewidth=0.5, color="k"
)
ax.add_collection(mutation_collection)


# Link between SNP and genome
links = [
    [
        (x_grid + i*width + width/4, y_grid - width/8),
        (x_grid + i*width + width/4, y_grid + gap/8),
        (pad_x + (snp+1.5)*mutation_step, y_gen-gap/30)
    ] 
    for i, snp in enumerate(snps)
]
link_collection = mpl.collections.LineCollection(
    links, linewidth=0.25, color="k"
)
ax.add_collection(link_collection)


# Adding colorbar
a = np.array([[0, 1]])
img = plt.imshow(a, cmap=cmap)
img.set_visible(False)
cax = plt.axes([0.05, 0.075, 0.2, 0.030])
cax.set_title("Linkage Disequilibrium (R$^2$)", fontsize=16, fontweight='bold')
cbar = plt.colorbar(orientation="horizontal", cax=cax)
cbar.set_ticks([0, 0.5, 1])
cbar.set_ticklabels([0, 0.5, 1])
cax.tick_params(axis='x', labelsize=14)
cbar.outline.set_linewidth(2) 

# Saving to file
plt.savefig(plot)
