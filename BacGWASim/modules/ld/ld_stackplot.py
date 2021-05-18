"""
ای که خواهان تولدی دیگری, نخست مرگ را پذیرا باش
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.clf()
plt.rc("font", family="serif")
matplotlib.use("Agg")

def format_kb(num):
    if num >= 1000 or num == 0:
        return str(int(num/1000))
    else:
        return str(num/1000)


def in_between(a, b, c):
    """ Returns True if b is a < b <= c """
    return a < b and b <= c


def LDstackplot(dist, r2, out, figsize, colors):
    boundaries = [0.2, 0.4, 0.6, 0.8, 1]

    # changing matplotlib the default style
    matplotlib.style.use("ggplot")

    # averaging all the r2 square for a specific distance and make a dictionary out of of it
    dist_dict = {}
    for _dist, _r2 in zip(dist, r2):
        if _dist not in dist_dict:
            dist_dict[_dist] = [_r2]
        else:
            dist_dict[_dist].append(_r2)
    for item in dist_dict:
        dist_dict[item] = np.mean(dist_dict[item])

    # summarizing the r2 linkage score for 1000bp window across the 500kb
    keylist = sorted(dist_dict.keys())

    bins = [50, 200, 500, 1000, 2000, 5000, 10000,
            20000, 50000, 100000, 200000, 500000, 1000000]
    inds = np.digitize(keylist, bins)
    max_bin = np.max(inds)

    dist_dict_sum = {}
    for Range in inds:
        dist_dict_sum[Range] = []
    for Value in zip(keylist, inds):
        # value[1] is equal to range value, here we accumulated all the r2 values relating to a bin
        dist_dict_sum[Value[1]].append(dist_dict[Value[0]])

    data = [list() for n in range(len(boundaries) - 1)]
    for item in dist_dict_sum:
        size_bin = len(dist_dict_sum[item])
        for i in range(len(boundaries) - 1):
            filt = filter(
                lambda x: in_between(boundaries[i], x, boundaries[i+1]), 
                dist_dict_sum[item]
            )
            data[i].append(len(list(filt)) / size_bin)
            
    # Data
    r = [item for item in range(max_bin + 1)]  # number of bins

    # plot
    plt.figure(figsize=(float(figsize.split(",")[0])/2.54, float(figsize.split(",")[1])/2.54))

    # Flipping some data to reorder the bars
    data = data[::-1]
    boundaries = boundaries[::-1]
    labels = list()
    for i in range(len(boundaries) - 1):
        bottom = [sum(x) for x in zip(*data[:i], [0]*len(data[0]))]
        label = str(boundaries[i+1]) + " < $R^2$ < " + str(boundaries[i])
        plt.bar(r, data[i], bottom=bottom, color=colors[i], label=label, edgecolor="white", width=1)
        labels.append(label)

    # Custom x axis
    names = list()
    named_bins = [0] + bins
    for i in range(max_bin):
        name = "[" + format_kb(named_bins[i]) + ", " + format_kb(named_bins[i+1]) + "["
        names.append(name)
    names.append("[" + format_kb(named_bins[max_bin]) + ", ")
    plt.xticks(r, names, rotation="vertical")
    plt.xlabel("Distance (kb)")
    plt.ylabel("$R^2$ score")

    # Add a legend
    plt.legend(reversed(plt.legend().legendHandles), reversed(labels))

    # Show graphic
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out, dpi=300)


# Running LDstackplot
ld = snakemake.input.ld
vcf = snakemake.input.vcf_subset
out = snakemake.output.ld_stackplot
colors = ["#b5ffb9", "#f9bc86", "#a3acff", "#FF0000"]

# Loading data
data = np.genfromtxt(ld)
N = data.shape[0]

# Loading VCF data
vcf = pd.read_csv(vcf, sep="\t", skiprows=13, usecols=["POS"]).to_dict()["POS"]

# Parsing data to get [(DIST, R2), ...]
dist = [
    np.abs(vcf[li] - vcf[col])
    for line in range(N) 
    for li, col in zip(range(line, N), range(0, N-line))
]

r2 = [
    data[li, col]
    for line in range(N) 
    for li, col in zip(range(line, N), range(0, N-line))
]

LDstackplot(dist, r2, out, "16,16", colors)
