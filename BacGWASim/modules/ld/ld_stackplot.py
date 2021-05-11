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


def LDstackplot(LD, out, figsize):
    boundaries = [0.2, 0.4, 0.6, 0.8, 1]

    # changing matplotlib the default style
    matplotlib.style.use("ggplot")
    # Getting the table of distances and r2 squares
    LDTable = pd.read_csv(LD, sep="\t")

    # getting the column of distance and r2 and converting them to series
    Dist = LDTable[["Dist"]].values.tolist()
    r2 = LDTable[["r^2"]].values.tolist()

    # averaging all the r2 square for a specific distance and make a dictionary out of of it
    DistDict = {}
    for item in zip(Dist, r2):
        if item[0][0] not in DistDict:
            DistDict[item[0][0]] = [item[1]]
        else:
            DistDict[item[0][0]].append(item[1])
    for item in DistDict:
        DistDict[item] = np.mean(DistDict[item])

    # summarizing the r2 linkage score for 1000bp window across the 500kb
    keylist = sorted(DistDict.keys())

    bins = [50, 200, 500, 1000, 2000, 5000, 10000,
            20000, 50000, 100000, 200000, 500000, 1000000]
    inds = np.digitize(keylist, bins)
    max_bin = np.max(inds)

    DistDictSum = {}
    for Range in inds:
        DistDictSum[Range] = []
    for Value in zip(keylist, inds):
        # value[1] is equal to range value, here we accumulated all the r2 values relating to a bin
        DistDictSum[Value[1]].append(DistDict[Value[0]])

    data = [list() for n in range(len(boundaries) - 1)]
    for item in DistDictSum:
        size_bin = len(DistDictSum[item])
        for i in range(len(boundaries) - 1):
            filt = filter(lambda x: in_between(boundaries[i], x, boundaries[i+1]), DistDictSum[item])
            data[i].append(len(list(filt)) / size_bin)
            
    # Data
    r = [item for item in range(max_bin + 1)]  # number of bins

    # plot
    plt.figure(figsize=(float(figsize.split(",")[0])/2.54, float(figsize.split(",")[1])/2.54))
    barWidth = 1

    colors = ["#b5ffb9", "#f9bc86", "#a3acff", "#FF0000"]
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
ld = snakemake.input.haploview_ld
out = snakemake.output.ld_stackplot
LDstackplot(ld, out, "16,16")
