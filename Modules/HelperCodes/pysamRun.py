#Checks whether the produced alignment by the bwa aligner is of proper quality
#The output of valid sam file and corrupted sam file is clearly visible by running this code
#in sudden shut down of program, bwa produces the corrupted result and since snakemake has not time
#to tag them as corrupted due to sudden termination, in the re-run the corrupted alignment causes error in variant-calling


import pysam
iter=1
samfile = pysam.AlignmentFile("/project/6006375/masih/LDLgtFinal4/simulations/4/Assembly/SE015.sam", "r")
reads=samfile.fetch()
for item in reads:
     while iter<5:
             print(item)
             iter+=1
