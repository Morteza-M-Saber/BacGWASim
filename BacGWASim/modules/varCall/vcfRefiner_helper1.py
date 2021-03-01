import numpy as np
from shutil import copyfile

vcf_input = snakemake.input.vcf
vcf_output = snakemake.output.vcf
varNumber = snakemake.config["num_var"]

if varNumber != -1:
    info_lines = []
    info_line_count = 0
    with open(vcf_input, 'r') as file:
        line = file.readline()
        while line[0] == '#':
            info_lines.append(line)
            info_line_count += 1
            line = file.readline()

    var_lines = open(vcf_input, 'r').readlines()[info_line_count:]

    if (len(var_lines) - info_line_count) <= varNumber:
        print('Number of simulated markers is less than requested!\nConsider increasing mutation rate!')
    else:
        rand_lines = np.random.choice(range(len(var_lines)), varNumber, replace=False)
        rand_indx = np.sort(rand_lines)
        with open(vcf_output, 'w') as txt:
            for info_ in info_lines:
                txt.write(info_)
            for var_choice in rand_indx:
                txt.write(var_lines[var_choice])

else:
    copyfile(vcf_input, vcf_output)