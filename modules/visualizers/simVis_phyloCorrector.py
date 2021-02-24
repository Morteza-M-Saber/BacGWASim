# Snakemake variables
input_file = snakemake.input.phylo
output_file = snakemake.output.phylo

def phylo_tip_corrector(phylo, old_id, new_id):
    pos1 = phylo.find(','+old_id+':')
    pos2 = phylo.find('('+old_id+':')
    if pos1 == -1 and pos2 == -1:
        print('Sample %s does not exist in tree!' % old_id)
        return (phylo)
    else:
        if pos1 != -1:
            return phylo[:pos1]+','+new_id+':'+phylo[pos1+len(old_id)+2:]
        else:
            return phylo[:pos2]+'('+new_id+':'+phylo[pos2+len(old_id)+2:]


txt = open(input_file, 'r')
phylo = txt.readline().strip()
old, new = '0', 'zero'
phylo_corrected = phylo_tip_corrector(phylo, old, new)
with open(output_file, 'w') as file:
    file.write(phylo_corrected)
