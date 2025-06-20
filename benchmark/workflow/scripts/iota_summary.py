# extract for each bin, which queries match
from Bio import SeqIO

dummy = snakemake.input[0]
bin_dic = {}
with open(dummy, "r") as dummy_f:
    for line in dummy_f:
        bin_path = line.rstrip('\n')
        bin_id = bin_path.split('.')[-2].split('_')[-1]
        bin_dic[bin_id] = []
        for record in SeqIO.parse(bin_path, "fasta"):
            query_id = record.id
            bin_dic[bin_id].append(query_id)

with open(snakemake.output[0], "w") as results_file:
    for bin_id in bin_dic.keys():
        results_file.write(bin_id + ', ' + ', '.join(bin_dic[bin_id]) + '\n')
