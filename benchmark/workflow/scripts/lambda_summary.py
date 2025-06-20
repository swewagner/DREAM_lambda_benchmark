bin_dic = {}

for f in snakemake.input:
    bin_id = int(f.split('.')[-2].split('_')[-1])
    with open(f) as lambda_out:
        for line in lambda_out:
            bin_dic.setdefault(bin_id, set()).add(line.split('\t')[0])

with open(snakemake.output[0], "w") as out_file:
    for bin_id in bin_dic:
        out_file.write(str(bin_id) + ", " + ",  ".join(list(bin_dic[bin_id])) + '\n')