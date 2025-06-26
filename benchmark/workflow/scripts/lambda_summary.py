bin_dic = {}

with open(snakemake.input.lambda_out) as f:
    for line in f:
        bin_dic.setdefault(line.split('\t')[1], set()).add(line.split('\t')[0])

with open(snakemake.output[0], "w") as out_file:
    for sbj_id in bin_dic:
        out_file.write(str(sbj_id) + ", " + ",  ".join(list(bin_dic[sbj_id])) + '\n')