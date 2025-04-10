
with open(snakemake.output[0], 'w') as f:
    for i in range(len(snakemake.input)):
        f.write(snakemake.input[i])
        if i < len(snakemake.input)-1:
            f.write('\n')