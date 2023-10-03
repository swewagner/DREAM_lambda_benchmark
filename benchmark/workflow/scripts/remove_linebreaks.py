import re

with open(snakemake.output[0],"w") as out_f:
    with open(snakemake.input[0]) as in_f:
        out_f.write(re.sub(r'(?<!\d)\n(?!>)', '', in_f.read()))
