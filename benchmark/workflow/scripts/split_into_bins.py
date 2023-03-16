from Bio import SeqIO

ref_file = snakemake.input[0]
bin_num = len(snakemake.output)
ref_num = int(snakemake.config["num_of_refseqs"])

def endOfBin():
    if bin_index < ref_num % bin_num:
        return index % ((ref_num // bin_num) + 1) == 0
    else:
        return (index - offset) % (ref_num // bin_num) == 0

current_records = []
index = 0
bin_index = 0
offset = (ref_num % bin_num) * ((ref_num // bin_num) + 1)
print("ref_num//bin_num: ", ref_num // bin_num)

with open(ref_file) as ref_f:
    for record in SeqIO.parse(ref_f, "fasta"):
        current_records.append(record)
        index += 1
        if endOfBin():
            SeqIO.write(current_records, snakemake.output[bin_index], "fasta")
            current_records = []
            bin_index += 1
