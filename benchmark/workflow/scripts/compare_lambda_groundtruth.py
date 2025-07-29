import csv
import math

gt_file = snakemake.input[0]
results_file = snakemake.input[1]
#gt_file = "../../results/er_0/prot/ground_truth.txt"
#results_file = "../../results/er_0/blastP/lambda/results.txt"

bin_num = int(snakemake.config["num_of_bins"])
ref_num = int(snakemake.config["num_of_refseqs"])
num_fuller_bins = ref_num % bin_num
num_refs_in_lesser_bin = (ref_num // bin_num)
num_refs_in_full_bin = (ref_num // bin_num)+1
num_refs_in_fuller_bins = num_fuller_bins * num_refs_in_full_bin

def ref_to_bin(ref_id):
    if (ref_id <= num_refs_in_fuller_bins):
        bin_id = math.ceil(ref_id / num_refs_in_full_bin) -1
    else:
        pos_in_lesser_bins = ref_id - num_refs_in_fuller_bins
        lesser_bin_id = math.ceil(pos_in_lesser_bins / num_refs_in_lesser_bin)-1
        bin_id = num_fuller_bins + lesser_bin_id
    return bin_id

gt_dic = {}
TP = 0
FP = 0
FN = 0
lambda_gt = []


with open(gt_file) as gt:
    csv_reader = csv.reader(gt, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            row = list(map(int, row))
            # row[1] = ref_id, row[2] = query_id, row[0] = bin_id
            gt_dic.setdefault(row[1], set()).add(row[2])

#print("gt: ", gt_dic)

with open(results_file) as res:
    for line in res:
        hit = line.rstrip('\n').split(',')
        hit = list(map(int, hit))
        for i in range(1,len(hit)):
            #hit[0] = ref_id
            #hit[i] = qry_id(s)
            #prepare lambda_as_gt: bin, ref, qry
            lambda_gt.append([ref_to_bin(hit[0]), hit[0], hit[i]])
            # check FP
            if hit[0] in gt_dic:
                if hit[i] in gt_dic[hit[0]]:
                    # true positive
                    TP += 1
                    # delete already found matches, so we can check FN later
                    gt_dic[hit[0]].remove(hit[i])
                    if len(gt_dic[hit[0]]) == 0:
                        del gt_dic[hit[0]]
                else:
                    # bin has hit but not this one
                    FP += 1
            else:
                # this bin has no hit (anymore), everything that maps to it is a false positive
                FP += 1

for elem in gt_dic:
    FN += len(gt_dic[elem])

#print("TP: " + str(TP) + " FP: " + str(FP) + " FN: " + str(FN))

with open(snakemake.output[0], 'w') as out_file:
    out_file.write("TP: " + str(TP) + " FP: " + str(FP) + " FN: " + str(FN))

# write lambda as ground truth
lambda_gt_sorted = sorted(lambda_gt, key=lambda x: x[0])

with open(snakemake.output[1], 'w') as lambda_gt_file:
    lambda_gt_file.write("bin_id, ref_id, qry_id\n")
    for elem in lambda_gt_sorted:
        lambda_gt_file.write(str(elem[0]) + ", " + str(elem[1]) + ", " + str(elem[2]) + "\n")