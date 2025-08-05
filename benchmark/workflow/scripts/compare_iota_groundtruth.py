import csv

gt_file = snakemake.input[0]
results_file = snakemake.input[1]
#gt_file = "../../results/er_0.025/ground_truth.txt"
#results_file = "../../results/er_0.025/blastN/iota/results.txt"

blast_to_frame = {"blastN": 2, "blastP": 1, "blastX": 6, "tBlastN":6}

qseq_num = int(snakemake.config["num_of_qseqs"])
bin_num = int(snakemake.config["num_of_bins"])
frame_num = blast_to_frame[snakemake.wildcards.blast_mode]

num_of_possible_matches = qseq_num * bin_num * frame_num


gt_dic = {}
num_actual_matches = 0
TP = 0
FP = 0
FN = 0
bins_with_hit = 0
bins_with_true_hit = 0

with open(gt_file) as gt:
    csv_reader = csv.reader(gt, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            row = list(map(int, row))
            gt_dic.setdefault(row[0], set()).add(row[2])

for elem in gt_dic:
    bins_with_true_hit +=1
    num_actual_matches += len(gt_dic[elem])

with open(results_file) as res:
    for line in res:
        bins_with_hit += 1
        hit = line.rstrip('\n').split(',')
        hit = list(map(int, hit))
        for i in range(1,len(hit)):
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

TN = num_of_possible_matches - num_actual_matches - FP


tpr = "undefined"
tnr = "undefined"
if (TP+FN > 0):
    tpr = TP / (TP+FN)
if (TN+FP > 0):
    tnr = TN / (TN+FP)

with open(snakemake.output[0], 'w') as out_file:
    out_file.write("TP: " + str(TP) + " FP: " + str(FP) + " FN: " + str(FN) + " TN: " + str(TN) + " BTM: "+ str(bins_with_true_hit) + " B: " + str(bins_with_hit) + "\n" + "Sens: " + str(tpr) + " Spec: " + str(tnr))