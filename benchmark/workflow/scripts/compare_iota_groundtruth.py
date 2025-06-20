import csv

gt_file = snakemake.input[0]
results_file = snakemake.input[1]
#gt_file = "../../results/er_0.025/ground_truth.txt"
#results_file = "../../results/er_0.025/blastN/iota/results.txt"


gt_dic = {}
TP = 0
FP = 0
FN = 0

with open(gt_file) as gt:
    csv_reader = csv.reader(gt, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            row = list(map(int, row))
            gt_dic.setdefault(row[0], set()).add(row[2])

print(gt_dic)

with open(results_file) as res:
    for line in res:
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

#print("TP: " + str(TP) + " FP: " + str(FP) + " FN: " + str(FN))

with open(snakemake.output[0], 'w') as out_file:
    out_file.write("TP: " + str(TP) + " FP: " + str(FP) + " FN: " + str(FN))