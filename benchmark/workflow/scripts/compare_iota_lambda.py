import csv

lambda_gt = snakemake.input[0]
iota_file = snakemake.input[1]

gt_dic = {}
TP = 0
FP = 0
FN = 0

with open(lambda_gt) as lf:
    csv_reader = csv.reader(lf, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            row = list(map(int, row))
            # row[0] = bin_id, row[1] = ref_id, row[2] = query_id
            gt_dic.setdefault(row[0], set()).add(row[2])

with open(iota_file) as res:
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

with open(snakemake.output[0], 'w') as out_file:
    out_file.write("TP: " + str(TP) + " FP: " + str(FP) + " FN: " + str(FN))