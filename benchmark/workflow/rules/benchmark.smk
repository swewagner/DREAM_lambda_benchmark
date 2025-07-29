def getGroundTruthDomain(wildcards):
    if wildcards.blast_mode == "blastN" or wildcards.blast_mode == "tBlastN":
        return "results/er_{er}/nuc/ground_truth.txt"
    else:
        return "results/er_{er}/prot/ground_truth.txt"

rule compare_iota_groundtruth:
    input:
        ground_truth = getGroundTruthDomain,
        iota_results = "results/er_{er}/{blast_mode}/iota/kmer{kmer}_error{error}_alphabet{alph}/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/kmer{kmer}_error{error}_alphabet{alph}/iota_groundtruth.txt"
    script:
        "../scripts/compare_iota_groundtruth.py"

rule compare_lambda_groundtruth:
    input:
        ground_truth = getGroundTruthDomain,
        lambda_results = "results/er_{er}/{blast_mode}/lambda/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/lambda_groundtruth.txt",
        l_gt = "results/er_{er}/{blast_mode}/benchmark/lambda_as_groundtruth.txt"
    script:
        "../scripts/compare_lambda_groundtruth.py"


rule extract_lambda_results:
    input:
        lambda_out = "results/er_{er}/{blast_mode}/lambda/search_out.m8"
    output:
        out = "results/er_{er}/{blast_mode}/lambda/results.txt"
    script:
        "../scripts/lambda_summary.py"


rule compare_iota_lambda:
    input:
        lambda_gt = "results/er_{er}/{blast_mode}/benchmark/lambda_as_groundtruth.txt",
        iota_res = "results/er_{er}/{blast_mode}/iota/kmer{kmer}_error{error}_alphabet{alph}/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/kmer{kmer}_error{error}_alphabet{alph}/iota_lambda.txt"
    script:
        "../scripts/compare_iota_lambda.py"



