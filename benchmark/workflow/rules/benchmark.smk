rule extract_iota_results:
    input:
        dummy = "results/er_{er}/{blast_mode}/iota/dummy.txt"
    output:
        results = "results/er_{er}/{blast_mode}/iota/results.txt"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/iota_summary.py"


rule compare_iota_groundtruth:
    input:
        ground_truth = "results/er_{er}/ground_truth.txt",
        iota_results = "results/er_{er}/{blast_mode}/iota/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/iota_groundtruth.txt"
    script:
        "../scripts/compare_iota_groundtruth.py"


rule extract_lambda_results:
    input:
        lambda_out = expand("results/er_{{er}}/{{blast_mode}}/lambda/out_{bin_id}.m8", bin_id = bin_ids)
    output:
        out = "results/er_{er}/{blast_mode}/lambda/results.txt"
    script:
        "../scripts/lambda_summary.py"


rule compare_iota_lambda:
    input:
        lambda_res = "results/er_{er}/{blast_mode}/lambda/results.txt",
        iota_res = "results/er_{er}/{blast_mode}/iota/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/iota_lambda.txt"
    script:
        "../scripts/compare_iota_lambda.py"