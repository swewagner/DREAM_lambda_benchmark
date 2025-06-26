def getGroundTruthDomain(wildcards):
    if wildcards.blast_mode == "blastN" or wildcards.blast_mode == "tblastN":
        return "results/er_{er}/nuc/ground_truth.txt"
    else:
        return "results/er_{er}/prot/ground_truth.txt"

rule compare_iota_groundtruth:
    input:
        ground_truth = getGroundTruthDomain,
        iota_results = "results/er_{er}/{blast_mode}/iota/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/iota_groundtruth.txt"
    script:
        "../scripts/compare_iota_groundtruth.py"

rule compare_lambda_groundtruth:
    input:
        ground_truth = getGroundTruthDomain,
        lambda_results = "results/er_{er}/{blast_mode}/lambda/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/lambda_groundtruth.txt"
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
        lambda_res = "results/er_{er}/{blast_mode}/lambda/results.txt",
        iota_res = "results/er_{er}/{blast_mode}/iota/results.txt"
    output:
        out = "results/er_{er}/{blast_mode}/benchmark/iota_lambda.txt"
    script:
        "../scripts/compare_iota_lambda.py"


rule setup_time_benchmarks:
    output:
        expand("benchmarks/{name}.time", name=benchmark_file_names)
    run:
        shell("mkdir -p benchmarks")
        for out in output:
            if "sparse" in out:
                with open(out, "w") as f:
                    f.write("time\tmem\terror-code\tcommand\terror-rate\tblast-mode\tbin-id\n")
            else:
                with open(out, "w") as f:
                    f.write("time\tmem\terror-code\tcommand\terror-rate\tblast-mode\n")

