rule run_iota_blastN:
    input:
        bin_paths = "results/er_{er}/all_bin_paths.txt",
        query = "results/queries.fasta"
    output:
        "results/er_{er}/iota/blastN/dummy.txt"
    log:
        "logs/iota/er_{er}/blastN.log"
    params:
        out_dir = "results/er_{er}/iota/blastN"
    shell:
        "./workflow/scripts/run_iota_blastN.sh {input.bin_paths} {input.query} {params.out_dir} {log} {output}"