rule run_iota_blastN:
    input:
        bin_paths = "results/er_{er}/all_bin_paths.txt",
        query = "results/queries.fasta"
    output:
        "results/er_{er}/blastN/iota/dummy.txt"
    log:
        "logs/iota/er_{er}/blastN.log"
    params:
        out_dir = "results/er_{er}/blastN/iota"
    shell:
        "./workflow/scripts/run_iota_blastN.sh {input.bin_paths} {input.query} {params.out_dir} {log} {output}"


rule run_iota_tBlastN:
    input:
        bin_paths = "results/er_{er}/all_bin_paths.txt",
        prot_query = "results/queries_prot.fasta"
    output:
        "results/er_{er}/tBlastN/iota/dummy.txt"
    log:
        "logs/iota/er_{er}/tBlastN.log"
    params:
        out_dir = "results/er_{er}/tBlastN/iota"
    shell:
        "./workflow/scripts/run_iota_tBlastN.sh {input.bin_paths} {input.prot_query} {params.out_dir} {log} {output}"