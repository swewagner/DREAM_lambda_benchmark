rule run_iota_blastN:
    input:
        bin_paths = "results/er_{er}/all_bin_paths.txt",
        query = "results/queries.fasta",
        bench = "benchmarks/iota.time"
    output:
        "results/er_{er}/blastN/iota/dummy.txt"
    log:
        "logs/iota/er_{er}/blastN.log"
    params:
        out_dir = "results/er_{er}/blastN/iota"
    shell:
        """
        /usr/bin/time -a -o {input.bench} -f "%e\t%M\t%x\tiota\tblastN\t{wildcards.er}" \
        ./workflow/scripts/run_iota_blastN.sh {input.bin_paths} {input.query} {params.out_dir} {log} {output}"""


rule run_iota_tBlastN:
    input:
        bin_paths = "results/er_{er}/all_bin_paths.txt",
        prot_query = "results/queries_prot.fasta",
        bench = "benchmarks/iota.time"
    output:
        "results/er_{er}/tBlastN/iota/dummy.txt"
    log:
        "logs/iota/er_{er}/tBlastN.log"
    params:
        out_dir = "results/er_{er}/tBlastN/iota"
    shell:
        """
        /usr/bin/time -a -o {input.bench} -f "%e\t%M\t%x\tiota\ttBlastN\t{wildcards.er}" \
        ./workflow/scripts/run_iota_tBlastN.sh {input.bin_paths} {input.prot_query} {params.out_dir} {log} {output}
        """


checkpoint extract_iota_results:
    input:
        dummy = "results/er_{er}/{blast_mode}/iota/dummy.txt"
    output:
        results = "results/er_{er}/{blast_mode}/iota/results.txt"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/iota_summary.py"

def get_sparse_ids(wildcards):
    co = checkpoints.extract_iota_results.get(**wildcards).output.results
    if wildcards.blast_mode == "blastN":
        sparse_ids = glob_wildcards(os.path.join(os.path.dirname(co), "queries_bin_{sparse_id}.fasta")).sparse_id
    else:
        sparse_ids = glob_wildcards(os.path.join(os.path.dirname(co), "queries_prot_bin_{sparse_id}.fasta")).sparse_id
    return expand(
        "results/er_{er}/{blast_mode}/lambda/sparse_out/out_{sparse_id}.m8",
        er = wildcards.er,
        blast_mode = wildcards.blast_mode,
        sparse_id = sparse_ids
    )