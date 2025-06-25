# lambda on whole database
rule make_index:
    input:
        db = "results/er_{er}/ref_seqs.fasta",
        bench = "benchmarks/lambda_index.time"
    output:
        index = "results/er_{er}/lambda/nuc_db.lba"
    params:
        binary_dir = "../build/lambda3/src/lambda3-build/bin"
    log:
        "logs/lambda/index/{er}/nuc.log"
    shell:
        """
        /usr/bin/time -a -o {input.bench} -f "%e\t%M\t%x\tlambda-index\t{wildcards.er}" \
        {params.binary_dir}/lambda3 mkindexn -d {input.db} -i {output.index} &> {log}
        """


rule search:
    input:
        query = "results/queries.fasta",
        index = "results/er_{er}/lambda/nuc_db.lba",
        bench = "benchmarks/lambda_search.time"
    output:
        "results/er_{er}/lambda/blastN_search.m8"
    params:
        binary_dir = "../build/lambda3/src/lambda3-build/bin"
    log:
        "logs/lambda/search/{er}/log_combined_blastN.log"
    shell:
        """
        /usr/bin/time -a -o {input.bench} -f "%e\t%M\t%x\tlambda-search\t{wildcards.er}" \
        {params.binary_dir}/lambda3 searchn -q {input.query} -i {input.index} -o {output} &> {log}
        """


# lambda without prefilter, but on each bin

rule make_index_bins:
    input:
        db = "results/er_{er}/bins/bin_{bin_id}.fasta"
    output:
        index = "results/er_{er}/{blast_mode}/lambda/dbs/db_{bin_id}.lba"
    log:
        "logs/lambda/index/{er}/{blast_mode}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_index.sh {input.db} {output.index} {wildcards.blast_mode} {log}"


rule search_bins:
    input:
        query = "results/queries.fasta",
        index = "results/er_{er}/{blast_mode}/lambda/dbs/db_{bin_id}.lba"
    output:
        "results/er_{er}/{blast_mode}/lambda/out_{bin_id}.m8"
    log:
        "logs/lambda/search/{er}/{blast_mode}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_search.sh {input.query} {input.index} {output} {wildcards.blast_mode} {log}"


rule make_index_sparse:
    input:
        db = "results/er_{er}/bins/bin_{sparse_id}.fasta",
        bench = "benchmarks/lambda_sparse_index.time"
    output:
        index = "results/er_{er}/{blast_mode}/lambda/sparse_dbs/db_{sparse_id}.lba"
    log:
        "logs/lambda/index_sparse/{er}/{blast_mode}/log_{sparse_id}.log"
    shell:
        """
        /usr/bin/time -a -o {input.bench} -f "%e\t%M\t%x\tlambda-sparse-index\t{wildcards.blast_mode}\t{wildcards.er}\t{wildcards.sparse_id}" \
        ./workflow/scripts/lambda_index.sh {input.db} {output.index} {wildcards.blast_mode} {log}
        """
        
def getQueryForSparseLambda(wildcards):
    if wildcards.blast_mode == "blastN":
        query = "results/er_{er}/{blast_mode}/iota/queries_bin_{sparse_id}.fasta"
    else:
        query = "results/er_{er}/{blast_mode}/iota/queries_prot_bin_{sparse_id}.fasta"
    return query

rule search_sparse:
    input:
        #query = "results/er_{er}/{blast_mode}/iota/queries_bin_{sparse_id}.fasta",
        query = getQueryForSparseLambda,
        index = "results/er_{er}/{blast_mode}/lambda/sparse_dbs/db_{sparse_id}.lba",
        bench = "benchmarks/lambda_sparse_search.time"
    output:
        "results/er_{er}/{blast_mode}/lambda/sparse_out/out_{sparse_id}.m8"
    log:
        "logs/lambda/search_sparse/{er}/{blast_mode}/log_{sparse_id}.log"
    shell:
        """
        /usr/bin/time -a -o {input.bench} -f "%e\t%M\t%x\tlambda-sparse-search\t{wildcards.blast_mode}\t{wildcards.er}\t{wildcards.sparse_id}" \
        ./workflow/scripts/lambda_search.sh {input.query} {input.index} {output} {wildcards.blast_mode} {log}"""


rule extract_lambda_results_sparse:
    input:
        get_sparse_ids
    output:
        out = "results/er_{er}/{blast_mode}/lambda/sparse_out/results.txt"
    script:
        "../scripts/lambda_summary.py"