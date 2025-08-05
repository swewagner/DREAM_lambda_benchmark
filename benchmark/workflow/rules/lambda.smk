# lambda on whole database
def getSbjDomain(wildcards):
    if wildcards.blast_mode == "blastN" or wildcards.blast_mode == "tBlastN":
        return "results/er_{er}/nuc/ref_seqs.fasta"
    else:
        return "results/er_{er}/prot/ref_seqs.fasta"

rule make_index:
    input:
        db = getSbjDomain
    output:
        index = "results/er_{er}/{blast_mode}/lambda/db.lba"
    params:
        binary_dir = "../build/lambda3/src/lambda3-build/bin",
        bench = "benchmarks/lambda_index.time"
    log:
        "logs/lambda/index/{er}/{blast_mode}.log"
    shell:
        """
        if [ {wildcards.blast_mode} = "blastN" ]; then
            /usr/bin/time -a -o {params.bench} -f "%e\t%M\t%x\tlambda-index\t{wildcards.er}\t{wildcards.blast_mode}" \
            {params.binary_dir}/lambda3 mkindexn -d {input.db} -i {output.index} -v 2 &> {log}
        else
            /usr/bin/time -a -o {params.bench} -f "%e\t%M\t%x\tlambda-index\t{wildcards.er}\t{wildcards.blast_mode}" \
            {params.binary_dir}/lambda3 mkindexp -d {input.db} -i {output.index} -v 2 &> {log}
        fi
        """

rule search:
    input:
        query = getQueryDomain,
        index = "results/er_{er}/{blast_mode}/lambda/db.lba"
    output:
        "results/er_{er}/{blast_mode}/lambda/search_out.m8"
    params:
        binary_dir = "../build/lambda3/src/lambda3-build/bin",
        bench = "benchmarks/lambda_search.time"
    log:
        "logs/lambda/search/{er}/combined_{blast_mode}.log"
    shell:
        """
        if [ {wildcards.blast_mode} = "blastN" ]; then
            /usr/bin/time -a -o {params.bench} -f "%e\t%M\t%x\tlambda-search\t{wildcards.er}\t{wildcards.blast_mode}" \
            {params.binary_dir}/lambda3 searchn -q {input.query} -i {input.index} -o {output} -v 2 -t 2 &> {log}
        else
            /usr/bin/time -a -o {params.bench} -f "%e\t%M\t%x\tlambda-search\t{wildcards.er}\t{wildcards.blast_mode}" \
            {params.binary_dir}/lambda3 searchp -q {input.query} -i {input.index} -o {output} -v 2 &> {log}
        fi
        """

def getSbjDomainSparse(wildcards):
    if wildcards.blast_mode == "blastN":
        return "results/er_{er}/nuc/bins/bin_{sparse_id}.fasta"
    else:
        return "results/er_{er}/prot/bins/bin_{sparse_id}.fasta"

rule make_index_sparse:
    input:
        db = getSbjDomainSparse
    output:
        index = "results/er_{er}/{blast_mode}/lambda/sparse_dbs/db_{sparse_id}.lba"
    params:
        bench = "benchmarks/lambda_sparse_index.time"
    log:
        "logs/lambda/index_sparse/{er}/{blast_mode}/log_{sparse_id}.log"
    shell:
        """
        /usr/bin/time -a -o {params.bench} -f "%e\t%M\t%x\tlambda-sparse-index\t{wildcards.er}\t{wildcards.blast_mode}\t{wildcards.sparse_id}" \
        ./workflow/scripts/lambda_index.sh {input.db} {output.index} {wildcards.blast_mode} {log}
        """
        
def getQueryForSparseLambda(wildcards):
    if wildcards.blast_mode == "blastN" or wildcards.blast_mode == "blastX":
        query = "results/er_{er}/{blast_mode}/iota/kmer{kmer}_error{error}_alphabet{alph}/queries_nuc_bin_{sparse_id}.fasta"
    else:
        query = "results/er_{er}/{blast_mode}/iota/kmer{kmer}_error{error}_alphabet{alph}/queries_prot_bin_{sparse_id}.fasta"
    return query

rule search_sparse:
    input:
        query = getQueryForSparseLambda,
        index = "results/er_{er}/{blast_mode}/lambda/sparse_dbs/db_{sparse_id}.lba"
    output:
        "results/er_{er}/{blast_mode}/lambda/sparse_out/kmer{kmer}_error{error}_alphabet{alph}/out_{sparse_id}.m8"
    params:
        bench = "benchmarks/lambda_sparse_search.time"
    log:
        "logs/lambda/search_sparse/{er}/{blast_mode}/kmer{kmer}_error{error}_alphabet{alph}/log_{sparse_id}.log"
    shell:
        """
        /usr/bin/time -a -o {params.bench} -f "%e\t%M\t%x\tlambda-sparse-search\t{wildcards.er}\t{wildcards.blast_mode}\t{wildcards.sparse_id}" \
        ./workflow/scripts/lambda_search.sh {input.query} {input.index} {output} {wildcards.blast_mode} {log}"""


rule extract_lambda_results_sparse:
    input:
        get_sparse_ids
    output:
        out = "results/er_{er}/{blast_mode}/lambda/sparse_out/kmer{kmer}_error{error}_alphabet{alph}/results.txt"
    script:
        "../scripts/lambda_summary_sparse.py"