# lambda without prefilter, but on each bin

rule make_index_nuc:
    input:
        db = "results/er_{er}/bins/bin_{bin_id}.fasta"
    output:
        index = "results/er_{er}/lambda_nuc_dbs/db_{bin_id}.lba"
    log:
        "logs/lambda/indexn_{er}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_indexn.sh {input.db} {output.index} {log}"


rule search_nuc:
    input:
        query = "results/queries.fasta",
        index = "results/er_{er}/lambda_nuc_dbs/db_{bin_id}.lba"
    output:
        "results/er_{er}/blastN/lambda/out_{bin_id}.m8"
    log:
        "logs/lambda/searchn_{er}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_searchn.sh {input.query} {input.index} {output} {log}"


rule make_index_prot:
    input:
        db = "results/er_{er}/bins/bin_{bin_id}.fasta"
    output:
        index = "results/er_{er}/lambda_prot_dbs/db_{bin_id}.lba"
    log:
        "logs/lambda/indexp_{er}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_indexp.sh {input.db} {output.index} {log}"


rule search_prot:
    input:
        query = "results/queries_prot.fasta",
        index = "results/er_{er}/lambda_prot_dbs/db_{bin_id}.lba"
    output:
        "results/er_{er}/tBlastN/lambda/out_{bin_id}.m8"
    log:
        "logs/lambda/searchp_{er}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_searchp.sh {input.query} {input.index} {output} {log}"


# lambda on bins found by iota

# rule make_index_sparse_nuc:
#     input:
#         db = "results/er_{er}/bins/bin_{sparse_bin_id}.fasta"
#     output:
#         index = "results/er_{er}/lambda_sparse_nuc_dbs/db_{sparse_bin_id}.lba"
#     log:
#         "logs/lambda/sparse/indexn_{er}/log_{sparse_bin_id}.log"
#     shell:
#         "./workflow/scripts/lambda_indexn.sh {input.db} {output.index} {log}"


# rule search_sparse_nuc:
#     input:
#         query = "results/er_{er}/blastN/iota/queries_bin_{sparse_bin_id}.fasta",
#         index = "results/er_{er}/lambda_sparse_nuc_dbs/db_{sparse_bin_id}.lba"
#     output:
#         "results/er_{er}/blastN/lambda_sparse/out_{sparse_bin_id}.m8"
#     log:
#         "logs/lambda/sparse/searchn_{er}/log_{sparse_bin_id}.log"
#     shell:
#         "./workflow/scripts/lambda_searchn.sh {input.query} {input.index} {output} {log}"
