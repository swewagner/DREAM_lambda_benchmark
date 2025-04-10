rule make_index:
    input:
        db = "results/er_{er}/bins/bin_{bin_id}.fasta"
    output:
        index = "results/er_{er}/dbs/db_{bin_id}.lba"
    log:
        "logs/lambda/indexp_{er}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_indexp.sh {input.db} {output.index} {log}"


rule search_prot:
    input:
        query = "results/queries_prot.fasta",
        index = "results/er_{er}/dbs/db_{bin_id}.lba"
    output:
        "results/er_{er}/lambda_out/lambda_out_{bin_id}.m8"
    log:
        "logs/lambda/searchp_{er}/log_{bin_id}.log"
    shell:
        "./workflow/scripts/lambda_searchp.sh {input.query} {input.index} {output} {log}"
