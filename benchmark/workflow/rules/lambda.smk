rule make_index:
    input:
        db = "results/ref_seqs_prot.fasta"
    output:
        index = "results/db.lba"
    log:
        "logs/lambda/indexp.log"
    shell:
        "./workflow/scripts/lambda_indexp.sh {input.db} {output.index} {log}"

rule search_prot:
    input:
        query = "results/er_{er}/query_seqs.fasta",
        index = "results/db.lba"
    output:
        "results/er_{er}/lambda_out.m8"
    log:
        "logs/lambda/searchp_{er}.log"
    shell:
        "./workflow/scripts/lambda_searchp.sh {input.query} {input.index} {output} {log}"
