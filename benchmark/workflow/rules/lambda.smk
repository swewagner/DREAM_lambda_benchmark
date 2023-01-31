rule make_index:
    input:
        db = "results/ref_seqs_prot.fasta"
    output:
        index = "results/db.lba"
    shell:
        "./workflow/scripts/lambda_indexp.sh {input.db} {output.index}"

rule search_prot:
    input:
        query = "results/query_seqs.fasta",
        index = "results/db.lba"
    output:
        "results/lambda_out.m8"
    shell:
        "./workflow/scripts/lambda_searchp.sh {input.query} {input.index} {output}"
