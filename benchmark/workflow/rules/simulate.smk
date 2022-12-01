rule simulate_ref_dna:
    output:
        ref = "results/ref_seqs_dna.fasta"
    threads:
        8
    shell:
        "./workflow/scripts/simulate_refs.sh {output.ref} {num_of_refseqs} {ref_seq_len}"


rule simulate_queries:
    input:
        ref = "results/ref_seqs_dna.fasta"
    output:
        queries = "results/query_seqs.fasta"
    shell:
        "./workflow/scripts/simulate_queries.sh {input.ref} {output.queries} {num_of_qseqs} {query_read_len}"
## TODO write script that takes random chunks of ref sequences
## simulate_queries should use a c++ binary
## write ground truth? should every query match? yes -> #queries = matches
#
#
#rule translate_ref:
#    input:
#        dna = "results/ref_seqs_dna.fasta"
#    output:
#        prot = "results/ref_seqs_prot.fasta"
#    shell:
#        "../scripts/translate_ref.sh {input.dna} {output.prot}"