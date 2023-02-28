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
        queries = "results/er_{er}/query_seqs.fasta",
        ground_truth = "results/er_{er}/ground_truth.txt"
    shell:
        "./workflow/scripts/simulate_queries.sh {input.ref} {output.queries} {output.ground_truth} {num_of_qseqs} {query_read_len} {num_of_refseqs} {wildcards.er}"


rule translate_ref:
    input:
        dna = "results/ref_seqs_dna.fasta"
    output:
        prot = "results/ref_seqs_prot.fasta"
    shell:
        "./workflow/scripts/translate_ref.sh {input.dna} {output.prot}"