rule simulate_ref_dna:
    output:
        ref = "results/ref_seqs_dna.fasta"
    shell:
        "./workflow/scripts/simulate_seq.sh {output.ref} {num_of_refseqs} {ref_seq_len}"


rule simulate_query_pre:
    output:
        queries = "results/pre_queries.fasta"
    shell:
        "./workflow/scripts/simulate_seq.sh {output.queries} {num_of_qseqs} {query_read_len}"


rule simulate_matches:
    input:
        ref = "results/ref_seqs_dna.fasta",
    output:
        matches = "results/er_{er}/matches.fasta",
        ground_truth = "results/er_{er}/ground_truth.txt"
    shell:
        "./workflow/scripts/simulate_matches.sh {input.ref} {output.matches} {output.ground_truth} {num_of_qseqs} {match_len} {num_of_refseqs} {wildcards.er}"


rule insert_matches:
    input:
        queries = "results/pre_queries.fasta",
        matches = "results/er_{er}/matches.fasta"
    output:
        "results/er_{er}/final_queries.fasta"
    shell:
        "./workflow/scripts/insert_matches.sh {input.queries} {input.matches} {output}"


rule translate_ref:
    input:
        dna = "results/ref_seqs_dna.fasta"
    output:
        prot = "results/ref_seqs_prot.fasta"
    shell:
        "./workflow/scripts/translate_ref.sh {input.dna} {output.prot}"


rule split_into_bins:
    input:
        ref = "results/ref_seqs_dna.fasta"
    output:
        bins = expand("results/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/split_into_bins.py"