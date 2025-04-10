rule simulate_ref_dna:
    output:
        ref = "results/pre_ref_seqs_dna.fasta"
    shell:
        "./workflow/scripts/simulate_seq.sh {output.ref} {num_of_refseqs} {ref_seq_len}"


rule simulate_query_dna:
    output:
        queries = "results/queries.fasta"
    shell:
        "./workflow/scripts/simulate_seq.sh {output.queries} {num_of_qseqs} {query_read_len}"


rule distribute_matches:
    input:
        queries = "results/queries.fasta",
        pre_refs = "results/pre_ref_seqs_dna.fasta"
    output:
        ground_truth = "results/er_{er}/ground_truth.txt",
        refs = "results/er_{er}/bins/bin_00.fasta"
    params:
        outdir = "results/er_{er}/bins"
    shell:
        "./workflow/scripts/distribute_matches.sh {input.queries} {input.pre_refs} {output.ground_truth} {params.outdir} {match_len} {num_of_refseqs} {num_of_bins} {wildcards.er}"


rule translate_qry:
    input:
        dna = "results/queries.fasta"
    output:
        prot = "results/queries_prot.fasta"
    shell:
        "./workflow/scripts/translate_ref.sh {input.dna} {output.prot}"


rule create_all_bin_paths:
    input:
        refs = "results/er_{er}/bins"
    output:
        bin_paths = "results/er_{er}/all_bin_paths.txt"
    params:
        ls_cmd = "bins/$output"
    shell:
        "for output in $(ls {input.refs}); do echo {params.ls_cmd} >> {output.bin_paths}; done"