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
        refs = expand("results/er_{{er}}/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
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
        refs = expand("results/er_{{er}}/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    output:
        bin_paths = "results/er_{er}/all_bin_paths.txt"
    params:
        max_bin = num_of_bins-1,
        bin_path = "bins/bin_",
        ext = ".fasta"
    shell:
        """
        for i in $(seq 0 {params.max_bin}); do
            echo {params.bin_path}$i{params.ext} >> {output.bin_paths}
        done
        """

rule combine_bins_into_one_file:
    input:
        expand("results/er_{{er}}/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    output:
        combined_bins = "results/er_{er}/ref_seqs.fasta"
    params:
        max_bin = num_of_bins-1,
        bin_path = "results/er_{er}/bins/bin_",
        ext = ".fasta"
    shell:
        """
        for i in $(seq 0 {params.max_bin}); do cat {params.bin_path}$i{params.ext} >> {output.combined_bins}; done
        """