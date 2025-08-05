# -------------------------------
# simulate dna sequences
# -------------------------------

rule simulate_ref_nuc:
    output:
        ref = "results/pre_ref_seqs_nuc.fasta"
    log:
        "logs/simulation/simulate_ref_nuc.log"
    shell:
        "./workflow/scripts/simulate_seq.sh {output.ref} {num_of_refseqs} {ref_seq_len} {log}"


rule simulate_query_nuc:
    output:
        queries = "results/queries_nuc.fasta"
    log:
        "logs/simulation/simulate_ref_nuc.log"
    shell:
        "./workflow/scripts/simulate_seq.sh {output.queries} {num_of_qseqs} {query_read_len} {log}"


rule distribute_matches_nuc:
    input:
        queries = "results/queries_nuc.fasta",
        pre_refs = "results/pre_ref_seqs_nuc.fasta"
    output:
        ground_truth = "results/er_{er}/nuc/ground_truth.txt",
        refs = expand("results/er_{{er}}/nuc/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    params:
        outdir = "results/er_{er}/nuc/bins"
    log:
        "logs/simulation/{er}_distribute_nuc.log"
    shell:
        "./workflow/scripts/distribute_matches.sh {input.queries} {input.pre_refs} {output.ground_truth} {params.outdir} {match_len} {num_of_refseqs} {num_of_bins} {wildcards.er} {log}"


# -------------------------------
# translate to prot sequences
# -------------------------------

rule translate_qry:
    input:
        dna = "results/queries_nuc.fasta"
    output:
        prot = "results/queries_prot.fasta"
    shell:
        "./workflow/scripts/translate_ref.sh {input.dna} {output.prot}"

rule translate_ref:
    input:
        dna = "results/pre_ref_seqs_nuc.fasta"
    output:
        prot = "results/pre_ref_seqs_prot.fasta"
    shell:
        "./workflow/scripts/translate_ref.sh {input.dna} {output.prot}"

rule distribute_matches_prot:
    input:
        queries = "results/queries_prot.fasta",
        pre_refs = "results/pre_ref_seqs_prot.fasta"
    output:
        ground_truth = "results/er_{er}/prot/ground_truth.txt",
        refs = expand("results/er_{{er}}/prot/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    log:
        "logs/simulation/{er}_distribute_prot.log"
    params:
        outdir = "results/er_{er}/prot/bins"
    shell:
        "./workflow/scripts/distribute_matches_prot.sh {input.queries} {input.pre_refs} {output.ground_truth} {params.outdir} {match_len} {num_of_refseqs} {num_of_bins} {wildcards.er} {log}"


rule create_all_bin_paths:
    input:
        refs = expand("results/er_{{er}}/{{domain}}/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    output:
        bin_paths = "results/er_{er}/{domain}/all_bin_paths.txt"
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
        expand("results/er_{{er}}/{{domain}}/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    output:
        combined_bins = "results/er_{er}/{domain}/ref_seqs.fasta"
    params:
        max_bin = num_of_bins-1,
        bin_path = "results/er_{er}/{domain}/bins/bin_",
        ext = ".fasta"
    shell:
        """
        for i in $(seq 0 {params.max_bin}); do cat {params.bin_path}$i{params.ext} >> {output.combined_bins}; done
        """