rule create_bin_path_list:
    input:
        bins_without_linebreaks = expand("results/bins/bin_{bin_id}.fasta", bin_id=bin_ids)
    output:
        binfile = "results/all_bin_paths.txt"
    log:
        "logs/valik/create_bin_path_list.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_bin_list.py"


rule build_index:
    input:
        binfile = "results/all_bin_paths.txt"
    output:
        index = "results/valik_index.ibf"
    log:
        "logs/valik/build.log"
    params:
        window = config["window_size"],
        kmer = config["kmer_size"],
        size = config["index_size"]
    shell:
        "./workflow/scripts/valik_build.sh {input.binfile} {output.index} {log} {params.window} {params.kmer} {params.size}"


rule search:
    input:
        query = "results/er_{er}/final_queries.fasta",
        index = "results/valik_index.ibf"
    output:
        "results/er_{er}/search.gff"
    log:
        "logs/valik/search_{er}.log"
    params:
        pattern = config["pattern_size"],
        error = lambda w: str(int(int(config["pattern_size"]) * float(w.er)))
    threads:
        4
    shell:
        "./workflow/scripts/valik_search.sh {input.query} {input.index} {output} {log} {params.pattern} {params.error} {threads}"