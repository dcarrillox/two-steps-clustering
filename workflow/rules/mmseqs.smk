configfile: "config/config.yaml"

results_dir = config["results_dir"]

rule run_mmseqs2:
    input:
        config["proteome"]
    output:
        results_dir + "/0_clustering/results/clustering.tsv"
    params:
        db = "proteome",
        db_path = results_dir + "/0_clustering/db/",
        cluster_db = "clustering",
        cluster_db_path = results_dir + "/0_clustering/results/",
        min_seq_id = config["mmseqs2_params"]["min_seq_id"],
        cov_mode = config["mmseqs2_params"]["cov_mode"],
        sensitivity = config["mmseqs2_params"]["sensitivity"],
        coverage = config["mmseqs2_params"]["coverage"]
    threads: config["mmseqs2_params"]["threads"]
    conda:
        "../envs/env.yaml"
    shell:
        '''
        mkdir -p tmp ;
        mkdir -p {params.db_path} ;
        mkdir -p {params.cluster_db_path} ;

        mmseqs createdb {input} {params.db_path}/{params.db} ;

        mmseqs cluster {params.db_path}/{params.db} {params.cluster_db_path}/{params.cluster_db} tmp \
        --min-seq-id {params.min_seq_id} \
        --cov-mode {params.cov_mode} \
        --threads {threads} \
        -s {params.sensitivity} \
        -c {params.coverage} ;
        rm -rf tmp/*
        mmseqs createtsv {params.db_path}/{params.db} {params.db_path}/{params.db} {params.cluster_db_path}/{params.cluster_db} {output}

        '''


checkpoint create_faa_clusters:
    input:
        rules.run_mmseqs2.output
    output:
        directory(results_dir + "/0_clustering/clusters_faa")
    params:
        script = "workflow/scripts/1_mmseqs2tsv_to_faa.py",
        faa_all = config["proteome"]
    conda:
        "../envs/env.yaml"
    shell:
        "mkdir -p {output}; python {params.script} -t {input} -f {params.faa_all} -o {output}"



rule msa_clusters:
    input:
        results_dir + "/0_clustering/clusters_faa/{cluster}.faa"
    output:
        results_dir + "/0_clustering/clusters_msa/{cluster}.mafft"
    threads: 40
    conda:
        "../envs/env.yaml"
    shell:
        "mafft-fftnsi --quiet --maxiterate 200 --thread {threads} {input} > {output}"
