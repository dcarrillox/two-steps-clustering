rule mmseqs2_subfamily:
    input:
        config["proteome"]
    output:
        "results/0_mmseqs/results/mmseqs_clustering.tsv"
    conda:
        "../envs/mmseqs.yml"
    threads:
        config["mmseqs2"].get("threads", 4)
    log:
        createdb = "log/mmseqs_createdb.stderr",
        cluster = "log/mmseqs_cluster.stderr",
        createtsv = "log/mmseqs_createtsv.stderr"

    params:
        db = "proteome",
        db_path = "results/0_mmseqs/db/",
        cluster_db = "clustering",
        cluster_db_path = "results/0_mmseqs/results/",
        min_seq_id = config["mmseqs2"].get("min_seq_id", 0),
        cov_mode = config["mmseqs2"].get("cov_mode", 0),
        sensitivity = config["mmseqs2"].get("sensitivity", 7.5),
        coverage = config["mmseqs2"].get("coverage", 0.5),
        evalue = config["mmseqs2"].get("evalue", 0.001),
        cluster_mode = config["mmseqs2"].get("cluster_mode", 0)
    shell:
        '''
        mkdir -p tmp ;
        mkdir -p {params.db_path} ;
        mkdir -p {params.cluster_db_path} ;

        mmseqs createdb {input} {params.db_path}/{params.db} > {log.createdb};

        mmseqs cluster {params.db_path}/{params.db} {params.cluster_db_path}/{params.cluster_db} tmp \
        --min-seq-id {params.min_seq_id} \
        --cov-mode {params.cov_mode} \
        --threads {threads} \
        -s {params.sensitivity} \
        -c {params.coverage} > {log.cluster};

        mmseqs createtsv {params.db_path}/{params.db} {params.db_path}/{params.db} {params.cluster_db_path}/{params.cluster_db} {output} > {log.createtsv}

        rm -rf tmp/*
        '''


checkpoint faa_subfamily:
    input:
        rules.mmseqs2_subfamily.output
    output:
        directory("results/0_mmseqs/subfamilies_faa"),
        "results/0_mmseqs/subfamilies_summary.tsv"
    conda:
        "../envs/mmseqs.yml"
    params:
        script = "workflow/scripts/mmseqstsv_to_faa.py",
        all_proteins = config["proteome"],
        out_faa_dir = "results/0_mmseqs/subfamilies_faa"
    shell:
        '''
        mkdir -p {output[0]}
        python {params.script} \
        -t {input} \
        -f {params.all_proteins} \
        -o {params.out_faa_dir} \
        -s {output[1]}
        '''


rule mafft_subfamily:
    input:
        "results/0_mmseqs/subfamilies_faa/{subfam}.faa"
    output:
        "results/0_mmseqs/subfamilies_msa/{subfam}.mafft"
    conda:
        "../envs/mmseqs.yml"
    threads:
        config["mafft"].get("threads", 4)
    shell:
        "mafft-fftnsi --quiet --maxiterate 200 --thread {threads} {input} > {output}"
