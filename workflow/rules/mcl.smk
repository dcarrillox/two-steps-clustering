rule create_mcl_network:
    input:
        gather_hhblits
    output:
        "results/1_hhblits/parsed_results.tsv",
        "results/2_mcl/mcl_subfamilies_network.txt"
    conda:
        "../envs/mcl.yml"
    params:
        script = "workflow/scripts/create_mcl_network.py",
        coverage = config["hhblits"].get("coverage", 0.5),
        probability = config["hhblits"].get("probability", 90),
        subfam_faa_dir = "results/0_mmseqs/subfamilies_faa",
        hhblits_dir = "results/1_hhblits/results"
    shell:
        '''
        python {params.script} \
        -f {params.subfam_faa_dir} \
        -hh {params.hhblits_dir}\
        -c {params.coverage} \
        -p {params.probability}
        '''

rule mcl_family:
    input:
        rules.create_mcl_network.output[1]
    output:
        "results/2_mcl/mcl_family_results.txt"
    conda:
        "../envs/mcl.yml"
    threads:
        config["mcl"].get("threads", 2)
    params:
        inflation = config["mcl"].get("inflation", 2)
    log:
        stderr = "log/mcl.stderr"
    shell:
        '''
        mcl {input} --abc \
        -I {params.inflation} \
        -te {threads} \
        -o {output} > {log.stderr}
        '''

rule faa_family:
    input:
        rules.mcl_family.output
    output:
        directory("results/2_mcl/families_faa"),
        "results/2_mcl/families_summary.tsv"
    conda:
        "../envs/mcl.yml"
    params:
        script = "workflow/scripts/create_families.py",
        out_faa_dir = "results/2_mcl/families_faa"
    shell:
        '''
        mkdir -p {params.out_faa_dir}
        python {params.script} -m {input} -o {params.out_faa_dir}
        '''


rule summarize_final_results:
    input:
        rules.faa_family.output[1],
        config["proteome"]
    output:
        "results/protein_clustering_results.tsv"
    conda:
        "../envs/mcl.yml"
    params:
        script = "workflow/scripts/summarize_final_results.py"
    shell:
        "python {params.script} -p {input[1]} -o {output}"
