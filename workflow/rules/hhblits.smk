configfile: "config/config.yaml"

results_dir = config["results_dir"]

rule build_hhdb:
    input:
        gather_msas
    output:
        results_dir + "/1_hhblits/db/hhblits_db.done"
    params:
        db = "hhblits_db",
        db_dir = results_dir + "/1_hhblits/db",
        msa_dir = results_dir + "/0_clustering/clusters_msa",
        outfile = "hhblits_db.done"
    threads: 2
    conda:
        "../envs/env.yaml"
    shell:
        '''
        mkdir -p {params.db_dir} ;
        ffindex_build -s {params.db_dir}/{params.db}.ffdata {params.db_dir}/{params.db}.ffindex {params.msa_dir} ;
        cd {params.db_dir} ;
        mpirun -np {threads} ffindex_apply {params.db}.ffdata {params.db}.ffindex -i {params.db}_a3m_wo_ss.ffindex -d {params.db}_a3m_wo_ss.ffdata -- hhconsensus -M 40 -maxres 65535 -i stdin -oa3m stdout -v 0 ;
        rm -f {params.db}.ff* ;
        mv {params.db}_a3m_wo_ss.ffdata {params.db}_a3m.ffdata ;
        mv {params.db}_a3m_wo_ss.ffindex {params.db}_a3m.ffindex ;
        mpirun -np {threads} ffindex_apply {params.db}_a3m.ffdata {params.db}_a3m.ffindex -i {params.db}_hhm.ffindex -d {params.db}_hhm.ffdata -- hhmake -i stdin -o stdout -v 0 ;
        mpirun -np {threads} cstranslate -f -x 0.3 -c 4 -I a3m -i {params.db}_a3m -o {params.db}_cs219 > cfstranslate.log 2>&1 ;
        sort -k3 -n -r {params.db}_cs219.ffindex | cut -f1 > sorting.dat ;
        ffindex_order sorting.dat {params.db}_hhm.ffdata {params.db}_hhm.ffindex {params.db}_hhm_ordered.ffdata {params.db}_hhm_ordered.ffindex ;
        ffindex_order sorting.dat {params.db}_a3m.ffdata {params.db}_a3m.ffindex {params.db}_a3m_ordered.ffdata {params.db}_a3m_ordered.ffindex ;
        mv {params.db}_a3m_ordered.ffdata {params.db}_a3m.ffdata ;
        mv {params.db}_a3m_ordered.ffindex {params.db}_a3m.ffindex ;
        mv {params.db}_hhm_ordered.ffdata {params.db}_hhm.ffdata ;
        mv {params.db}_hhm_ordered.ffindex {params.db}_hhm.ffindex ;

        touch {params.outfile}
        '''

rule run_hhblits:
    input:
        rules.build_hhdb.output,
        cluster_msa = results_dir + "/0_clustering/clusters_msa/{cluster}.mafft"
    output:
        results_dir + "/1_hhblits/results/{cluster}.hhr"
    params:
        db = results_dir + "/1_hhblits/db/hhblits_db"
    threads: 30
    conda:
        "../envs/env.yaml"
    shell:
        '''
        hhblits -i {input.cluster_msa} -M 50 -o {output} -cpu {threads} \
        -d {params.db} -v 0 -aliw 6000 -z 300 -Z 300 -b 300 -B 300
        '''


rule aggregate_results:
    input:
        aggregate_final
    output:
        results_dir + "/done"
    shell:
        "touch {output}"
