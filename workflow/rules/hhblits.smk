rule build_subfamily_hhdb:
    input:
        gather_msas
    output:
        touch("results/1_hhblits/db/hhblits_db.done")
    conda:
        "../envs/hhsuite.yml"
    params:
        db = "hhblits_db",
        db_dir = "results/1_hhblits/db",
        msa_dir = "results/0_mmseqs/subfamilies_msa"
    threads: 2

    shell:
        '''
        mkdir -p {params.db_dir} ;
        ffindex_build -s {params.db_dir}/{params.db}.ffdata {params.db_dir}/{params.db}.ffindex {params.msa_dir} ;
        cd {params.db_dir} ;
        ffindex_apply {params.db}.ffdata {params.db}.ffindex -i {params.db}_a3m_wo_ss.ffindex -d {params.db}_a3m_wo_ss.ffdata -- hhconsensus -M 40 -maxres 65535 -i stdin -oa3m stdout -v 0 ;
        rm -f {params.db}.ff* ;
        mv {params.db}_a3m_wo_ss.ffdata {params.db}_a3m.ffdata ;
        mv {params.db}_a3m_wo_ss.ffindex {params.db}_a3m.ffindex ;
        ffindex_apply {params.db}_a3m.ffdata {params.db}_a3m.ffindex -i {params.db}_hhm.ffindex -d {params.db}_hhm.ffdata -- hhmake -i stdin -o stdout -v 0 ;
        cstranslate -f -x 0.3 -c 4 -I a3m -i {params.db}_a3m -o {params.db}_cs219 > cfstranslate.log 2>&1 ;
        sort -k3 -n -r {params.db}_cs219.ffindex | cut -f1 > sorting.dat ;
        ffindex_order sorting.dat {params.db}_hhm.ffdata {params.db}_hhm.ffindex {params.db}_hhm_ordered.ffdata {params.db}_hhm_ordered.ffindex ;
        ffindex_order sorting.dat {params.db}_a3m.ffdata {params.db}_a3m.ffindex {params.db}_a3m_ordered.ffdata {params.db}_a3m_ordered.ffindex ;
        mv {params.db}_a3m_ordered.ffdata {params.db}_a3m.ffdata ;
        mv {params.db}_a3m_ordered.ffindex {params.db}_a3m.ffindex ;
        mv {params.db}_hhm_ordered.ffdata {params.db}_hhm.ffdata ;
        mv {params.db}_hhm_ordered.ffindex {params.db}_hhm.ffindex ;
        '''

rule hhblits_family:
    input:
        rules.build_subfamily_hhdb.output,
        subfamily_msa = "results/0_mmseqs/subfamilies_msa/{subfam}.mafft"
    output:
        "results/1_hhblits/results/{subfam}.hhr"
    conda:
        "../envs/hhsuite.yml"
    params:
        db = "results/1_hhblits/db/hhblits_db"
    threads:
        config["hhblits"].get("threads", 4)
    shell:
        '''
        hhblits -i {input.subfamily_msa} -M 50 -o {output} -cpu {threads} \
        -d {params.db} -v  0 -aliw 6000 -z 300 -Z 300 -b 300 -B 300
        '''
