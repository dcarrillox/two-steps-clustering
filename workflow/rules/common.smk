configfile: "config/config.yaml"


def gather_msas(wildcards):
    checkpoint_output = checkpoints.faa_subfamily.get(**wildcards).output[0]

    faa_subfamilies = expand("results/0_mmseqs/subfamilies_msa/{subfam}.mafft",
                    subfam=glob_wildcards(f"{checkpoint_output}/{{subfam}}.faa").subfam,
                    )
    return faa_subfamilies

def gather_hhblits(wildcards):
    checkpoint_output = checkpoints.faa_subfamily.get(**wildcards).output[0]

    hhblits_family = expand("results/1_hhblits/results/{subfam}.hhr",
                    subfam=glob_wildcards(f"{checkpoint_output}/{{subfam}}.faa").subfam,
                    )
    return hhblits_family
