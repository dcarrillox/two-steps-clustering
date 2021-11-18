configfile: "config/config.yaml"


def gather_msas(wildcards):
    checkpoint_output = checkpoints.create_faa_clusters.get(**wildcards).output[0]

    print(checkpoint_output)
    results_dir = config["results_dir"]

    clusters_msas = expand(results_dir + "/0_clustering/clusters_msa/{cluster}.mafft",
                    cluster=glob_wildcards(f"{checkpoint_output}/{{cluster}}.faa").cluster,
                    )
    return clusters_msas

def aggregate_final(wildcards):
    checkpoint_output = checkpoints.create_faa_clusters.get(**wildcards).output[0]
    print(checkpoint_output)

    results_dir = config["results_dir"]

    clusters_hhr = expand(results_dir + "/1_hhblits/results/{cluster}.hhr",
                    cluster=glob_wildcards(f"{checkpoint_output}/{{cluster}}.faa").cluster,
                    )
    return clusters_hhr
