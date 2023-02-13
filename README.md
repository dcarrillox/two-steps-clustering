# two-steps-clustering
Snakemake pipeline to cluster proteomes in two steps with [MMseqs2](https://github.com/soedinglab/MMseqs2),
, [MAFFT](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html),
[HHsuite](https://github.com/soedinglab/hh-suite) and [MCL](http://micans.org/mcl/#:~:text=The%20MCL%20algorithm%20is%20short,in%20bioinformatics%20and%20other%20disciplines.):

1) Proteins are clustered into subfamilies in a first round using MMseq2.
2) For each subfamily, a MSA is generated using MAFTT.
3) MSA profiles are compared all-vs-all using hhblits.
4) Subfamilies are clustered using MCL based on probability and coverage obtained
with hhblits.


## Instalation

Pipeline is written to run as an Snakemake workflow. After cloning the repo,
create a conda environment `spc` containing Snakemake:

~~~
# clone the repo
$ git clone https://github.com/dcarrillox/two-steps-clustering.git
$ cd two-steps-clustering

# create minimum environment for execution
$ conda create -n spc --file=conda-linux-64.lock

# activate environment
$ conda activate spc

# test Snakemake
$ (spc) snakemake --version
6.6.1
~~~


## Execution

Pipeline's input is a FASTA file with all the proteins to be clustered. The path
to this file needs to be provided in the config file `config/config.yaml` under the
`proteome` section. Parameters for execution can be tuned in the config file as
well.

To run the pipeline:

~~~
# adjust number of threads with the -j argument
$ snakemake -j 16  --use-conda --conda-frontend mamba
~~~

Final results with the proteins clustered into subfamilies and families should be
under `results/protein_clustering_results.tsv`.
