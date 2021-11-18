import glob, os, argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-t', '--mmseqs_tsv',
                               dest='mmseqs_tsv',
                               required=True,
                               help='tsv from "createtsv"'
                               )
    requiredArgs.add_argument('-f', '--all_seqs_faa',
                               dest='all_seqs_faa',
                               required=True,
                               help='faa file that underwent the clustering'
                               )

    requiredArgs.add_argument('-o', '--out_faa_dir',
                               dest='out_faa_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def main():

    args = parse_args()

    # keys are the representative sequences
    repr_cluster = {line.strip().split("\t")[0]:list() for line in open(args.mmseqs_tsv).readlines()}

    # count how many proteins are associated with the representative sequences
    # cluster identifiers will be sorted based on this
    lines = [line.strip().split("\t") for line in open(args.mmseqs_tsv).readlines()]
    for line in lines:
        repr_cluster[line[0]] += [line[1]]

    sorted_repr_cluster = {repr:seqs for repr, seqs in sorted(repr_cluster.items(), key=lambda item: len(item[1]), reverse=True)}


    # read input faa file with all the sequences
    records = {record.id:record for record in SeqIO.parse(args.all_seqs_faa, "fasta")}

    # write clusters to .faa file. Only clusters with at least 2 sequences.
    # on the way, store the pairs that will be in the two columns table later
    to_write_table = list()
    cluster_n = 1
    for reprs, seqs in sorted_repr_cluster.items():
        if len(seqs) > 1:
            # create the cluster id
            n = str(cluster_n).zfill(4)
            cl_id = f"cl_{n}"

            to_write_table += [[cl_id, seq] for seq in seqs]

            to_write = [records[seq] for seq in seqs]
            with open(f"{args.out_faa_dir}/{cl_id}.faa", "w") as fout:
                SeqIO.write(to_write, fout, "fasta")

            cluster_n += 1



    # write the two columns table:  [0] cl_xx, [1] protein_id
    with open(f"{args.out_faa_dir}/clusters_table.txt", "w") as fout:
        for protein in to_write_table:
            fout.write("\t".join(protein) + "\n")

if __name__ == "__main__":
    main()
