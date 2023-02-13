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
    requiredArgs.add_argument('-s', '--out_summary',
                               dest='out_summary',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def main():

    args = parse_args()

    lines = [line.strip().split("\t") for line in open(args.mmseqs_tsv, "r").readlines()]
    subfamilies = dict()
    for line in lines:
        if len(line) == 2:
            if line[0] not in subfamilies:
                subfamilies[line[0]] = [line[1]]
            else:
                subfamilies[line[0]] += [line[1]]

    n_prots_list = [(repr, len(prots)) for repr, prots in subfamilies.items() if len(prots) > 1]
    n_prots_sorted = sorted(n_prots_list, key=lambda x: x[1], reverse=True)
    nzfill = len(str((len(n_prots_sorted))))

    records = {record.id:record for record in SeqIO.parse(args.all_seqs_faa, "fasta")}

    n = 1
    summary = list()
    for cluster in n_prots_sorted:
        repr = cluster[0]
        cluster_id = "subfamily_" + str(n).zfill(nzfill)

        with open(f"{args.out_faa_dir}/{cluster_id}.faa", "w") as fout:
            to_write = [records[prot] for prot in subfamilies[repr]]
            SeqIO.write(to_write, fout, "fasta")

        summary.append([cluster_id, repr, str(cluster[1])])
        n += 1


    with open(args.out_summary, "w") as fout:
        fout.write("protein\tsubfamily\tn\n")
        for line in summary:
            fout.write("\t".join(line) + "\n")


if __name__ == "__main__":
    main()
