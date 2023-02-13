import pandas as pd
from Bio import SeqIO
import os, glob
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-m', '--mcl',
                               dest='mcl',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-o', '--out_faa_dir',
                               dest='out_faa_dir',
                               required=True,
                               help=''
                               )

    return parser.parse_args()



def main():

    args = parse_args()

    # store subfamilies secuencies in a dict
    subfs_faas = glob.glob("results/0_mmseqs/subfamilies_faa/*.faa")
    subf_seqs = {os.path.basename(subf).replace(".faa", ""):[record for record in SeqIO.parse(subf, "fasta")] for subf in subfs_faas}


    # read MCL clustering
    lines =[line.strip().split() for line in open(args.mcl).readlines()]
    nfill = len(str(len(lines)))


    # store number of prots in a family and their index line in the file
    fams_len = list()

    for i, fam in enumerate(lines):
        nprots = 0
        for subfam in fam:
            nprots += len(subf_seqs[subfam])

        fams_len.append([nprots, i])

    sorted_fams_len = sorted(fams_len, key=lambda x: x[0], reverse=True)


    summary = list()

    # iterate the sorted by index list
    n = 1
    for fam in sorted_fams_len:
        family = f"family_{str(n).zfill(nfill)}"
        n += 1

        to_write = list()
        for subfam in lines[fam[1]]:
            to_write += subf_seqs[subfam]

        with open(f"{args.out_faa_dir}/{family}.faa", "w") as fout:
            SeqIO.write(to_write, fout, "fasta")

        summary.append([family, str(len(to_write))])

    with open("results/2_mcl/families_summary.tsv", "w") as fout:
        fout.write("protein\tfamily\tn\n")
        for line in summary:
            fout.write("\t".join(line) + "\n")


if __name__ == "__main__":
    main()
