import glob, os, argparse
from Bio import SeqIO
import pandas as pd
import glob
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-p', '--proteome',
                               dest='proteome',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-o', '--out_file',
                               dest='out_file',
                               required=True,
                               help=''
                               )


    return parser.parse_args()



def main():

    args = parse_args()

    proteome = [record for record in SeqIO.parse(args.proteome, "fasta")]

    to_df = list()
    for record in proteome:
        id = record.id
        description = record.description.split(id, maxsplit=1)[1].strip()
        to_df.append([id, description, "", ""])


    df = pd.DataFrame(to_df, columns=["protein", "description", "subfamily", "family"])
    df.set_index("protein", inplace=True)


    subfamilies = glob.glob("results/0_mmseqs/subfamilies_faa/*.faa")
    for faa in subfamilies:
        subfamily = os.path.basename(faa).replace(".faa", "")
        records = [record.id for record in SeqIO.parse(faa, "fasta")]
        for record in records:
            df.loc[record, "subfamily"] = subfamily

    families = glob.glob("results/2_mcl/families_faa/*.faa")
    for faa in families:
        family = os.path.basename(faa).replace(".faa", "")
        records = [record.id for record in SeqIO.parse(faa, "fasta")]
        for record in records:
            df.loc[record, "family"] = family


    df.to_csv(args.out_file, sep="\t", header=True, index=True)


if __name__ == "__main__":
    main()
