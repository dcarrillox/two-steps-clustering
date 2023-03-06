import glob, os
import os,sys,re
from collections import defaultdict
import argparse
import logging
from datetime import date, time, datetime
import shutil
import json
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description='')

    requiredArgs = parser.add_argument_group("Required Arguments")

    requiredArgs.add_argument('-f', '--subfam_faa_dir',
                               dest='subfam_faa_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-hh', '--hhblits_dir',
                               dest='hhblits_dir',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-c', '--coverage',
                               dest='coverage',
                               required=True,
                               help=''
                               )
    requiredArgs.add_argument('-p', '--probability',
                               dest='probability',
                               required=True,
                               help=''
                               )

    return parser.parse_args()

def readingHhrFile(hhr_filename, coverage_threshold, probs_threshold, seqs_subfam) :
    hits = list()
    summary = list()
    lines = open(hhr_filename,'r').readlines()
    query = seqs_subfam[lines[0].split()[1].split(" ")[0]]
    qlen  = int(lines[1].strip().split()[1])

    # find last line of hits
    for i, line in enumerate(lines):
        if line.startswith("No 1\n"):
            last = i -1

    # get target complete ids
    target_ids = dict()
    for i in range(last, len(lines)-1):
        if lines[i].startswith("No ") and lines[i+1].startswith(">"):
            n_id = lines[i].strip().split(" ")[1]
            complete_id = lines[i+1].strip().replace(">", "")
            target_ids[n_id] = complete_id

    for line in lines[9:last]:
        fields = line.strip().split()
        target = seqs_subfam[target_ids[fields[0]]]
        if "-" not in fields[-1]:
            values = fields[-9:]
            values[-1] = values[-1].replace("(", "").replace(")", "")
        else:
            values = fields[-8:]
            toadd = values[-1].split("(")[1].replace(")", "")
            values[-1] = values[-1].split("(")[0]
            values.append(toadd)

        values.insert(0, target)

        prob = float(values[1])
        qstart = int(values[7].split("-")[0])
        qend   = int(values[7].split("-")[1])
        qcover = (qend - qstart + 1.0) / qlen
        slen   = int(values[-1])
        sstart = int(values[8].split("-")[0])
        send   = int(values[8].split("-")[1])
        scover = (send - sstart + 1.0) / slen

        if prob > probs_threshold and (scover > coverage_threshold or qcover > coverage_threshold) :
            if query != target : # removing self hit
                if query < target :
                    hits.append([query,target,prob,qcover,scover])
                else:
                    hits.append([target,query,prob,scover,qcover])

        summary.append([query, target, str(prob), str(qcover), str(scover), str(qlen), str(slen)])
    return hits, summary

def main():

    args = parse_args()

    subfam_faas = glob.glob(f"{args.subfam_faa_dir}/*.faa")
    seqs_subfam = {record.id:os.path.basename(file).replace(".faa", "") for file in subfam_faas for record in SeqIO.parse(file, "fasta") }

    hhr_files = sorted(glob.glob(f"{args.hhblits_dir}/*.hhr"))
    edges_weights = dict()
    summary_results = list()
    for file in hhr_files:
        hits, summary = readingHhrFile(file, float(args.coverage), float(args.probability), seqs_subfam)
        summary_results += summary
        for hit in hits:
            query  = hit[0]
            target = hit[1]
            prob = float(hit[2]) / 100
            qcov = float(hit[3])
            scov = float(hit[4])
            weight = prob * max([qcov, scov])
            edge = f"{query}\t{target}"
            if edge in edges_weights:
                if weight > edges_weights[edge]:
                    edges_weights[edge] = weight
            else:
                edges_weights[edge] = weight


    with open("results/2_mcl/mcl_subfamilies_network.txt", "w") as fout:
        for edge, weight in edges_weights.items():
            fout.write(f"{edge}\t{weight}\n")

    with open("results/1_hhblits/parsed_results.tsv", "w") as fout:
        header = ["query", "target", "prob", "qcover", "scover", "qlen", "slen"]
        fout.write("\t".join(header) + "\n")
        for line in summary_results:
            fout.write("\t".join(line) + "\n")

if __name__ == "__main__":
    main()
