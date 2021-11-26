#!/usr/bin/env python3
import os
import operator
import argparse
import subprocess
from Bio import SeqIO
from Bio import codonalign
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

import gffutils

from alfpy.utils import seqrecords
from alfpy import word_distance
from alfpy import word_pattern
from alfpy import word_vector
from alfpy.utils import distmatrix
from Bio.Seq import Seq


def safe_div(x, y):
    if y == 0:
        return "NA"
    return x / y


def return_dups(filename):
    dup_dict = dict()
    filehandle = open(filename + "/full_table.tsv", 'r')
    for line in filehandle:
        chomp = line.rstrip("\n")
        line_list = chomp.split("\t")
        # print(line_list)
        if not line.startswith("#"):
            if line_list[1] == "Duplicated":
                if line_list[0] not in dup_dict.keys():
                    dup_dict[line_list[0]] = line_list[0:5]
                else:
                    dup_dict[line_list[0]].extend(line_list[0:5])
    return (dup_dict)

def return_full(filename):
    dup_dict = dict()
    filehandle = open(filename + "/full_table.tsv", 'r')
    for line in filehandle:
        chomp = line.rstrip("\n")
        line_list = chomp.split("\t")
        # print(line_list)
        if not line.startswith("#"):
            if line_list[0] not in dup_dict.keys():
                dup_dict[line_list[0]] = line_list[0:5]
            else:
                dup_dict[line_list[0]].extend(line_list[0:5])
    return (dup_dict)


def main(dup_dict, in_fasta, busco_dir):
    wd = os.getcwd()
    # #make a singular gff database of augustus preds
    # subprocess.call("cat "+ busco_dir + "gff/*.gff > ./tmp/all.gff",shell=True)
    # fn = gffutils.example_filename(wd+'/tmp/all.gff')
    # db = gffutils.create_db(fn, dbfn='test.db', force=True, keep_order=True, merge_strategy = 'merge', sort_attribute_values = True)

    # Load ref genome
    with open(in_fasta) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))
    # Load multicopy Busco DNA seqs
    # with open(wd+"/tmp/all.fna") as in_handle:
    #    fna_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    # prot alignment
    out_dict = dict()
    for key in dup_dict.keys():

        subprocess.call(
            "mafft --auto " + busco_dir + "busco_sequences/multi_copy_busco_sequences/" + key + ".faa > " + wd + "/tmp/tmp_aligned.faa",
            shell=True, stderr=subprocess.DEVNULL)
        subprocess.call(
            "pal2nal.pl " + wd + "/tmp/tmp_aligned.faa " + busco_dir + "busco_sequences/multi_copy_busco_sequences/" + key + ".fna -nogap -nomismatch -output fasta > ./tmp/nt_aligned.fna",
            shell=True)
        # print("pal2nal.pl " + wd +"/tmp/tmp_aligned.faa " + busco_dir + "busco_sequences/multi_copy_busco_sequences/" + key + ".fna -nogap -nomismatch -output fasta > ./tmp/nt_aligned.fna")
        with open(wd + "/tmp/nt_aligned.fna") as in_handle:
            nt_align_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

        # This snippet allows us to go through the dictionary of dups regardless of the number of dups
        # get number of dups per busco (should be mostly 2)
        recs = len(dup_dict[key]) / 5
        chunk_size = 5
        name_list = []
        seq_list = []
        len_list = []
        busco_cds_list = []
        busco_cds_len = []
        busco_gDNA = dict()
        for i in range(0, len(dup_dict[key]), chunk_size):
            chunk = dup_dict[key][i:i + chunk_size]
            tig_name = chunk[2].split(":")[0]
            name_list.append(tig_name + ":" + chunk[3] + "-" + chunk[4])
            seq_temp = ref_recs[tig_name].seq[int(chunk[3]):int(chunk[4])]
            # seq_temp = subprocess.check_output(
            #    "samtools faidx " + in_fasta + " " + tig_name + ":" + chunk[3] + "-" + chunk[4] + " | grep -v '>' ",
            #    shell=True).decode('utf-8').replace("\n", "")
            seq_list.append(seq_temp)
            len_list.append(len(seq_temp))
            busco_cds_list.append(str(nt_align_recs[tig_name + ":" + chunk[3] + "-" + chunk[4]].seq))
            busco_cds_len.append(str(len(nt_align_recs[tig_name + ":" + chunk[3] + "-" + chunk[4]].seq)))
            busco_gDNA[tig_name + ":" + chunk[3] + "-" + chunk[4]] = str(seq_temp)
        bp = word_pattern.create(busco_cds_list, word_size=2)

        bcounts = word_vector.Counts(busco_cds_len, bp)
        bdist = word_distance.Distance(bcounts, 'canberra')
        # print(dist)
        bmatrix = distmatrix.create(name_list, bdist)

        # BUSCO whole locus DNA canberra dist
        p = word_pattern.create(seq_list, word_size=2)

        counts = word_vector.Counts(len_list, p)
        dist = word_distance.Distance(counts, 'canberra')
        # print(dist)
        matrix = distmatrix.create(name_list, dist)

        # for i, j in zip(matrix, pr_matrix):
        #    print(key, int(recs), i[2:], j[2:])

        # prot alignment

        for i, j in zip(matrix, bmatrix):
            # CodonSeq("AAATTTGGGCCAAATTT", rf_table=(0,3,6,8,11,14))
            cdn1 = codonalign.codonseq.CodonSeq(str(nt_align_recs[i[2]].seq))
            cdn2 = codonalign.codonseq.CodonSeq(str(nt_align_recs[i[3]].seq))

            gDNA1 = busco_gDNA[i[2]]
            gDNA2 = busco_gDNA[i[3]]

            dn_ds = codonalign.codonseq.cal_dn_ds(cdn1, cdn2, method='NG86')
            if dn_ds[1] == -1:
                ds = 0
            else:
                ds = dn_ds[1]

            out_dict[key] = [int(recs), abs(len(busco_gDNA[i[2]]) - len(busco_gDNA[i[3]])), *i[2:], j[4], dn_ds[0], ds,
                  safe_div(dn_ds[0], ds)]
            #print(key, int(recs), abs(len(busco_gDNA[i[2]]) - len(busco_gDNA[i[3]])), *i[2:], j[4], dn_ds[0], ds,
            #      safe_div(dn_ds[0], ds), sep=",")

            # print(key, len_list, int(recs), *i[2:], j[4] , *dn_ds, dn_ds[0]/dn_ds[1], j[4]/i[4] , sep ="\t")
            # print(key, int(recs), i[2:])
    return out_dict


if __name__ == '__main__':
    my_parser = argparse.ArgumentParser(description='Return a bed file of motifs from a sequence and fasta file')

    my_parser.add_argument('-t', '--threads', action='store', help='threads', type=str)
    my_parser.add_argument('-f', '--input_fasta', action='store', help='Input fasta file', type=str, required=True)
    my_parser.add_argument('-i', '--input_busco_dir', action='store', help='Input busco dir', type=str, required=True)
    my_parser.add_argument('-p', '--input_purged_busco_dir', action='store', help='Input busco dir for purged genome', type=str, required=True)

    # Execute parse_args()

    subprocess.call("mkdir tmp", shell=True)
    args = my_parser.parse_args()


    dup_dict = return_dups(args.input_busco_dir)
    full_purged_dict = return_full(args.input_purged_busco_dir)
   

    final_dict = main(dup_dict, args.input_fasta, args.input_busco_dir)

    for key in final_dict:
        print(key, *final_dict[key],*full_purged_dict[key],sep=",")
        
