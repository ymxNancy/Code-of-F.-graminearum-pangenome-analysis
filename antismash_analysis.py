#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from multiprocessing.pool import Pool
from sys import argv 

'''
conda activate antismash

'''

def parse_strain_name(strain_name_file):
    with open(strain_name_file, 'r') as f:
        name_list = []
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

def generate_no_locus_gff(strain_name):
    if strain_name == 'PH1':
        original_file = 'PH1_nomt.gff3'
        cmd = 'sed "s/;logic_name[^;]*$//" /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation/{0} > /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation_nolocus/{1}_nolocus.gff3'.format(original_file, strain_name) 
        subprocess.check_call(cmd, shell=True)

    else:
        original_file = strain_name + '_renamed_EVM.gff3'
        cmd = 'sed "s/;locus_tag[^;]*$//" /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation/{0} > /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation_nolocus/{1}_nolocus.gff3'.format(original_file, strain_name)
        subprocess.check_call(cmd, shell=True)

def parse_antismash(strain_name):
    if strain_name == 'PH1':
        fasta_file = strain_name + '_nomt.fasta'
    else:
        fasta_file = strain_name + '_ordered_soft.fas'

    subprocess.check_call("~/bin/run_antismash {0} functional_annotation/antismash_annotation --taxon fungi --cassis --fullhmmer --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --cc-mibig --cpus 2 --genefinding-gff3 /input/gene_prediction/renamed_EVM_annotation_nolocus/{1}_nolocus.gff3".format(fasta_file, strain_name), shell=True)

if __name__ == '__main__':
    name_list_file = argv[1]
    getted_name_list = parse_strain_name(name_list_file)

    pool = Pool(42)
    generated_nolocus_gff = pool.map(generate_no_locus_gff, getted_name_list)
    antismash_analysis = pool.map(parse_antismash, getted_name_list)


# work dictionary: /home/zhanghao/graminearum/scripts
# command line: screen -L ./antismash_analysis.py isolates_244.txt

