#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess, re
from multiprocessing import Pool

def parse_strain_names(name_file):
    name_list = []
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list


def generate_protein_seqs(strain_name):
    if strain_name == 'PH1':
        cmd = 'cp /home/zhanghao/graminearum/rothamsted_research/Fusarium_graminearum.RR1.pep.all.fa /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation_proteins/PH1_proteins.fas'
        subprocess.check_call(cmd, shell=True)
    else:
        cmd = 'gffread /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation/{0}_renamed_EVM.gff3 -g /home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas -y /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation_proteins/{0}_proteins.fas'.format(strain_name)
        subprocess.check_call(cmd, shell=True)

if __name__ == '__main__':
    name_file = '/home/zhanghao/graminearum/scripts/isolates_244.txt'

    got_name_list = parse_strain_names(name_file)

    pool = Pool(125)
    pool.map(generate_protein_seqs, got_name_list)



