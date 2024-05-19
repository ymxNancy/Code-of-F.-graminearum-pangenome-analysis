#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from multiprocessing import Pool

def parse_strain_names(name_file):
    name_list = []
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list


def perform_assembly(strain_name):
    '''
    cmd = 'spades.py -o /home/zhanghao/graminearum/spades_analyses/spades_assembly_files/{0} -t 8 -1 /home/zhanghao/graminearum/qc_data/{0}_1_qc.fastq.gz -2 /home/zhanghao/graminearum/qc_data/{0}_2_qc.fastq.gz -k 31,61,91,121'.format(strain_name)
    res = subprocess.check_call(cmd, shell=True)
    if res == 0:
        print(strain_name + ' spades run successfully')
    else:
        print(strain_name + ' spades run failed')
    '''

    cmd_2 = 'cp /home/zhanghao/graminearum/spades_analyses/spades_assembly_files/{0}/contigs.fasta /home/zhanghao/graminearum/spades_analyses/spades_genomes/{0}_spades.fasta'.format(strain_name)
    res_2 = subprocess.check_call(cmd_2, shell=True)
    if res_2 == 0:
        print(strain_name + ' moved successfully')
    else:
        print(strain_name + ' moved failed')

if __name__ == "__main__":

    # wur isolates
    #getted_name_list = parse_strain_names('wur_isolates.txt')

    # ipp isolates
    #getted_name_list = parse_strain_names('ipp_isolates.txt')

    # ncbi isolates
    getted_name_list = parse_strain_names('ncbi_isolates.txt')


    # Create a multiprocessing Pool
    pool = Pool(15)
    # proces strains_list iterable with pool
    pool.map(perform_assembly, getted_name_list)
