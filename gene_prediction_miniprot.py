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


def perform_miniprot(strain_name):
    cmd = 'mkdir /home/zhanghao/graminearum/gene_prediction/homo_based_miniprot/miniprot_tmp/{0} && cd /home/zhanghao/graminearum/gene_prediction/homo_based_miniprot/miniprot_tmp/{0} && ln -s /home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas ./ && miniprot -t8 -d {0}.mpi {0}_ordered_soft.fas && miniprot -t8 -G5k --gff {0}.mpi /home/zhanghao/graminearum/gene_prediction/homo_based_miniprot/PH1_RR1_tri7_tri13.pep.fas > /home/zhanghao/graminearum/gene_prediction/homo_based_miniprot/miniprot_result/{0}_miniprot.gff'.format(strain_name)
    res = subprocess.check_call(cmd, shell=True)
    if res == 0:
        print(strain_name + ' run miniprot successfully.')
    else:
        print(strain_name + ' run miniprot failed.')


if __name__ == "__main__":
    getted_name_list = parse_strain_names('isolates_271.txt')
    # Create a multiprocessing Pool
    pool = Pool(13)
    # proces strains_list iterable with pool
    pool.map(perform_miniprot, getted_name_list)

# one strain need 1min using 8 threads.