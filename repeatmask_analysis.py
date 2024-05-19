#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from sys import argv
from multiprocessing import Pool

def get_genomes_name(genome_name_file):
    name_list = []
    with open(genome_name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

def conduct_repeat_analysis(strain_name):
    cmd = 'mkdir /home/zhanghao/graminearum/repeatmask_analyses/repeatmask_tmp_files/{0} && cp /home/zhanghao/graminearum/order_contigs/ordered_genomes/{0}_ordered.fas /home/zhanghao/graminearum/repeatmask_analyses/repeatmask_tmp_files/{0}/ && cd /home/zhanghao/graminearum/repeatmask_analyses/repeatmask_tmp_files/{0} && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 BuildDatabase -name {0} -engine ncbi {0}_ordered.fas && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 RepeatModeler -database {0} -pa 2 -LTRStruct && docker run -v $(pwd):/work --workdir /work --rm -it --user $(id -u):$(id -u) dfam/tetools:1.2 RepeatMasker {0}_ordered.fas -pa 2 -e ncbi -lib RM*/consensi.fa.classified -dir result_dir -xsmall -gff -html && cp result_dir/*.masked /home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas'.format(strain_name)

    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
        print(strain_name + ' run successfully')
    else:
        print(strain_name + ' run unsuccessfully')


def main(strain_name):
    repeat_analysis = conduct_repeat_analysis(strain_name)

if __name__ == '__main__':
    genome_name_file = argv[1]
    getted_genomes_name = get_genomes_name(genome_name_file)
    
    pool = Pool(62)
    pool.map(main, getted_genomes_name)

    # command line: screen -L ./repeat_analysis.py all_isolates.txt
