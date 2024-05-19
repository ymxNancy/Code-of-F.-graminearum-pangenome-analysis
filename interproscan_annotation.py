#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from multiprocessing import Pool, pool
import subprocess
from sys import argv

def get_name_list(name_file):
    with open(name_file, 'r') as f:
        name_list = []
        for line in f:
            line = line.strip()
            name_list.append(line)
        return name_list

def parse_interproscan(strain_num):
    cmd = 'interproscan.sh -i /home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation_proteins/{0}_proteins.fas -f TSV,GFF3,XML -cpu 8 -goterms -iprlookup -d /home/zhanghao/graminearum/functional_annotation/interproscan_annotation'.format(strain_num)
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
        print(strain_num + ' run successfully.')
    else:
        print(strain_num + ' run unsuccessfully.')
    
if __name__ == '__main__':
    name_file = '/home/zhanghao/graminearum/scripts/isolates_122_1_2.txt'
    #name_file = '/home/zhanghao/graminearum/scripts/isolates_test.txt'
    #name_file = '/home/zhanghao/graminearum/scripts/test.txt'
    getted_name_list = get_name_list(name_file)
    pool = Pool(8)
    pool.map(parse_interproscan, getted_name_list)


# command line: ./interproscan_annotation.py isolates_244.txt











