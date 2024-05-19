#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from sys import argv
import os 


def parse_ragtag(strain_name):
    cmd = 'ragtag.py scaffold -t 40 /home/zhanghao/graminearum/rothamsted_research/PH1_nomt_renamed_ordered.fasta /home/zhanghao/graminearum/filter_much_contigs/filter_much_contigs_tmp/{0}/{0}_left_contigs.fas -o {0}_ragtag'.format(strain_name)
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
        print(strain_name + ' ragtag.py was run successfully')
    else:
        print(strain_name + 'ragtag.py was not run successfully')

def reverse_compliment(s):
    rc = s[::-1].upper().replace('A','t').replace('G','c').replace('T','a').replace('C','g').upper()
    return rc

def get_original_contig_dir(strain_name):
    original_contig_dir = {}
    with open('/home/zhanghao/graminearum/filter_much_contigs/filter_much_contigs_tmp/{0}/{0}_left_contigs.fas'.format(strain_name), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                key = line.replace('>', '')
                original_contig_dir[key] = []
            else:
                original_contig_dir[key].append(line)
        return original_contig_dir

def get_ordered_contig_file(strain_name, original_contig_dir):
    with open(strain_name + '_ragtag/ragtag.scaffold.agp', 'r') as f1, open('/home/zhanghao/graminearum/filter_much_contigs/filter_much_contigs_genomes/' + strain_name + '_ordered_left_contigs.fas', 'w') as f2:
        n = 1
        for line in f1:
            if not line.startswith('#'):
                line = line.strip().split()
                if line[4] == "W":
                    if n < 10:
                        ctg = strain_name + '_V1_ctg00' + str(n)
                    elif 10 <= n < 100:
                        ctg = strain_name + '_V1_ctg0' + str(n)
                    else:
                        ctg = strain_name + '_V1_ctg' + str(n)
                    n = n + 1

                    if line[8] == '+':
                        seq = ''.join(original_contig_dir[line[5]])
                    else:
                        seq = reverse_compliment(''.join(original_contig_dir[line[5]]))

                    f2.write('>Fgra_' + ctg + '\n' + seq + '\n')


def get_strain_name(strain_name_file):
    strain_name_list = []
    with open(strain_name_file, 'r') as f:
        for line in f:
            line = line.strip()
            strain_name_list.append(line)
        return strain_name_list

if __name__ == '__main__':
    strain_name_file = argv[1]
    getted_strain_name = get_strain_name(strain_name_file)
    for strain_name in getted_strain_name:
        os.chdir('/home/zhanghao/graminearum/filter_much_contigs/filter_much_contigs_tmp/{0}'.format(strain_name))
        getted_ragtag_result = parse_ragtag(strain_name)
        getted_original_contig_dir = get_original_contig_dir(strain_name)
        getted_ordered_contig_file = get_ordered_contig_file(strain_name, getted_original_contig_dir)


# Execution directory: /home/zhanghao/graminearum/scripts
# Execution command: screen -L ./rename_reorientation_much_contigs_keep.py much_contigs_keep_strains.txt


