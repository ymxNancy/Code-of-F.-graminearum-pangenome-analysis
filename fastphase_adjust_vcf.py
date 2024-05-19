#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


def parse_strain_name(strain_name_file):
    with open(strain_name_file, 'r') as f:
        strain_dir = {}
        for line in f:
            line = line.strip().split()
            strain = line[0]
            strain_dir[strain] = line[1]
        return strain_dir

def adjust_name(strain_dir, original_vcf_file, renamed_vcf_file):
    with open(original_vcf_file, 'r') as f1, open(renamed_vcf_file, 'w') as f2:
        n = 0
        for line in f1:
            n += 1
            if n != 38:
                line = line.strip()
                f2.write(line + '\n')
            else:
                line = line.strip().split()
                for strain, combined in strain_dir.items():
                    line = [combined if i == strain else i for i in line]
                line = '\t'.join(line)
                f2.write(line + '\n')

def get_strain_order(renamed_vcf_file):
    with open(renamed_vcf_file, 'r') as f:
        n = 0
        for line in f:
            n += 1
            line = line.strip().split()
            if n == 38:
                header = list(range(0,9))
                strains_info = line[9:]
                order_strain = sorted(enumerate(strains_info), key=lambda x:x[1])
                ordered_list = header + [(strains_info[0] + 9) for strains_info in order_strain]
        return ordered_list

def adjust_order(renamed_vcf_file, reordered_vcf_file, ordered_list):
    with open(renamed_vcf_file, 'r') as f1, open(reordered_vcf_file, 'w') as f2:
        n = 0
        for line in f1:
            n += 1
            line = line.strip()
            if n < 38:
                f2.write(line + '\n')
            else:
                line = line.strip().split()

                new_line_list = []
                
                for position in ordered_list:
                    new_line_list.append(line[position])
                new_line = '\t'.join(new_line_list)
                f2.write(new_line + '\n')


if __name__ == '__main__':
    strain_name_file = '/home/zhanghao/graminearum/call_snp_289/output/fg_pop6_210/isolates_pop_match_name.txt'
    original_vcf_file = '/home/zhanghao/graminearum/call_snp_289/output/fg_pop6_210/fg_pop6_210_com_snp_filtered_gatk_vcfR_noasterisk.vcf'
    renamed_vcf_file = '/home/zhanghao/graminearum/call_snp_289/output/fg_pop6_210/fg_pop6_210_com_snp_filtered_gatk_vcfR_noasterisk_renamed.vcf'
    reordered_vcf_file = '/home/zhanghao/graminearum/call_snp_289/output/fg_pop6_210/fg_pop6_210_com_snp_filtered_gatk_vcfR_noasterisk_renamed_reordered.vcf'

    strain_name_dir = parse_strain_name(strain_name_file)
    getted_renamed_vcf = adjust_name(strain_name_dir, original_vcf_file, renamed_vcf_file)
    getted_strain_order = get_strain_order(renamed_vcf_file)
    getted_ordered_vcf = adjust_order(renamed_vcf_file, reordered_vcf_file, getted_strain_order)

