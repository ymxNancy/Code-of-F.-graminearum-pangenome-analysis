#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os,subprocess
from sys import argv
from multiprocessing import Pool

def get_genomes_name(genome_name_file):
    name_list = []
    with open(genome_name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

def extract_matching_contig_name(strain_name):
    cmd_1 = 'mkdir /home/zhanghao/graminearum/filter_much_contigs/filter_much_contigs_tmp/{0}'.format(strain_name)
    subprocess.check_call(cmd_1, shell=True)
    os.chdir('/home/zhanghao/graminearum/filter_much_contigs/filter_much_contigs_tmp/{0}'.format(strain_name))
    cmd_2 = "ln -s /home/zhanghao/graminearum/order_contigs/ordered_genomes/{0}_ordered.fas ./ && nucmer /home/zhanghao/graminearum/rothamsted_research/PH1_nomt_renamed_ordered.fasta {0}_ordered.fas -p ph1_{0} && show-coords -d -c -r ph1_{0}.delta > ph1_{0}.delta.matching && cut -f 2 ph1_{0}.delta.matching | tail -n +6 | sort -n | uniq > filtered_contigs_name.txt && cp filtered_contigs_name.txt filtered_contigs_name_adjust_contigs.txt".format(strain_name)
    subprocess.check_call(cmd_2, shell=True)
    os.chdir('/home/zhanghao/graminearum/scripts')

# extract_extracted_contigs.sh
## After then, manually check each name file, based on the .delta file to decide the left contigs (mummerplot ph1_SRR5948960.delta -p ph1_SRR5948960)

def extract_left_contigs(strain_name):
    os.chdir('/home/zhanghao/graminearum/filter_much_contigs/filter_much_contigs_tmp/{0}'.format(strain_name))
    cmd_1 = "tr '\n' '\$' < filtered_contigs_name_adjust_contigs.txt > filtered_contigs_name_adjust_contigs_wrap.txt && sed -i 's/\$/\$ /g' filtered_contigs_name_adjust_contigs_wrap.txt"
    subprocess.check_call(cmd_1, shell=True)
    cmd_2 = 'echo -e "fasta_extract $(cat filtered_contigs_name_adjust_contigs_wrap.txt) {0}_ordered.fas > {0}_left_contigs.fas" > extract_extracted_contigs.sh && chmod +x extract_extracted_contigs.sh && ./extract_extracted_contigs.sh'.format(strain_name)
    subprocess.check_call(cmd_2, shell=True)
    os.chdir('/home/zhanghao/graminearum/scripts')


if __name__ == '__main__':
    pool = Pool(10)
    # Divide into 2 steps:
    # 1st 
    #genome_name_file = '/home/zhanghao/graminearum/scripts/much_contigs_name.txt'
    #getted_genomes_name = get_genomes_name(genome_name_file)
    #pool.map(extract_matching_contig_name, getted_genomes_name)

    # 2nd
    keep_strains_name = '/home/zhanghao/graminearum/scripts/much_contigs_keep_strains.txt'
    getted_keep_contigs_name = get_genomes_name(keep_strains_name)
    pool.map(extract_left_contigs, getted_keep_contigs_name)



