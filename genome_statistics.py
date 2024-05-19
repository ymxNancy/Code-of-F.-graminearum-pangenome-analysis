#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
server: ipp

workdic: /home/zhanghao/graminearum/scripts

command line: screen -L ./genome_statistics.py all_isolates.txt

qc_result: /home/zhanghao/graminearum/qc_data/qc_statistics/180201.json

after spades: /home/zhanghao/graminearum/spades_analyses/spades_genomes/180201_spades.fasta

after decontamination: /home/zhanghao/graminearum/sidr_analyses/sidr_genomes/180201_sidr.fasta

after filtering contigs: /home/zhanghao/graminearum/filter_contigs/filtered_genomes/180201_filtered.fasta

after reorder contigs: /home/zhanghao/graminearum/order_contigs/ordered_genomes/180201_ordered.fas

'''


import subprocess, re
from sys import argv


def parse_strain_name(strain_name_file):
    strain_name_list = []
    with open(strain_name_file, 'r') as f:
        for line in f:
            line = line.strip()
            strain_name_list.append(line)
        return strain_name_list

def genome_statistics(strain_name_list):
    for strain_name in strain_name_list:
        cmd_spades = 'fasta_assembly_statistics /home/zhanghao/graminearum/spades_analyses/spades_genomes/{0}_spades.fasta > /home/zhanghao/graminearum/genomes_statistics/after_spades/{0}_spades.txt'.format(strain_name)
        subprocess.check_call(cmd_spades, shell=True)

        cmd_decontamination = 'fasta_assembly_statistics /home/zhanghao/graminearum/sidr_analyses/sidr_genomes/{0}_sidr.fasta > /home/zhanghao/graminearum/genomes_statistics/after_decontamination/{0}_decontamination.txt'.format(strain_name)
        subprocess.check_call(cmd_decontamination, shell=True)

        cmd_filtering = 'fasta_assembly_statistics /home/zhanghao/graminearum/filter_contigs/filtered_genomes/{0}_filtered.fasta > /home/zhanghao/graminearum/genomes_statistics/after_filtering/{0}_filtering.txt'.format(strain_name)
        subprocess.check_call(cmd_filtering, shell=True)
        print(strain_name + ': 3 fasta_assembly_statistics run successfully')


def parse_qc(strain_name):
    qc_list = []
    with open('/home/zhanghao/graminearum/qc_data/qc_statistics/' + strain_name + '.json', 'r') as f:
        n = 0
        for line in f:
            n += 1 
            line = line.strip()
            if n < 26:
                if line.startswith('"total_reads"'):
                    re_total_reads = re.match(r'"total_reads":(.+),', line)
                    total_reads = re_total_reads.group(1)
                    qc_list.append(total_reads)
                elif line.startswith('"total_bases"'):
                    re_total_bases = re.match(r'"total_bases":(.+),', line)
                    total_bases = re_total_bases.group(1)
                    qc_list.append(total_bases)
                elif line.startswith('"q20_rate"'):
                    re_q20_rate = re.match(r'"q20_rate":(.+),', line)
                    q20_rate = re_q20_rate.group(1)
                    qc_list.append(format((float(q20_rate)), '.2f'))
                elif line.startswith('"q30_rate"'):
                    re_q30_rate = re.match(r'"q30_rate":(.+),', line)
                    q30_rate = re_q30_rate.group(1)
                    qc_list.append(format(float(q30_rate), '.2f'))
                elif line.startswith('"gc_content"'):
                    re_gc_content = re.match(r'"gc_content":(.+)', line)
                    gc_content = re_gc_content.group(1)
                    qc_list.append(format(float(gc_content), '.2f'))
        return qc_list

def parse_spades(strain_name):
    spades_list = []
    with open('/home/zhanghao/graminearum/genomes_statistics/after_spades/' + strain_name + '_spades.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Number of contigs'):
                line = line.split(':\t')
                spades_list.append(line[1])
            elif line.startswith('Total size'):
                line = line.split(':\t')
                spades_list.append(line[1])
            elif line.startswith('Mean contig size'):
                line = line.split(':\t')
                spades_list.append(format(float(line[1]), '.2f'))
            elif line.startswith('Shortest contig'):
                line = line.split(':\t')
                spades_list.append(line[1])
        return spades_list

def parse_decontamination(strain_name):
    decontamination_list = []
    with open('/home/zhanghao/graminearum/genomes_statistics/after_decontamination/' + strain_name + '_decontamination.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Number of contigs'):
                line = line.split(':\t')
                decontamination_list.append(line[1])
            elif line.startswith('Total size'):
                line = line.split(':\t')
                decontamination_list.append(line[1])
            elif line.startswith('Mean contig size'):
                line = line.split(':\t')
                decontamination_list.append(format(float(line[1]), '.2f'))
            elif line.startswith('Shortest contig'):
                line = line.split(':\t')
                decontamination_list.append(line[1])
        return decontamination_list

def parse_filtering(strain_name):
    filtering_list = []
    with open('/home/zhanghao/graminearum/genomes_statistics/after_filtering/' + strain_name + '_filtering.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('Number of contigs'):
                line = line.split(':\t')
                filtering_list.append(line[1])
            elif line.startswith('Total size'):
                line = line.split(':\t')
                filtering_list.append(line[1])
            elif line.startswith('Mean contig size'):
                line = line.split(':\t')
                filtering_list.append(format(float(line[1]), '.2f'))
            elif line.startswith('Shortest contig'):
                line = line.split(':\t')
                filtering_list.append(line[1])
        return filtering_list

def main(strain_name_list):

    with open('/home/zhanghao/graminearum/genomes_statistics/genome_statistics_result_271.csv', 'w') as f:
        f.write('strain_name\tRaw_Total_reads\tRaw_Total_bases\tRaw_Q20_bases\tRaw_Q30_bases\tRaw_GC_content\t-\tQC_Total_reads\tQC_Total_bases\tQC_Q20_bases\tQC_Q30_bases\tQC_GC_content\t-\tSpades_Number_of_contigs\tSpades_Total_size(bp)\tSpades_Mean_contig_size(bp)\tSpades_Shortest_contig(bp)\t-\tDecon_Number_of_contigs\tDecon_Total_size(bp)\tDecon_Mean_contig_size(bp)\tDecon_Shortest_contig(bp)\t-\tFilter_Number_of_contigs\tFilter_Total_size(bp)\tFilter_Mean_contig_size(bp)\tFilter_Shortest_contig(bp)\t-\tRemoved_length_percentage(%)\n')
        for strain_name in strain_name_list:
            qc = parse_qc(strain_name)
            spades = parse_spades(strain_name)
            decontamination = parse_decontamination(strain_name)
            filtering = parse_filtering(strain_name)
            removed_percentage = format((int(decontamination[1]) - int(filtering[1]))/int(decontamination[1]) * 100, '.2f')
            f.write(strain_name + '\t' + qc[0] + '\t' + qc[1] + '\t' + qc[2] + '\t' + qc[3] + '\t' + qc[4] + '\t-\t' + qc[5] + '\t' + qc[6] + '\t' + qc[7] + '\t' + qc[8] + '\t' + qc[9] + '\t-\t' + spades[0] + '\t' + spades[1] + '\t' + spades[2] + '\t' + spades[3] + '\t-\t' + decontamination[0] + '\t' + decontamination[1] + '\t' + decontamination[2] + '\t' + decontamination[3] + '\t-\t' + filtering[0] + '\t' + filtering[1] + '\t' + filtering[2] + '\t' + filtering[3] + '\t-\t' + removed_percentage + '\n')

if __name__ == '__main__':
    strain_names_file = argv[1]
    getted_strain_name_list = parse_strain_name(strain_names_file)
    assembly_result = genome_statistics(getted_strain_name_list)

    main(getted_strain_name_list)

# screen -L ./genome_statistics.py all_isolates.txt