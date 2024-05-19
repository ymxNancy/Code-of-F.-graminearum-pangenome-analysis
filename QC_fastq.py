#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess,json
from multiprocessing import Pool

def parse_strain_names(name_file):
    name_list = []
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

# wur isolates:
def parse_qc_wur(name_list):
    for name in name_list:
        cmd = 'fastp -i /home/zhanghao/graminearum/raw_data/{0}_R1_001.fastq.gz -I /home/zhanghao/graminearum/raw_data/{0}_R2_001.fastq.gz -o /home/zhanghao/graminearum/qc_data/{0}_1_qc.fastq.gz -O /home/zhanghao/graminearum/qc_data/{0}_2_qc.fastq.gz -w 100 -j /home/zhanghao/graminearum/qc_data/qc_statistics/{0}.json -h /home/zhanghao/graminearum/qc_data/qc_statistics/{0}.html'.format(name)
        check = subprocess.check_call(cmd, shell=True)
        if check == 0:
                print(name + ' run successfully')
        else:
            print(name + ' run failed')


# IPP isolates:
def parse_qc_ipp(strain_name):
    cmd = 'fastp -i /home/zhanghao/graminearum/raw_data/{0}_1.clean.fq.gz -I /home/zhanghao/graminearum/raw_data/{0}_2.clean.fq.gz -o /home/zhanghao/graminearum/qc_data/{0}_1_qc.fastq.gz -O /home/zhanghao/graminearum/qc_data/{0}_2_qc.fastq.gz -w 16 -j /home/zhanghao/graminearum/qc_data/qc_statistics/{0}.json -h /home/zhanghao/graminearum/qc_data/qc_statistics/{0}.html'.format(strain_name)
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
            print(strain_name + ' run successfully')
    else:
        print(strain_name + ' run failed')

# NCBI isolates:
def parse_qc_ncbi(strain_name):
    cmd = 'fastp -i /home/zhanghao/graminearum/raw_data/{0}_1.fastq.gz -I /home/zhanghao/graminearum/raw_data/{0}_2.fastq.gz -o /home/zhanghao/graminearum/qc_data/{0}_1_qc.fastq.gz -O /home/zhanghao/graminearum/qc_data/{0}_2_qc.fastq.gz -w 16 -j /home/zhanghao/graminearum/qc_data/qc_statistics/{0}.json -h /home/zhanghao/graminearum/qc_data/qc_statistics/{0}.html'.format(strain_name)
    check = subprocess.check_call(cmd, shell=True)
    if check == 0:
            print(strain_name + ' run successfully')
    else:
        print(strain_name + ' run failed')


def qc_statistics(name_list, qc_statistics_file):
    with open(qc_statistics_file, 'w') as f:
        f.write('strain_name\tbefore_total_bases\tbefore_q30\tbefore_gc\tafter_total_bases\tafter_q30\tafter_gc\n')
        for name in name_list:
            json_file = '/home/zhanghao/graminearum/qc_data/qc_statistics/{0}.json'.format(name)
            with open(json_file, 'r') as f1:
                parsed_json = json.load(f1)
                before_total_bases = parsed_json['summary']['before_filtering']['total_bases']
                before_q30 = parsed_json['summary']['before_filtering']['q30_rate']
                before_gc = parsed_json['summary']['before_filtering']['gc_content']

                after_total_bases = parsed_json['summary']['after_filtering']['total_bases']
                after_q30 = parsed_json['summary']['after_filtering']['q30_rate']
                after_gc = parsed_json['summary']['after_filtering']['gc_content']
                
                f.write(str(name) + '\t' + str(before_total_bases) + '\t' + str(before_q30) + '\t' + str(before_gc) + '\t' + str(after_total_bases) + '\t' + str(after_q30) + '\t' + str(after_gc) + '\n')


if __name__ == '__main__':

    
    # wur isolates 
    '''
    wur_name_file = '/home/zhanghao/graminearum/scripts/wur_isolates.txt'
    wur_qc_statistics_file = '/home/zhanghao/graminearum/qc_data/qc_statistics/wur_qc_statistics.txt'
    wur_getted_name_list = parse_strain_names(wur_name_file)
    wur_perform_qc = parse_qc_wur(wur_getted_name_list)
    wur_got_qc_statistics = qc_statistics(wur_getted_name_list, wur_qc_statistics_file)
    

    # ipp isolates
    ipp_name_file = '/home/zhanghao/graminearum/scripts/ipp_isolates.txt'    
    ipp_qc_statistics_file = '/home/zhanghao/graminearum/qc_data/qc_statistics/ipp_qc_statistics.txt'
    ipp_getted_name_list = parse_strain_names(ipp_name_file)
    pool = Pool(8)
    pool.map(parse_qc_ipp, ipp_getted_name_list)
    ipp_got_qc_statistics = qc_statistics(ipp_getted_name_list, ipp_qc_statistics_file)
    '''

    # NCBI isolates
    ncbi_name_file = '/home/zhanghao/graminearum/scripts/ncbi_isolates.txt'    
    ncbi_qc_statistics_file = '/home/zhanghao/graminearum/qc_data/qc_statistics/ncbi_qc_statistics.txt'
    ncbi_getted_name_list = parse_strain_names(ncbi_name_file)
    pool = Pool(8)
    pool.map(parse_qc_ncbi, ncbi_getted_name_list)
    ncbi_got_qc_statistics = qc_statistics(ncbi_getted_name_list, ncbi_qc_statistics_file)

