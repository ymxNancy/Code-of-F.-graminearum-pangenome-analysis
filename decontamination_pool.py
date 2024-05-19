#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess, os
from multiprocessing import Pool

'''
nt (blast database name): nt 

'''

from sys import argv

def get_strain_list(strain_name_file):
    with open(strain_name_file, 'r') as f:
        strain_name_list = []
        for line in f:
            line = line.strip()
            strain_name_list.append(line)
    return strain_name_list

def build_index(strain_name):
    '''
    input_fn: fas file of the sequence, such as '140001.fas'
    generate five files: .amb, .ann, .bwt, .pac, .sa
    '''
    subprocess.check_call('cp /home/zhanghao/graminearum/spades_analyses/spades_genomes/{0}_spades.fasta /home/zhanghao/graminearum/sidr_analyses/sidr_tmp_files/{0}/{0}_spades.fasta'.format(strain_name), shell=True)
    cmd = 'bwa-mem2 index {0}_spades.fasta'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def mem_alignment(strain_name):
    #cmd = 'bwa-mem2 mem {0}_spades.fasta /home/zhanghao/graminearum/qc_data/{0}_1_qc.fastq.gz /home/zhanghao/graminearum/qc_data/{0}_2_qc.fastq.gz -t 5 > {0}.sam'.format(strain_name)
    cmd = 'bwa-mem2 mem {0}_spades.fasta /home/zhanghao/graminearum/qc_data/{0}_1_qc.fastq.gz /home/zhanghao/graminearum/qc_data/{0}_2_qc.fastq.gz -t 95 > {0}.sam'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def transfer_bam(strain_name):
    #cmd = 'samtools view -b -S -@ 5 {0}.sam | samtools sort -@ 5 -m 8G -o {0}.sort.bam'.format(strain_name)
    cmd = 'samtools view -b -S -@ 95 {0}.sam | samtools sort -@ 25 -m 20G -o {0}.sort.bam'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def build_bam_index(strain_name):
    cmd = 'samtools index {0}.sort.bam'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def bulid_blast(strain_name):
    cmd = 'blastn -task megablast -query {0}_spades.fasta -db nt -culling_limit 5 -evalue 1e-25 -num_threads 5 -outfmt "6 std staxids" -out {0}_nt.blast'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def build_blast2(strain_name):
    '''
    cut: print selected parts of lines
    -f: select pointed colcums to print
    '''
    cmd = 'cut -f 1,13 {0}_nt.blast > {0}_nt.blast2'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def remove_contaminated(strain_name):
    cmd = 'sidr default -d /home/zhanghao/software/ncbi_db/taxdump -b {0}.sort.bam -f {0}_spades.fasta -r {0}_nt.blast2 -k {0}_tokeep.contigids -x {0}_toremove.contigids -t ascomycota'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def get_seq_dic(strain_name):
    with open(strain_name + '_spades.fasta', 'r') as f:
        total_contigs_dic = {}
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
                total_contigs_dic[key] = []
            else:
                total_contigs_dic[key].append(line)
        return total_contigs_dic

def get_keep_contigs(strain_name):
    with open('{0}_tokeep.contigids'.format(strain_name), 'r') as f:
        keep_contigs_list = []
        for line in f:
            line = line.strip()
            keep_contigs_list.append(line)
        return keep_contigs_list

def get_filtered_contig_seq(strain_name, total_contigs_dic, keep_contigs_list):
    with open('/home/zhanghao/graminearum/sidr_analyses/sidr_genomes/{0}_sidr.fasta'.format(strain_name), 'w') as f:
        for contigs_name in keep_contigs_list:
            seq = ''.join(total_contigs_dic[contigs_name])
            f.write('>' + contigs_name + '\n' + seq + '\n')

def main(strain_name):

    
    path = os.path.join('/home/zhanghao/graminearum/sidr_analyses/sidr_tmp_files', strain_name)
    os.mkdir(path)
    os.chdir(path)
    build_index(strain_name)
    mem_alignment(strain_name)
    transfer_bam(strain_name)
    build_bam_index(strain_name)
    bulid_blast(strain_name)
    build_blast2(strain_name)
    remove_contaminated(strain_name)



    # save sample fas in dictionary {strain_name:seq}
    fas_seq_dic = get_seq_dic(strain_name)
    keep_contigs_list = get_keep_contigs(strain_name)
    clean_contigs = get_filtered_contig_seq(strain_name, fas_seq_dic, keep_contigs_list)
    os.chdir('/home/zhanghao/graminearum/scripts')

if __name__ == "__main__":
    # wur strain name file:
    #strain_name_file = '/home/zhanghao/graminearum/scripts/wur_isolates.txt'

    # ipp strain name file:
    #strain_name_file = '/home/zhanghao/graminearum/scripts/ipp_isolates.txt'

    # ncbi strain name file:
    #strain_name_file = '/home/zhanghao/graminearum/scripts/ncbi_isolates.txt'

    #left_sidr name file:
    strain_name_file = '/home/zhanghao/graminearum/scripts/isolates_the_missing_3.txt'

    got_strain_list = get_strain_list(strain_name_file)

    pool = Pool(25)
    pool.map(main, got_strain_list)

