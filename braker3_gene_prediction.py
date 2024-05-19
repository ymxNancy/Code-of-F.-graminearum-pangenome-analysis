#!/usr/bin/env python3
# -*- coding: utf-8 -*-


'''
docker run braker3
'''

import subprocess, os
from sys import argv
from multiprocessing import Pool

def get_name_list(name_file):
    with open(name_file, 'r') as f:
        name_list = []
        for line in f:
            line = line.strip()
            name_list.append(line)
        return name_list

def conduct_annotation(strain_name):
    cmd_1 = 'mkdir /home/zhanghao/graminearum/gene_prediction/de_novo_based/braker3_tmp/{0}'.format(strain_name)
    subprocess.check_call(cmd_1, shell=True)

    cmd_2 = 'docker run -u root -v /home/zhanghao/graminearum:/work --rm -it braker3_genmark-etp:v2 sh -c "cp -r /work/gene_prediction/de_novo_based/braker3_model/ph1_braker3/species/graminearum_ph_model /usr/share/augustus/config/species/ && braker.pl --skipAllTraining --genome=repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas --species=graminearum_ph_model --bam=gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_1.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_2.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_3.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_4.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_5.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_6.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_7.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_8.sort.bam,gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_9.sort.bam --gff3 --threads 2 --workingdir=gene_prediction/de_novo_based/braker3_tmp/{0}" && cp /home/zhanghao/graminearum/gene_prediction/de_novo_based/braker3_tmp/{0}/braker.gff3 /home/zhanghao/graminearum/gene_prediction/de_novo_based/braker3_result/{0}_braker3.gff3'.format(strain_name)

    check = subprocess.check_call(cmd_2, shell=True)
    if check == 0:
        print(strain_name + ' run braker3 successfully')
    else:
        print(strain_name + ' run braker3 failed')


if __name__ == '__main__':
    pool = Pool(55)
    name_file = '/home/zhanghao/graminearum/scripts/isolates_270.txt'

    list_name = get_name_list(name_file)
    pool.map(conduct_annotation, list_name)

