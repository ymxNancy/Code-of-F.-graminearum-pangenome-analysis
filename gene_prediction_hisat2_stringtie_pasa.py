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

def perform_hisat2_stringtie(strain_name):
    cmd_1 = 'mkdir /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0} && cd /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0} && ln -s /home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas ./ && hisat2-build -p9 {0}_ordered_soft.fas {0}'.format(strain_name)
    subprocess.check_call(cmd_1, shell=True)


    for i in range(1,10):
        cmd_2 = 'cd /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0} && hisat2 -x {0} --max-intronlen 5000 --dta -p 9 -1 /home/zhanghao/graminearum/northwest_af_university/RNAseq_data_qc/group_{1}_1.clean.fq.gz -2 /home/zhanghao/graminearum/northwest_af_university/RNAseq_data_qc/group_{1}_2.clean.fq.gz -S {0}_group_{1}.sam && samtools view -buS -@ 9 {0}_group_{1}.sam | samtools sort -m 7G -@ 9 -o {0}_group_{1}.sort.bam'.format(strain_name, i)
        res = subprocess.check_call(cmd_2, shell=True)
        if res == 0:
            print(strain_name + ': group-' + str(i) + ' run hisat2 successfully')
        else:
            print(strain_name + ': group-' + str(i) + ' run hisat2 failed')

        cmd_3 = 'cd /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0} && stringtie {0}_group_{1}.sort.bam -l {0}_group_{1} -o {0}_group_{1}.gtf -p 9'.format(strain_name, i)
        res = subprocess.check_call(cmd_3, shell=True)
        if res == 0:
            print(strain_name + ': group-' + str(i) + ' run stringtie successfully')
        else:
            print(strain_name + ': group-' + str(i) + ' run stringtie failed')

    cmd_4 = 'cd /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0} && stringtie --merge {0}_group_1.gtf {0}_group_2.gtf {0}_group_3.gtf {0}_group_4.gtf {0}_group_5.gtf {0}_group_6.gtf {0}_group_7.gtf {0}_group_8.gtf {0}_group_9.gtf -o /home/zhanghao/graminearum/gene_prediction/rna_based/stringtie_result/{0}_combined_stringtie.gtf'.format(strain_name)
    res = subprocess.check_call(cmd_4, shell=True)
    if res == 0:
        print(strain_name + ' run stringtie merge successfully')
    else:
        print(strain_name + ' run stringtie merge failed')
    
    cmd_5 = 'rm /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/*sam'
    subprocess.check_call(cmd_5, shell=True)

def perform_pasa(strain_name):
    cmd_1 = 'mkdir /home/zhanghao/graminearum/gene_prediction/rna_based/pasa_tmp/{0} && cd /home/zhanghao/graminearum/gene_prediction/rna_based/pasa_tmp/{0} && cp /home/zhanghao/graminearum/gene_prediction/rna_based/pasa_tmp/alignAssembly.config ./ && sed -i "6s/$/_{0}/" alignAssembly.config && cp /home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas ./ && gffread /home/zhanghao/graminearum/gene_prediction/rna_based/stringtie_result/{0}_combined_stringtie.gtf -g {0}_ordered_soft.fas -w {0}_stringtie_transcripts.fasta'.format(strain_name)
    subprocess.check_call(cmd_1, shell=True)

    cmd_2 = 'docker run -it -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/gene_prediction/rna_based/:/data -w /data/pasa_tmp/{0} pasapipeline/pasapipeline:latest /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g /data/pasa_tmp/{0}/{0}_ordered_soft.fas -t /data/pasa_tmp/{0}/{0}_stringtie_transcripts.fasta --ALIGNERS blat,gmap,minimap2 --CPU 10 && cp /home/zhanghao/graminearum/gene_prediction/rna_based/pasa_tmp/{0}/sample_mydb_pasa.sqlite_{0}.valid_minimap2_alignments.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/pasa_result/{0}_pasa_alignment.gff3'.format(strain_name)
    res = subprocess.check_call(cmd_2, shell=True)
    if res == 0:
        print(strain_name + ' run pasa successfully')
    else:
        print(strain_name + ' run pasa failed')


if __name__ == "__main__":
    getted_name_list = parse_strain_names('isolates_271.txt')
    #getted_name_list = parse_strain_names('test.txt')
    # Create a multiprocessing Pool
    pool = Pool(13)
    # proces strains_list iterable with pool
    #pool.map(perform_hisat2_stringtie, getted_name_list)
    pool.map(perform_pasa, getted_name_list)

# only one strain excute perform_pasa, spending 25:38.16 using 20 threads.