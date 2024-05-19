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

def perform_transdecoder(strain_name):

    cmd_1 = 'mkdir /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}'.format(strain_name)
    subprocess.check_call(cmd_1, shell=True)

    for i in range(1,10):
        cmd_2 = 'cd /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0} && gtf_genome_to_cdna_fasta.pl /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_{1}.gtf /home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas > {0}_group_{1}_transcripts.fasta && gtf_to_alignment_gff3.pl /home/zhanghao/graminearum/gene_prediction/rna_based/hisat2_stringtie_tmp/{0}/{0}_group_{1}.gtf > {0}_group_{1}.gff3 && TransDecoder.LongOrfs -t {0}_group_{1}_transcripts.fasta && TransDecoder.Predict -t {0}_group_{1}_transcripts.fasta && cdna_alignment_orf_to_genome_orf.pl {0}_group_{1}_transcripts.fasta.transdecoder.gff3 {0}_group_{1}.gff3 {0}_group_{1}_transcripts.fasta > {0}_group_{1}_transcripts.fasta.transdecoder.genome.gff3 && sed $"s/transdecoder\t/transdecoder_{1}\t/g" {0}_group_{1}_transcripts.fasta.transdecoder.genome.gff3 > {0}_group_{1}_transcripts.fasta.transdecoder.genome_rename.gff3 && sed $"s/Cufflinks\t/Cufflinks_{1}\t/g" {0}_group_{1}.gff3 > {0}_group_{1}_rename.gff3'.format(strain_name, i)
        subprocess.check_call(cmd_2, shell=True)
    
    cmd_3 = 'cd /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0} && cat {0}_group_1_rename.gff3 {0}_group_2_rename.gff3 {0}_group_3_rename.gff3 {0}_group_4_rename.gff3 {0}_group_5_rename.gff3 {0}_group_6_rename.gff3 {0}_group_7_rename.gff3 {0}_group_8_rename.gff3 {0}_group_9_rename.gff3 > {0}_group_1_9_rename.gff3'.format(strain_name)
    subprocess.check_call(cmd_3, shell=True)

def perform_evm(strain_name):
#    cmd_1 = "sed 's/Parent=/ID=/g' /home/zhanghao/graminearum/gene_prediction/homo_based_miniprot/miniprot_result/{0}_miniprot.gff > /home/zhanghao/graminearum/gene_prediction/homo_based_miniprot/miniprot_result/{0}_miniprot.gff3".format(strain_name)
#    subprocess.check_call(cmd_1, shell=True)

    cmd_2 = 'mkdir /home/zhanghao/graminearum/gene_prediction/merged_EVidenceModeler_tmp/{0} && cat /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_1_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_2_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_3_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_4_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_5_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_6_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_7_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_8_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_9_transcripts.fasta.transdecoder.genome_rename.gff3 /home/zhanghao/graminearum/gene_prediction/de_novo_based/braker3_result/{0}_braker3.gff3 > /home/zhanghao/graminearum/gene_prediction/merged_EVidenceModeler_tmp/{0}/{0}_combined_braker3_transdecoder_1_9.gff3 && docker run --rm  -v /home/zhanghao/graminearum:/data -it brianjohnhaas/evidencemodeler:latest \
   bash -c "cd /data/gene_prediction/merged_EVidenceModeler_tmp/{0} && EVidenceModeler \
               --sample_id {0} \
               --genome /data/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas \
               --weights /data/scripts/weights_EVM.txt \
               --gene_predictions /data/gene_prediction/merged_EVidenceModeler_tmp/{0}/{0}_combined_braker3_transdecoder_1_9.gff3 \
               --protein_alignments /data/gene_prediction/homo_based_miniprot/miniprot_result/{0}_miniprot.gff3 \
               --transcript_alignments /data/gene_prediction/rna_based/transdecoder_tmp/{0}/{0}_group_1_9_rename.gff3 \
               --segmentSize 100000 \
               --overlapSize 10000 \
               --CPU 1"'.format(strain_name)
    res = subprocess.check_call(cmd_2, shell=True)
    if res == 0:
        print(strain_name + ' run EVM successfully')
    else:
        print(strain_name + ' run EVM failed')

    cmd_3 = 'cp /home/zhanghao/graminearum/gene_prediction/merged_EVidenceModeler_tmp/{0}/{0}.EVM.gff3 /home/zhanghao/graminearum/gene_prediction/merged_EVidenceModeler_result/'.format(strain_name)
    subprocess.check_call(cmd_3, shell=True)

if __name__ == "__main__":
    #getted_name_list = parse_strain_names('test.txt')
    getted_name_list = parse_strain_names('isolates_271.txt')
    
    # Create a multiprocessing Pool
    pool = Pool(91)
    pool.map(perform_transdecoder, getted_name_list)
    # proces strains_list iterable with pool
    pool.map(perform_evm, getted_name_list)







