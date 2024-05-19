#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# /home/zhanghao/graminearum/functional_annotation/interproscan_annotation/PH1_proteins.fas.tsv

def extract_go(interpro_gff):
    go_dic = {}
    with open(interpro_gff, 'r') as f:
        for line in f:
            line = line.strip()
            if 'GO:' in line:
                line = line.strip().split('\t')
                gene_name = line[0]
                go_dic[gene_name] = line[13]
    return go_dic

def generate_ph1_go_file(go_dic, ph1_go_file):
    with open(ph1_go_file, 'w') as f:
        for gene_name, go_info in go_dic.items():
            if '|' in go_info:
                for each_go in go_info.split('|'):
                    f.write(gene_name + '\t' + each_go + '\n')
            else:
                f.write(gene_name + '\t' + go_info + '\n')

if __name__ == '__main__':
    interpro_gff = '/home/zhanghao/graminearum/functional_annotation/interproscan_annotation/PH1_proteins.fas.tsv'
    ph1_go_file = '/home/zhanghao/graminearum/functional_annotation/interproscan_annotation/PH1_GO_info.txt'
    got_go_dic = extract_go(interpro_gff)
    generated_ph1_go_file = generate_ph1_go_file(got_go_dic, ph1_go_file)












