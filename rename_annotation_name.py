#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

def parse_strain_names(name_file):
    name_list = []
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

def parse_gff_file(origin_gff_file):
    content_list = []
    gene_dir = {}
    with open(origin_gff_file, 'r') as f:
        paragraphs = f.read().split('\n\n')
        for paragraph in paragraphs:
            gene_dir = {}
            gene_line = paragraph.strip().split('\n')[0]
            if gene_line != '':
                gene_dir[gene_line] = paragraph
                content_list.append(gene_dir)
    return content_list

# 定义排序函数，按照键的第一个元素和第四个元素排序
def custom_sort(item):
    for key in item.keys():
        return item[key].split('\t')[0], int(item[key].split('\t')[3])

def parse_content_list(generated_content_list):
    # 使用sorted函数对列表进行排序，指定key参数为custom_sort
    sorted_content_list = sorted(generated_content_list, key=custom_sort)
    return sorted_content_list

def new_annotation_file(sorted_content_list, new_file, strain_name):
    with open(new_file, 'w') as f:
        gene_num = 0
        for gene in sorted_content_list:
            for gene_content in gene.values():
                gene_content_lines = gene_content.split('\n')
                one_gene_num = 0
                gene_id = ''
                for line in gene_content_lines:
                    line = line.strip().split('\t')
                    if line[2] == 'gene':
                        one_gene_num += 1
                        one_mrna_num = 0
                        exon_num = 0
                        cds_num = 0

                        if one_gene_num == 1:
                            gene_num += 1
                            gene_id = 'Fgra' + strain_name +'_G' + str(gene_num).zfill(5)
                            gene_line_8 = 'ID=' + gene_id + ';locus_tag=' + gene_id
                            f.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + gene_line_8 + '\n')
                        else:
                            print('More times "gene line" present in one gene')
                    elif line[2] == 'mRNA':
                        one_mrna_num += 1
                        if one_mrna_num == 1:
                            mrna_line_8 = 'ID=' + gene_id + '.t' + str(one_mrna_num) + ';Parent=' + gene_id + ';locus_tag=' + gene_id
                            f.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + mrna_line_8 + '\n')
                        else:
                            print('More times "mRNA line" present in one gene')
                    elif line[2] == 'exon':
                        exon_num += 1
                        exon_line_8 = 'ID=' + gene_id + '.t' + str(one_mrna_num)  + '.exon' + str(exon_num) + ';Parent=' + gene_id  + '.t' + str(one_mrna_num) + ';locus_tag=' + gene_id
                        f.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + exon_line_8 + '\n')
                    elif line[2] == 'CDS':
                        cds_num += 1
                        cds_line_8 = 'ID=' + gene_id + '.t' + str(one_mrna_num)  + '.CDS' + str(cds_num) + ';Parent=' + gene_id + '.t' + str(one_mrna_num)  + ';locus_tag=' + gene_id
                        f.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t' + line[7] + '\t' + cds_line_8 + '\n')
                    else:
                        print(line)

def compare_files(origin_file, new_annotation, strain_name):
    cmd = '[ $(grep -v "^$" {0} | wc -l ) -eq $(less -S {1} | wc -l) ] && echo "{2}: Same" || echo "{2}: Not the same"'.format(origin_file, new_annotation, strain_name)
    subprocess.check_call(cmd, shell=True)

if __name__ == '__main__':
    name_file = '/home/zhanghao/graminearum/scripts/isolates_271.txt'

    got_name_list = parse_strain_names(name_file)
    for strain_name in got_name_list:
        origin_gff_file = '/home/zhanghao/graminearum/gene_prediction/merged_EVidenceModeler_result/{0}.EVM.gff3'.format(strain_name)
        renamed_gff_file = '/home/zhanghao/graminearum/gene_prediction/renamed_EVM_annotation/{0}_renamed_EVM.gff3'.format(strain_name)

        got_content_list = parse_gff_file(origin_gff_file)
        got_sorted_content_list = parse_content_list(got_content_list)
        got_new_annotation = new_annotation_file(got_sorted_content_list, renamed_gff_file, strain_name)
        compare_files(origin_gff_file, renamed_gff_file, strain_name)



