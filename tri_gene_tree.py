#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess

def generate_tri_fas():
    for i in 1,3,12:
        with open('/home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/tri{0}_tree/TRI{0}_graminearum_271_match.blastn'.format(i), 'r') as f1, open('/home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/tri{0}_tree/TRI{0}_graminearum_271.fas'.format(i), 'w') as f2:
            for line in f1:
                line = line.strip().split()
                strain_name = line[1].split('_')[1]
                sequence = line[13]
                f2.write('>' + str(strain_name) + '\n' + sequence + '\n')

def generate_tree():
    for i in 1,3,12:
        if i == 1:
            outgroup_name = 'XM_009265144.1'
        elif i == 3:
            outgroup_name = 'XM_009265280.1'
        elif i == 12:
            outgroup_name = 'XM_009265274.1'
        cmd = 'cd /home/zhanghao/graminearum/repeatmask_analyses/masked_genomes/tri{0}_tree/ && mafft --maxiterate 1000 --localpair --thread 30 TRI{0}_outgroup_272.fas > TRI{0}_outgroup_272_mafft.fas && iqtree2 -s TRI{0}_outgroup_272_mafft.fas -m MFP -B 1000 -T AUTO -o {1}'.format(i, outgroup_name)
        subprocess.check_call(cmd, shell=True)

if __name__ == "__main__":
    generated_tri_fas = generate_tri_fas()
    generated_tree = generate_tree()


