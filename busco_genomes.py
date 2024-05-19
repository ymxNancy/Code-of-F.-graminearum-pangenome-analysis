#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess, re
from multiprocessing import Pool

def parse_strain_names(name_file):
    name_list = []
    with open(name_file, 'r') as f:
        for line in f:
            line = line.strip()
            name_list.append(line)
    return name_list

def parse_busco(strain_name):
    cmd = 'docker run -it --rm -u $(id -u):$(id -u) -v /home/zhanghao:/busco_wd ezlabgva/busco:v5.5.0_cv1 busco -i graminearum/repeatmask_analyses/masked_genomes/{0}_ordered_soft.fas --out {0}_busco_genomes --out_path graminearum/repeatmask_analyses/busco_genomes/ -m genome --cpu 10 -l graminearum/sordariomycetes_odb10'.format(strain_name)
    subprocess.check_call(cmd, shell=True)

def get_statistics(strain_name_list, statistics_file):
    with open(statistics_file, 'w') as f1:
        f1.write('Strains\tComplete\tFragmented\tMissing\n')
        for strain_name in strain_name_list:
            with open('/home/zhanghao/graminearum/repeatmask_analyses/busco_genomes/{0}_busco_genomes/run_sordariomycetes_odb10/short_summary.txt'.format(strain_name), 'r') as f2:
                line_n = 0
                for line in f2:
                    line_n += 1
                    if line_n == 9:
                        line = line.strip()
                        searched = re.search(r'C:(.+)\[S.+,F:(.+),M:(.+),n.+', line)
                        com = searched.group(1)
                        fra = searched.group(2)
                        mis = searched.group(3)
                        f1.write('{0}\t{1}\t{2}\t{3}\n'.format(strain_name, com, fra, mis))

if __name__ == '__main__':
    name_file = '/home/zhanghao/graminearum/scripts/isolates_271.txt'
    statistics_file = '/home/zhanghao/graminearum/repeatmask_analyses/busco_genomes/busco_genomes_statistics.txt'

    got_name_list = parse_strain_names(name_file)

    pool = Pool(10)
    pool.map(parse_busco, got_name_list)

    got_statistics = get_statistics(got_name_list, statistics_file)









