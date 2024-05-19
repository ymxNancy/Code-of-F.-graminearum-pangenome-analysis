#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from multiprocessing import Pool
import os

'''
work dictionary: /home/zhanghao/graminearum/scripts

before running, motify ulimit -u 31308 or larger
command line: screen -L python3 bin/call_snp_289.py

use docker run gatk: 
    1. docker run -it -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd broadinstitute/gatk:latest bash
    2. docker run -it -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd broadinstitute/gatk:latest bash

    docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' GenomicsDBImport --genomicsdb-workspace-path my_database -R 180197_pilon4.fas -L intervals.list --sample-name-map output/input.list.test

PH-1 Illumina data download from: https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR15032575&display=download

'''

def prepare_ref(ref_fa):
    subprocess.call("bwa-mem2 index {}".format(ref_fa), shell=True)
    print('generate reference index successfully')
    subprocess.call("samtools faidx {}".format(ref_fa), shell=True)
    print('generate reference .fai file successfully')

    subprocess.call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk CreateSequenceDictionary -R {} -O PH1_nomt_renamed_ordered.dict".format(ref_fa), shell=True)
    print('generate reference .dict file successfully')

def generate_strain_list(list_file):
    strains_list = []
    with open(list_file, 'r') as f:
        for line in f:
            line = line.strip()
            strains_list.append(line)
    return strains_list

def get_vcf(strain_names): 
    line = strain_names.strip()
    if line == '180197':
        subprocess.call("bwa-mem2 mem -t 50 -M -R '@RG\\tID:{0}\\tLB:{0}\\tPL:illumina\\tSM:{0}\\tPU:{0}' PH1_nomt_renamed_ordered.fasta /home/zhanghao/asiaticum/filtered/R180_filtered_1.fq.gz /home/zhanghao/asiaticum/filtered/R180_filtered_2.fq.gz > output/bam/{0}.sam".format(line), shell=True)
    elif line == 'SRR15255548':
        subprocess.call("bwa-mem2 mem -t 50 -M -R '@RG\\tID:{0}\\tLB:{0}\\tPL:illumina\\tSM:{0}\\tPU:{0}' PH1_nomt_renamed_ordered.fasta /home/zhanghao/culmorum/ncbi/SRR15255548/qc_SRR15255548/{0}_1_qc.fastq.gz /home/zhanghao/culmorum/ncbi/SRR15255548/qc_SRR15255548/{0}_2_qc.fastq.gz > output/bam/{0}.sam".format(line), shell=True)

    subprocess.call("samtools view -buS output/bam/{0}.sam | samtools sort -m 300G -o output/bam/{0}.sort.bam".format(line), shell=True)

    subprocess.call("samtools index output/bam/{}.sort.bam".format(line), shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' FixMateInformation -I output/bam/{0}.sort.bam -O output/bam/{0}.fix.bam -MC true".format(line), shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' MarkDuplicates -I output/bam/{0}.fix.bam -O output/bam/{0}.mkd.bam -ASO coordinate -M output/bam/{0}_marked_dup_metrics.txt".format(line), shell=True)

    subprocess.check_call("samtools index output/bam/{}.mkd.bam".format(line), shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx4g' HaplotypeCaller -R PH1_nomt_renamed_ordered.fasta -I output/bam/{0}.mkd.bam -ERC GVCF -ploidy 1 -O output/gvcf/{0}.g.vcf.gz".format(line), shell=True)

    #   subprocess.call("rm output/bam/{0}.sam output/bam/{0}.sort.bam output/bam/{0}.fix.bam".format(line), shell=True)
    print(line + ' runs gvcf successfully')

def generate_vcf():
# 1st time use fm and hn9-1 as outgroups:
#     subprocess.check_call("find output/gvcf/ -name '*.g.vcf.gz' > output/add_outgroup/vcf.list", shell=True)

# 2nd time use fm and PH-1 as outgroups:
#    subprocess.check_call("find output/gvcf/ -name '*.g.vcf.gz' > output/add_outgroup_ph1_180197/vcf.list", shell=True)

# only fg, no outgroup
#    subprocess.check_call("find output/gvcf/ -name '*.g.vcf.gz' > output/fg_fa_265_273/vcf.list", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' CombineGVCFs -R PH1_nomt_renamed_ordered.fasta --variant output/fg_fa_265/vcf.list -O output/fg_fa_265/cohort.g.vcf", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' GenotypeGVCFs -R PH1_nomt_renamed_ordered.fasta -V output/fg_fa_265/cohort.g.vcf -O output/fg_fa_265/fg_fa_265_com_var.vcf", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' VariantFiltration -V output/fg_fa_265/fg_fa_265_com_var.vcf --filter-expression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filter-name 'my_filter' -O output/fg_fa_265/fg_fa_265_com_var_marked.vcf", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v /home/zhanghao/graminearum/call_snp_289:/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R PH1_nomt_renamed_ordered.fasta -V output/fg_fa_265/fg_fa_265_com_var_marked.vcf -O output/fg_fa_265/fg_fa_265_com_var_filtered_gatk.vcf -select 'vc.isNotFiltered()'", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R PH1_nomt_renamed_ordered.fasta -V output/fg_fa_265/fg_fa_265_com_var_filtered_gatk.vcf -O output/fg_fa_265/fg_fa_265_com_snp_filtered_gatk.vcf --select-type-to-include SNP", shell=True)

    subprocess.check_call("docker run --rm -u $(id -u):$(id -u) -v $(pwd):/gatk_wd -w /gatk_wd broadinstitute/gatk:latest gatk --java-options '-Xmx500g' SelectVariants -R PH1_nomt_renamed_ordered.fasta -V output/fg_fa_265/fg_fa_265_com_var_filtered_gatk.vcf -O output/fg_fa_265/fg_fa_265_com_indel_filtered_gatk.vcf --select-type-to-include INDEL", shell=True)

    print('finished')

if __name__ == '__main__':
    os.chdir('/home/zhanghao/graminearum/call_snp_289/')
    subprocess.check_call("cp /home/zhanghao/graminearum/rothamsted_research/PH1_nomt_renamed_ordered.fasta /home/zhanghao/graminearum/call_snp_289", shell=True)
    ref_fg_file = 'PH1_nomt_renamed_ordered.fasta'
# add outgroup do not need to re do this step:
    #prepared_ref = prepare_ref(ref_fg_file)

    #strain_name_file = '/home/zhanghao/graminearum/scripts/outgroup_name.txt'
    #generated_strain_list = generate_strain_list(strain_name_file)
    pool = Pool(2)
    #pool.map(get_vcf, generated_strain_list)

    generated_vcf = generate_vcf()


'''
conda activate R4.1

Rscript vcf_filter.R
# generate file: fg_fa_265_com_snp_filtered_gatk_vcfR.vcf.gz

# Then using vcftools filter the max missing:
screen -L vcftools --gzvcf fg_fa_265_com_snp_filtered_gatk_vcfR.vcf.gz --max-missing 1 --recode --recode-INFO-all --out fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools
# generate file: fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.vcf

sed '/\*/d' fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.vcf > fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.noasterisk.vcf

screen -L vcf2phylip.py -i fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.noasterisk.vcf -f -n

screen -L iqtree2 -s fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.noasterisk.min4.phy -m GTR_ASC -o 180197 -T 100 -B 1000
# ERROR: Invalid use of +ASC because of 30 invariant sites in the alignment
#The above command line generate a new file:`fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.noasterisk.min4.phy.varsites.phy`. Removing the invariant sites, Then using the new file to generate the tree.

screen -L iqtree2 -s fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.noasterisk.min4.phy.varsites.phy -m GTR+ASC -o 180197 -T 100 -B 1000

nw_reroot -l fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.noasterisk.min4.phy.varsites.treefile 180197 > reroot_fg_fa_265_com_snp_filtered_gatk_vcfR_vcftools.recode.noasterisk.min4.phy.varsites.treefile
'''

# CombineGVCFs --variant file name must be vcf.list


