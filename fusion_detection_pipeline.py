'''
David van IJzendoorn, 20151210

Fusion detection pipeline running:
* STAR-fusion
* Tophat-fusion
* Defuse
* FusionCatcher
* Chimerascan
'''

from __future__ import print_function

from ruffus import *
from subprocess import call
import os
import time

#read files for pipeline
read_files = [('L1079-N_R1.fastq.gz', 'L1079-N_R2.fastq.gz'),
                ('L1079-T_R1.fastq.gz', 'L1079-T_R2.fastq.gz'),
                ('L1357-N_R1.fastq.gz', 'L1357-N_R2.fastq.gz'),
                ('L1357-T_R1.fastq.gz', 'L1357-T_R2.fastq.gz'),
                ('L2039-T1-T_R1.fastq.gz', 'L2039-T1-T_R2.fastq.gz')]

#locations for tools and references
star_reference = '/data/DIV3/ngs-sarcoma/Reference/star-fusion/Hg19_CTAT_resource_lib'
bowtie_reference = '/data/DIV3/ngs-sarcoma/Reference/hg19/bowtie/reference'
tophat_annotation = '/home/dgpvanijzendoorn/ngs-sarcoma/Reference/hg19/annotation'
tophat_blast = '/home/dgpvanijzendoorn/ngs-sarcoma/Reference/blast_human'
fusioncatcher_bin = '/home/dgpvanijzendoorn/ngs-sarcoma/Tools/fusioncatcher/bin'
fusioncatcher_data = '/home/dgpvanijzendoorn/ngs-sarcoma/Tools/fusioncatcher/data/current'
defuse_scripts = '/home/dgpvanijzendoorn/ngs-sarcoma/Tools/defuse/scripts'
chimerascan_bin = '/home/dgpvanijzendoorn/ngs-sarcoma/Tools/chimerascan/bin'
chimerascan_reference = '/home/dgpvanijzendoorn/ngs-sarcoma/Reference/chimerascan'
cpu = '10'

@transform(read_files, suffix('.fastq.gz'), r'star_fusion_\1')
def star_fusion(input_file, output_file):
    '''Run star_fusion.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    run_time = time.time()
    call(['STAR-Fusion',
    '--genome_lib_dir', star_reference,
    '--left_fq', ii1,
    '--right_fq', ii2,
    '--output_dir', output_file,
    '--CPU', cpu])
    run_time = time.time() - run_time
    os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    os.system('rm -rf '+output_file+'/Aligned.sortedByCoord.out.bam')
    return None

@transform(read_files, suffix('.fastq.gz'), r'tophat_\1')
def tophat_fusion(input_file, output_file):
    '''Run tophat_fusion.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    run_time = time.time()
    os.system('mkdir '+output_file)
    call(['tophat',
    '-o', output_file+'/tophat_out',
    '-p', cpu,
    '--fusion-search',
    '--keep-fasta-order',
    '--bowtie1',
    '--no-coverage-search',
    '-r', '0',
    '--mate-std-dev', '80',
    '--max-intron-length', '100000',
    '--fusion-min-dist', '100000',
    '--fusion-anchor-length', '13',
    '--fusion-ignore-chromosomes', 'chrM',
    bowtie_reference,
    ii1,
    ii2])
    run_time = time.time() - run_time
    os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

@transform(tophat_fusion, suffix(''), r'\1/tophatfusion_out')
def tophat_fusion_post(input_file, output_file):
    '''Run tophat_fusion_post.'''
    run_time = time.time()
    os.system('ln -s '+tophat_annotation+'/* '+input_file+'/')
    os.system('ln -s '+tophat_blast+' '+input_file+'/blast')
    os.system('(cd '+input_file+' && exec tophat-fusion-post '+bowtie_reference+')')
    os.system('rm '+input_file+'/blast '+input_file+'/ensGene.txt '+input_file+'/Homo_sapiens.GRCh37.75.gtf '+input_file+'/refGene.txt')
    run_time = time.time() - run_time
    os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    os.system('rm -rf '+output_file+'/../tophat_out')
    return None

@transform(read_files, suffix('.fastq.gz'), r'fusioncatcher_\1')
def fusioncatcher(input_file, output_file):
    '''Run FusionCatcher.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    run_time = time.time()
    call([fusioncatcher_bin+'/fusioncatcher',
    '-p', cpu,
    '-d', fusioncatcher_data,
    '-i', ii1+','+ii2,
    '-o', output_file])
    run_time = time.time() - run_time
    os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

@transform(read_files, suffix('.fastq.gz'), '.fastq')
def gunzip_files(input_file, output_file):
    '''Gunzip files for Defuse.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    out1 = ii1.split('.')[0]
    out2 = ii2.split('.')[0]
    os.system('gunzip -c '+ii1+' >'+out1+'.fastq')
    os.system('gunzip -c '+ii2+' >'+out2+'.fastq')
    return None

@transform(gunzip_files, suffix('.fastq'), r'defuse_\1')
def defuse(input_file, output_file):
    '''Run Defuse on extracted file.'''
    base_name = input_file.split('_')[0]
    run_time = time.time()
    call([defuse_scripts+'/defuse.pl',
    '-c', defuse_scripts+'/config.txt',
    '-1', base_name+'_R1.fq',
    '-2', base_name+'_R2.fq',
    '-o', output_file,
    '-p', cpu])
    run_time = time.time() - run_time
    os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    os.system('rm -rf '+output_file+'/cdna.pair.sam '
                +output_file+'/reads*')
    return None

@transform(read_files, suffix('.fastq.gz'), r'chimerascan_\1')
def chimerascan(input_file, output_file):
    '''Run Chimerascan.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    run_time = time.time()
    call(['python', chimerascan_bin+'/chimerascan_run.py',
    '-p', cpu,
    chimerascan_reference,
    ii1,
    ii2,
    output_file])
    run_time = time.time() - run_time
    os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    os.system('rm -rf '+output_file+'/tmp '+output_file+'/aligned_reads.bam '
                +output_file+'/sorted_aligned_reads.bam*')
    return None

pipeline_run(multiprocess=1)
#pipeline_printout() multiprocess=2
#pipeline_printout_graph('flowchart.jpg', 'jpg')
