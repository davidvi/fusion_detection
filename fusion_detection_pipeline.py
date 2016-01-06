'''
David van IJzendoorn, 20151210

Fusion detection pipeline running:
* STAR-fusion
* Tophat-fusion
* Defuse
* FusionCatcher
* Chimerascan

TODO:
* Tophat-fusion-post implementation
'''

from __future__ import print_function

from ruffus import *
from subprocess import call
import os
import time

read_files = [('01LLUMC_L4027_1.fq.gz', '01LLUMC_L4027_2.fq.gz'),
            ('02LLUMC_L872_1.fq.gz', '02LLUMC_L872_2.fq.gz'),
            ('03LLUMC_L1062_1.fq.gz', '03LLUMC_L1062_2.fq.gz'),
            ('04LLUMC_RDES_1.fq.gz', '04LLUMC_RDES_2.fq.gz'),
            ('05LLUMC_EW3_1.fq.gz', '05LLUMC_EW3_2.fq.gz'),
            ('06LLUMC_EW7_1.fq.gz', '06LLUMC_EW7_2.fq.gz'),
            ('07LLUMC_SKES1_1.fq.gz', '07LLUMC_SKES1_2.fq.gz'),
            ('08LLUMC_EW3mutp53_1.fq.gz', '08LLUMC_EW3mutp53_2.fq.gz')]

@transform(read_files, suffix('.fq.gz'), r'star_fusion_\1')
def star_fusion(input_file, output_file):
    '''Run star_fusion.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    if not os.path.isfile(output_file+'/star-fusion.fusion_candidates.final'):
        run_time = time.time()
        call(['STAR-Fusion',
        '--genome_lib_dir', '/home/dijzendoorn/reference/hg19/star-fusion/Hg19_CTAT_resource_lib',
        '--left_fq', ii1,
        '--right_fq', ii2,
        '--output_dir', output_file])
        run_time = time.time() - run_time
        os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

@transform(read_files, suffix('.fq.gz'), r'tophat_\1')
def tophat_fusion(input_file, output_file):
    '''Run tophat_fusion.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    if not os.path.isfile(output_file+'/tophat_out/fusions.out'):
        run_time = time.time()
        os.system('mkdir '+output_file)
        call(['tophat',
        '-o', output_file+'/tophat_out',
        '-p', '10',
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
        '/home/dijzendoorn/reference/hg19/bowtie/reference',
        ii1,
        ii2])
        run_time = time.time() - run_time
        os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

@transform(tophat_fusion, suffix(''), r'\1/tophatfusion_out')
def tophat_fusion_post(input_file, output_file):
    '''Run tophat_fusion_post.'''
    if not os.path.isfile(output_file+'/final-list_candidate-fusion-genes.GRCh37.txt'):
        run_time = time.time()
        os.system('ln -s /home/dijzendoorn/reference/hg19/annotation/* '+input_file+'/')
        os.system('ln -s /home/dijzendoorn/blast_human '+input_file+'/blast')
        os.system('(cd '+input_file+' && exec tophat-fusion-post /home/dijzendoorn/reference/hg19/bowtie/reference)')
        os.system('rm '+input_file+'/blast '+input_file+'/ensGene.txt '+input_file+'/Homo_sapiens.GRCh37.75.gtf '+input_file+'/refGene.txt')
        run_time = time.time() - run_time
        os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

@transform(read_files, suffix('.fq.gz'), r'fusioncatcher_\1')
def fusioncatcher(input_file, output_file):
    '''Run FusionCatcher.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    if not os.path.isfile(output_file+'/final-list_candidate-fusion-genes.GRCh37.txt'):
        run_time = time.time()
        call(['/home/dijzendoorn/tools/fusioncatcher/bin/fusioncatcher',
        '-p', '10',
        '-d', '/home/dijzendoorn/tools/fusioncatcher/data/current/',
        '-i', ii1+','+ii2,
        '-o', output_file])
        run_time = time.time() - run_time
        os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

@transform(read_files, suffix('.fq.gz'), '.fq')
def gunzip_files(input_file, output_file):
    '''Gunzip files for Defuse.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    out1 = ii1.split('.')[0]
    out2 = ii2.split('.')[0]
    if not os.path.isfile(output_file):
        os.system('gunzip -c '+ii1+' >'+out1+'.fq')
        os.system('gunzip -c '+ii2+' >'+out2+'.fq')
    return None

@transform(gunzip_files, suffix('.fq'), r'defuse_\1')
def defuse(input_file, output_file):
    print(output_file)
    '''Run Defuse on extracted file.'''
    base_name = input_file.split('_')[0] +'_'+ input_file.split('_')[1]
    print(base_name)
    if not os.path.isfile(output_file+'/results.tsv'):
        run_time = time.time()
        call(['/home/dijzendoorn/tools/defuse/scripts/defuse.pl',
        '-c', '/home/dijzendoorn/tools/defuse/scripts/config.txt',
        '-1', base_name+'_1.fq',
        '-2', base_name+'_2.fq',
        '-o', output_file,
        '-p', '10'])
        run_time = time.time() - run_time
        os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

@transform(read_files, suffix('.fq.gz'), r'chimerascan_\1')
def chimerascan(input_file, output_file):
    '''Run Chimerascan.'''
    ii1 = input_file[0]
    ii2 = input_file[1]
    if not os.path.isfile(output_file):
        run_time = time.time()
        call(['python', '/home/dijzendoorn/tools/chimerascan/bin/chimerascan_run.py',
        '-p', '10',
        '/home/dijzendoorn/reference/chimerascan',
        ii1,
        ii2,
        output_file])
        run_time = time.time() - run_time
        os.system('echo '+str(run_time)+ '> '+output_file+'/run_time')
    return None

pipeline_run(multiprocess = 3)
#pipeline_printout() multiprocess=2
#pipeline_printout_graph('flowchart.jpg', 'jpg')
