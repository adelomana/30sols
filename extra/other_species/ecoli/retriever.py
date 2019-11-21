import os,sys
import multiprocessing,multiprocessing.pool

def sample_retriever(sample):

    print('')
    cmd1='time prefetch {} -v --log-level 6 --output-directory {}'.format(sample,output_dir)
    print('\t {}'.format(cmd1))
    os.system(cmd1)

    print('')
    cmd2='time fastq-dump -v -W {}{}.sra --outdir {}'.format(output_dir,sample,output_dir)
    print('\t {}'.format(cmd2))
    os.system(cmd2)
        
    return None

###
### MAIN
###

# 0. user defined variables
accession_list_file='list.txt'
output_dir='/Users/alomana/scratch/third_run/'
threads=6

# 1. read samples to download
samples=[]
with open(accession_list_file,'r') as f:
    for line in f:
        v=line.split()
        samples.append(v[0])
print(samples)

# 2. iterate
hydra=multiprocessing.pool.Pool(threads)
hydra.map(sample_retriever,samples)

# 3. last message
print('')
print('work completed.')
