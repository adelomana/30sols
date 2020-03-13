import os,numpy,sys

def caller(sample):

    '''
    this function calls kallisto
    '''

    
    fastq_files=' '.join([fastqDir+element+'_clean.fastq' for element in samples[sample]])

    strandFlag='--fr-stranded'
    
    cmd='time kallisto quant -i {} -o {}{} --bias --single -l 180 -s 20 -t 8 -b {} {} {}'.format(transcriptomeIndex,quantDir,sample,boots,strandFlag,fastq_files)

    print()
    print(cmd)
    print()

    os.system(cmd)
    
    return None

### MAIN

# 0. user defined variables
fastqDir='/Users/alomana/scratch/ecoli_GSE53767/clean/'
transcriptomeIndex='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/transcriptome/511145.transcriptomes.fasta.index'
boots=int(1e3)
quantDir='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli_GSE53767/kallisto.1e{}/'.format(int(numpy.log10(boots)))

if os.path.exists(quantDir) == False:
    os.mkdir(quantDir)

# 1. reading files
print('reading files...')

# 1.1. defining fastq files
samples={}
samples['mRNA']=['SRR1067773','SRR1067774']
samples['footprints']=['SRR1067765','SRR1067766','SRR1067767','SRR1067768']

# 2. processing
print('processing files...')
for sample in samples:
    caller(sample)
