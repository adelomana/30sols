import os,numpy,sys

def caller(element):

    '''
    this function calls kallisto
    '''

    tag=element.split('_clean.fastq')[0]

    strandFlag='--fr-stranded' # --rf-stranded /// --fr-stranded
    
    cmd='time kallisto quant -i {} -o {}{} --bias --single -l 180 -s 20 -t 8 -b {} {} {}{}'.format(transcriptomeIndex,quantDir,tag,boots,strandFlag,fastqDir,element)

    print()
    print(cmd)
    print()

    os.system(cmd)
    
    return None

### MAIN

# 0. user defined variables
fastqDir='/Users/alomana/scratch/not_clean/'
transcriptomeIndex='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/transcriptome/511145.transcriptomes.fasta.index'

boots=int(1e1)

quantDir='/Volumes/omics4tb2/alomana/projects/TLR/data/ecoli/kallisto1e{}/'.format(int(numpy.log10(boots)))

if os.path.exists(quantDir) == False:
    os.mkdir(quantDir)

# 1. reading files
print('reading files...')

# 1.1. defining fastq files
items=os.listdir(fastqDir)
files=[element for element in items if '.fastq' in element]

print(files)

# 2. processing
print('processing files...')
for element in files:
    caller(element)
