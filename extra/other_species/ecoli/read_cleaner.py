import sys,os

def cleaner(sample_label):

    input_file=fastq_dir+sample_label+'.fastq'
    output_file=clean_fastq_dir+sample_label+'_clean.fastq'
    
    blocks=['time java -jar /Users/alomana/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads {}'.format(threads),
                '{} {} ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(input_file,output_file,adapters_file)
        ]
    cmd=' '.join(blocks)

    print('')
    print(cmd)
    print('')

    os.system(cmd)

    return None

###
### MAIN
###

# 0. user defined variables
adapters_file='/Users/alomana/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa'
fastq_dir='/Users/alomana/scratch/'
clean_fastq_dir='/Users/alomana/scratch/clean_fastq/'
threads=8

# 1. read sample labels
files=os.listdir(fastq_dir)
sample_labels=[element.split('.fastq')[0] for element in files if '.fastq' in element]

# 2. iterate sample labels
for sample_label in sample_labels:
    cleaner(sample_label)
