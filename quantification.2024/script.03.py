import os

#
# user-defined variables
#
raw_fastq_dir = '/Users/adrian/research/tlr/data/fastq/'
clean_fastq_dir = '/Users/adrian/research/tlr/data/clean_fastq/'
trimmomatic_path = '/Users/adrian/software/Trimmomatic-0.39/trimmomatic-0.39.jar'
adapter_file = '/Users/adrian/software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa'

number_threads = 4

#
# read information
#
samples = os.listdir(raw_fastq_dir)
samples.sort()
print(samples)

#
# call trimmomatic
#
for sample in samples:

    print(sample)

    executable = 'time java -jar {} SE -threads {} -phred33 '.format(trimmomatic_path, number_threads)
    options=' ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'.format(adapter_file)
    input_file = '{}{}/{}.fastq'.format(raw_fastq_dir, sample, sample)
    output_file = '{}/{}.fastq'.format(clean_fastq_dir, sample, sample)

    command = executable + input_file + ' ' + output_file + options

    print('')
    print(command)
    print('')

    os.system(command)
