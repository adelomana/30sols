import os

#
# user-defined variables
#
sra_dir = '/Users/adrian/research/tlr/data/sra/'
fastq_dir = '/Users/adrian/research/tlr/data/fastq/'
executable = '/Users/adrian/software/sratoolkit.3.0.10-mac-x86_64/bin/fastq-dump'

#
# read folders
#
all_files = os.listdir(sra_dir)
all_dirs = [element for element in all_files if 'SRR6' in element]
all_dirs.sort()
print(all_dirs)

#
# execute command
#
for input_dir in all_dirs:

    output_dir = fastq_dir + input_dir

    cmd = 'time {} {}{}/{}.sra --outdir {} --split-3 --v'.format(executable, sra_dir, input_dir, input_dir, output_dir)

    print()
    print(cmd)
    print()

    os.system(cmd)
