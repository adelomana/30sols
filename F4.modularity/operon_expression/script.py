import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':6,'font.family':'Arial','xtick.labelsize':6,'ytick.labelsize':6})
#matplotlib.rcParams['pdf.fonttype']=42

###
### FUNCTIONS
###

def membership_reader():

    operon_membership={}
    operon_genes=[]
    with open(ribosomal_operons_file,'r') as f:
        next(f)
        for line in f:
            v=line.split()
            name=v[0]
            operon_membership[name]=[]
            elements=v[1:]
            for element in elements:
                operon_membership[name].append(element)
                operon_genes.append(element)

    return operon_membership,operon_genes

###
### MAIN
###

# 0. preamble

# 0.1. user defined variables
figure_file='figure.pdf'
ribosomal_operons_file='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/riboPtOperons.txt'
ribosomal_genes_file='/Volumes/omics4tb/alomana/projects/TLR/data/rp.transcription.groups/ribo.groupings.csv'
annotation_file='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/NC_002607.1.cs.NC_001869.cs.NC_002608.1.fasta'

# 0.2. read data

# 0.2.1. read operon memberships and positions
operon_membership,operon_genes=membership_reader()

# 0.2.

# 1. find operon memberships







# 2. iterate over operon memberships and build two figures, one with operon architecture and another one with gene expression dynamics

fig = matplotlib.pyplot.figure(constrained_layout=True,figsize=(2,8))
gs = fig.add_gridspec(8,2)

ax_big_operon = fig.add_subplot(gs[0:2,0:2])
ax_big_operon.set_title('gs[0, :]')

operon00 = fig.add_subplot(gs[2,0])
operon00.set_title('b')

operon01 = fig.add_subplot(gs[2,1])
operon01.set_title('c')

operon10 = fig.add_subplot(gs[3,0])
operon10.set_title('d')

operon11 = fig.add_subplot(gs[3,1])
operon11.set_title('e')


f3_ax5 = fig.add_subplot(gs[-2:,-2:])
f3_ax5.set_title('gs[-1, -2]')




matplotlib.pyplot.savefig(figure_file)

for operon in operon_membership.keys():
    print(operon,len(operon_membership[operon]))

# 3. make another figure with all ribosomal protein genes that are not in operons

