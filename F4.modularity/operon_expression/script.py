import sys, operator, numpy
import matplotlib ,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':6,'font.family':'Arial','xtick.labelsize':6,'ytick.labelsize':6})
matplotlib.rcParams['pdf.fonttype']=42

###
### FUNCTIONS
###

def gene_positions_reader():

    gene_positions={}

    with open(annotation_file,'r') as f:
        for line in f:
            if line[0] == '>':
                name='gene-'+line.split('locus_tag=')[1].split(']')[0]
                positions=line.split('location=')[1].split(']')[0].split('..')
                clean_positions=[]
                direction='forward'
                for position in positions:
                    if 'complement' in position:
                        direction='reverse'
                    position=position.replace('complement(','')
                    position=position.replace(')','')
                    position=position.replace('>','')
                    position=position.replace('<','')
                    position=int(position)
                    clean_positions.append(position)
                gene_positions[name]=[clean_positions,direction]

    return gene_positions

def inset_plotter(working_genes,operon_name,hand):

    cmap=matplotlib.cm.get_cmap('tab10',10)
    color_count=0
    for gene in working_genes:

        # f.1. select color
         
        if (gene not in ribosomal_genes_in_operons) and (gene not in ribosomal_genes_not_in_operons):
            theColor='black'; theAlpha=0.8
        else:
            theColor=cmap(color_count); theAlpha=0.8
            color_count=color_count+1
            if color_count == 10:
                color_count=0

        # f.2. compute trajectory
        gene_name=gene.replace('gene-VNG_RS','VNGRS')
        x=[1,2,3,4]
        y=[]
        for timepoint in timepoints:
            average=numpy.mean(numpy.log10([expression['trna'][replicate][timepoint][gene_name]+1 for replicate in replicates]))
            y.append(average)
        top=numpy.max(y)
        if top == 0:
            y=[0,0,0,0]
        else:
            y=numpy.log2(numpy.array(y)/top)

        # f.3. plot
        matplotlib.pyplot.plot(x,y,'o-',lw=2/3,alpha=theAlpha,color=theColor,ms=2,mew=0)

    # limits and ticks
    if hand == 'right':
        matplotlib.pyplot.yticks((-1,0),[])
    else:
        matplotlib.pyplot.yticks((-1,0),('-1','0'))
        
    matplotlib.pyplot.xticks((1,2,3,4),('TP1','TP2','TP3','TP4'))
    matplotlib.pyplot.xlim([0.85,4.15])
    matplotlib.pyplot.ylim([-1.5,0.2])
    matplotlib.pyplot.grid(ls=':',alpha=1/6)
    matplotlib.pyplot.title(operon_name,pad=3)

    return None

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

    # sort operons based on size
    operon_sizes={}
    for operon in operon_membership.keys():
        operon_sizes[operon]=len(operon_membership[operon])
    sorted_operon_sizes=sorted(operon_sizes.items(),key=operator.itemgetter(1),reverse=True)
    
    return operon_membership,operon_genes,sorted_operon_sizes

def operon_structure_builder(ax,working_genes):


    print(working_genes)
    all_values=[]
    for name in working_genes:
        print(name,gene_positions[name])
        all_values.append(gene_positions[name][0][0])
        all_values.append(gene_positions[name][0][1])
        
    five_end=numpy.min(all_values)
    three_end=numpy.max(all_values)
    span=three_end-five_end
    print(five_end,three_end)

    cmap = matplotlib.cm.get_cmap('tab10',10)
    color_count=0

    direction=None
    
    for gene_name in working_genes:
        a=((gene_positions[gene_name][0][0]-five_end)/span)*2
        b=((gene_positions[gene_name][0][1]-five_end)/span)*2
        direction=gene_positions[gene_name][1]
        
        shift=1
        low_left=(shift+a,-1.25)
        w=b-a
        h=0.1

        # define the color of the box
        if (gene_name not in ribosomal_genes_in_operons) and (gene_name not in ribosomal_genes_not_in_operons):
            theColor='black'; theAlpha=1
        else:
            theColor=cmap(color_count); theAlpha=1
            color_count=color_count+1
            if color_count == 10:
                color_count=0

        rec=matplotlib.patches.Rectangle(low_left,w,h,color=theColor,alpha=theAlpha,zorder=100)
        ax.add_patch(rec)

    # draw an arrow
    print(direction)
    if direction == 'forward':
        x2, y2 = 0.92, -1.4
        x1, y1 = 1.5, -0.9
    else:
        x1, y1 = 2.5, -0.9
        x2, y2 = 3.08, -1.4
    ax.annotate("",xy=(x1,y1),xycoords='data',xytext=(x2,y2),textcoords='data',arrowprops=dict(arrowstyle="->",color="0.5",patchA=None,patchB=None,connectionstyle='angle,angleA=-90,angleB=180,rad=0'))

    print()
        
    return None

def ribosomal_gene_names_reader():

    ribosomal_gene_names=[]
    with open(ribosomal_genes_file,'r') as f:
        next(f)
        for line in f:
            v=line.split(',')
            ribosomal_gene_names.append(v[1])

    return ribosomal_gene_names

def transcriptomicsReader():

    '''
    this function reads transcriptomics data as in
    transcriptomics[trna/rbf][replicate][timepoint][gene]
    '''

    data={}
    geneNames=[]; timepoints=[]; replicates=[]
    
    with open(expression_file,'r') as f:
        header=f.readline()
        labels=header.split('\t')[1:-1]
        for label in labels:
            crumbles=label.split('.')

            fraction=crumbles[0]
            replicate='br'+crumbles[2]
            timepoint='tp.'+crumbles[4]

            if replicate not in replicates:
                replicates.append(replicate)
            if timepoint not in timepoints:
                timepoints.append(timepoint)

            if fraction not in data.keys():
                data[fraction]={}
            if replicate not in data[fraction].keys():
                data[fraction][replicate]={}
            if timepoint not in data[fraction][replicate].keys():
                data[fraction][replicate][timepoint]={}
            
        for line in f:
            vector=line.split('\t')[:-1]
            values=[float(element) for element in vector[1:]]
            geneName=vector[0].replace('_','')
            if geneName not in geneNames:
                geneNames.append(geneName)
            for i in range(len(values)):
                crumbles=labels[i].split('.')
                fraction=crumbles[0]
                replicate='br'+crumbles[2]
                timepoint='tp.'+crumbles[4]
                data[fraction][replicate][timepoint][geneName]=values[i]
    
    return data,geneNames,timepoints,replicates

###
### MAIN
###

# 0. preamble

# 0.1. user defined variables
figure_file='figure.pdf'
ribosomal_operons_file='/Volumes/omics4tb/alomana/projects/TLR/data/microbesOnline/riboPtOperons.txt'
ribosomal_genes_file='/Volumes/omics4tb/alomana/projects/TLR/data/rp.transcription.groups/ribo.groupings.csv'
annotation_file='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/NC_002607.1.cs.NC_001869.1.cs.NC_002608.1.fasta'
expression_file='/Volumes/omics4tb/alomana/projects/TLR/data/expression1e3/expressionMatrix.kallisto.txt'

# 0.2. read data
# 0.2.1. read operon memberships and positions
operon_membership,operon_genes,sorted_operon_sizes=membership_reader()
# 0.2.2. read expression
expression,tempo,timepoints,replicates=transcriptomicsReader()
# 0.2.3. read ribosomal genes
ribosomal_gene_names=ribosomal_gene_names_reader()
# 0.2.4 read gene positions
gene_positions=gene_positions_reader()

# 1. find operon memberships
ribosomal_genes_in_operons=[element for element in operon_genes if element in ribosomal_gene_names]
ribosomal_genes_not_in_operons=[element for element in ribosomal_gene_names if element not in ribosomal_genes_in_operons]
#! print(list(set(ribosomal_genes_in_operons) & set(ribosomal_genes_not_in_operons))) # this should be an empty list


# 2. iterate over operon memberships and build figures
fig=matplotlib.pyplot.figure(figsize=(2,6))
fig.set_tight_layout(True) #constrained_layout=True
gs=fig.add_gridspec(11,2)

# big operon
big_operon=fig.add_subplot(gs[0:2,0:2])
working_operon=sorted_operon_sizes[0][0]
working_genes=operon_membership[working_operon]
operon_structure_builder(big_operon,working_genes)
inset_plotter(working_genes,working_operon,None)
print('** for table ** \t {}\t{}'.format(working_operon, ', '.join([element.replace('gene-', '') for element in working_genes])))

# small operons
for i in range(7):
    # left panel
    operon_left=fig.add_subplot(gs[i+2,0])
    working_operon=sorted_operon_sizes[i+1][0]
    working_genes=operon_membership[working_operon]
    operon_structure_builder(operon_left,working_genes)
    inset_plotter(working_genes,working_operon,'left')
    print('** for table ** \t {}\t{}'.format(working_operon, ', '.join([element.replace('gene-', '') for element in working_genes])))

    # right panel
    operon_right=fig.add_subplot(gs[i+2,1])
    working_operon=sorted_operon_sizes[i+2][0]
    working_genes=operon_membership[working_operon]
    operon_structure_builder(operon_right,working_genes)
    inset_plotter(working_genes,working_operon,'right')
    print('** for table ** \t {}\t{}'.format(working_operon, ', '.join([element.replace('gene-', '') for element in working_genes])))

# final row with all ribosomal protein genes that are not in operons
final=fig.add_subplot(gs[-2:,-2:])
inset_plotter(ribosomal_genes_not_in_operons,'RP genes not in operons','None')
print('** for table ** \t {}\t{}'.format(working_operon, ', '.join([element.replace('gene-', '') for element in working_genes])))

# close figure
matplotlib.pyplot.savefig(figure_file)
matplotlib.pyplot.clf()
