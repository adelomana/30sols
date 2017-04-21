import os,sys,numpy
import sklearn,sklearn.decomposition,sklearn.manifold
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

def caller(element):

    '''
    this function calls kallisto
    '''

    tag=element.split('.fastq')[0]
    
    cmd='time kallisto quant -i {} -o {}{} --bias --plaintext --single -l 180 -s 20 -t 4 -b 100 {}{}'.format(transcriptomeIndex,quantDir,tag,fastqDir,element)

    print()
    print(cmd)
    print()

    sys.exit()
    os.system(cmd)
    
    return None

### MAIN

# 0. user defined variables
fastqDir='/Volumes/omics4tb/alomana/projects/TLR/data/cleanFASTQ/'
quantDir='/Volumes/omics4tb/alomana/projects/TLR/data/kallisto/'
transcriptomeIndex='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/hsa.ASM680v1.20170116.idx'
transcriptomeFastaFile='/Volumes/omics4tb/alomana/projects/TLR/data/transcriptome/hsa.ASM680v1.fasta'
resultsDir='/Volumes/omics4tb/alomana/projects/TLR/data/expression/'

# 1. reading files
print('reading files...')

# 1.1. defining fastq files
items=os.listdir(fastqDir)
files=[element for element in items if '.fastq' in element]

# 1.2. reading gene IDs and gene names
synonyms={}
with open(transcriptomeFastaFile,'r') as f:
    for line in f:
        if line[0] == '>':

            id=line.split('>')[1].split(' ')[0]
            if 'locus_tag' in line:
                name=line.split('[locus_tag=')[1].split(' ')[0].replace(']','')
            elif 'gene' in line:
                name=line.split('[gene=')[1].split(' ')[0].replace(']','')
            else:
                print('error while parsing transcriptome fasta file...')
                sys.exit()

            synonyms[id]=name

# 2. processing 
### print('processing files...')
### for element in files:
###     caller(element)

# 3. generating full expression matrix
print('generating expression matrix file...')

# 3.1. reading the genes
genes=[]
oneFile=quantDir+files[0].split('.fastq')[0]+'/abundance.tsv'
f=open(oneFile,'r')
next(f)
for line in f:
    vector=line.split()
    geneName=vector[0]
    genes.append(geneName)
f.close()

# 3.2. reading expression
expression={}
for element in files:
    tag=element.split('.clean.fastq')[0]
    if tag not in expression.keys():
        expression[tag]={}

    workingFile=quantDir+element.split('.fastq')[0]+'/abundance.tsv'
    f=open(workingFile,'r')
    next(f)
    for line in f:
        vector=line.split()
        geneName=vector[0]
        abundance=float(vector[-1])
        expression[tag][geneName]=abundance
    f.close()

# 3.3. writing expression matrix
expressionFile=resultsDir+'expressionMatrix.kallisto.txt'
conditionNames=list(expression.keys())

rbfConditions=[element for element in conditionNames if 'rbf' in element]
inverse=[element[::-1] for element in rbfConditions]
inverse.sort()
revertedRBF=[element[::-1] for element in inverse]

trnaConditions=[element for element in conditionNames if 'trna' in element]
inverse=[element[::-1] for element in trnaConditions]
inverse.sort()
revertedTRNA=[element[::-1] for element in inverse]

reverted=revertedTRNA+revertedRBF

x=[]
theEdgeColors=[]
theFaceColors=[]
theMarkers=[]

g=open(expressionFile,'w')

g.write('\t')
for i in range(len(reverted)):
    g.write('{}\t'.format(reverted[i]))
    
    x.append([])

    if 'tp.1' in reverted[i]:
        theEdgeColors.append('blue')
    elif 'tp.2' in reverted[i]:
        theEdgeColors.append('green')
    elif 'tp.3' in reverted[i]:
        theEdgeColors.append('orange')
    else:
        theEdgeColors.append('red')

    if 'trna' in reverted[i]:
        theFaceColors.append('w')
    else:
        theFaceColors.append(theEdgeColors[i])

    if 'rep.1' in reverted[i]:
        theMarkers.append('o')
    elif 'rep.2' in reverted[i]:
        theMarkers.append('s')
    else:
        theMarkers.append('^')
    
g.write('\n')

for i in range(len(genes)):
    g.write('{}\t'.format(synonyms[genes[i]]))
    for j in range(len(reverted)):
        value=expression[reverted[j]][genes[i]]
        g.write('{}\t'.format(value))

        x[j].append(value)
        
    g.write('\n')
    
g.close()

# 4. exploring the data
print('visualizing the data...')
original=numpy.array(x)

# 4.1. PCA of samples
print('running PCA...')
pcaMethod=sklearn.decomposition.PCA(n_components=5)
pcaObject=pcaMethod.fit(original)
new=pcaObject.transform(original)
explainedVar=pcaObject.explained_variance_ratio_
print('cumsum explained variance...')
print(numpy.cumsum(explainedVar))

for i in range(len(new)):
    matplotlib.pyplot.scatter(new[i,0],new[i,1],c=theFaceColors[i],marker=theMarkers[i],s=60,edgecolors=theEdgeColors[i])

matplotlib.pyplot.xlabel('PCA 1 ({0:.2f} var)'.format(explainedVar[0]))
matplotlib.pyplot.ylabel('PCA 2 ({0:.2f} var)'.format(explainedVar[1]))
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('figure.pca.png')
matplotlib.pyplot.clf()
print()

# 4.2. t-SNE of samples
print('running t-SNE...')
tSNE_Method=sklearn.manifold.TSNE(method='exact',verbose=1,init='pca')
tSNE_Object=tSNE_Method.fit(original)
new=tSNE_Object.fit_transform(original)

for i in range(len(new)):
    matplotlib.pyplot.scatter(new[i,0],new[i,1],c=theFaceColors[i],marker=theMarkers[i],s=60,edgecolors=theEdgeColors[i])
matplotlib.pyplot.savefig('figure.tSNE.png')
matplotlib.pyplot.clf()
print()
