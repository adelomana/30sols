import os,sys,numpy
import sklearn,sklearn.decomposition,sklearn.manifold
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})

# 0. user defined variables
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/cufflinks/allSamples/isoforms.fpkm_table'

# 1. reading the data
x=[]
theEdgeColors=[]
theFaceColors=[]
theMarkers=[]

with open(expressionDataFile,'r') as f:

    header=f.readline()
    vectorHeader=header.split('\t')[1:]
    for sampleName in vectorHeader:
        x.append([])

        if 'tp.1' in sampleName:
            theEdgeColors.append('blue')
        elif 'tp.2' in sampleName:
            theEdgeColors.append('green')
        elif 'tp.3' in sampleName:
            theEdgeColors.append('orange')
        else:
            theEdgeColors.append('red')

        if 'trna' in sampleName:
            theFaceColors.append('w')
        else:
            theFaceColors.append(theEdgeColors[-1])

        if 'rep.1' in sampleName:
            theMarkers.append('o')
        elif 'rep.2' in sampleName:
            theMarkers.append('s')
        else:
            theMarkers.append('^')

    for line in f:
        vector=line.split('\t')[1:]
        for i in range(len(vector)):
            value=float(vector[i])
            x[i].append(value)

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
