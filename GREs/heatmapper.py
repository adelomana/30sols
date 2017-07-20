import sys,numpy,scipy
import matplotlib,matplotlib.pyplot
import sklearn,sklearn.cluster,sklearn.datasets

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42


# 0. define user variables
print('initializing...')
dataFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.txt'
jarFile=''
initializingSetSize=int(1e5)
scoreThreshold=0.026
goodResults=3

# 1. reading data
print('reading data...')
data={} # data[coremName]=[(GRE,frequency),(GRE,frequency),...]

allCoremNames=[]
allGREnames=[]

with open(dataFile,'r') as f:
    next(f)
    for line in f:
        vector=line.split('\t')

        coremName=vector[0].replace('"','')
        motifName=vector[1].replace('"','')
        frequencyString=vector[2].replace('"','')
        frequency=int(frequencyString)
        
        if coremName not in data.keys():
            data[coremName]={}

        if motifName not in data[coremName].keys():
            data[coremName][motifName]=frequency
        
        if coremName not in allCoremNames:
            allCoremNames.append(coremName)
        if motifName not in allGREnames:
            allGREnames.append(motifName)

# 2. arranging data
print('preparing data...')
coremInt=[int(element) for element in allCoremNames]
coremInt.sort()
sortedCoremNames=[str(element) for element in coremInt]

GREnamesInt=[int(element.split('_')[1]) for element in allGREnames]
GREnamesInt.sort()
sortedGREnames=['MOTC_{}'.format(element) for element in GREnamesInt]

x=[]
for corem in sortedCoremNames:
    y=[]
    for gre in sortedGREnames:
        try:
            value=data[corem][gre]
        except:
            value=0
        if value > 50:
            value=50
        y.append(value)
    x.append(y)
Z=numpy.array(x)

###
numpy.savetxt("foo.csv",Z,delimiter=",")
sys.exit() 
###

# 3. clustering
print('performing biclustering...')

# generating a random shuffled comparable data set
data, rows, columns = sklearn.datasets.make_biclusters(shape=Z.shape,n_clusters=3,shuffle=True)
data, row_idx, col_idx = sklearn.datasets.samples_generator._shuffle(data)

currentBestScore=1
iteration=0
stopping=False
successes=0

while stopping == False:

    iteration=iteration+1
    model=sklearn.cluster.bicluster.SpectralBiclustering(svd_method='arpack',n_init=initializingSetSize,n_jobs=-1) 
    model.fit(Z)
    clusteredData = Z[numpy.argsort(model.row_labels_)]
    clusteredData = Z[:, numpy.argsort(model.column_labels_)]

    score = sklearn.metrics.consensus_score(model.biclusters_,(rows[:, row_idx], columns[:, col_idx]))

    print('\t iteration {}, score {}...'.format(iteration,score))

    if score < scoreThreshold:
        successes=successes+1
        print('\t plotting good result ({}/{})...'.format(successes,goodResults))
        # plotting figure
        label='i.{}.s.{:.5f}'.format(iteration,score)
        figureName='figures/{}.png'.format(label)
        matplotlib.pyplot.matshow(clusteredData,cmap=matplotlib.pyplot.cm.Blues)
        matplotlib.pyplot.xlabel('GREs')
        matplotlib.pyplot.ylabel('corems')
        matplotlib.pyplot.title(label)
        matplotlib.pyplot.savefig(figureName)
        matplotlib.pyplot.clf()        

    if successes == goodResults:
        print('... all done.')
        stopping=True


#get at least 10 <0.3 in niint of 10 also other 10 for hiher levels.
### for the entire set, rows will be corems, columns will be GREs
