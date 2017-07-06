import sys,numpy
import matplotlib,matplotlib.pyplot

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42


# 0. define user variables
dataFile='/Volumes/omics4tb/alomana/projects/TLR/data/GREs/ribosomal_corems_motif_counts.txt'



# 1. reading data
data={} # data[coremName]=[(GRE,frequency),(GRE,frequency),...]

allCoremNames=[]
allGREnames=[]

with open(dataFile,'r') as f:
    next(f)
    for line in f:
        print(line)
        vector=line.split('\t')

        coremName=vector[0].replace('"','')
        motifName=vector[1].replace('"','')
        frequencyString=vector[2].replace('"','')
        frequency=int(frequencyString)

        print(coremName,motifName,frequency)

        if coremName not in data.keys():
            data[coremName]={}

        if motifName not in data[coremName].keys():
            data[coremName][motifName]=frequency
        
        if coremName not in allCoremNames:
            allCoremNames.append(coremName)
        if motifName not in allGREnames:
            allGREnames.append(motifName)

# 2. plotting heatmap of gres and corems
coremInt=[int(element) for element in allCoremNames]
coremInt.sort()
sortedCoremNames=[str(element) for element in coremInt]

GREnamesInt=[int(element.split('_')[1]) for element in allGREnames]
GREnamesInt.sort()
sortedGREnames=['MOTC_{}'.format(element) for element in GREnamesInt]

print(sortedGREnames)

x=[]
for corem in sortedCoremNames:
    y=[]
    for gre in sortedGREnames:
        try:
            value=data[corem][gre]
        except:
            value=0
        y.append(value)
    x.append(y)
Z=numpy.array(x)

# plotting figure
figureName='figure.pdf'
    
matplotlib.pyplot.imshow(Z,interpolation='none',cmap='viridis')

cb=matplotlib.pyplot.colorbar(label='similarity',orientation='horizontal',fraction=0.025) 
cb.ax.tick_params(labelsize=10)
    
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()    



### for the entire set, rows will be corems, columns will be GREs
