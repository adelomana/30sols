###
### This script investigates if protein changes can be explained by ribosomal footprints, i.e., translational regulation.
###

import numpy,sys,os,pandas,seaborn
import scipy,scipy.stats

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def proteomicsDataReader():

    '''
    This function reads data available and outputs the defined dictionary.
    '''

    data={} # data[lysate/rbf][rep1/rep2/rep3][tp2vs1/tp3vs1/tp4vs1][geneName]=log2FC

    conditions=[]
    geneNames=[]
    timepoints=[]
    replicates=[]
    
    allFiles=os.listdir(proteomicsDataFolder)
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element]

    lostGenes=[]

    for csvFile in csvFiles:
        path=proteomicsDataFolder+csvFile

        brokenName=csvFile.split('.')
        condition=brokenName[0]
        replicate=brokenName[1]

        if condition not in conditions:
            conditions.append(condition)
        if replicate not in replicates:
            replicates.append(replicate)

        if condition not in data.keys():
            data[condition]={}
        if replicate not in data[condition].keys():
            data[condition][replicate]={}

        timepoints=['tp2vs1','tp3vs1','tp4vs1']
        for timepoint in timepoints:
            if timepoint not in data[condition][replicate].keys():
                data[condition][replicate][timepoint]={}

        with open(path,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')

                try:
                    geneName=synonyms[vector[0]]

                    if geneName not in geneNames:
                        geneNames.append(geneName)

                    a=float(vector[2])
                    b=float(vector[6])
                    c=float(vector[10])

                    data[condition][replicate]['tp2vs1'][geneName]=a
                    data[condition][replicate]['tp3vs1'][geneName]=b
                    data[condition][replicate]['tp4vs1'][geneName]=c
                    
                except:
                    if vector[0] not in lostGenes:
                        lostGenes.append(vector[0])
                
    # sort
    geneNames.sort()
    conditions.sort()
    replicates.sort()

    # print lost genes
    print('\t {} protein lost because of new annotation.'.format(len(lostGenes)))
    
    return data,geneNames,conditions,replicates,timepoints

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[0])
            
    return riboPtNames

def synonymsReader():

    '''
    This function reads the GFF3 file and returns a dictionary with synonyms between old and new locus names.
    '''

    synonyms={}
    with open(annotationFile,'r') as f:
        for line in f:
            vector=line.split('\t')
            if vector[0][0] != '#':
                info=vector[-1].replace('\n','')
                if 'old_locus_tag=' in info:
                    old=info.split('old_locus_tag=')[1].split(';')[0]
                    new=info.split('ID=')[1].split(';')[0]

                    if '%' in old:
                        olds=old.split('%2C')
                        for element in olds:
                            synonyms[element]=new
                    else:
                        synonyms[old]=new

    # reverting mapping
    synonymsReverseMapping={v: k for k, v in synonyms.items()}
                
    return synonyms,synonymsReverseMapping

def violinAnalysis():

    '''
    This function builds violin plots for protein fold-changes.
    '''

    foldChangesFCL=[]; timeStampsFCL=[]; foldChangesREF=[]; timeStampsREF=[]
    violinStructure={}; violinNames={}
    
    for fraction in proteinConditions:
        timeStamp=0
        for timepoint in proteinTimepoints:
            timeStamp=timeStamp+1

            label=fraction+timepoint
            if label not in violinStructure:
                violinStructure[label]=[]; violinNames[label]=[]
                
            for name in riboPtNames:
                values=[]
                for replicate in proteinReplicates:
                    value=None
                    try:
                        value=proteinAbundance[fraction][replicate][timepoint][name]
                    except:
                        pass
                    if value != None:
                        values.append(value)
                if len(values) >= 3 :
                    average=numpy.median(values)
                    sem=numpy.std(values)/numpy.sqrt(len(values))
                    rsem=sem/numpy.mean(values)
                    if rsem < 0.3:
                        
                        violinStructure[label].append(value); violinNames[label]=name

                        if fraction == 'lysate':
                            foldChangesFCL.append(value); timeStampsFCL.append(timeStamp)
                        else:
                            foldChangesREF.append(value); timeStampsREF.append(timeStamp)
                        if value > 0:
                            print('FC > 0: ',fraction,timepoint,name,value)
                    else:
                        print('\t\t loosing {} {} {} {} {} for low precision'.format(fraction,timepoint,synonymsReverseMapping[name],values,rsem))
                else:
                    print('\t\t no appropriate data for {} {} {}({}): {}'.format(timepoint,fraction,name,synonymsReverseMapping[name],values))
            localList=violinStructure[label]
            cv=numpy.std(localList)/numpy.mean(localList)
            print('{}; n = {}; median = {}; cv = {}'.format(label,len(localList),numpy.median(localList),cv))

            # define outliers
            #outliersDefiner(violinStructure[label],violinNames[label])

    matplotlib.pyplot.subplot(121)
    yLimits=[-4.5,1.5]
    
    # create a dataframe for plotting with seaborn
    foldChangeData=list(zip(timeStampsFCL,foldChangesFCL))
    df=pandas.DataFrame(data=foldChangeData,columns=['Time points','Fold change'])
    
    # plot violin and swarm plots with seaborn
    ax=seaborn.violinplot(x='Time points',y='Fold change',data=df,inner=None,linewidth=0,palette=['orange','green','blue'])
    matplotlib.pyplot.setp(ax.collections, alpha=0.5)
    ax=seaborn.swarmplot(x='Time points',y='Fold change',data=df,size=3,zorder=1,palette=['orange','green','blue'])

    # aesthetics
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Protein rel. abundance (log$_2$ FC)')
    matplotlib.pyplot.ylim(yLimits)

    matplotlib.pyplot.subplot(122)

    # create a dataframe for plotting with seaborn
    foldChangeData=list(zip(timeStampsREF,foldChangesREF))
    df=pandas.DataFrame(data=foldChangeData,columns=['Time points','Fold change'])
    
    # plot violin and swarm plots with seaborn
    ax=seaborn.violinplot(x='Time points',y='Fold change',data=df,inner=None,linewidth=0,palette=['orange','green','blue'])
    matplotlib.pyplot.setp(ax.collections, alpha=0.5)
    ax=seaborn.swarmplot(x='Time points',y='Fold change',data=df,size=3,zorder=1,palette=['orange','green','blue'])

    # final figure closing
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Protein rel. abundance (log$_2$ FC)')
    matplotlib.pyplot.ylim(yLimits)
    generalTickLabels=['t2.FCL','t2.REF','t3.FCL','t3.REF','t4.FCL','t4.REF']

    #specificTickLabels=[generalTickLabels[i]+'\nn={}'.format(countsDict[i+1]) for i in range(len(generalTickLabels))]
    #matplotlib.pyplot.xticks([0,1,2,3,4,5],specificTickLabels)

    figureName='figure.violin.ribo.pdf'
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()    

    return None

###
### M A I N
###

# 0. user defined variables.
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'

# 1. read data
print('reading data...')

# 1.1. define synonyms
synonyms,synonymsReverseMapping=synonymsReader()
riboPtNames=riboPtNamesReader()

# 1.3. define proteomics data
print('\t reading proteomics data...')
proteinAbundance,proteinNames,proteinConditions,proteinReplicates,proteinTimepoints=proteomicsDataReader()

# 2. analysis
print('running analysis...')

# 2.1. violin plot of ribosomal proteins along the timepoints
print('\t building general trends figure...')
violinAnalysis()
        
