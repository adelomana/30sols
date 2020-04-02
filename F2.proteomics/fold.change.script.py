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
    csvFiles=[element for element in allFiles if '.csv' in element and '._' not in element and 'rbf' in element]

    for csvFile in csvFiles:
        print(csvFile)
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
                geneName = vector[0]

                if geneName in riboPtNames:
                    
                    if geneName not in geneNames:
                        geneNames.append(geneName)
                        
                    a=float(vector[2])
                    b=float(vector[6])
                    c=float(vector[10])
                    data[condition][replicate]['tp2vs1'][geneName]=a
                    data[condition][replicate]['tp3vs1'][geneName]=b
                    data[condition][replicate]['tp4vs1'][geneName]=c
                                
    # sort
    conditions.sort()
    replicates.sort()
    geneNames.sort()

    # missing ones
    for element in riboPtNames:
        if element not in geneNames:
            print('lost {}'.format(element))
    print(geneNames, len(geneNames))
    
    return data, geneNames, conditions, replicates, timepoints

def riboPtNamesReader():

    '''
    This function reads the ribosomal protein names.
    '''

    riboPtNames=[]
    with open(ribosomalProteinsFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            riboPtNames.append(vector[1])

    print(riboPtNames)
    riboPtNames=['VNG1139Gm' if element == 'VNG1139G' else element for element in riboPtNames]
    print(riboPtNames)
            
    return riboPtNames

def violinAnalysis():

    '''
    This function builds violin plots for protein fold-changes.
    '''

    foldChangesFCL=[]; timeStampsFCL=[]; foldChangesREF=[]; timeStampsREF=[]
    violinStructure={}; violinNames={}

    L14pREF=[]

    # define a file to save plotting value
    f=open(plotValuesFile,'w')
    
    for fraction in proteinConditions:
        f.write('{}\n'.format(fraction))
        timeStamp=0
        for timepoint in proteinTimepoints:
            f.write('{}\n'.format(timepoint))
            timeStamp=timeStamp+1
            
            label=fraction+timepoint
            if label not in violinStructure:
                violinStructure[label]=[]; violinNames[label]=[]
                
            for name in riboPtNames:
                f.write('{}\t'.format(name))
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
                        f.write('{}\n'.format(average))
                        violinStructure[label].append(average); violinNames[label]=name

                        if fraction == 'lysate':
                            foldChangesFCL.append(average); timeStampsFCL.append(timeStamp)
                        else:
                            ###
                            conversion = {}
                            conversion['VNG1158G'] = 'S28E'
                            conversion['VNG0551G'] = 'L44E'
                            conversion['VNG1159G'] = 'L24E'
                            conversion['VNG1132G'] = 'S13'   
                            conversion['VNG0433C'] = 'S10-like'
                            if name in conversion.keys():
                                pt = conversion[name]
                                locx = timeStamp - 1 + 0.1
                                locy = average
                                matplotlib.pyplot.text(locx, locy, pt, fontsize=8)
                            ###
                            foldChangesREF.append(average); timeStampsREF.append(timeStamp)
                            if name == 'VNG1701G':
                               L14pREF.append([average,timeStamp])
                        if value > 0:
                            print('FC > 0: ',fraction,timepoint,name,value)
                    else:
                        f.write('HCV\n')
                        print('\t\t loosing {} {} {} {} {} for low precision'.format(fraction,timepoint,name,values,rsem))
                else:
                    f.write('ND\n')
                    print('\t\t no appropriate data for {} {} {}({}): {}'.format(timepoint,fraction,name,name,values))
            localList=violinStructure[label]
            cv=numpy.std(localList)/numpy.mean(localList)
            print('{}; n = {}; median = {}; cv = {}'.format(label,len(localList),numpy.median(localList),cv))

    f.close()
    
    # create a dataframe for plotting with seaborn
    foldChangeData=list(zip(timeStampsREF,foldChangesREF))
    df=pandas.DataFrame(data=foldChangeData,columns=['Time points','Fold change'])
    
    # plot violin and swarm plots with seaborn
    ax=seaborn.violinplot(x='Time points',y='Fold change',data=df,inner=None,linewidth=0,palette=['orange','green','blue'])
    matplotlib.pyplot.setp(ax.collections, alpha=0.5)
    ax=seaborn.swarmplot(x='Time points',y='Fold change',data=df,size=3,zorder=1,palette=['orange','green','blue'])

    # aesthetics
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Protein rel. abundance (log$_2$ FC)')

    # final figure closing
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Protein rel. abundance (log$_2$ FC)')

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
proteomicsDataFolder='/Users/alomana/scratch/tempo/'

plotValuesFile='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/violin/foldChanges.txt'

# 1. read data
print('reading data...')
riboPtNames=riboPtNamesReader()
print(riboPtNames, len(riboPtNames))

# 1.3. define proteomics data
print('\t reading proteomics data...')
proteinAbundance,proteinNames,proteinConditions,proteinReplicates,proteinTimepoints=proteomicsDataReader()

# 2. analysis
print('running analysis...')

# 2.1. violin plot of ribosomal proteins along the timepoints
print('\t building general trends figure...')
violinAnalysis()
