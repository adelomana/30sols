###
### This script investigates if protein changes can be explained by ribosomal footprints, i.e., translational regulation.
###

import numpy,sys,os,pandas,seaborn
import scipy,scipy.stats

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    This function creates a dictionary for expression values as
    expression[trna/rbf][ribo-pt gene name][timepoint][replicate]=value
    '''

    expression={}
    
    sampleTypes=[]
    geneNames=[]
    timepoints=[]
    replicates=[]

    with open(expressionDataFile,'r') as f:

        firstLine=f.readline()
        header=firstLine.split(',')
        sampleNames=header[1:]
        sampleNames[-1]=sampleNames[-1].replace('\n','')

        for line in f:
            vector=line.split(',')

            # geneName
            geneName=vector[0]
            if geneName not in geneNames:
                geneNames.append(geneName)

            for i in range(len(sampleNames)):

                # sampleType
                sampleType=sampleNames[i].split('.')[0]
                if sampleType not in sampleTypes:
                    sampleTypes.append(sampleType)

                # timepoint
                timepoint='tp.{}'.format(int(sampleNames[i].split('.')[-1]))
                if timepoint not in timepoints:
                    timepoints.append(timepoint)

                # replicate
                replicate='rep.{}'.format(int(sampleNames[i].split('rep.')[1][0]))
                if replicate not in replicates:
                    replicates.append(replicate)

                # value
                value=float(vector[i+1])

                # make sure keys exist
                if sampleType not in expression.keys():
                    expression[sampleType]={}
                if geneName not in expression[sampleType].keys():
                    expression[sampleType][geneName]={}
                if timepoint not in expression[sampleType][geneName].keys():
                    expression[sampleType][geneName][timepoint]={}

                expression[sampleType][geneName][timepoint][replicate]=value

    # sort variables
    sampleTypes.sort()
    geneNames.sort()
    timepoints.sort()
    replicates.sort()

    return expression,sampleTypes,geneNames,timepoints,replicates

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

    countsDict={}
    noInfoSet=[]
    foldChanges=[]; timeStamps=[]
    timeStamp=1
    for fraction in proteinConditions:
        for timepoint in proteinTimepoints:
            count=0
            for name in riboPtNames:
            #for name in proteinNames:
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
                        foldChanges.append(value); timeStamps.append(timeStamp)
                        if value > 0:
                            print('values',fraction,timepoint,name,value)
                        count=count+1
                        print(synonymsReverseMapping[name],values,average,sem,count)
                    else:
                        print('\t\t loosing {} {} {} {} {} for low precision'.format(fraction,timepoint,synonymsReverseMapping[name],values,rsem))
                else:
                    print('\t\t no appropriate data for {} {} {}({}): {}'.format(timepoint,fraction,name,synonymsReverseMapping[name],values))
                    if name not in noInfoSet:
                        noInfoSet.append(name)
                        
            print('\t found {} proteins in {} {}'.format(count,fraction,timepoint))
            countsDict[timeStamp]=count
            timeStamp=timeStamp+2
        timeStamp=2
        
    # create a dataframe for plotting with seaborn
    foldChangeData=list(zip(timeStamps,foldChanges))
    df=pandas.DataFrame(data=foldChangeData,columns=['Time points','Fold change'])

    print(countsDict)
    
    # plot violin and swarm plots with seaborn
    ax=seaborn.violinplot(x='Time points',y='Fold change',data=df,inner=None,linewidth=0,palette=['orange','orange','green','green','blue','blue'])
    matplotlib.pyplot.setp(ax.collections, alpha=0.5)
    #matplotlib.pyplot.setp(ax.collections, alpha=0.75)
    ax=seaborn.swarmplot(x='Time points',y='Fold change',data=df,size=3,zorder=1,palette=['orange','orange','green','green','blue','blue'])

    # final figure closing
    matplotlib.pyplot.grid(alpha=0.5, ls=':')
    matplotlib.pyplot.xlabel('Time point')
    matplotlib.pyplot.ylabel('Protein rel. abundance (log$_2$ FC)')
    matplotlib.pyplot.ylim([-6,6])
    generalTickLabels=['t2.FCL','t2.REF','t3.FCL','t3.REF','t4.FCL','t4.REF']
    specificTickLabels=[generalTickLabels[i]+'\nn={}'.format(countsDict[i+1]) for i in range(len(generalTickLabels))]
    matplotlib.pyplot.xticks([0,1,2,3,4,5],specificTickLabels)

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
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv'
proteomicsDataFolder='/Volumes/omics4tb/alomana/projects/TLR/data/proteomics/all/'
annotationFile='/Volumes/omics4tb/alomana/projects/TLR/data/genome/alo.build.NC002607.NC001869.NC002608.gff3'
scratchDir='/Volumes/omics4tb/alomana/scratch/'

# 1. read data
print('reading data...')

# 1.1. define synonyms
synonyms,synonymsReverseMapping=synonymsReader()

# 1.2. define transcriptome data
riboPtNames=riboPtNamesReader()
print('\t reading transcriptome data...')
expression,expressionSampleTypes,expressionGeneNames,expressionTimepoints,expressionReplicates=expressionReader()

# 1.3. define proteomics data
print('\t reading proteomics data...')
proteinAbundance,proteinNames,proteinConditions,proteinReplicates,proteinTimepoints=proteomicsDataReader()

# 3. analysis
print('running analysis...')

# 3.1. violin plot of ribosomal pt along the timepoints
print('\t building general trends figure...')
violinAnalysis()
        
# 3.2. scatter plot of FC_pt vs FC_mRNA. This figure would reveal how well mRNA explains pt changes.
print('\t building pt vs mRNA...')

# 2.1. empty figure calling to maintain sizes
matplotlib.pyplot.plot([0,0],[1,1],'ok')
matplotlib.pyplot.savefig('{}temp.pdf'.format(scratchDir))
matplotlib.pyplot.clf()

# make sure that you do rsem filtering here too for pt and transcripts.

timepointLate=expressionTimepoints[-1]
timepointEarly=expressionTimepoints[0]
ptTimepoint='tp4vs1'

for fraction in proteinConditions:
    for sampleType in expressionSampleTypes:

        print('\t\t working with ',fraction,sampleType)
        figureFileName='figure.all.{}.{}.pdf'.format(fraction,sampleType)
        x=[]; y=[]

        #for name in riboPtNames:
        for name in proteinNames:

            log2fcx=None; log2fcy=None
            
            try:
                # compute averages for transcript late
                a=numpy.mean([expression[sampleType][name][timepointLate][replicate] for replicate in expressionReplicates])
                # compute averages for transcript late
                b=numpy.mean([expression[sampleType][name][timepointEarly][replicate] for replicate in expressionReplicates])
                # compute log2 fold-change
                log2fcx=a-b
            except:
                print('\t\t\t no transcript data for {}'.format(name))

            values=[]
            for replicate in proteinReplicates:
                value=None
                try:
                    value=proteinAbundance[fraction][replicate]['tp4vs1'][name]
                except:
                    pass
                if value != None:
                    values.append(value)
                if len(values) >= 3:
                    average=numpy.median(values)
                    sem=numpy.std(values)/numpy.sqrt(len(values))
                    rsem=sem/numpy.mean(values)
                    if rsem < 0.3:
                        log2fcy=average
                    else:
                        print('\t\t\t loosing {} {} {} {} {} for low precision.'.format(fraction,sampleType,name,values,rsem))
                else:
                    print('\t\t\t not enough protein data for {} {}'.format(name,values))
            
            # append if valid values are found
            if log2fcy != None and log2fcx != None:
                x.append(log2fcx); y.append(log2fcy)

        # plot data
        matplotlib.pyplot.plot(x,y,'o',alpha=0.5,mew=0,ms=8,color='black')

        # compute regression
        slope,intercept,r_value,p_value,std_err=scipy.stats.linregress(x,y)
        print('linear regression')
        print('slope',slope)
        print('intercept',intercept)
        print('r_value',r_value)
        print('pvalue',p_value)
        print('std_err',std_err)
        print()

        # compute and plot linear regression line
        resolution=0.1
        newx=numpy.arange(min(x),max(x),resolution)
        newy=slope*newx+intercept
        #idx=numpy.where(newy>0)

        matplotlib.pyplot.plot(newx,newy,lw=4,color='red')

        description='R$^2$={:.2f}\np={:.2e}\na={:.2f}'.format(r_value**2,p_value,slope)
        matplotlib.pyplot.text(-8,-8,description)
        
        # labels
        matplotlib.pyplot.xlabel('Transcript ({}), log$_2$ FC'.format(sampleType))
        matplotlib.pyplot.ylabel('Protein ({}), log$_2$ FC'.format(fraction))
        matplotlib.pyplot.title('n={}'.format(len(x)))
        
        matplotlib.pyplot.xlim([-8.05,8.05])
        matplotlib.pyplot.ylim([-8.05,8.05])

        tickPositions=numpy.arange(-8,8+2,2)
        matplotlib.pyplot.xticks(tickPositions)
        matplotlib.pyplot.yticks(tickPositions)

        matplotlib.pyplot.grid(alpha=0.5, ls=':')
            
        # closing figure
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.axes().set_aspect('equal')
        matplotlib.pyplot.savefig(figureFileName)
        matplotlib.pyplot.clf()
 
sys.exit()
# x.2. scatter plot of FC_pt vs FC_RF. This figure would reveal how well RF explains pt changes.
print('\t building pt vs RF...')

# x.3 scatter plot of FC_pt vs FC_mRNA + FC_RF. This figure would reveal how well an equal linear model would explain pt changes.
print('\t building pt vs equal model...')

# x.4. scatter plot of model: FC_pt = w_1 FC_mRNA + w_2 (FC_RF - FC_RF_0).
# For this section I need to obtain FC_RF_0 from linear regression on FC_mRNA vs FC_RF, from all genes.
# w_1 is defined within the range [0,1].
# w_2 is defined within the range [-1,1].
print('\t building weighted model...')

# x.5. visualization of w_1 vs w_2 to distinguish TCR from TLR.
# I should have one for each time point - t0
print('\t visualizing regulatory weights...')

# x. a sanity check for the script is the capture of RF causing pausing on rps10-like.
