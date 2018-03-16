###
### This script computes outliers using prediction intervals.
### More info, http://markthegraph.blogspot.com/2015/05/using-python-statsmodels-for-ols-linear.html
###

#
#make sure you can color points based on PI border
#then run all data
#get lines and plot ribo-pt genes with color accordingly


import numpy,sys
import statsmodels,statsmodels.api,statsmodels.sandbox,statsmodels.sandbox.regression,statsmodels.sandbox.regression.predstd

import scipy,scipy.stats

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    this function creates a dictionary for expression values as
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
            geneName=vector[0].replace('_','')
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

def regressionAnalysis(x,y):

    '''
    This function performs regression analysis based on:
    http://markthegraph.blogspot.com/2015/05/using-python-statsmodels-for-ols-linear.html
    '''

    # f.1. build regression model
    xc=statsmodels.api.add_constant(x) # constant intercept term
    model=statsmodels.api.OLS(y,xc)
    fitted=model.fit()

    #print(fitted.params)     # the estimated parameters for the regression line
    #print(fitted.summary())  # summary statistics for the regression


    # f.2. interpolate model
    a=x.min()
    b=x.max()
    x_pred=numpy.linspace(a,b,100)
    x_pred2=statsmodels.api.add_constant(x_pred)
    y_pred=fitted.predict(x_pred2)
    regressionLine=[x_pred,y_pred]

    # f.3. compute CI
    y_hat=fitted.predict(xc)
    y_err=y-y_hat
    mean_x=xc.T[1].mean()
    n=len(xc)
    dof=n-fitted.df_model-1
    t=scipy.stats.t.ppf(1-0.025,df=dof)
    s_err=numpy.sum(numpy.power(y_err, 2))
    conf = t * numpy.sqrt((s_err/(n-2))*(1.0/n + (numpy.power((x_pred-mean_x),2) / ((numpy.sum(numpy.power(x_pred,2))) - n*(numpy.power(mean_x,2))))))
    upper=y_pred+abs(conf)
    lower=y_pred-abs(conf)
    CI=[upper,lower]

    # f.4. compute PI
    sdevP,lowerP,upperP=statsmodels.sandbox.regression.predstd.wls_prediction_std(fitted,exog=x_pred2,alpha=0.05)
    PI=[upperP,lowerP]

    return regressionLine,CI,PI

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

### MAIN

# 0. user defined variables
expressionDataFile='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/normalizedCounts.all.csv'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'

# 1. read data
print('reading data...')
riboPtNames=riboPtNamesReader()
expression,sampleTypes,geneNames,timepoints,replicates=expressionReader()

# 2. compute and plot regression and prediction intervals
print('computing prediction intervals...')

# 2.1. defining x and y data points
fcx=[]; fcy=[]
ribox=[]; riboy=[]
orderedRiboNames=[]

timepointLate=timepoints[-1]
timepointEarly=timepoints[0]

for geneName in geneNames:
        
        # compute averages for RNA-seq late
        x=numpy.mean([expression['trna'][geneName][timepointLate][replicate] for replicate in replicates])

        # compute averages for RNA-seq early
        y=numpy.mean([expression['trna'][geneName][timepointEarly][replicate] for replicate in replicates])
        
        # compute averages for Ribo-seq late
        z=numpy.mean([expression['rbf'][geneName][timepointLate][replicate] for replicate in replicates])

        # compute averages for Ribo-seq early
        w=numpy.mean([expression['rbf'][geneName][timepointEarly][replicate] for replicate in replicates])
        
        # compute fold-changes
        fcx.append(x-y)
        fcy.append(z-w)

        # add ribo-pt genes
        if geneName in riboPtNames:
            ribox.append(x-y)
            riboy.append(z-w)
            orderedRiboNames.append(geneName)

        # just checking some hard-coded values
        if geneName in ['VNG0433C','VNG1701G','VNG1133G']:
            print(geneName,x-y,z-w)        

# 2.2. compute regression line and intervals
regressionLine,CI,PI=regressionAnalysis(numpy.array(ribox),numpy.array(riboy))

matplotlib.pyplot.plot(regressionLine[0],regressionLine[1],color='black',lw=2)

matplotlib.pyplot.fill_between(regressionLine[0],CI[1],CI[0],color='black',alpha=0.4,lw=0)
matplotlib.pyplot.fill_between(regressionLine[0],PI[1],PI[0],color='black',alpha=0.1,lw=0)

# 3. define colors of scatter plot ribo-pt genes based on prediction intervals
matplotlib.pyplot.plot(ribox,riboy,'o',alpha=0.1,mew=0,ms=8,color='black')

matplotlib.pyplot.plot(ribox[3],riboy[3],'o',alpha=1,mew=0,ms=8,color='red')
matplotlib.pyplot.plot(ribox[35],riboy[35],'o',alpha=1,mew=0,ms=8,color='red')
matplotlib.pyplot.plot(ribox[13],riboy[13],'o',alpha=1,mew=0,ms=8,color='blue')

for outlier in [3,35,13]:
    print(orderedRiboNames[outlier],ribox[outlier],riboy[outlier])

# 3.1. close figure
matplotlib.pyplot.xlim([-5,0.5])
matplotlib.pyplot.ylim([-5,0.5])

matplotlib.pyplot.grid(alpha=0.25, ls=':')

matplotlib.pyplot.plot([-5,0.5],[0,0],':',color='black',alpha=0.5)
matplotlib.pyplot.plot([0,0],[-5,0.5],':',color='black',alpha=0.5)
matplotlib.pyplot.plot([-5,0.5],[-5,0.5],':',color='black',alpha=0.5)
    
matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')

figureName='stats.pdf'
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()
