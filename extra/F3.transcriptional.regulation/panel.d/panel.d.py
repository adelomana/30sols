###
### This script computes outliers using prediction intervals.
### More info, http://markthegraph.blogspot.com/2015/05/using-python-statsmodels-for-ols-linear.html
###

import numpy,sys
import statsmodels,statsmodels.api,statsmodels.sandbox,statsmodels.sandbox.regression,statsmodels.sandbox.regression.predstd

import scipy,scipy.stats

import matplotlib,matplotlib.pyplot
matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def expressionReader():

    '''
    this function creates a dictionary for expression values as
    expression[trna/rbf][ribo-pt gene name]=foldchange
    '''

    expression={}

    for sampleType in sampleTypes:
        expressionDataFile=expressionDataDir+'significance.{}.condition_tp.4_vs_tp.1.csv'.format(sampleType)
        expression[sampleType]={}

        with open(expressionDataFile,'r') as f:
            next(f)
            for line in f:
                vector=line.split(',')
                
                geneName=vector[0].replace('"','')
                log2FC=float(vector[2])

                if geneName in riboPtNames:
                    expression[sampleType][geneName]=log2FC

    return expression

def regressionAnalysis(x,y):

    '''
    This function performs regression analysis based on:
    http://markthegraph.blogspot.com/2015/05/using-python-statsmodels-for-ols-linear.html
    '''

    # f.1. build regression model
    xc=statsmodels.api.add_constant(x) # constant intercept term
    model=statsmodels.api.OLS(y,xc)
    fitted=model.fit()

    print(fitted.params)     # the estimated parameters for the regression line
    print(fitted.summary())  # summary statistics for the regression

    # f.2. interpolate model
    a=x.min()
    b=x.max()
    x_pred=numpy.linspace(a,b,100)
    x_pred2=statsmodels.api.add_constant(x_pred)
    y_pred=fitted.predict(x_pred2)
    regressionLine=[x_pred,y_pred]

    print()

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

    riboPtNames.sort()
            
    return riboPtNames

### MAIN

# 0. user defined variables
expressionDataDir='/Volumes/omics4tb/alomana/projects/TLR/data/DESeq2/'
ribosomalProteinsFile='/Volumes/omics4tb/alomana/projects/TLR/data/ribosomalGeneNames.txt'
sampleTypes=['trna','rbf']

# 1. read data
print('reading data...')
riboPtNames=riboPtNamesReader()
expression=expressionReader()

# 2. compute and plot regression and prediction intervals
print('computing prediction intervals...')

# 2.1. defining x and y data points
ribox=[]; riboy=[]
for geneName in riboPtNames:

        # compute fold-changes
        ribox.append(expression['trna'][geneName])
        riboy.append(expression['rbf'][geneName])
        
# check
if len(ribox) != len(riboPtNames):
    print(len(ribox),len(riboPtNames))
    print('Mismatch. Exiting...')
    sys.exit()

# 2.2. compute regression line and intervals
regressionLine,CI,PI=regressionAnalysis(numpy.array(ribox),numpy.array(riboy))

matplotlib.pyplot.plot(regressionLine[0],regressionLine[1],color='black',lw=2)

matplotlib.pyplot.fill_between(regressionLine[0],CI[1],CI[0],color='black',alpha=0.4,lw=0)
matplotlib.pyplot.fill_between(regressionLine[0],PI[1],PI[0],color='black',alpha=0.1,lw=0)

# 3. define colors of scatter plot ribo-pt genes based on prediction intervals
for i in range(len(riboPtNames)):

    # find the regression point closer to ribox[i] and define the approximate limit from PI
    distances=[abs(position-ribox[i]) for position in regressionLine[0]]
    index=distances.index(min(distances))
    limitTop=PI[0][index]
    limitBottom=PI[1][index]

    # check if values are above or below PI
    if riboy[i] > limitTop:
        theColor='red'; theAlpha=1
        expected=
        print('\t {} detected as upper outlier at x={}; y={}.'.format(riboPtNames[i],ribox[i],riboy[i]))
        
    elif riboy[i] < limitBottom:
        theColor='blue'; theAlpha=1
        print('\t {} detected as bottom outlier at x={}; y={}.'.format(riboPtNames[i],ribox[i],riboy[i]))
    else:
        theColor='black'; theAlpha=0.1
        print('black\t{}\t{}\t{}\t{}'.format(riboPtNames[i],ribox[i],riboy[i],abs(riboy[i]-limitBottom)))

    # plot the point
    matplotlib.pyplot.plot(ribox[i],riboy[i],'o',alpha=theAlpha,mew=0,ms=8,color=theColor)

# 3.1. close figure
matplotlib.pyplot.xlim([-6.5,0.5])
matplotlib.pyplot.ylim([-6.5,0.5])

matplotlib.pyplot.grid(alpha=0.25, ls=':')

matplotlib.pyplot.plot([-6,0.5],[-6,0.5],':',color='black',alpha=0.5)
    
matplotlib.pyplot.xlabel('RNA-seq, log$_2$ FC')
matplotlib.pyplot.ylabel('Ribo-seq, log$_2$ FC')

matplotlib.pyplot.xticks([-6,-5,-4,-3,-2,-1,0])

figureName='stats.pdf'
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig(figureName)
matplotlib.pyplot.clf()
