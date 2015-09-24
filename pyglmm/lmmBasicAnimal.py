__author__ = 'Igor Siveroni'

import math

import numpy as np
from pandas import Series,DataFrame
import pandas as pd

from scipy.optimize import minimize, minimize_scalar

import phylotree
#reload(phylotree)
from phylotree import PhyloTree

# The basic Animal model uses one random effect: the phylogeny.
# Number of response variables = 1 (even though we are passing a list as argument)
class LMMBasicAnimal:
    # Constructor
    # Only one trait allowed for basic model
    def __init__(self,ylist=[],plist=[]):
        if isinstance(ylist,(list)):
            self.num_traits = len(ylist)
            if (self.num_traits > 1):
                print "Only one trait allowed - taking first in list"
            self.ylist = ylist
        else:
            if isinstance(ylist,(basestring)):
                self.num_traits = 1
                self.ylist = [ ylist ]
            else:
                raise Exception("Incorrect response variable name(s)")
        self.num_fixed = len(plist)
        self.plist = plist
        # design matrices and rest of data to be provided from a file or 'by hand'
        # initialize to None
        print 'Model: ' + str(self.ylist) + " ~ " + str(self.plist)

    # sep: 'TAB' , 'SPACE'
    def loadFromFile(self,dataFile,treeFile,taxa_column='taxa',sep='TAB'):
        # read from csv
        # url = 'http://vincentarelbundock.github.io/Rdatasets/csv/HistData/Guerry.csv'
        # dat = pd.read_csv(url)
        # read from table
        if (sep not in ['TAB','SPACE']):
             raise Exception("unknown data file separator type: "+sep)
        try:
            if (sep=='TAB'):
                data = pd.read_table(dataFile)
            else:
                if (sep=='SPACE'):
                    data = pd.read_table(dataFile,delim_whitespace=True)
        except Exception as exc:
            raise Exception("Couldn't parse data file")

        data = data.dropna()  # subset = ['col1'col2']
        # May need to duplicate column if it's needed for data manipulation
        if (taxa_column not in data.columns):
            raise Exception("Taxa column "+taxa_column+" not present in dataframe")
        print "Reindexing dataframe with column: "+taxa_column
        data = data.set_index(taxa_column)
        dataIndex = data.index

        tree = PhyloTree(treeFile)
        taxaIndex  = pd.Index(tree.tipNames)
        # Sanity check
        if (tree.num_tips != taxaIndex.size):
            raise Exception("size of array of tip names and number of tips disagree")
        # Compare taxa from data with tip names - there must be enough data
        missing = taxaIndex.difference(dataIndex)
        if (missing.size > 0):
            raise Exception("Taxa mismatch: Not enough data matching tree tips. Missing: "+str(missing.values))
        else:
            if ((dataIndex.difference(taxaIndex)).size > 0):
                print "Extra rows (taxa) in data: deleting extra rows"
        # Readjust dataframe so it matches tip names
        data = data.reindex(taxaIndex)

        # ********* Check model variables ********
        if (self.ylist[0] not in data.columns):
            raise Exception("Variable " + str(self.ylist[0]) + " not present in dataframe")
        for p in self.plist:
            if (p not in data.columns):
                raise Exception ("Predictor "+p+ " not a column in dataframe")
        self.tree = tree
        self.lambdaV = None
        self.oldLambda = None
        # LMM related
        self.G = None
        self.R = None
        self.Z = None
        self.b = None
        self.sigmag = None
        self.sigmar = None

        # prepare general design matrix
        self.prepareMatrices(data)


    def prepareMatrices(self,all_data):
        # Reindex dataframe so it only contains model data
        data = all_data.reindex(columns=self.ylist+self.plist)
        yvar = self.ylist[0]
        self.num_predictors = len(self.plist)
        self.ntips = self.tree.num_tips
        self.y = np.matrix(data[yvar]).reshape(self.ntips,1)
        self.X = np.matrix(np.zeros((self.ntips,self.num_predictors+1)))
        self.X[:,0] = 1
        for idx,cname in enumerate(self.plist):
            self.X[:,1+idx] = data[cname].values.reshape(self.ntips,1)
        self.tree.normalize()
        self.G = self.tree.getV()
        self.R = np.asmatrix(np.identity(self.ntips))
        self.Z = np.asmatrix(np.identity(self.ntips))

    def loadFromMatrices(self,y,X,Z,G,R):
        self.y = y
        (ny,m) = np.shape(y)
        if (m != 1):
            raise Exception("y must be a column vector")
        (n,m) = np.shape(X)
        if (n == ny):
            self.ntips = ny
            self.num_predictors = m
            self.X = X
        else:
            raise Exception("Incorrect X dimensions")
        self.Z = Z
        self.G = G
        self.R = R

    def fitAll(self):
        #print self.X
        #print self.y
        self.lambda_list = []
        self.lh_list = []
        self.fitGLS()
        print "-----Lambda OPtimization ----"
        self.fitGLSLambda()


    def _fitGLSWorker(self,lamb=None,sigma=None):
        # Need to estimate beta and sigma
        # beta = np.matrix(np.zeros((self.num_predictors+1,1)))   -- no need to allocate
        self.lamb=lamb
        msg = "lambda: "+str(lamb)
        V = self.G
        print "V " +str(V)
        if (lamb is not None):
            # for some reason, lamb turns into a matrix
            lamb = self.fixLambda(lamb)
            Vtmp = self.G
            V = lamb * Vtmp
            np.fill_diagonal(V,Vtmp.diagonal())
        (s,logdetV) = np.linalg.slogdet(V)
        invV = np.linalg.inv(V)
        # estimator beta = (X' invV X)^{-1} (X' invV y)
        left = np.linalg.inv(self.X.T * invV * self.X)
        right = self.X.T * invV * self.y
        self.beta =  left * right
        e = self.y - self.X * self.beta
        if (sigma is None): # Calculate sigma that optimizes likelihood
            sigma = (e.T * invV  *e)/(self.ntips-1)
            if (sigma < 0):
                raise Exception("Sigma negative")
        self.sigma = sigma
        # print "Cholesky"
        #np.linalg.cholesky(invV)
        # print "finished cholesky"
        self.lhval = self.lh(self.beta,sigma,invV,logdetV)
        msg += " newlambda: "+str(lamb)+" lh:"+str(self.lhval)
        print msg
        # Below, an example on how to optimize for sigma - it gives the same result as the formula
        # lh1 = lambda(s): -self.lh(beta,s,invV,logdetV)
        # res = minimize_scalar(lh1,method='Golden',bracket=(0.5,1.5))
        # newsigma = res.x
        return self.lhval

    def fixLambda(self,lamb):
        lamb = float(lamb)
        if (lamb < 0.0):
            lamb = -lamb
        if (lamb > 1.0):
            intlamb = int(lamb)
            if ((intlamb%2)==1):
                lamb = 1-(lamb-intlamb)
            else:
                lamb = lamb-intlamb
        return lamb

    # Computes log-likelihood of MVN
    # This shouldn't be part of the class(?)
    # Move this somewhere else - this is mvn lh
    def lh(self,beta,sigma,invV,logdetV):
        # (-1/2)[log(2PI) + log(det(V)) + (y-Xb)' invV (y-Xb)]
        # recall realV = sigma * V and realinvV = invV / sigma
        #print 'sigma: '+str(sigma)+ ' logdet: '+ str(logdetV)
        l = self.ntips*(math.log(2*math.pi*sigma)) + logdetV
        print "first part lh = "+str(l);
        e = self.y - self.X * beta
        l += (e.T * invV *e) / sigma
        l /= -2
        return l

    def fitGLS(self,sigma=None):
        self._fitGLSWorker(lamb=None,sigma=sigma)
        print "beta: " + str(self.beta)
        print "sigma: " + str(self.sigma)
        print "Likelihood:" + str(self.lhval)


    def fitGLSLambda(self,sigma=None):
        lhLambda = lambda(l): -self._fitGLSWorker(lamb=l,sigma=sigma)
        #res = minimize_scalar(lhLambda,method='Golden',bracket=(0.05,0.1))
        res = minimize_scalar(lhLambda,method='Bounded',bounds=(0,1))
        optlambda = res.x
        self.sigmag = self.sigma*optlambda
        self.sigmae = self.sigma*(1-optlambda)
        #for l in [ x/200.0 for x in range(0,205,5)]:
        #    optlh = self.fitGLSWorker(lamb=l)
        #   #print "Lambda: "+str(l) + " lh: "+ str(optlh)
        print "beta: "+str(self.beta)
        print "sigma_b: " + str(self.sigmag)
        print "sigma_e: " + str(self.sigmae)
        print "Last lambda: "+str(self.lamb)
        print "optlambda: "+str(optlambda)
        print "Likelihood: "+str(self.lhval)

    def fitLMMBasic(self):
        n = self.ntips
        # we are assuming that Z is the identity
        lhLMM = lambda(h2): -self._fitLMMBasicWorker(h2)
        res = minimize_scalar(lhLMM,method='Bounded',bounds=(0,1))
        h2 = res.x
        # Compute small sigmas and b
        self.h2 =  h2
        self.sigmag = self.sigma*h2
        self.sigmar = self.sigma*(1-h2)
        self.b = float(h2)*(self.G*self.invH*self.e)
        # output
        print "H2 = "+str(self.h2)
        print "beta ="+ str(self.beta)
        print "Sigma Random = " + str(self.sigmag)
        print "Sigma Residual = " + str(self.sigmar)
        print "lh = "+ str(self.lhval)

    # Basic animal model using inheritability
    # V = sigma2(h2*Z'*G*Z + (1-h2)*R)
    def _fitLMMBasicWorker(self,h2):
        h2 = float(h2)
        H = h2*self.G + (1-h2)*self.R
        self.invH = np.linalg.inv(H)
        (s,logdetH) = np.linalg.slogdet(H)
        # estimator beta = (X' invV X)^{-1} (X' invV y)
        left = np.linalg.inv(self.X.T * self.invH * self.X)
        right = self.X.T * self.invH * self.y
        self.beta =  left * right
        self.e = self.y - self.X * self.beta
        self.sigma = (self.e.T * self.invH  * self.e)/(self.ntips-1)
        if (self.sigma < 0):
            raise Exception("Sigma negative")
        self.lhval = self.lh(self.beta,self.sigma,self.invH,logdetH)
        return self.lhval
