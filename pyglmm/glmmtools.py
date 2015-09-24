__author__ = 'Igor Siveroni'

import dendropy
import numpy as np
import sys
from scipy.stats import invgamma
from scipy.stats import chi2
from scipy.linalg import cho_solve
from scipy.linalg import solve_triangular

# schemas supported = "nexus", "newick"
# Coverts between file representation of Phylogenetic trees
def convertTreeFile(infile, schema1, schema2,outfile=""):
    assert(schema1 in ['nexus','newick'])
    assert(schema2 in ['nexus','newick'])
    t = dendropy.Tree(stream=open(infile), schema=schema1)
    if (outfile==""):
        outfile = infile+"."+schema2
    t.write(stream=open(outfile,'w'),schema=schema2)

def applyLambda(V,l):
    d = V.diagonal().copy()
    V *= l
    np.fill_diagonal(V,d)

# column_names: a list of strings denoting column/header names
# column_vectors: a list of one dimensional arrays (or list) of values associated to each column
# float arrays are formated to 4 decimal positions
def createDataFile(fileName,column_names,column_vectors,sep=' '):
    f = open(fileName,'w')
    if (len(column_names) !=  len(column_vectors)):
        print "Number of Column Headers and Data do not match"
        return -1
    n = min([np.size(d) for d in column_vectors])
    # write headers
    for header in column_names[:-1]:
        f.write(header)
        f.write(sep)
    f.write(column_names[-1])
    f.write('\n')
    # write rows
    for i in range(0,n):
        for d in column_vectors[:-1]:
            v = d[i]
            if isinstance(v,float):
                v = '%.4f' % v
            f.write(str(v))
            f.write(sep)
        v = column_vectors[-1][i]
        if isinstance(v,float):
            v = '%.4f' % v
        f.write(str(v))
        f.write('\n')
    f.close()

# ************************** Prior Distributions ******************************************


# Scaled inverted chi-square. Tuple (v,s2) = (degrees freedom,scale)
# Implemented as a dictionary
def create_chisquare_si(v,s2):
    return {'dist': 'chisquare-si','df':v,'scale':s2 }


# *************************** Random Sampling from distributions ***********************

# Our first assumption is that mu, V and L are matrices.
# mu is a column vector (n,1)
# V = L * L.T = U.T * U --> U=L.T
def rmvn(mu,V=None,cholL=None,inverse=False,sigma=None):
    if (cholL is None):
        if (V is None):
            raise Exception('rmvn:Either V matrix or its Cholesky factorization must be provided')
        cholL = np.linalg.cholesky(V)
        # cholLU = np.matrix(scipy.linalg.cholesky(V,lower=True),copy=False)
        (n,m) = np.shape(V)
    else:
        (n,m) = np.shape(cholL)
    # Cholesky: V = L * L.T
    Z = np.matrix(np.random.normal(size=n),copy=False).reshape(n,1)
    if (inverse):
        #Linv = np.linalg.inv(L)
        #VM = Linv.T * Z
        VM = np.asmatrix(solve_triangular(cholL.T,Z,lower=False))
    else:
        VM = cholL * Z
    if (sigma is not None):
        VM *= np.sqrt(sigma)
    return mu + VM

# X ~ scaled-inv-chisquare(v,s) <=> X ~ inv-gamma(v/2,v*s/2)
# X ~ inv-gamma(a,b) <=> KX ~ inv-gamma(a,Kb)
# X ~ inv-gamma(a,1) <=> bX ~ inv-gamma(a,b)
def rchisquare_inv_scaled(df,scale,size=None):
    a = df/2.0
    b = (df*scale)/2.0
    # size=None returns a scalar
    return invgamma.rvs(a,scale=b,size=size)

def rchisquare_inv_scaled2(df,scale,size=None):
    x = chi2.rvs(df,size=size)
    if (size is None):
        return (1/x)*df*scale
    else:
        return [(1/e)*df*scale for e in x]

