__author__ = 'Igor Siveroni'

import sys
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import sparse
import random
import pandas as pd

import phylotree
reload(phylotree)
from phylotree import PhyloTree

import lmmBasicAnimal
reload(lmmBasicAnimal)
from lmmBasicAnimal import LMMBasicAnimal

import gibbsAnimal
reload(gibbsAnimal)
from gibbsAnimal import GibbsAnimal
import mc2lmm
reload(mc2lmm)
from mc2lmm import MC2Lmm

import phylmm
reload(phylmm)
from phylmm import PhyLmm

import glmmtools
reload(glmmtools)
import glmmtools as gt


import time

# beta = list of fixed locations
# X = design matrix or range [l,h]
# Ginvs = list of variance inverses
# First one corresponds to Phylo
def lmm_sample(beta,Gs,sigma2=[],X=[5,30]):
    p = len(beta)
    beta = np.matrix(beta).reshape(p,1)
    c = len(Gs)
    qs = [np.shape(g)[0] for g in Gs]
    q = sum(qs)
    n =  qs[0]
    if (isinstance(X,list)):
        # X contains list of ranges
        rangex0 = X[0]
        rangex1 = X[1]
        X = np.matrix(np.zeros((n,p)))
        X[:,0] = 1
        X[:,1] = rangex0 + (rangex1-rangex0)*np.random.rand(n).reshape(n,1)
    else:
        assert(isinstance(X,np.matrix))
        (n1,p1) = np.shape(X)
        assert((n1==n)and(p1==n))
    # Create Z: n x q  --- Sparse
    # sparse(y)
    zs = [sparse.eye(n)]
    data = [1 for u in range(0,n)]
    for qi in qs[1:]:
        # col = [random.randint(0,qi-1) for i in range(0,n)]
        col = [i%qi for i in range(0,n)]
        zi = sparse.csc_matrix((data,(range(0,n),col)),shape=(n,qi))
        zs.append(zi)
    Z = sparse.hstack(zs)
    assert(len(sigma2)==(c+1))
    # Sample random effects
    us = []
    cholLs=[]
    for g,s2 in zip(Gs,sigma2[:-1]):
        cholL = np.linalg.cholesky(g)
        cholLs.append(cholL)
        ui = gt.rmvn(0,cholL=cholL,sigma=s2)
        us.append(ui)
    u = np.vstack(us)
    # Sample residual
    sigma2_e = sigma2[-1]
    R = np.asmatrix(np.eye(n))
    cholL = np.linalg.cholesky(R)
    e = gt.rmvn(0,cholL=cholL,sigma=sigma2_e)
    # Calculate y
    y = X*beta + Z*u + e
    return (y,X,Z,cholLs,u,e)

def lmmsample(tree):
    # The script

    #tree = PhyloTree(treeFile, schema='newick')
    #tree.normalize()
    G = tree.getV()

    (n,n) = np.shape(G)
    num_fixed = 2
    # Create Beta
    beta = np.matrix([15,2.5]).reshape(num_fixed,1)

    # Create X
    rangex0 = 5
    rangex1 = 30
    X = np.matrix(np.zeros((n,num_fixed)))
    X[:,0] = 1
    X[:,1] = rangex0 + (rangex1-rangex0)*np.random.rand(n).reshape(n,1)

    # Create Z
    Z = np.asmatrix(np.eye(n))

    # Sample b ~ N(0,sigma2_b * V)
    sigma2_b = 5
    cholL = np.linalg.cholesky(G)
    b = gt.rmvn(0,cholL=cholL,sigma=sigma2_b)

    # Sample residual e ~ N(0, sigma2_e * R)
    sigma2_e = 1.0
    R = np.asmatrix(np.eye(n))
    cholL = np.linalg.cholesky(R)
    e = gt.rmvn(0,cholL=cholL,sigma=sigma2_e)

    # Calculate y
    y = X*beta + Z*b + e

    #yl = y.ravel().tolist()[0]
    #xl =  X[:,1].ravel().tolist()[0]
    #plt.plot(xl,yl,'ro')
    #ymax = max(yl)
    #plt.axis([0,rangex1,0, 10*(int(ymax/10)+1)])
    #plt.show()
    return (y,X,beta,Z,b,G,cholL,sigma2_b,R,sigma2_e)

# fixed_names: fixed effects associated to columns in matrix X
# tree: random effect 'species'
# we should get a list of trees and random effect names instead
def saveModelData(filename,y_name,fixed_names,y,X,tree):
    headers =  [y_name] + fixed_names + ['species']
    (n,p) = np.shape(X)
    fixed = [np.asarray(X[:,i]).reshape(-1) for i in range(1,p)]
    data_columns = [y] + fixed + [tree.tipNames]
    #data_columns = [y] + [ X[,i] for i in range(1,w)]
    gt.createDataFile(filename,headers,data_columns)

# Z is sparse
def saveModelData2(filename,y_name,fixed_names,random_names,y,X,Z,trees):
    headers =  [y_name] + fixed_names + random_names
    (n,p) = np.shape(X)
    # = [np.asarray(d).reshape(np.size(d)) for d in column_arrays]
    fixed = [np.asarray(X[:,i]).reshape(-1) for i in range(1,p)]
    # Match Z against each tree in trees
    random_columns = convertZ(Z,n,[tree.tipNames for tree in trees])
    # concatenate all data columns
    data_columns = [np.asarray(y).reshape(-1)] + fixed + random_columns
    #data_columns = [y] + [ X[,i] for i in range(1,w)]
    gt.createDataFile(filename,headers,data_columns)

# Maps incidence matrix to vector of names
def convertZ(Z,n,names_list):
    Z = sparse.csr_matrix(Z)
    (nn,q) = np.shape(Z)
    qis = [np.size(l) for l in names_list]
    c = len(qis)
    assert(nn==n)
    last = np.cumsum(qis).tolist()
    first = [0] +last
    first = first[:-1]
    columns = [[] for qi in qis]
    for row in range (0,n):
        r = Z.getrow(row).toarray().reshape(-1)
        # same as np.where(z == 1)
        hits =  [i for i in range(q) if r[i] != 0]
        for (column_idx,f,l,names) in zip(range(c),first,last,names_list):
            hit_idx = hits[0]
            if ((hit_idx >= f) and (hit_idx < l)):
                data = names[hit_idx - f]
                hits = hits[1:]
            else:
                data = 'NA'
            columns[column_idx].append(data)
    return columns


def simulate_two_trees(file=True):
    if (file):
        f1 = 'data/MorphoSpeciesTree.Trees.newick'
        t1 = PhyloTree(f1, schema='newick')
        t1.normalize()
        a = t1.getV()
        f2 = 'data/LineageTree.trees.newick'
        t2 = PhyloTree(f2, schema='newick')
        t2.normalize()
        b = t2.getV()
        trees = [t1,t2]
    else:
        # Generate two random positive definite matrices (variances)
        a = np.matrix([random.uniform(-3,3) for i in range(0,25)]).reshape(5,5)
        a = a*a.T
        b = np.matrix([random.uniform(-2,2) for i in range(0,9)]).reshape(3,3)
        b = b*b.T
        trees = []
    beta = [15,2.5]
    Gs = [a,b]
    (y,X,Z,Ls,u,e) =lmm_sample(beta,Gs=Gs,sigma2=[5,2,0.5])
    saveModelData2('data/bodybraindata.txt','Body',['Brain'],['r1','r2'],y,X,Z,trees)
   # ---- The current MLE for LMM does not accept more then one variance component ---
    # Compute inverses
    Ginvs = [np.linalg.inv(G)  for G in Gs]
    mcmc = MC2Lmm(y,X,Z,Ginvs=Ginvs)
    mcmc.set_init_params(init_beta=beta,init_sigma2_u=0.5,init_sigma2_e=0.5)
    # Time this
    t1 = time.time()
    ret = mcmc.simulate(num_samples=500,burnin=1000,thin=20)
    t2 = time.time()
    print 'time='+str(t2-t1)
    print "**** Finished sampling ****"
    print "NUM SAMPLES = " + str(mcmc.num_samples)
    if (ret != 0):
        print "Horror!"
        if (mcmc.num_samples < 1):
            return -1
        mcmc.plot()
    else:
        print 'Finished successfully'
        mcmc.plot()

def simulate_one_tree():
    #treeFile = 'data/test.newick'
    treeFile = "data/bz.trees.newick"
    tree = PhyloTree(treeFile, schema='newick')
    tree.normalize()
    G = tree.getV()
    Gs = [G]
    #(y,X,beta,Z,b,G,cholL,sigma2_b,R,sigma2_e) = lmmsample(tree)
    beta = [15,2.5]
    (y,X,Z,Ls,u,e) =lmm_sample(beta,Gs=Gs,sigma2=[5,1])
    (n,xx) = np.shape(y)
    # saveModelData('data/bodybraindata.txt','Body',['Brain'],y,X,tree)
    model =  LMMBasicAnimal(ylist = ['Body'],plist=['Brain'])
    #return
    #model = LMMBasicAnimal(ylist = ['length'],plist=['size'])
    R = np.asmatrix(np.eye(n))
    model.loadFromMatrices(y,X,Z,G,R)
    model.fitLMMBasic()
    #model.fitGLSLambda()
    #return 0
    invG = np.linalg.inv(model.G)
    #cholG = np.asmatrix(scipy.linalg.cholesky(G,lower=True))
    print ("*** initializing LmmMcmc object")
    # gibbs = GibbsAnimal(y,X,Z,invG,cholG)
    gibbs = MC2Lmm(y,X,Z,Ginvs=[invG])
    gibbs.set_init_params(init_beta=[10,1.5],init_sigma2_u=0.5,init_sigma2_e=0.5)
    # Time this
    t1 = time.time()
    ret = gibbs.simulate(num_samples=100,burnin=2000,thin=20)
    t2 = time.time()
    print 'time='+str(t2-t1)
    print "**** Finished sampling ****"
    print "NUM SAMPLES = " + str(gibbs.num_samples)
    if (ret != 0):
        print "Horror!"
        if (gibbs.num_samples < 1):
            return -1
        gibbs.plot()
    else:
        print 'Finished successfully'
        gibbs.plot()
    #gibbs.save_to_file('samples2.log',filetype='xls')
    #print samples
    #print np.shape(samples)
    #(num_params,num_samples) = np.shape(samples)

    return 0


def fit_simulate():
    #treeFile = 'data/test.newick'
    treeFile = "data/bz.trees.newick"
    tree = PhyloTree(treeFile, schema='newick')
    tree.normalize()
    (y,X,beta,Z,b,G,cholL,sigma2_b,R,sigma2_e) = lmmsample(tree)
    #saveModelData2('data/bodybraindata.txt','Body',['Brain'],y,X,tree)
    model =  LMMBasicAnimal(ylist = ['Body'],plist=['Brain'])
    #return
    #model = LMMBasicAnimal(ylist = ['length'],plist=['size'])
    model.loadFromMatrices(y,X,Z,G,R)
    model.fitLMMBasic()
    #model.fitGLSLambda()
    #return 0
    invG = np.linalg.inv(model.G)
    #cholG = np.asmatrix(scipy.linalg.cholesky(G,lower=True))
    print ("*** initializing LmmMcmc object")
    # gibbs = GibbsAnimal(y,X,Z,invG,cholG)
    gibbs = MC2Lmm(y,X,Z,Ginvs=[invG])
    gibbs.set_init_params(init_beta=[10,1.5],init_sigma2_u=0.5,init_sigma2_e=0.5)
    # Time this
    t1 = time.time()
    ret = gibbs.simulate(num_samples=100,burnin=2000,thin=20)
    t2 = time.time()
    print 'time='+str(t2-t1)
    print "**** Finished sampling ****"
    print "NUM SAMPLES = " + str(gibbs.num_samples)
    if (ret != 0):
        print "Horror!"
        if (gibbs.num_samples < 1):
            return -1
        gibbs.plot()
    else:
        print 'Finished successfully'
        gibbs.plot()
    #gibbs.save_to_file('samples2.log',filetype='xls')
    #print samples
    #print np.shape(samples)
    #(num_params,num_samples) = np.shape(samples)

    return 0

def fit_simulate2():
    treeFile = 'data/test.newick'
    dataFile = "data/test.data"
    #treeFile = "data/bz.trees.newick"
    tree = PhyloTree(treeFile, schema='newick')
    tree.normalize()
    #(y,X,beta,Z,b,G,cholL,sigma2_b,R,sigma2_e) = lmmsample(tree)
    #saveModelData2('data/bodybraindata.txt','Body',['Brain'],y,X,tree)
    #model =  LMMBasicAnimal(ylist = ['Body'],plist=['Brain'])
    #return
    model = LMMBasicAnimal(ylist = ['length'],plist=['size'])
    model.loadFromFile(dataFile,treeFile,sep='SPACE')

    invG = np.linalg.inv(model.G)
    cholG = np.asmatrix(scipy.linalg.cholesky(model.G,lower=True))
    print ("*** initializing LmmMcmc object")
    gibbs = GibbsAnimal(model.y,model.X,model.Z,invG,cholG)
    #gibbs = MC2Lmm(y,X,Z,Ginvs=[invG])
    #gibbs.set_init_params(init_beta=[10,1.5],init_sigma2_u=0.5,init_sigma2_e=0.5)
    # Time this
    t1 = time.time()
    ret = gibbs.simulate(num_samples=100,burnin=2000,thin=20)
    t2 = time.time()
    print 'time='+str(t2-t1)
    print "**** Finished sampling ****"
    print "NUM SAMPLES = " + str(gibbs.num_samples)
    if (ret != 0):
        print "Horror!"
        if (gibbs.num_samples < 1):
            return -1
        gibbs.plot()
    else:
        print 'Finished successfully'
        gibbs.plot()
    #gibbs.save_to_file('samples2.log',filetype='xls')
    #print samples
    #print np.shape(samples)
    #(num_params,num_samples) = np.shape(samples)

    return 0



def fit_test():
    print "Testing GLS"
    model = LMMBasicAnimal(ylist = ['length'],plist=['size'])
    dataFile = "data/test.data"
    treeFile = "data/test.newick"
    # loadFromFile(self,dataFile,treeFile,taxa_column='taxa',sep='TAB'):
    model.loadFromFile(dataFile,treeFile,sep='SPACE')
    model.fitGLS()
    # model.fitLMMBasic()
    print "--------------------"
    # model.fitGLSLambda()

def fit_mammals():
    model = LMMBasicAnimal(ylist = ['Body'],plist=['Brain'])
    dataFile = "data/bzdata.txt"
    treeFile = "data/bz.trees.newick"
    # loadFromFile(self,dataFile,treeFile,taxa_column='taxa',sep='TAB'):
    model.loadFromFile(dataFile,treeFile,taxa_column='Name',sep='SPACE')
    model.fitGLS()
    # model.fitLMMBasic()
    print "--------------------"
    model.fitGLSLambda()

def test_phylmm():
    dataFile = "data/test.data"
    df = data = pd.read_table(dataFile,delim_whitespace=True)
    print df
    return
    model = PhyLmm(ynames=['Body'], fnames=['Brain'], rnames=['morpho'])


def main_exec(par):
    #fit_test()
    #fit_mammals()
    fit_simulate2()
    #simulate_two_trees()
    #simulate_one_tree()
    #test_phylmm()


if __name__ == '__main__':
    # Called from command line.
    # execfile does not pass arguments
    main_exec(sys.argv[1:])


