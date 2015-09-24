__author__ = 'Igor Siveroni'

__author__ = 'Igor Siveroni'

import numpy as np
import glmmtools as gt
import matplotlib.pyplot as plt
import scipy
from scipy import sparse

# Implementation of a Gibbs Sampler for Linear Mixed Models
# Single trait
# Multiple random effects.
# The model: y = XB + Zu + e
# Vector u = [u_1 u_2 ... u_c ]
# u_i ~ N(0, G_i * sigma2_u_i
# e ~ N(0, R_1 * sigma2_e
# I am not considering priors for B (beta) - don't understand that part.
# Maybe for glm makes more sense?


# Input:
# y = Vector of observations (response variable/trait)
# X = Design matrix for fixed effects
# Z = Design matrix for random effects (sparse)
# t = number of traits (t=1)
# u_len = list with the size of each random effect vector.
# Calculated c = number of random effects vectors. c = len(u_len)
# Ginvs = list of covariance matrix inverses associated to each random effect vector
# Rinvs = list of covariance matrix inverses associated to the residuals of each trait.
# Gpriors = list of priors for each variance component
# Rpriors = list of priors
# Note: Since we are dealing with a single trait model, there's only 1 element
#       in Rinvs and Rpriors.
# Vectors are all of class Matrix
class MC2Lmm:
    def __init__(self,y,X,Z,t=1,u_len=None,Ginvs=None,Rinvs=None):
        # ********* Check dimensions first ****************
        if (t > 1):
            raise Exception("Number of traits must be 1")
        # Response data column vector - single trait
        # A bit flexible here - flatten to column vector
        self.y_size = np.size(y)
        self.n = self.y_size  # y_size / t
        self.y = np.asmatrix(y).reshape(self.y_size,1)
        # ***************** X Fixed effects design matrix
        dim = np.shape(X)
        if (len(dim) != 2):
            raise Exception('X must be di-dimensional')
        if (dim[0] != self.n):
            raise Exception("X: height doesn't match group size")
        self.X = np.asmatrix(X)
        self.p = dim[1]
        # ***************** Z Random effects design matrix
        dim = np.shape(Z)
        if (len(dim) != 2):
            raise Exception('Z must be di-dimensional')
        if (dim[0] != self.n):
            raise Exception("Z: height doesn't match group size")
        #self.Z = np.asmatrix(Z)  - changes dimension of matrix
        self.Z = Z
        self.q = dim[1]
        self.pq = self.p + self.q
        # *********** Random vectors and Ginvs - sizes
        if ((u_len is not None) and (Ginvs is not None)):
            # Check that group sizes coincide
            if (u_len != [np.shape(ginv)[0] for ginv in Ginvs]):
                raise Exception("Random effect vector and matrix inverse sizes don't agree")
        else:
            if ((u_len is None) and (Ginvs is None)):
                raise Exception("Can't figure out number of random effects vectors. Must enter Ginvs or u_len")
            else:
                if (u_len is None):
                    u_len = [np.shape(ginv)[0] for ginv in Ginvs]
                else:
                    Ginvs = [sparse.eye(l) for l in u_len]
        # u_len and Ginvs agree
        # Now check if it agrees with Z
        if (self.q !=  sum(u_len)):
            raise Exception('Size mismatch: Width of Z and total random vectors')
        self.c = len(u_len) # number of random effect vectors
        self.qi = u_len
        self.Ginvs = Ginvs
        # Check Rinvs
        if (Rinvs is None):
            Rinvs = [sparse.eye(self.n) ]
        else:
            if (np.shape(Rinvs[0])[0] != self.n):
                raise Exception('R matrix size must match size of y')
        self.Rinvs = Rinvs
        self.t = 1 # number of traits
        # Initialize intermediate structures.
        self._initialize()
        # Initialize priors to defaults
        self.set_default_priors()
        # Default initial values for sampling
        self.set_init_params(init_beta=0,init_u=0,init_sigma2_u=0,init_sigma2_e=0)

    def _initialize(self):
        # Create W)
        # self.W = np.hstack((self.X,self.Z)) -- stores W as a sparse matrix
        # Location parameters = theta
        self.theta = np.asmatrix(np.zeros(self.pq).reshape(self.pq,1))
        self.theta_hat = np.asmatrix(np.zeros(self.pq).reshape(self.pq,1))
        ## Aliases - sub-arrays (does not create a copy)
        # Fixed Parameters
        self.beta = np.matrix(self.theta[0:self.p,0],copy=False)
        # All random vectors
        self.u = np.matrix(self.theta[self.p:,0],copy=False)
        # Might be useful to have each random vector as sub-matrices
        # todo: implement if needed
        # Variance Components - note these are vectors (ndarrays)
        self.sigma2_u = np.asmatrix(np.zeros(self.c).reshape(self.c,1))
        self.sigma2_e = np.asmatrix(np.zeros(self.t).reshape(self.t,1))

        # other initializations
        self._compute_mm_equations()
        # Rest
        self.store_random_effects=False

    def set_default_priors(self):
        default_gpriors = [gt.create_chisquare_si(1.0,0.5)  for i in range(0,self.c)]
        default_rpriors = [gt.create_chisquare_si(1.0,0.5) ]
        # Set priors to defaults scaled inverse chi-square. User can change them later.
        assert(self.set_priors(Gpriors=default_gpriors, Rpriors=default_rpriors))


    # For the time being, we will use scale inverse chi-square priors
    # This affects the sampling of individual variance components
    # Return True if input is correct. False otherwise.
    def set_priors(self,Gpriors=None,Rpriors=None):
        # Check they are all chi-square
        if (Gpriors is not None):
            if (any([ (p['dist'] != 'chisquare-si') for p in Gpriors ])):
                print('All Gpriors must be scaled inverse chi-square')
                print('Ignoring all new Gpriors')
                return False
            # One prior for each random effect variance
            if (len(Gpriors) != self.c):
                print("Incorrect number of random effect variance component priors (Gpriors)")
                print('Ignoring all new Gpriors')
                return False
            # A list of random effect variance component priors
            self.priors_sigma2u = Gpriors
        if (Rpriors is not None):
            if (any([ (p['dist'] != 'chisquare-si') for p in Rpriors ])):
                print('All Rpriors must be scaled inverse chi-square')
                print('Ignoring new Rprior')
                return False
            if (len(Rpriors) != 1):
                print('The model uses ONE Residual variance component')
                print('Ignoring new Rpriors')
            self.priors_sigma2e = Rpriors
        return True

    def set_init_params(self,init_beta=None,init_u=None,init_sigma2_u=None,init_sigma2_e=None):
        if (init_beta is not None):
            self.init_beta = init_beta
        if (init_u is not None):
            self.init_u = init_u
        if (init_sigma2_u is not None):
            self.init_sigma2_u = init_sigma2_u
        if (init_sigma2_e is not None):
            self.init_sigma2_e = init_sigma2_e


       # Mixed Model Equations
    def _compute_mm_equations(self):
        # C beta_hat = r
        # C dimension = (p+q)^2
        self.mmeq_C = np.asmatrix(np.zeros(self.pq*self.pq).reshape(self.pq,self.pq))
        self.mmeq_rhs = np.asmatrix(np.zeros(self.pq).reshape(self.pq,1))
        self.mmeq_C22 = np.matrix(self.mmeq_C[self.p:,self.p:],copy=False)
        Rinv = self.Rinvs[0] # One Residual structure
        self.mmeq_rhs[0:self.p,0] = self.X.T * Rinv * self.y
        self.mmeq_rhs[self.p:,0] = self.Z.T * Rinv *self.y
        # Coefficient matrix
        self.mmeq_C[0:self.p,0:self.p] = self.X.T * Rinv *self.X
        self.mmeq_C[0:self.p,self.p:] = self.X.T * Rinv * self.Z
        self.mmeq_C[self.p:,0:self.p] = self.Z.T * Rinv * self.X
        self.ZTRZ = self.Z.T * Rinv * self.Z
        #self._update_mm_equations() - executed before each variance sampling

    # Updates/Computes bottom-right corner of the coefficient matrix of the model equations
    def _update_mm_equations(self):
        se = self.sigma2_e[0,0]
        # self.mmeq_C22[:,:] = self.ZTRZ
        self.mmeq_C22[:,:] = self.ZTRZ.todense()
        l=0
        r=0
        for i in range(0,self.c):
            l = r
            r += self.qi[i]
            k = se / self.sigma2_u[i,0]
            self.mmeq_C22[l:r,l:r] += k*self.Ginvs[i]


    # It should return an MCMC-run object (class to be implemented)
    def simulate(self,num_samples=10000,burnin=1000,thin=100):
        # MCMC runnning parameters
        if (burnin < 0):
            print "Incorrect burn-in value"
            return None
        if (num_samples < 1):
            print "Incorrect number of samples value"
            return None
        # Preparations
        self._initialize_params()
        # Output - Samples
        ct = self.c + self.t # t=1 for first model
        self.beta_samples = np.asmatrix(np.zeros(self.p*num_samples).reshape(self.p,num_samples))
        if (self.store_random_effects):
            self.u_samples = np.asmatrix(np.zeros(self.q*num_samples).reshape(self.q,num_samples))
        self.sigma2_samples = np.asmatrix(np.zeros(ct*num_samples).reshape(ct,num_samples))
        #print 'Initial Value'
        #print self.theta.reshape(1,self.pq)
        print '*** Burn-in ***'
        self.num_samples = 0
        try:
            for i in range(0,burnin):
                self._update_mm_equations()
                self._sample_model()
        except Exception, e:
            print 'Error during burn-in'
            print e
            raise e
            return -1
        print '*** Sampling ***'
        try:
            for i in range(0,num_samples):
                # Blocked updates
                #print "Sample:"+str(i)
                for j in range(0,thin):
                    self._update_mm_equations()
                    self._sample_model()
                self.beta_samples[:,i] = self.theta[:self.p,0]
                if (self.store_random_effects):
                    self.u_samples[:,i] = self.theta[self.p:,0]
                self.sigma2_samples[:self.c,i] = self.sigma2_u
                self.sigma2_samples[self.c:,i] = self.sigma2_e
                #print '--> sampled: '+str(i)
        except Exception, e:
            print 'Error while Sampling'
            print e
            self.num_samples = i+1
            return -1
        self.num_samples = num_samples
        return 0

    def set_vector_matrix(self,m,l):
        if (isinstance(l,list)):
            for i in range(0,len(l)):
                m[i,0] = l[i]
        elif (isinstance(l,(int,long,float))): # l is a number
            m[:,0] = l
        else:
            raise Exception('Unknown type of initial value to matrix')

    def _initialize_params(self):
        self.set_vector_matrix(self.theta[0:self.p,0], self.init_beta)
        self.set_vector_matrix(self.theta[self.p:,0], self.init_u)
        self.set_vector_matrix(self.sigma2_u[:,0], self.init_sigma2_u)
        self.set_vector_matrix(self.sigma2_e[:,0], self.init_sigma2_e)


    # Samples location parameters one at a time using uni-variate normal random distributions
    def _sample_model(self):
        # Sample location parameters
        self._sample_effects_block_chol()
        #self._sample_effects_block()
        # Sample variance components
        self._sample_variance_components()
        # check for nan or inf
        if (np.any(map((lambda x: np.isnan(x) or np.isinf(x)),self.theta))):
            print self.theta
            raise Exception("NAN/INF")

    # Gibbs sampling - uses full conditional distribution of a single location (i)
    def _sample_effects(self):
        #self._update_mm_equations()
        C = self.mmeq_C
        sigma2_e0 = self.sigma2_e[0,0]
        for i in range(0,self.pq):
            self.theta[i,0] = 0
            #print np.shape(C[i,:])
            #print np.shape(self.theta[0:self.num_locations,:])
            mu_theta = (C[i,:] * self.theta)[0,0]
            mu_theta = (self.mmeq_rhs[i] - mu_theta)/C[i,i]
            s2_theta = self.sigma2_e0/C[i,i]
            self.theta[i,0] = np.random.normal(mu_theta,s2_theta)

    def _sample_effects_block(self):
        #self._update_mm_equations()
        C = self.mmeq_C
        invC = np.linalg.inv(C)
        # Compute location estimators - expected values
        # Equation C * loc = rhs -> loc = C^{-1}*rhs
        self.theta_hat = invC * self.mmeq_rhs
        # Sample from multivariate normal
        sigma2_e0 = self.sigma2_e[0,0]
        self.theta[:,0] = gt.rmvn(self.theta_hat,V=invC,sigma=sigma2_e0)

    # Cholesky factorization
    def _sample_effects_block_chol(self):
        #self._update_mm_equations()
        C = self.mmeq_C
        L = np.linalg.cholesky(C)
        #L = np.asmatrix(scipy.linalg.cholesky(C,lower=True))
        #print 'before solve'
        self.theta_hat = np.asmatrix(scipy.linalg.cho_solve((L,True),self.mmeq_rhs))
        #print 'after solve'
        # invChol = np.linalg.inv(cholC)
        sigma2_e0 = self.sigma2_e[0,0]
        self.theta[:,0] = gt.rmvn(self.theta_hat,cholL=L,inverse=True,sigma=sigma2_e0)
        #print 'after rvn'


    def _sample_variance_components(self):
        ## Sample sigma2_u
        first = self.p
        for i in range(0,self.c):
            prior = self.priors_sigma2u[i]
            qi = self.qi[i]  # size of random vector ui
            df = prior['df'] + qi
            #print 'num_fixed: '+str(self.num_fixed)+' num_locations: '+ str(self.num_locations)
            #print 'theta:'+str(np.shape(self.theta))
            ui = self.theta[first:first+qi,0]
            first += qi
            #print np.shape(u)
            #print np.shape(self.invG)
            try:
                scale1 = (ui.T * self.Ginvs[i] * ui)[0,0]
                #print 'u.T * invG * u = ' + str(scale0)
                scale = (scale1 + prior['df'] * prior['scale']) / df
                #print 'scale = ' + str(scale) + ' scalebig = ' + str(scalebig)
                self.sigma2_u[i,0] = gt.rchisquare_inv_scaled(df,scale)
                #print 'sigma_b: '+str(x)+ ' df='+str(df)+' scale='+str(scale)
            except Exception, e:
                print "Error while sampling sigma2_u"
                print "df="+str(df)+" scale="+str(scale)
                print e
                raise e
        ## Sample sigma2_e
        #print 'sigma_e  index: '+ str(self.sigma2_e_idx)
        prior = self.priors_sigma2e[0]
        try:
            df = prior['df'] + self.n
            e = self.X*self.theta[0:self.p,0] + self.Z*self.theta[self.p:]
            e = self.y - e
            scale = ((e.T * e)[0,0] + prior['df']*prior['scale']) / df
            #print 'df '+ str(df) + ' scale ' + str(scale)
            self.sigma2_e[0,0] = gt.rchisquare_inv_scaled(df,scale)
            #print 'sigma_e: '+str(x)+ ' sample_e'+str(self.sample)
        except Exception, e:
            print "Error while sampling sigma2_e"
            print "df="+str(df)+" scale="+str(scale)
            print e
            raise e

    def plot(self):
        #n = len(params)
        num_samples = self.num_samples
        plt.figure(1)
        #print "NUM SAMPLES="+str(num_samples)

        x = range(0,num_samples)

        beta = [self.beta_samples[i,:] for i in range(0,self.p)]
        s2u = [self.sigma2_samples[i,:] for i in range(0,self.c)]
        s2e = self.sigma2_samples[self.c,:]
        print 'MEANS:'
        for i in range(0,self.p):
            print np.mean(beta[i])
        for i in range(0,self.c):
            print np.mean(s2u[i])
        print np.mean(s2e)

        plt.subplot(411)
        y = beta[0].tolist()[0]
        plt.plot(x,y)
        plt.subplot(412)
        y =  beta[1].tolist()[0]
        plt.plot(x,y)

        plt.subplot(413)
        y = s2u[0].tolist()[0]
        plt.plot(x,y)
        if (self.c > 1):
            plt.subplot(414)
            y = s2u[1].tolist()[0]
            plt.plot(x,y)
        else:
            # Residuals
            plt.subplot(414)
            y = s2e.tolist()[0]
            plt.plot(x,y)
        # show
        plt.show()
