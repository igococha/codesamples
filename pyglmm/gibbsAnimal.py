__author__ = 'Igor Siveroni'

import numpy as np
import glmmtools as gt
import matplotlib.pyplot as plt
import scipy

# Implementation of a Gibbs Sampler for Animal Models
# Single trait
# one random effect vector u associated with variance A (phylogeny)
# (not difficult to extend to more than one)

# The model:
# y = XB + Zu + e
# E(y) = XB, E(u) = 0, E(e) = 0
# y|B,u ~ N(XB+Zu,R)
# y ~ N(XB, ZAZ'+R), e ~ N(0,R), u ~ N(0,A)
# Simple version with one random effect, variance A
# R = I*sigma2_e , A = G*sigma2_b
# One trait (one response variable)
# Must link this to function that reads data and builds the design matrices
class GibbsAnimal:
    def __init__(self,y,X,Z,invG,cholG):
        # ********* Check dimensions first ****************
        # single trait extract dimension
        self.num_y,n = np.shape(y)
        if (n != 1):
            raise Exception("Response data should be stored as a column vector")
        self.y = y
        n,self.num_fixed = np.shape(X)
        if (n != self.num_y):
            raise Exception("X: height doesn't match number of species")
        self.X = X
        n,self.num_random = np.shape(Z)
        if (n != self.num_y):
            raise Exception("Z: height doesn't match number of species")
        self.Z = Z
        # Create W
        self.W = np.hstack((X,Z))
        # Check invG size
        n,m = np.shape(invG)
        if ((n != self.num_random) or (m != self.num_random)):
            raise Exception("invG: Incorrect dimensions")
        self.invG = invG
        self.cholG  = cholG
        # Other variables
        self.num_locations = self.num_fixed + self.num_random
        self.sigma2_b_idx = self.num_locations
        self.sigma2_e_idx = self.num_locations+1
        self.num_parameters = self.num_locations + 2
        # Sampling stats
        self.samples = None
        self.num_samples = 0
        # Restore option
        self.restore_size = 0
        self.restore_buffer = None
        self.restore_idx = 0
        self.restore_used = 0
        # debug
        self.num_debug = 0
        self.utinvu_idx = 0
        self.u_df_idx = 1
        self.u_scale_idx = 2
        self.eeT_idx = 3
        self.e_df_idx = 4
        self.e_s_idx = 5

    def _prepare_simulate(self):
        # Priors
        # scaled inverted chi-square. Tuple (v,s2) = (degrees freedom,scale)
        self.prior_b = gt.create_chisquare_si(1.0,0.5)
        self.prior_e = gt.create_chisquare_si(1.0,0.5)

        # Parameters - theta
        self.theta = np.asmatrix(np.zeros(self.num_parameters).reshape(self.num_parameters,1))
        self.loc_hat = np.asmatrix(np.zeros(self.num_locations).reshape(self.num_locations,1))
        # used for Garcia-Cortes & Sorensen method
        self.loc_star = np.asmatrix(np.zeros(self.num_locations).reshape(self.num_locations,1))
        self.y_star = np.asmatrix(np.zeros(self.num_y).reshape(self.num_y,1))
        self.loc_tilde = np.asmatrix(np.zeros(self.num_locations).reshape(self.num_locations,1))
        ## Aliases
        self.B = np.matrix(self.theta[0:self.num_fixed,0],copy=False)
        self.B_hat = np.matrix(self.loc_hat[0:self.num_fixed,0],copy=False)
        # self.B = np.asmatrix(np.zeros(self.num_fixed).reshape(self.num_fixed,1))
        self.b = np.matrix(self.theta[self.num_fixed:self.num_locations,0],copy=False)
        self.b_hat = np.matrix(self.loc_hat[self.num_fixed:,0],copy=False)
        # self.b = np.asmatrix(np.zeros(self.num_random).reshape(self.num_random,1))
        ## Initial Values
        # B and b initialised to 0 - is this correct for B?
        self.B[0,0] = 10
        self.B[1,0] = 1.5
        ### Variance components
        # one sigma for each random effect
        self.theta[self.sigma2_b_idx,0] = 0.5
        # one sigma for each random effect
        self.theta[self.sigma2_e_idx,0] = 0.5
        # Restore buffer
        self.restore_initialize()
        # other initializations
        self._compute_mm_equations()


    # It should return an MCMC-run object (class to be implemented)
    def simulate(self,num_samples=10000,burnin=1000,thin=100,restore=20):
        # MCMC runnning parameters
        if (burnin < 0):
            print "Incorrect burn-in value"
            return None
        if (num_samples < 1):
            print "Incorrect number of samples value"
            return None
        # Preparations
        self.restore_size = restore
        self._prepare_simulate()
        # Output
        self.samples = np.zeros(num_samples*self.num_parameters).reshape(self.num_parameters,num_samples)
        self.samples = np.matrix(self.samples,copy=False)
        # Debugging
        if (self.num_debug>0):
            self.debugm = np.zeros(num_samples*self.num_debug).reshape(self.num_debug,num_samples)
            self.debugm = np.matrix(self.debugm,copy=False)
        #print 'Initial Value'
        #print self.theta.reshape(1,self.num_parameters)
        print '*** Burn-in ***'
        self.num_samples = 0
        try:
            for i in range(0,burnin):
                self._update_mm_equations()
                self._sample_gibbs_animal()
        except Exception, e:
            print 'Error during burn-in'
            print e
            return -1
        self.samples[:,0] = self.theta[:,0]
        print '*** Sampling ***'
        try:
            for i in range(1,num_samples):
                # Blocked updates
                #print "Sample:"+str(i)
                for j in range(0,thin):
                    self._update_mm_equations()
                    self._sample_gibbs_animal()
                self.samples[:,i] = self.theta[:,0]
                #print '--> sampled: '+str(i)
                if (self.num_debug>0):
                    self.debugm[self.utinvu_idx,i] = self.utinvu
                    self.debugm[self.u_df_idx,i] = self.u_df
                    self.debugm[self.u_scale_idx,i] = self.u_scale
                    self.debugm[self.eeT_idx,i] = self.eeT
                    self.debugm[self.e_df_idx,i] = self.e_df
                    self.debugm[self.e_s_idx,i] = self.e_scale
        except Exception, e:
            print 'Error while Sampling'
            print e
            self.num_samples = i
            return -1
        self.num_samples = num_samples
        return 0


    # Samples location parameters one at a time using univariate normal random distributions
    def _sample_gibbs_animal(self):
        # Sample location parameters
        self._sample_effects_block_chol()
        #self._sample_effects_block()

        print self.theta[0:self.num_locations,0]
        exit(0)

        # Sample variance components
        self._sample_variance_components()

        # -- admin --
        #print self.theta.reshape(1,self.num_parameters)
        # check for nan or inf
        if (np.any(map((lambda x: np.isnan(x) or np.isinf(x)),self.theta))):
            print self.theta
            raise Exception("NAN/INF")

    # Variance components are sampled using scaled inverse chi-square distributions
    def _sample_effects(self):
        self._update_mm_equations()
        C = self.mmeq_C
        #print 'num_locations='+str(self.num_locations)
        for i in range(0,self.num_locations):
            self.theta[i,0] = 0
            #print np.shape(C[i,:])
            #print np.shape(self.theta[0:self.num_locations,:])
            mu_theta = (C[i,:] * self.theta[0:self.num_locations,:])[0,0]
            mu_theta = (self.mmeq_rhs[i] - mu_theta)/C[i,i]
            s2_theta = self.theta[self.sigma2_e_idx,0]/C[i,i]
            self.theta[i,0] = np.random.normal(mu_theta,s2_theta)

    def _sample_effects_block(self):
        self._update_mm_equations()
        C = self.mmeq_C
        invC = np.linalg.inv(C)
        # Compute location estimators - expected values
        # Equation C * loc = rhs -> loc = C^{-1}*rhs
        self.loc_hat = invC * self.mmeq_rhs

        # Sample from multivariate normal
        self.theta[0:self.num_locations,0] = gt.rmvn(self.loc_hat,V=invC,sigma=self.theta[self.sigma2_e_idx])

    # Cholesky factorization
    def _sample_effects_block_chol(self):
        self._update_mm_equations()
        C = self.mmeq_C
        L = np.linalg.cholesky(C)
        #L = np.asmatrix(scipy.linalg.cholesky(C,lower=True))
        print 'before solve'
        print L
        print self.mmeq_rhs
        self.loc_hat = np.asmatrix(scipy.linalg.cho_solve((L,True),self.mmeq_rhs))
        print "first loc hat"
        print self.loc_hat
        #print 'after solve'
        # invChol = np.linalg.inv(cholC)
        sigma2e = self.theta[self.sigma2_e_idx]
        self.theta[0:self.num_locations,0] = gt.rmvn(self.loc_hat,cholL=L,inverse=True,sigma=sigma2e)
        print 'sample'
        print self.theta[0:self.num_locations,0]
        exit(0)
        #print 'after rvn'

    def _sample_effects_block_garcia(self):
        self._update_mm_equations()
        C = self.mmeq_C
        # sample loc_star
        # -- fixed locations are left as 0 since they don't affect result (GC&S claim)
        sigma2b = self.theta[self.sigma2_b_idx]
        self.loc_star[self.num_fixed:,0] = gt.rmvn(0,cholL=self.cholG,sigma=sigma2b)
        # Sample y*
        sigma2e = self.theta[self.sigma2_e_idx]
        Ie = np.matrix(np.identity(self.num_y), copy=False)
        e_star = gt.rmvn(0,cholLU=Ie,lower=True,sigma=sigma2e)
        self.y_star = self.W*self.loc_star + e_star
        cholC = np.asmatrix(scipy.linalg.cholesky(C,lower=True))
        rhs = self.W.T * (self.y - self.y_star)
        self.loc_tilde =  np.asmatrix(scipy.linalg.cho_solve((cholC,True),rhs))
        # Sampling result
        self.theta[0:self.num_locations,0] = self.loc_star -  self.loc_tilde

    def _sample_variance_components(self):
        ## Sample sigma2_b
        df = self.prior_b['df'] + self.num_random
        #print 'num_fixed: '+str(self.num_fixed)+' num_locations: '+ str(self.num_locations)
        #print 'theta:'+str(np.shape(self.theta))
        u = self.theta[self.num_fixed:self.num_locations,0]
        #print np.shape(u)
        #print np.shape(self.invG)
        try:
            scale1 = (u.T * self.invG * u)[0,0]
            #print 'u.T * invG * u = ' + str(scale0)
            scale = (scale1 + self.prior_b['df']*self.prior_b['scale']) / df
            #print 'scale = ' + str(scale) + ' scalebig = ' + str(scalebig)
            self.utinvu = scale1
            self.u_df = df
            self.u_scale = scale
            self.theta[self.sigma2_b_idx,0] = gt.rchisquare_inv_scaled(df,scale)
            #print 'sigma_b: '+str(x)+ ' df='+str(df)+' scale='+str(scale)
        except Exception, e:
            print "Error while sampling sigma2_b"
            print "df="+str(df)+" scale="+str(scale)
            print e
            raise e
        ## Sample sigma2_e
        #print 'sigma_e  index: '+ str(self.sigma2_e_idx)
        try:
            df = self.prior_e['df'] + self.num_y
            e = self.X*self.theta[0:self.num_fixed,0] + self.Z*self.theta[self.num_fixed:self.num_locations]
            e = self.y - e
            scale = ((e.T * e)[0,0] + self.prior_e['df']*self.prior_e['scale']) / df
            self.eeT = (e.T * e)[0,0]
            self.e_df = df
            self.e_scale = scale
            #print 'df '+ str(df) + ' scale ' + str(scale)
            self.theta[self.sigma2_e_idx] = gt.rchisquare_inv_scaled(df,scale)
            #print 'sigma_e: '+str(x)+ ' sample_e'+str(self.sample)
        except Exception, e:
            print "Error while sampling sigma2_e"
            print "df="+str(df)+" scale="+scale
            print e
            raise e

    # Mixed Model Equations
    def _compute_mm_equations(self):
        # C beta_hat = r
        # C dimension = (self.num_fixed+self.num_random)^2
        # we trust in C
        p = self.num_fixed
        pq = self.num_locations
        self.mmeq_C = np.asmatrix(np.zeros(pq*pq).reshape(pq,pq))
        self.mmeq_rhs = np.asmatrix(np.zeros(pq).reshape(pq,1))
        self.mmeq_C22 = np.matrix(self.mmeq_C[p:,p:],copy=False)
        # Right-hand side - assuming R is identity matrix
        self.mmeq_rhs[0:p,0] = self.X.T * self.y
        self.mmeq_rhs[p:,0] = self.Z.T * self.y
        # Coefficient matrix
        self.mmeq_C[0:p,0:p] = self.X.T * self.X
        self.mmeq_C[0:p,p:] = self.X.T * self.Z
        self.mmeq_C[p:,0:p] = self.Z.T * self.X
        self.ZZT = self.Z.T * self.Z
        #self._update_mm_equations() - executed before each variance sampling



    # Updates/Computes bottom-right corner of the coefficient matrix of the model equations
    def _update_mm_equations(self):
        c = self.theta[self.sigma2_e_idx,0] / self.theta[self.sigma2_b_idx,0]
        self.mmeq_C[self.num_fixed:,self.num_fixed:] = self.ZZT + c*self.invG

    def plot(self,last=None):
        #n = len(params)
        num_samples = self.num_samples
        plt.figure(1)
        #print "NUM SAMPLES="+str(num_samples)
        if (last):
            x = range(0,5)
            samples = self.samples[:,num_samples-5:num_samples]
            #print x
            #print samples
        else:
            x = range(0,num_samples)
            samples = self.samples[:,0:num_samples]

        print 'MEANS:'
        print np.mean(self.samples[0,:])
        print np.mean(self.samples[1,:])
        print np.mean(self.samples[self.sigma2_b_idx])
        print np.mean(self.samples[self.sigma2_e_idx])

        plt.subplot(411)
        y = samples[0,:].tolist()[0]
        plt.plot(x,y)
        plt.subplot(412)
        y = samples[1,:].tolist()[0]
        plt.plot(x,y)

        plt.subplot(413)
        y = samples[self.sigma2_b_idx,:].tolist()[0]
        plt.plot(x,y)
        plt.subplot(414)
        y = samples[self.sigma2_e_idx,:].tolist()[0]
        plt.plot(x,y)
        # show
        plt.show()

    def save_to_file(self,filename,filetype='cvs'):
        if (self.num_debug > 0):
            s = np.vstack((self.samples[:,0:self.num_samples],self.debugm[:,0:self.num_samples]))
        else:
            s = self.samples[:,0:self.num_samples]
        if (filetype=='cvs'):
            delim=','
        elif (filetype=='xls'):
            delim='\t'
        else:
            raise Exception('Incorrect filetype(cvs,xls)')
        np.savetxt(filename,s.T,delimiter=delim,fmt='%.4f')

    # **************** Restore structures ************************
    # ***** probably not useful *****

    def restore_initialize(self):
        self.restore = np.zeros(self.restore_size*self.num_parameters)
        self.restore.reshape(self.num_parameters,self.restore_size)
        self.restore_idx = 0
        self.restore_used = 0

    # Save theta into restore buffer
    def restore_save(self):
        self.restore_buffer[:,self.restore_idx] = self.theta[:,0]
        self.restore_dx += 1
        # Check if wrap-around applies
        if (self.restore_idx >= self.restore_size):
            self.restore_idx = 0
        # Keep track of number of buckets used
        if (self.restore_used < self.restore_size):
            self.restore_used += 1

    # Restore theta to past state
    def restore_do(self):
        if (self.restore_used == 0):
            return -1
        if (self.restore_used < self.restore_size):
            idx = 0
        else:
            idx = self.restore_idx - 1
            if (idx < 0):
                idx = self.restore_size -1
        self.theta = self.restore_buffer[:,idx]
        self.restore_idx = 0
        self.restore_used = 0
