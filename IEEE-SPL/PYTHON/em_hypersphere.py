import numpy as np
from scipy.special import iv

# From the paper "Hypersphere Fitting from Noisy Data Using an EM Algorithm" by Julien Lesouple, Barbara Pilastre,
# Yoann Altmann, and Jean-Yves Tourneret, accepted in IEEE Signal Processing Letters

# By Julien Lesouple (https://perso.tesa.prd.fr/jlesouple)
# 8 Jan. 2021


class em_hypersphere(object):
    
    def __init__(self, n_iter = 100):
        self.n_iter = n_iter
        
    def fit_hypersphere(self, Z, kappa_prior = 0, mu_prior = None):
        """Fits the hypersphere using observations Z and von Mises-Fisher prior with parmeters 
        kappa_prior (concentration parameter) and mu_prior (mean direction). If kappa_prior = 0, it is equivalent to a uniform prior."""

        [n, d] = Z.shape
        
        if mu_prior is None :
            if kappa_prior == 0:
                mu_prior = np.zeros([d,1])
                print("Uniform prior is used")
            else:
                raise ValueError("Please provide also a prior for mu")  
        else:
            print("Von-Mises Fisher prior used with kappa =", kappa_prior,"and mu =",mu_prior)
        
        
        # Initialization
        c_em = np.mean(Z,axis=0)
        r_em = np.linalg.norm(Z-c_em,axis = 1).mean()
        sigma2_em = (1/d)*((np.linalg.norm(Z-c_em,axis = 1)**2).mean()-r_em**2)
        r_mem = np.zeros([self.n_iter+1])
        c_mem = np.zeros([d,self.n_iter+1])
        sigma2_mem = np.zeros([self.n_iter+1])
        
        r_mem[0] = r_em
        c_mem[:,0] = c_em
        sigma2_mem[0] = sigma2_em
        
        for k in range(self.n_iter):
            # E-step
            delta_i = r_em*(Z-c_em).transpose()+sigma2_em*kappa_prior*mu_prior
            mu_i = delta_i/np.linalg.norm(delta_i,axis = 0)
            kappa_i = np.linalg.norm(delta_i,axis = 0)/sigma2_em

            alpha_i = iv(d/2,kappa_i)/iv(d/2-1,kappa_i)*mu_i

            f = np.r_[np.sum(alpha_i.transpose()*Z),np.sum(Z,axis=0)]

            alpha = np.sum(alpha_i, axis = 1)
            
            # M-step
            H_inv = 1/(n**2-np.inner(alpha,alpha))*np.r_[np.r_[n,-alpha].reshape(1,d+1),\
                np.c_[-alpha,((n**2-np.inner(alpha,alpha))/n)*np.eye(d)+(1/n)*np.outer(alpha,alpha)]]
            
            sol_em = np.dot(H_inv,f)
            r_em = sol_em[0]    
            c_em = sol_em[1:]

            M = np.sum(np.linalg.norm(Z,axis=1)**2)-np.dot(f,sol_em)
            sigma2_em = M/(n*d)
            
            r_mem[k] = r_em
            c_mem[:,k] = c_em
            sigma2_mem[k] = sigma2_em
            
        return r_em, c_em, sigma2_em, r_mem, c_mem, sigma2_mem
            
            
            