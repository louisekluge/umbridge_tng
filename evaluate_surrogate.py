import numpy as np
import matplotlib.pyplot as plt
import umbridge
import tinyDA as tda

import scipy.stats as stats
import pickle

link = "http://0.0.0.0:4242"

np.random.seed(987)

model = umbridge.HTTPModel(link, "forward_level2")
my_model = tda.UmBridgeModel(model)

a1, a2 = np.log(0.1), np.log(10)
b1, b2 = np.log(0.74), np.log(74)
c1, c2 = np.log(0.000001), np.log(0.001)

prior_WindEn = stats.uniform(loc=a1, scale=a2-a1) 
prior_WindVel = stats.uniform(loc=b1,scale=b2-b1)
prior_SeedBH = stats.uniform(loc=c1, scale=c2-c1)

my_prior = tda.CompositePrior([prior_WindEn, prior_WindVel, prior_SeedBH])

Omega_obs = np.array([0.000480, 0.000713, 0.000899])
cov_likelihood = (10**(-8))*np.eye(len(Omega_obs))

my_loglike_fine = tda.GaussianLogLike(Omega_obs, cov_likelihood)
my_loglike_coarse = tda.AdaptiveGaussianLogLike(Omega_obs, cov_likelihood)
#my_loglike_coarse = tda.GaussianLogLike(Omega_obs, cov_likelihood)

fiducial_param = np.array([np.log(3.6), np.log(7.4), np.log(0.00008)])

# Adaptive Metropolis
factor = 7
a = np.array([(a2-a1)/factor, (b2-b1)/factor, (c2-c1)/factor])
am_cov = np.diag(a)
am_t0 = 90
am_sd = None
am_epsilon = 1e-6
am_adaptive = True
my_proposal_ad = tda.AdaptiveMetropolis(C0=am_cov, t0=am_t0, sd=am_sd, epsilon=am_epsilon)

gp_multi = pickle.load(open('/home/hd/hd_hd/hd_uv175/gpfs/hd_uv175-uq_proj/gp_surrogate/gp_multioutput_unnorm.pkl', 'rb'))
class GPMulti:
    def __init__(self,gp):
        self.gp = gp
    
    def __call__(self, parameters):
        norm = 10*1.274*25736.791845552525
        out = self.gp.predict([parameters])
        return(np.array([out[0][0]/norm, out[0][1]/norm, out[0][2]/norm]))
    
my_surrogate = GPMulti(gp_multi)
my_posterior_coarse = tda.Posterior(my_prior, my_loglike_coarse, my_surrogate)
my_posterior_fine = tda.Posterior(my_prior, my_loglike_fine, my_model)

my_posteriors = [my_posterior_coarse, my_posterior_fine]

my_chain = tda.sample(my_posteriors, my_proposal_ad, iterations=50, n_chains=1, initial_parameters=fiducial_param, subsampling_rate=20)#, adaptive_error_model='state-dependant')

with open('test_surrogate_no_errormodel_50steps_3.pkl', 'wb') as f:
    pickle.dump(my_chain, f)
