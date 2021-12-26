import numpy as np
import pandas as pd
import os, sys, time
import warnings
import matplotlib.pyplot as plt 
warnings.filterwarnings("ignore")
 
def sigma_clip(orbit_number, t, y_err, nsigma_clip = 5):
    niter_sigmaclip = 10

    outlier = np.full(t.shape[0], False, dtype=bool)
    idxs = []

    for i in range(niter_sigmaclip):
        n_outliers_before = np.sum(outlier)

        # learn this is a magical function - it makes exactly what we want for this design matrix
        X = np.vander(orbit_number, N=2, increasing=True)
        # OR:
        # X = np.vstack((np.ones_like(x), x)).T

        Cov = np.diag(y_err**2)
        Cinv = np.linalg.inv(Cov) # we need the inverse covariance matrix

        # get the parameter covariance matrix ...
        theta_Cov = np.linalg.inv(X.T @ Cinv @ X)

        # and the best parameters using the new Python matrix operator
        theta_best = theta_Cov @ (X.T @ Cinv @ t)
        #print ("a, b = {:.3f}, {:.3f}".format(*theta_best))

        # add MLE estimate
        X_ = np.vander(orbit_number, N=2, increasing=True)
        t_model = X_ @ theta_best

        t0, p = theta_best


        t_model = p*orbit_number + t0
        t_res = t-t_model
        scatter = np.std(t_res[~outlier])
        outlier = (np.abs(t_res/scatter) > nsigma_clip)
        idx = np.where(outlier)[0]
        idxs.append(idx)

        t0_unc = np.sqrt(np.diag(theta_Cov))[0]
        per_unc = np.sqrt(np.diag(theta_Cov))[1]

        if (np.sum(outlier) == n_outliers_before): break
    idxs = np.concatenate(idxs).ravel()
    idxs = np.unique(idxs)
 
    return idxs, p, t0, t0_unc, per_unc

 
def sigma_clip_quadratic(orbit_number, t, y_err, nsigma_clip = 5):
    niter_sigmaclip = 10

    outlier = np.full(t.shape[0], False, dtype=bool)
    idxs = []

    for i in range(niter_sigmaclip):
        n_outliers_before = np.sum(outlier)



        # learn this is a magical function - it makes exactly what we want for this design matrix
        X = np.vander(orbit_number, N=3, increasing=True)
        # OR:
        # X = np.vstack((np.ones_like(x), x)).T

        Cov = np.diag(y_err**2)
        Cinv = np.linalg.inv(Cov) # we need the inverse covariance matrix

        # get the parameter covariance matrix ...
        theta_Cov = np.linalg.inv(X.T @ Cinv @ X)

        # and the best parameters using the new Python matrix operator
        theta_best = theta_Cov @ (X.T @ Cinv @ t)
        #print ("a, b = {:.3f}, {:.3f}".format(*theta_best))

        # add MLE estimate
        X_ = np.vander(orbit_number, N=3, increasing=True)
        t_model = X_ @ theta_best

        t0, p, dpde_over_2 = theta_best


        t_model = t0 + p * orbit_number + dpde_over_2 * orbit_number**2
        t_res = t - t_model
        scatter = np.std(t_res[~outlier])
        outlier = (np.abs(t_res/scatter) > nsigma_clip)
        idx = np.where(outlier)[0]
        idxs.append(idx)

        t0_unc = np.sqrt(np.diag(theta_Cov))[0]
        per_unc = np.sqrt(np.diag(theta_Cov))[1]
        dpde_over_2_unc = np.sqrt(np.diag(theta_Cov))[2]


        #print('   Sigma-clipping iteration ', i, ': total number of outliers clipped = ', np.sum(outlier))
        if (np.sum(outlier) == n_outliers_before): break
    idxs = np.concatenate(idxs).ravel()
    idxs = np.unique(idxs)
 
    return idxs, p, dpde_over_2, t0, t0_unc, per_unc, dpde_over_2_unc


def fold_time(t, t0, period):
    # returns folded time in between -period/2 and +period/2
    # as well as an orbit number, assigning t0 = orbit zero

    orbit_number, t_fold = np.divmod(t-(t0-0.5*period), period)
    t_fold = t_fold - 0.5*period
    
    return orbit_number.astype(np.int)
 
def f(sigma0, t_after_sigma_clip, model, t_err_after_sigma_clip):
  n_dof = t_after_sigma_clip.shape[0] - 2
  return np.sum((t_after_sigma_clip-model)**2/(t_err_after_sigma_clip**2+sigma0**2)) - n_dof 



def find_root(t_after_sigma_clip, model, t_err_after_sigma_clip):
  a = 0                     # beginning of initial bracket
  b = 10                 # end of initial bracket
  fctr = 1.6                # factor for increasing bracket
  itrtnmax = 10000           # maximum number of iterations

  itrtn = 0                 # initialize iteration counter
  fa = f(a, t_after_sigma_clip, model, t_err_after_sigma_clip)                 # compute f(x=a)
  fb = f(b, t_after_sigma_clip, model, t_err_after_sigma_clip)                 # compute f(x=b)

  while (itrtn < itrtnmax): # start bracketing procedure
    if fa*fb < 0:           # (a,b) already brackets root
      break
    if abs(fa) < abs(fb):
      a += fctr*(a-b)       # move a uphill
      fa = f(a, t_after_sigma_clip, model, t_err_after_sigma_clip)             # recompute f(x=a)
    else:
      b += fctr*(b-a)       # move b downhill
      fb = f(b, t_after_sigma_clip, model, t_err_after_sigma_clip)             # recompute f(x=b)
    itrtn += 1
    print(itrtn,a,b)

  eps = 1e-10               # tolerance for bisection
  notdone = 1
  itrtn = 0                 # initialize iteration counter
  fstart = f(a, t_after_sigma_clip, model, t_err_after_sigma_clip)             # compute f at x=a

  while (notdone==1 and itrtn < itrtnmax):
    midp = 0.5*(a+b)        # compute midpoint in [a,b]
    fmid = f(midp, t_after_sigma_clip, model, t_err_after_sigma_clip)          # evaluate f at midpoint
    if (fmid==0):           # if f(midpoint)=0, success!
      root = midp
      break
    elif (fstart*fmid < 0): # if f changes sign between a
      b = midp              # and midpoint, set b = midpt
    else:
      a = midp              # otherwise, set a = midpt

    if (abs(b-a) > eps*a):  # if error is bigger than tolerance
      fstart = f(a, t_after_sigma_clip, model, t_err_after_sigma_clip)         # rinse and repeat
      itrtn += 1
    else:
      root = a              # otherwise accept as "root"
      notdone = 0

  #print(root,itrtn)
  return root

 