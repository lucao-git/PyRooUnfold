#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Bias_study
    
This provides tools to compare the unfolded results from different methods
and optimize the corresponding paramters.

For each method, a 'method_fom' function is defined to get all figures of merit (fom) as individual array for later plotting.

Five types of fom defined in 'study_complex_errors'. More definitions can be added/replaced there if needed.

To do:
    1) maybe better to merge all 'method_fom' into one function
    2) add fom for Bayes method
"""

import numpy as np
from pyroounfold.utils.roo_convertor import th1_to_arr, tvec_to_arr
from pyroounfold.utils.unfold_methods import do_unfold
from scipy.stats import norm
import ROOT


def study_complex_errors(unf, true, unf_cov):
    """ Define and calculate figures of merit
        
        Args:
        unf : array of unfolded central value
        true : array of true central value
        unf_cov : covariance matrix (TMatrixT) after unfolding, directly passed from RooUnfold
        
        
        Returns:
        err_a : sum of absolute biases
        err_b : sum of normalised absolute biases
        err_c : sum of biases
        err_d : total stat error after unfolding taking into account bin-to-bin correlations
        err_e : ratio of sumed absoluted biases and total error after unfolding
        err_f : sum of two terms in type e
        err_g : sum of bin-wise ratio between absolute bias and unfolding error
        
        """
   
    nbins = len(unf)
    bias_abs = np.abs((unf - true))
    bias = (unf - true)
    
    su_cov = 0
    bias_abs_nor = np.zeros(nbins)
    bias_ratio = np.zeros(nbins)
    
    for x in range(nbins):
        if true[x]>0 : ## physical region
            bias_abs_nor[x] = np.abs((unf[x] - true[x]))/true[x]
        for y in range(nbins):
            su_cov += unf_cov[x][y]
            if x==y:
                bias_ratio[x]=bias_abs[x]/np.sqrt(unf_cov[x][y])
    
    err_a = np.sum(bias_abs)
    err_b = np.sum(bias_abs_nor)
    err_c = np.sum(bias)
    err_d = np.sqrt(su_cov)
    err_e = err_a/np.sqrt(su_cov)
    err_f = np.sqrt(su_cov + err_a**2)
    err_g = np.sum(bias_ratio)


    return err_a, err_b, err_c, err_d, err_e, err_f, err_g




def svd_fom(hist_reco, hist_true, hist_respon):
    """ Interally run through all regularisation parameters (k=2..nbins, k=0 using default nbins/2, k=1 not good to use) and do unfolding with SVD method for each k. Calculate figures of merit for comparasion.
        
        Args:
        hist_reco : (TH1D) measured distrbution to be unfolded
        hist_true : (TH1D) true distrbution
        hist_respon : (RooUnfoldResponse) response matrix (x=reco, y=true)
        
        Returns:
        unf_cen_all : 2D-array of central values after unfolding, i.e. unf_cen_all[k, bin]
        unf_err_all : 2D-array of stat. error after unfolding, i.e. unf_err_all[k, bin]
        unf_coverage_all : 2D-array of bin-wise converage
        err_a : sum of absolute biases
        err_b : sum of normalised absolute biases
        err_c : sum of biases
        err_d : total stat error after unfolding
        err_e : ratio of sumed absoluted biases and total stat error after unfolding
        err_f : sum of two terms in type e
        err_g : sum of bin-wise ratio between absolute bias and unfolding error
        err_h : averaged coverage
        k_arr : array of tested k
        """
    
    nbins = hist_reco.GetNbinsX()
    unf_cen_all = np.zeros([nbins-1, nbins])
    unf_err_all = np.zeros([nbins-1, nbins])
    unf_coverage_all = np.zeros([nbins-1, nbins])
    err_a=[]
    err_b=[]
    err_c=[]
    err_d=[]
    err_e=[]
    err_f=[]
    err_g=[]
    err_h=[]
    k_arr=[]
    
    for x in range( 2, nbins+1):
        
        unfsvd = ROOT.RooUnfoldSvd(hist_respon, hist_reco, x)
        witherror = ROOT.RooUnfold.kCovariance
        unfsvd.IncludeSystematics(3)
        unfres = unfsvd.Hunfold(witherror)
        unfcov = unfsvd.Eunfold(witherror)
        
        unf_cen, unf_err = th1_to_arr(unfres)
        mc_cen, _ = th1_to_arr(hist_true)
        coverage_perbin = tvec_to_arr(unfsvd.CoverageProbV(1))
        coverage_avr = sum(coverage_perbin)/len(coverage_perbin)
        
        svd_err_a, svd_err_b, svd_err_c, svd_err_d, svd_err_e, svd_err_f, svd_err_g = study_complex_errors(unf_cen, mc_cen, unfcov)
        err_a.append(svd_err_a)
        err_b.append(svd_err_b)
        err_c.append(svd_err_c)
        err_d.append(svd_err_d)
        err_e.append(svd_err_e)
        err_f.append(svd_err_f)
        err_g.append(svd_err_g)
        err_g.append(svd_err_g)
        err_h.append(coverage_avr)
        k_arr.append(x)
        
        for y in range(nbins):
            unf_cen_all[x-2, y] = unf_cen[y]
            unf_err_all[x-2, y] = unf_err[y]
            unf_coverage_all[x-2, y] = coverage_perbin[y]

    return unf_cen_all, unf_err_all, unf_coverage_all, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h, k_arr


def ids_fom(hist_reco, hist_true, hist_respon, ntry=None):
    """ Interally run through a list iteration parameters (n=1..ntry, default ntry=nbins) and do unfolding with IDS method for each n. Calculate figures of merit for comparasion.
        
        Args:
        hist_reco : (TH1D) measured distrbution to be unfolded
        hist_true : (TH1D) true distrbution
        hist_respon : (RooUnfoldResponse) response matrix (x=reco, y=true)
        ntry : the max iteration parameter, default is #bins
        
        Returns:
        unf_cen_all : 2D-array of central values after unfolding, i.e. unf_cen_all[n, bin]
        unf_err_all : 2D-array of stat. error after unfolding, i.e. unf_err_all[n, bin]
        unf_coverage_all : 2D-array of bin-wise converage
        err_a : sum of absolute biases
        err_b : sum of normalised absolute biases
        err_c : sum of biases
        err_d : total stat error after unfolding
        err_e : ratio of sumed absoluted biases and total stat error after unfolding
        err_f : sum of two terms in type e
        err_g : sum of bin-wise ratio between absolute bias and unfolding error
        err_h : averaged coverage
        n_arr : array of tested n
        """
    nbins = hist_reco.GetNbinsX()
    if ntry is None: ntry = nbins
    unf_cen_all = np.zeros([ntry, nbins])
    unf_err_all = np.zeros([ntry, nbins])
    unf_coverage_all = np.zeros([ntry, nbins])
    err_a=[]
    err_b=[]
    err_c=[]
    err_d=[]
    err_e=[]
    err_f=[]
    err_g=[]
    err_h=[]
    n_arr=[]
    
    
    for x in range( 1, ntry+1):
        
        unfids = ROOT.RooUnfoldIds(hist_respon, hist_reco, x)
        witherror = ROOT.RooUnfold.kCovariance
        unfids.IncludeSystematics(3)
        unfres = unfids.Hunfold(witherror)
        unfcov = unfids.Eunfold(witherror)
        
        unf_cen, unf_err = th1_to_arr(unfres)
        mc_cen, _ = th1_to_arr(hist_true)
        coverage_perbin = tvec_to_arr(unfids.CoverageProbV(1))
        coverage_avr = sum(coverage_perbin)/len(coverage_perbin)
        
        ids_err_a, ids_err_b, ids_err_c, ids_err_d, ids_err_e, ids_err_f, ids_err_g = study_complex_errors(unf_cen, mc_cen, unfcov)
        err_a.append(ids_err_a)
        err_b.append(ids_err_b)
        err_c.append(ids_err_c)
        err_d.append(ids_err_d)
        err_e.append(ids_err_e)
        err_f.append(ids_err_f)
        err_g.append(ids_err_g)
        err_h.append(coverage_avr)
        n_arr.append(x)
        
        for y in range(nbins):
            unf_cen_all[x-1, y] = unf_cen[y]
            unf_err_all[x-1, y] = unf_err[y]
            unf_coverage_all[x-1, y] = coverage_perbin[y]

    return unf_cen_all, unf_err_all, unf_coverage_all, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h, n_arr

def bay_fom(hist_reco, hist_true, hist_respon, ntry=None):

    """ Interally run through a list iteration parameters (n=1..ntry, default ntry=nbins) and do unfolding with Bayes method for each n. Calculate figures of merit for comparasion.
    
        Args:
        hist_reco : (TH1D) measured distrbution to be unfolded
        hist_true : (TH1D) true distrbution
        hist_respon : (RooUnfoldResponse) response matrix (x=reco, y=true)
        ntry : the max iteration parameter, default is #bins
        
        Returns:
        unf_cen_all : 2D-array of central values after unfolding, i.e. unf_cen_all[n, bin]
        unf_err_all : 2D-array of stat. error after unfolding, i.e. unf_err_all[n, bin]
        unf_coverage_all : 2D-array of bin-wise converage
        err_a : sum of absolute biases
        err_b : sum of normalised absolute biases
        err_c : sum of biases
        err_d : total stat error after unfolding
        err_e : ratio of sumed absoluted biases and total stat error after unfolding
        err_f : sum of two terms in type e
        err_g : sum of bin-wise ratio between absolute bias and unfolding error
        err_h : averaged coverage
        n_arr : array of tested n
        """
        
    nbins = hist_reco.GetNbinsX()
    if ntry is None: ntry = nbins
    unf_cen_all = np.zeros([ntry, nbins])
    unf_err_all = np.zeros([ntry, nbins])
    unf_coverage_all = np.zeros([ntry, nbins])
    err_a=[]
    err_b=[]
    err_c=[]
    err_d=[]
    err_e=[]
    err_f=[]
    err_g=[]
    err_h=[]
    n_arr=[]
    
    
    for x in range( 1, ntry+1):
        
        unfbay = ROOT.RooUnfoldBayes(hist_respon, hist_reco, x)
        witherror = ROOT.RooUnfold.kCovariance
        unfbay.IncludeSystematics(3)
        unfres = unfbay.Hunfold(witherror)
        unfcov = unfbay.Eunfold(witherror)
        
        unf_cen, unf_err = th1_to_arr(unfres)
        mc_cen, _ = th1_to_arr(hist_true)
        coverage_perbin = tvec_to_arr(unfbay.CoverageProbV(1))
        coverage_avr = sum(coverage_perbin)/len(coverage_perbin)
        
        
        bay_err_a, bay_err_b, bay_err_c, bay_err_d, bay_err_e, bay_err_f, bay_err_g = study_complex_errors(unf_cen, mc_cen, unfcov)
        err_a.append(bay_err_a)
        err_b.append(bay_err_b)
        err_c.append(bay_err_c)
        err_d.append(bay_err_d)
        err_e.append(bay_err_e)
        err_f.append(bay_err_f)
        err_g.append(bay_err_g)
        err_h.append(coverage_avr)
        n_arr.append(x)
        
        for y in range(nbins):
            unf_cen_all[x-1, y] = unf_cen[y]
            unf_err_all[x-1, y] = unf_err[y]
            unf_coverage_all[x-1, y] = coverage_perbin[y]

    return unf_cen_all, unf_err_all, unf_coverage_all, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h, n_arr




def tuf_fom(hist_reco, hist_true, hist_respon):
    """ Calculate figures of merit for comparasion with TUnfold method.
        
        Args:
        hist_reco : (TH1D) measured distrbution to be unfolded
        hist_true : (TH1D) true distrbution
        hist_respon : (RooUnfoldResponse) response matrix (x=reco, y=true)
        
        Returns:
        unf_cen : 1D-array of central values after unfolding, i.e. unf_cen[bin]
        unf_err : 1D-array of stat. error after unfolding, i.e. unf_err[bin]
        unf_coverage_all : 1D-array of bin-wise converage
        err_a : sum of absolute biases
        err_b : sum of normalised absolute biases
        err_c : sum of biases
        err_d : total stat error after unfolding
        err_e : ratio of sumed absoluted biases and total stat error after unfolding
        err_f : sum of two terms in type e
        err_g : sum of bin-wise ratio between absolute bias and unfolding error
        err_h : averaged coverage
        """
    
    unftuf = ROOT.RooUnfoldTUnfold(hist_respon, hist_reco)
    witherror = ROOT.RooUnfold.kCovariance
    unftuf.IncludeSystematics(3)
    unfres = unftuf.Hunfold(witherror)
    unfcov = unftuf.Eunfold(witherror)
    
    unf_cen, unf_err = th1_to_arr(unfres)
    mc_cen, _ = th1_to_arr(hist_true)
    coverage_perbin = tvec_to_arr(unftuf.CoverageProbV(1))
    coverage_avr = sum(coverage_perbin)/len(coverage_perbin)
    
    tuf_err_a, tuf_err_b, tuf_err_c, tuf_err_d, tuf_err_e, tuf_err_f, tuf_err_g = study_complex_errors(unf_cen, mc_cen, unfcov)
    tuf_err_h = coverage_avr
    
    return unf_cen, unf_err, coverage_perbin, tuf_err_a, tuf_err_b, tuf_err_c, tuf_err_d, tuf_err_e, tuf_err_f, tuf_err_g, tuf_err_h


def inv_fom(hist_reco, hist_true, hist_respon):
    """ Calculate figures of merit for comparasion with matrix invert method.
        
        Args:
        hist_reco : (TH1D) measured distrbution to be unfolded
        hist_true : (TH1D) true distrbution
        hist_respon : (RooUnfoldResponse) response matrix (x=reco, y=true)
        
        Returns:
        unf_cen : 1D-array of central values after unfolding, i.e. unf_cen[bin]
        unf_err : 1D-array of stat. error after unfolding, i.e. unf_err[bin]
        unf_coverage_all : 1D-array of bin-wise converage
        err_a : sum of absolute biases
        err_b : sum of normalised absolute biases
        err_c : sum of biases
        err_d : total stat error after unfolding
        err_e : ratio of sumed absoluted biases and total stat error after unfolding
        err_f : sum of two terms in type e
        err_g : sum of bin-wise ratio between absolute bias and unfolding error
        err_h : averaged coverage
        """
    
    unfinv = ROOT.RooUnfoldInvert(hist_respon, hist_reco)
    witherror = ROOT.RooUnfold.kCovariance
    unfinv.IncludeSystematics(3)
    unfres = unfinv.Hunfold(witherror)
    unfcov = unfinv.Eunfold(witherror)
    
    unf_cen, unf_err = th1_to_arr(unfres)
    mc_cen, _ = th1_to_arr(hist_true)
    coverage_perbin = tvec_to_arr(unfinv.CoverageProbV(1))
    coverage_avr = sum(coverage_perbin)/len(coverage_perbin)
    
    inv_err_a, inv_err_b, inv_err_c, inv_err_d, inv_err_e, inv_err_f, inv_err_g = study_complex_errors(unf_cen, mc_cen, unfcov)
    inv_err_h = coverage_avr
    return unf_cen, unf_err, coverage_perbin, inv_err_a, inv_err_b, inv_err_c, inv_err_d, inv_err_e, inv_err_f, inv_err_g, inv_err_h


def byb_fom(hist_reco, hist_true, hist_respon):
    """ Calculate figures of merit for comparasion with bin-by-bin method.
        
        Args:
        hist_reco : (TH1D) measured distrbution to be unfolded
        hist_true : (TH1D) true distrbution
        hist_respon : (RooUnfoldResponse) response matrix (x=reco, y=true)
        
        Returns:
        unf_cen : 1D-array of central values after unfolding, i.e. unf_cen[bin]
        unf_err : 1D-array of stat. error after unfolding, i.e. unf_err[bin]
        unf_coverage_all : 1D-array of bin-wise converage
        err_a : sum of absolute biases
        err_b : sum of normalised absolute biases
        err_c : sum of biases
        err_d : total stat error after unfolding
        err_e : ratio of sumed absoluted biases and total stat error after unfolding
        err_f : sum of two terms in type e
        err_g : sum of bin-wise ratio between absolute bias and unfolding error
        err_h : averaged coverage
        """
    
    unfbyb = ROOT.RooUnfoldBinByBin(hist_respon, hist_reco)
    witherror = ROOT.RooUnfold.kCovariance
    unfbyb.IncludeSystematics(3)
    unfres = unfbyb.Hunfold(witherror)
    unfcov = unfbyb.Eunfold(witherror)
    
    unf_cen, unf_err = th1_to_arr(unfres)
    mc_cen, _ = th1_to_arr(hist_true)
    coverage_perbin = tvec_to_arr(unfbyb.CoverageProbV(1))
    coverage_avr = sum(coverage_perbin)/len(coverage_perbin)
    
    byb_err_a, byb_err_b, byb_err_c, byb_err_d, byb_err_e, byb_err_f, byb_err_g = study_complex_errors(unf_cen, mc_cen, unfcov)
    byb_err_h = coverage_avr
    return unf_cen, unf_err, coverage_perbin, byb_err_a, byb_err_b, byb_err_c, byb_err_d, byb_err_e, byb_err_f, byb_err_g, byb_err_h



def get_fom(hist_reco, hist_true, hist_respon, method=None, ite_ntry=None):
    
    #unfold_cen, unfold_err, err_a, err_b, err_c, err_d, err_e, err_f, err_g = None
    if method is None : print('Please indicate one method for unfolding: \'Ids\', \'Svd\', \'Bayes\', \'TUnfold\', \'Invert\', \'BinByBin\'.')
    
    elif method=='Ids':
        print('Use IDS method for error study.')
        if ite_ntry is None: print('No input max iteration parameter. Error study compares default n = 1..nbins')
        else: print('Error study compares iteration parameter n = 1..' + str(ite_ntry) + '.')
        unfold_cen, unfold_err, coverage_perbin, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h, k_arr = ids_fom(hist_reco, hist_true, hist_respon, ite_ntry)
        return unfold_cen, unfold_err, err_a, err_b, err_c, err_d, err_e, err_f, err_g, k_arr

    elif method=='Svd':
        print('Use SVD method for error study to compare regularisation parameter k = 2..nins.')
        unfold_cen, unfold_err, coverage_perbin, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h, n_arr = svd_fom(hist_reco, hist_true, hist_respon)
        return unfold_cen, unfold_err, err_a, err_b, err_c, err_d, err_e, err_f, err_g, n_arr
                              
    elif method=='Bayes':
        print('Use Bayes method for error study.')
        if ite_ntry is None: print('No input max iteration parameter. Error study compares default n = 1..nbins')
        else: print('Error study compares iteration parameter n = 1..' + str(ite_ntry) + '.')
        unfold_cen, unfold_err, coverage_perbin, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_g, err_h, n_arr = bay_fom(hist_reco, hist_true, hist_respon, ite_ntry)
        return unfold_cen, unfold_err, err_a, err_b, err_c, err_d, err_e, err_f, err_g, n_arr

    elif method=='Invert':
        print('Use matrix invert method for error study.')
        unfold_cen, unfold_err, coverage_perbin, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h = inv_fom(hist_reco, hist_true, hist_respon)
        return unfold_cen, unfold_err, err_a, err_b, err_c, err_d, err_e, err_f, err_g
    
    elif method=='TUnfold':
        print('Use TUnfold method for error study.')
        unfold_cen, unfold_err, coverage_perbin, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h = tuf_fom(hist_reco, hist_true, hist_respon)
        return unfold_cen, unfold_err, err_a, err_b, err_c, err_d, err_e, err_f, err_g

    elif method=='BinByBin':
        print('Use bin-by-bin correction method for error study.')
        unfold_cen, unfold_err, coverage_perbin, err_a, err_b, err_c, err_d, err_e, err_f, err_g, err_h = byb_fom(hist_reco, hist_true, hist_respon)
        return unfold_cen, unfold_err, err_a, err_b, err_c, err_d, err_e, err_f, err_g
                              


def cov2corr(cov):
    """Calculates the correlation matrix from a given
    covariance matrix.

    Arguments
    ---------
    cov : np.ndarray
        Covariance matrix. Shape is (n,n).

    Return
    ------
    out : np.ndarray
        Correlation matrix. Shape is (n,n).
    """
    Dinv = np.nan_to_num(np.diag(1 / np.sqrt(np.diag(cov))))
    return np.matmul(Dinv, np.matmul(cov, Dinv))


def check_coverage(bias, err, sigma=1):
    return norm.cdf(bias/err + sigma) - norm.cdf(bias/err - sigma)
