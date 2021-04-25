#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Unfold_methodes
    
    This code wrapps all the methods provided by RooUnfold with key words:
    
    1) "Ids": Iterative, Dynamically Stabilized(IDS) --- RooUnfoldIds
    2) "Svd": Singular Value Decomposition(SVD) --- RooUnfoldSvd
    3) "Bayes": Iterative Bayes --- RooUnfoldBayes
    4) "TUnfold": Regularised matrix inversion --- RooUnfoldTUnfold
    5) "Invert": Unregularised matrix inversion --- RooUnfoldInvert
    6) "BinByBin": Bin-by-bin --- RooUnfoldBinByBin

    
    
    
    Note:
    1) By default, stat. error on reponse matrix is included, by
        unf.IncludeSystematics()
    2) 'ROOUNFOLD_PATH' is the path for libRooUnfold.so
    
    To do:
    1) add toy flag to use toys generated internally by RooUnfold
    2) add regularised least squart fit method (RooUnofldTUnfold)
    
"""
import ROOT
import os
ROOT.gSystem.Load(os.getenv("ROOUNFOLD_PATH"))
from pyroounfold.utils.roo_convertor import th1_to_arr, ndarr_to_tmatrix
import pandas as pd
import numpy as np


def do_unfold(hist_true, hist_measure, hist_respon, method=None, para=None, mea_cov='False', kcovtoy=False, mc_stat_err=0):
    """ do unfolding on a measured distribution with a response matrix.
        
        Args:
        hist_respon (RooUnfoldResponse) : response matrix (x=reco, y=true)
        hist_measure (TH1D):  measured distrbution and to be unfolded
        method  : string for unfold method: 'Ids', 'Svd', 'Bayes', 'TUnfold', 'Invert', 'BinByBin'
        para    : parameters for 'Ids', 'Svd' and 'Bayes' methods.
        mea_cov (optional) : measured covariance matrix, default is statistical covariance
        kcovtoy (optional) : flag provided by ROOUNFOLD. Default is False and the full covariance matrix 'reco_cov' propagated through unfolding. If True, the error propagation is based on toys generated internally by RooUnfold.
        mc_stat_err (optional) : relate to ROOUNFOLD::includeSystematics().  Default "0" is to leave out the effect of statistical uncertainties on the migration matrix. "1" is to include the effect. "2" is for only counting the statistical uncertainties of measured distribtuon and migration matrix. The effect is valueated by internal toys.
        
        
        Returns:
        df_unf :  dataframe including bin_index, truth, measured and unfolded result
        cov_array :  2D array, covariance matrix after unfolding
        """
    unfres = None
    unfcov = None
    
    if(kcovtoy):
        witherror = ROOT.RooUnfold.kCovToy  #  error propagation based on toys generated internally by RooUnfold
    else:
        witherror = ROOT.RooUnfold.kCovariance   #  error propagation based on full covariance matrix 'mea_cov'
    
        
    if method is None : print('Please indicate one method for unfolding: \'Ids\', \'Svd\', \'Bayes\', \'TUnfold\', \'Invert\', \'BinByBin\'.'
                              + "\n" + 'e.g. do_unfold(hist_respon, hist_measure, \'Svd\', 5)')
        
    elif method=='Ids':
        #print('Use IDS method with iteration number = '+ str(para) + '.')
        if para is None or para <0 : print('Ids method requires a iteration number (>=0).')
        elif para>=0 :
                    unf = ROOT.RooUnfoldIds(hist_respon, hist_measure, para)
                    if(incl_mig_err):
                        unf.IncludeSystematics()
                    if(mea_cov!='False'):
                        unf.SetMeasuredCov(ndarr_to_tmatrix(mea_cov))
                    unfres = unf.Hreco()
                    unfcov = unf.Ereco(witherror)


    elif method=='Svd':
        #print('Use SVD method with regularisation number = '+ str(para) + '.')
        if para is None or para <0 : print('Svd method require a regularisation number (<= nbins, 0 using default nbins/2).')
        elif para > hist_measure.GetNbinsX(): print('Svd method do not work when regularisation number > nbins.')
        elif para>= 0 & para <= hist_measure.GetNbinsX():
            unf = ROOT.RooUnfoldSvd(hist_respon, hist_measure, para)
            if(mc_stat_err > 0):
                unf.IncludeSystematics(mc_stat_err)
            if(mea_cov!='False'):
                unf.SetMeasuredCov(ndarr_to_tmatrix(mea_cov))
            unfres = unf.Hreco()
            unfcov = unf.Ereco(witherror)

    
    elif method=='Bayes':
        #print('Use iterative Bayes method with iteration number = '+ str(para) + '.')
        if para is None :
            print('Bayes method requires a iteration number (default 4).')
            para=4
        unf = ROOT.RooUnfoldBayes(hist_respon, hist_measure, para)
        if(mc_stat_err > 0):
            unf.IncludeSystematics(mc_stat_err)
        if(mea_cov!='False'):
            unf.SetMeasuredCov(ndarr_to_tmatrix(mea_cov))
        unfres = unf.Hreco()
        unfcov = unf.Ereco(witherror)

    elif method=='Invert':
        #print('Use matrix invert method.')
        if para is not None: print('Unregularised matrix inerson method does not need parameter. Input parameter was ignored.')
        unf = ROOT.RooUnfoldInvert(hist_respon, hist_measure)
        if(mc_stat_err > 0):
            unf.IncludeSystematics(mc_stat_err)
        if(mea_cov!='False'):
            unf.SetMeasuredCov(ndarr_to_tmatrix(mea_cov))
        unfres = unf.Hreco()
        unfcov = unf.Ereco(witherror)

    
    elif method=='TUnfold':
        #print('Use TUnfold method.')
        if para is not None: print('TUfold method does not need input parameter. Input parameter was ignored.')
        unf = ROOT.RooUnfoldTUnfold(hist_respon, hist_measure)
        if(mc_stat_err > 0):
            unf.IncludeSystematics(mc_stat_err)
        if(mea_cov!='False'):
            unf.SetMeasuredCov(ndarr_to_tmatrix(mea_cov))
        unfres = unf.Hreco()
        unfcov = unf.Ereco(witherror)
    
    
    elif method=='BinByBin':
        #print('Use bin-by-bin correction method.')
        if para is not None: print('Bin-by-bin correction method does not need parameter. Input parameter was ignored.')
        unf = ROOT.RooUnfoldBinByBin(hist_respon, hist_measure)
        if(mc_stat_err > 0):
            unf.IncludeSystematics(mc_stat_err)
        if(mea_cov!='False'):
            unf.SetMeasuredCov(ndarr_to_tmatrix(mea_cov))
        unfres = unf.Hreco()
        unfcov = unf.Ereco(witherror)

    else : print('Method not found! Please select one from \'Ids\', \'Svd\', \'Bayes\', \'Invert\', \'BinByBin\'.')

    # convert results to numpy format
    bins=[]
    nbins = hist_measure.GetNbinsX()
    cov_array = np.zeros((nbins, nbins))
    
    for x in range(nbins):
        for y in range(nbins):
            cov_array[x][y]= unfcov[x][y]

    df_unf = pd.DataFrame(range(0, nbins), columns=['bin_index'])
    df_unf['true_central'], df_unf['true_error'] =th1_to_arr(hist_true)
    df_unf['measured_central'], df_unf['measured_error'] = th1_to_arr(hist_measure)
    df_unf['unfolded_central'], df_unf['unfolded_error'] = th1_to_arr(unfres)

  
    return df_unf, cov_array
