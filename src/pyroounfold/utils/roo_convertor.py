#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Roo_convertor
    
    This convert numpy array into ROOT class: TH1, TH2, TMatrixD and RooUnfoldResponse
    or vice versa.

"""
import matplotlib.pyplot as plt
import numpy as np
import ROOT

def get_stat_err(df_sig, weight_sig, df_bkg, weight_bkg, var, bin_arr):
    """ Extract expected statistical uncertainty by combining signal and background
        
        Args:
        df_sig : signal dataframe
        weight_sig : array of event-weights for signal
        df_bkg : background dataframe
        weight_bkg : array of event-weights for background
        var : string of index in the dataframe for a varibale, e.g. 'true_q2'
        bin_arr: an array of binning (length = number of bins + 1)
        
        Returns:
        err : an array of expected statistical uncertainty for each bin
        """
    n_exp, bins_exp, _ = plt.hist([df_bkg[var], df_sig[var]],
                                  bins = bin_arr, normed=False, stacked=True,
                                  weights = [weight_bkg, weight_sig])
    plt.close()
    err = np.sqrt(n_exp[1])
    return err

def df_to_roounf(train, weight_train, test, weight_test, true_name, reco_name, bins, bkg=None, weight_bkg=None):
    """ Convert dataframe input to ROOT hist and RooUnfoldResponse
    
        Args:
        train : train dataframe
        weight_train : array of event-weights for train sample
        test : test dataframe
        weight_test : array of event-weights for test sample
        true_name : string of index in the dataframe for a true varibale, e.g. 'true_q2'
        reco_name : string of index in the dataframe for a reconstructed varibale, e.g. 'reco_q2'
        bins: an array of binning (length = number of bins + 1)
        bkg (optional) : background dataframe
        weight_bkg (optional) : array of event-weights for background sample
        
        Returns:
        hist_train_true : TH1D for true distrbution (event weighted) in train sample
        hist_train_reco : TH1D for reconstructed distrbution (event weighted) in train sample
        hist_respon : response matrix as RooUnfoldResponse, using train sample
        hist_test_true : TH1D for true distrbution (event weighted) in test sample
        hist_test_reco : TH1D for reconstructed distrbution (event weighted) in test sample
        
        To do:
        1) add function to use RooUnfoldResponse.Miss() for scaling response matrix by the efficiency of each truth bin
        """
        
    ROOT.TH1.AddDirectory(False)
    
    hist_train_true = ROOT.TH1D("train_"+true_name, "train_"+true_name, bins.size-1, bins)
    hist_train_reco = ROOT.TH1D("train_"+reco_name, "train_"+reco_name, bins.size-1, bins)
    hist_train_Adet = ROOT.TH2D("train_Adet_"+reco_name, "train_Adet_"+reco_name, bins.size-1, bins, bins.size-1, bins)
    
    hist_test_true = ROOT.TH1D("test_"+true_name, "test_"+true_name, bins.size-1, bins)
    hist_test_reco = ROOT.TH1D("data_"+reco_name, "test_"+reco_name, bins.size-1, bins)
    
    for x in range(len(train)):
        hist_train_true.Fill(train[true_name].values[x], weight_train.values[x])
        hist_train_reco.Fill(train[reco_name].values[x], weight_train.values[x])
        hist_train_Adet.Fill(train[reco_name].values[x], train[true_name].values[x], weight_train.values[x])
        
    for x in range(len(test)):
        hist_test_true.Fill(test[true_name].values[x], weight_test.values[x])
        hist_test_reco.Fill(test[reco_name].values[x], weight_test.values[x])
        
    # replace sqrt(#sig) to expected stat. error: sqrt(#sig + #bkg)
    if bkg is not None:
        for x in range(bins.size-1):
            err = get_stat_err(test, weight_test, bkg, weight_bkg, reco_name, bins )[x]
            hist_test_reco.SetBinError(x+1, err)

    hist_respon = ROOT.RooUnfoldResponse(hist_train_reco, hist_train_true, hist_train_Adet)

    return hist_train_true, hist_train_reco, hist_respon, hist_test_true, hist_test_reco


def th1_to_arr(hist):
    """ Convert TH1D histogram's central value and error to numpy arraies, which are easier for later plotting.
        
        Args:
        hist : TH1D histogram
        
        Returns:
        central : an array of central value for each bin
        error : an array of error for each bin
        
        To do:
        add a similar function for 2D hist
        
        """
    central = []
    error = []
    nbins = hist.GetNbinsX()
    
    for x in range(nbins):
        central.append(hist.GetBinContent(x+1))
        error.append(hist.GetBinError(x+1))

    return np.asarray(central), np.asarray(error)
    
    
def arr_to_th1(bins, cen, err='False'):
    """ Convert central value and error in the format of numpy arraies to TH1D histogram.
        Args:
        bins : binning information
        central : an array of central value for each bin
        error (optional) : an array of error for each bin
        
        Returns:
        hist : TH1D histogram
        
        """
    ROOT.TH1.AddDirectory(False)
    hist = ROOT.TH1D("h1", "h1", bins.size-1, bins)
    if err=='False': # default is statistical error
        err=np.sqrt(cen)
    for x in range(bins.size-1):
            hist.SetBinContent(x+1, cen[x])
            if err.size==bins.size-1:
                hist.SetBinError(x+1, err[x])
    return hist
    
    
def ndarr_to_tmatrix(cov):
    """ Convert numpy covariance matrix to TMatrixD.
        Args:
        cov : numpy ndarray (n,n)
        
        Returns:
        m : TMatrixD (n,n)
        
        """
    ndim = len(cov)
    m = ROOT.TMatrixD (ndim, ndim)
    
    for x in range(ndim):
        for y in range(ndim):
            m[x][y] = cov[x][y]
    return m


def ndarr_to_th2(cov):
    """ Convert numpy 2D array to TH2D histogram.
        Args:
        cov : numpy ndarray (n,n)
        
        Returns:
        m : TH2D (n,n)
        
        """
    ROOT.TH1.AddDirectory(False)
    
    ndim = len(cov)
    hist = ROOT.TH2D("h2", "h2", bins.size-1, bins, bins.size-1, bins)
    
    for x in range(ndim):
        for y in range(ndim):
            hist.SetBinContent(x+1, y+1, cov[x][y])
    return hist
