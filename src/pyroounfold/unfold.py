#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" unfold
    
This provides main unfolding class.

    The methods provided by RooUnfold:
    
    1) "Ids": Iterative, Dynamically Stabilized(IDS) --- RooUnfoldIds
    2) "Svd": Singular Value Decomposition(SVD) --- RooUnfoldSvd
    3) "Bayes": Iterative Bayes --- RooUnfoldBayes
    4) "TUnfold": Regularised matrix inversion --- RooUnfoldTUnfold
    5) "Invert": Unregularised matrix inversion --- RooUnfoldInvert
    6) "BinByBin": Bin-by-bin --- RooUnfoldBinByBin

"""

import ROOT
import os
ROOT.gSystem.Load(os.getenv("ROOUNFOLD_PATH"))
from pyroounfold.utils.roo_convertor import df_to_roounf, th1_to_arr, arr_to_th1
from pyroounfold.utils.bias_study import study_complex_errors
from pyroounfold.utils.generate_toys import *
import pyroounfold.plotstyle as ps
from pyroounfold.plotstyle.colors import PaperColors as p
from pyroounfold.plotting import get_migration

import pandas as pd
import numpy as np

class unfold:

    def __init__(self, df_train, weight_train, df_test, weight_test, name_var_true, name_var_reco, show_var, bins, reco_bin_error='False', reco_cov='False'):
    
        self.witherror = ROOT.RooUnfold.kCovariance
        
        
        self.hist_train_true, self.hist_train_measure, self.hist_respon, self.hist_test_true,  self.hist_test_measure = df_to_roounf(
        train = df_train,
        weight_train = weight_train,
        test = df_test,
        weight_test = weight_test,
        true_name = name_var_true,
        reco_name = name_var_reco,
        bins = bins
        )
        
        self.reco_bin_error = reco_bin_error
        self.reco_cov = reco_cov
        
        if reco_cov != 'False':
            self.reco_bin_error = np.sqrt(reco_cov.diagonal())
            self.reco_cor = cov2corr(reco_cov)
            
        if (reco_cov == 'False')&(reco_bin_error == 'False'):
            self.reco_bin_error = np.sqrt(self.hist_test_measure)
            
        self.bins = bins
        
        self.show_var = show_var
        
        self.nbins = self.hist_test_measure.GetNbinsX()
        self.result_df = pd.DataFrame(np.array(range(0,self.nbins)), columns=['bin_index'])
        
        self.result_cov =  np.zeros((self.nbins, self.nbins))
        
        if len(self.reco_bin_error)==self.nbins:
            for x in range(self.nbins):
                self.hist_test_measure.SetBinError(x+1, self.reco_bin_error[x])
        
        
        
    def do_Ids(self, para):
        if para is None or para <0 : print('ERROR: Ids method requires a iteration number (>=0).')
        elif para>=0 :
                    unf = ROOT.RooUnfoldIds(self.hist_respon, self.hist_test_measure, para)
                    unf.IncludeSystematics()
                    self.unfres = unf.Hreco()
                    self.unfres_cov = unf.Ereco(self.witherror)
                    self.unf_result_cov()
                    self.unf_result_df()
                    
    
    def do_Svd(self, para):
        if para is None or para <0 : print('ERROR: Svd method require a regularisation number (<= nbins, 0 using default nbins/2).')
        elif para > self.hist_test_measure.GetNbinsX(): print('ERROR: Svd method do not work when regularisation number > nbins.')
        elif para>= 0 & para <= self.hist_test_measure.GetNbinsX():
            unf = ROOT.RooUnfoldSvd(self.hist_respon, self.hist_test_measure, para)
            unf.IncludeSystematics()
            self.unfres = unf.Hreco()
            self.unfres_cov = unf.Ereco(self.witherror)
            self.unf_result_cov()
            self.unf_result_df()
            
    def do_Bayes(self, para):
        if para is None :
            print('ERROR: Bayes method requires a iteration number.')
        else:
            unf = ROOT.RooUnfoldBayes(self.hist_respon, self.hist_test_measure, para)
            unf.IncludeSystematics()
            self.unfres = unf.Hreco()
            self.unfres_cov = unf.Ereco(self.witherror)
            self.unf_result_cov()
            self.unf_result_df()
            
    def do_Invert(self):
        unf = ROOT.RooUnfoldInvert(self.hist_respon, self.hist_test_measure)
        unf.IncludeSystematics()
        self.unfres = unf.Hreco()
        self.unfres_cov = unf.Ereco(self.witherror)
        self.unf_result_cov()
        self.unf_result_df()
        
    def do_TUnfold(self):
        unf = ROOT.RooUnfoldTUnfold(self.hist_respon, self.hist_test_measure)
        unf.IncludeSystematics()
        self.unfres = unf.Hreco()
        self.unfres_cov = unf.Ereco(self.witherror)
        self.unf_result_cov()
        self.unf_result_df()
    
    def do_BinByBin(self):
        unf = ROOT.RooUnfoldBinByBin(self.hist_respon, self.hist_test_measure)
        unf.IncludeSystematics()
        self.unfres = unf.Hreco()
        self.unfres_cov = unf.Ereco(self.witherror)
        self.unf_result_cov()
        self.unf_result_df()
        
    def unf_result_df(self):
        self.result_df['truth_central'], self.result_df['truth_stat_error'] = th1_to_arr(self.hist_test_true)
        self.result_df['measured_central'], self.result_df['measured_error'] = th1_to_arr(self.hist_test_measure)
        self.result_df['unfolded_central'], self.result_df['unfolded_error'] = th1_to_arr(self.unfres)
        
    def unf_result_cov(self):
        """  Convert TMatrixT to ndarray
        """
        for x in range(self.nbins):
            for y in range(self.nbins):
                self.result_cov[x][y] = self.unfres_cov[x][y]
                
    
    def set_hist_measure_cov(self, reco_cov):
        self.reco_cov = reco_cov
        
    def set_hist_measure(self, central_array, error_array='False'):
        self.hist_test_measure = arr_to_th1(self.bins, central_array, error_array)
            

            
    def plot_migration(self, txt_offset=0, txt_fontsize=12, fname_mig=None, fname_n=None):
        mig_matrix, n_matrix, mig_fig, n_fig = get_migration(
                   th1_to_arr(self.hist_train_true), th1_to_arr(self.hist_train_measure),
                   np.ones(len(th1_to_arr(self.hist_train_true))), self.bins, self.show_var, txt_offset=txt_offset, txt_fontsize=txt_fontsize, fname_mig=fname_mig, fname_n=fname_n)
                   
                   
    
    def check_bias(self):
        """ Define and calculate figures of merit
        Returns:
        a : sum of absolute biases
        b : sum of normalised absolute biases
        c : sum of biases
        d : total stat error after unfolding taking into account bin-to-bin correlations
        e : ratio of sumed absoluted biases and total stat error after unfolding
        f : sum of two terms in type e
        g : sum of bin-wise ratio between biases and unfolding error
        
        """
        self.bias_a, self.bias_b, self.bias_c, self.bias_d, self.bias_e, self.bias_f, self.bias_g = study_complex_errors(th1_to_arr(self.unfres)[0], th1_to_arr(self.hist_test_true)[0], self.result_cov)
    
    
       
