#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" toy_unfold
    
This provides main unfolding class for toy study.

"""

import ROOT
import os
ROOT.gSystem.Load(os.getenv("ROOUNFOLD_PATH"))
from pyroounfold.utils.roo_convertor import *
from pyroounfold.utils.bias_study import *
from pyroounfold.utils.generate_toys import *
from pyroounfold.utils.unfold_methods import do_unfold as uf
import pyroounfold.plotstyle as ps
from pyroounfold.plotstyle.colors import PaperColors as p
from pyroounfold.plotting import get_migration
from pyroounfold.unfold import unfold

import pandas as pd
import numpy as np
    
        
class toy_unfold:

    def __init__(self, df_train, weight_train, df_test, weight_test, name_var_true, name_var_reco, show_var, bins, reco_bin_error='False', reco_cov='False', toy_size=1000, poisson=True, kcovtoy=False, mc_stat_err=0):
    
        """
        
        Args:
        df_train : train dataframe used for migration matrix, e.g. MC
        weight_train : array of event-weights for train sample
        df_test : test dataframe used for unfolding target, e.g. data
        weight_test : array of event-weights for test sample
        name_var_true : string of index in the dataframe for a true varibale, e.g. 'true_q2'
        name_var_reco : string of index in the dataframe for a reconstructed varibale, e.g. 'reco_q2'
        show_var : string of variable to be shown in plot, e.g. r'$q^{2}$'
        bins: an array of binning
        reco_bin_error (optional) : measured bin-wiese uncertainty, and used to build toys (Gaussian smearing) if reco_cov is not provided.
        reco_cov (optional) : measured covariance matrix, and used to build toys based on multivariate Gaussian smearing. If none of reco_bin_error or reco_cov are provided, toys are produced based on stat. error with Gaussian ('poisson=False') or Poisson ('poisson=True') smearing.

        toy_size (optional) : number of toys, default is 1000
        poisson (optional) : flag to use Poisson smearing for statistical uncertainty
        kcovtoy (optional) : flag provided by ROOUNFOLD. Default is False and the full covariance matrix 'reco_cov' propagated through unfolding. If True, the error propagation is based on toys generated internally by RooUnfold.
        mc_stat_err (optional) : relate to ROOUNFOLD::includeSystematics().  Default "0" is to leave out the effect of statistical uncertainties on the migration matrix. "1" is to include the effect. "2" is for only counting the statistical uncertainties of measured distribtuon and migration matrix. The effect is valueated by internal toys.
        
        
        """
        
        if(kcovtoy):
            self.witherror = ROOT.RooUnfold.kCovToy  #  error propagation based on toys generated internally by RooUnfold
        
        else:
            self.witherror = ROOT.RooUnfold.kCovariance  #  error propagation based on full covariance matrix
        
        self.mc_stat_err = mc_stat_err
        self.hist_train_true, self.hist_train_measure, self.hist_respon, self.hist_test_true, self.hist_test_measure = df_to_roounf(
        train = df_train,
        weight_train = weight_train,
        test = df_test,
        weight_test = weight_test,
        true_name = name_var_true,
        reco_name = name_var_reco,
        bins = bins
        )
        
        self.size = toy_size
        self.bins = bins
        
        self.show_var = show_var

        self.nbins = self.hist_test_measure.GetNbinsX()
        
        self.reco_bin_error = reco_bin_error
        
        self.reco_cov = reco_cov
        
        self.reco_cor = np.diagonal(np.ones(self.nbins))
        
        self.kcovtoy = kcovtoy
        
        if reco_cov != 'False':
            self.reco_bin_error = np.sqrt(reco_cov.diagonal())
            self.reco_cor = cov2corr(reco_cov)
            
        if (reco_cov == 'False')&(reco_bin_error == 'False'):
            self.reco_bin_error = np.sqrt(self.hist_test_measure)
            
        if len(self.reco_bin_error)==self.nbins:
            for x in range(self.nbins):
                self.hist_test_measure.SetBinError(x+1, self.reco_bin_error[x])
        
        
        self.result_cov =  np.zeros((self.nbins, self.nbins))
        
        self.toys_df = get_toys_fromDataFrame_toBinnedData(df_test, name_var_reco, bins, weight_test, toy_size, self.reco_cov, self.reco_bin_error, poisson)
        


    def do_toyUnfold(self,  method=None, para=None, get_fom=False):
        if method is None : print('Please indicate one method for unfolding: \'Ids\', \'Svd\', \'Bayes\', \'TUnfold\', \'Invert\', \'BinByBin\'.')
        
        result_df = pd.DataFrame()
        for i in range(self.size):
            hist_test_measure_toy = arr_to_th1(self.bins, self.toys_df.iloc[i], self.reco_bin_error )
            df_unf, _ = uf(self.hist_test_true, hist_test_measure_toy, self.hist_respon, method, para, self.reco_cov, self.kcovtoy, self.mc_stat_err)
            result_df = result_df.append(df_unf)
            
        result_cen_mean = [result_df.loc[result_df.bin_index==i,'unfolded_central'].median() for i in range(0, len(self.bins)-1)]
        
        result_cen_err = [result_df.loc[result_df.bin_index==i,'unfolded_error'].std() + result_df.loc[result_df.bin_index==i,'unfolded_error'].median()
                   for i in range(0, len(self.bins)-1)]
        result_cov = np.outer(result_cen_err, result_cen_err) * self.reco_cor
        
        
        self.result_df = result_df
        self.result_cen_mean = np.asarray(result_cen_mean)
        self.result_cen_err = np.asarray(result_cen_err)
        self.result_cov = result_cov
        if get_fom==True:
            err_a, err_b, err_c, err_d, err_e, err_f, err_g = study_complex_errors(result_cen_mean, th1_to_arr(self.hist_test_true)[0], result_cov)
            dict_fom ={'fom_a': err_a, 'fom_b': err_b, 'fom_c': err_c, 'fom_d': err_d, 'fom_e': err_e, 'fom_f': err_f, 'fom_g': err_g }
            self.dict_fom = dict_fom
        
    
        
    def do_toyUnfold_scan(self, method=None, para_arr=None, get_fom=True):
        if method is None :
            print('Please indicate one method for unfolding: \'Ids\', \'Svd\', \'Bayes\'.')
        elif method in ['TUnfold', 'Invert', 'BinByBin']:
            print(method + " has no parameter for scanning.")
        else:
            print("Loop in given parameters......")
            
            unf_cen_all = np.array([])
            unf_err_all = np.array([])
            err_a=np.array([])
            err_b=np.array([])
            err_c=np.array([])
            err_d=np.array([])
            err_e=np.array([])
            err_f=np.array([])
            err_g=np.array([])
            
            for x in para_arr:
                print("para = " + str(x))
                self.do_toyUnfold(method, x, True)
                unf_cen_all = np.append(unf_cen_all, self.result_cen_mean)
                unf_err_all = np.append(unf_err_all, self.result_cen_err)
                err_a = np.append(err_a, self.dict_fom['fom_a'])
                err_b = np.append(err_b, self.dict_fom['fom_b'])
                err_c = np.append(err_c, self.dict_fom['fom_c'])
                err_d = np.append(err_d, self.dict_fom['fom_d'])
                err_e = np.append(err_e, self.dict_fom['fom_e'])
                err_f = np.append(err_f, self.dict_fom['fom_f'])
                err_g = np.append(err_g, self.dict_fom['fom_g'])
        
            unf_cen_all = unf_cen_all.reshape(len(para_arr), self.nbins)
            unf_err_all = unf_err_all.reshape(len(para_arr), self.nbins)
        if get_fom==True:
            dict_fom_all ={'fom_a': err_a, 'fom_b': err_b, 'fom_c': err_c, 'fom_d': err_d, 'fom_e': err_e, 'fom_f': err_f, 'fom_g': err_g }
            return unf_cen_all, unf_err_all, dict_fom_all
            
        return unf_cen_all, unf_err_all
        
        
        
        
        
    
        

