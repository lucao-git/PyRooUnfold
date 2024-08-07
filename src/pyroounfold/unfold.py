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
#ROOT.gSystem.Load(os.getenv("ROOUNFOLD_PATH"))
ROOT.gSystem.Load("libRooUnfold.so")
from pyroounfold.utils.roo_convertor import df_to_roounf, th1_to_arr, arr_to_th1_withErr, ndarr_to_tmatrix, tvec_to_arr, arr_to_tvec
from pyroounfold.utils.bias_study import study_complex_errors, cov2corr
from pyroounfold.utils.generate_toys import *
import pyroounfold.plotstyle as ps
from pyroounfold.plotstyle.colors import PaperColors as p


import pandas as pd
import numpy as np

class unfold:

    def __init__(self, df_train, weight_train, df_test, weight_test, name_var_true, name_var_reco, show_var, bins, reco_bin_error='False', reco_cov='False', kcovtoy=False, mc_stat_err=3):
        """
        
        Args:
        df_train : train dataframe used for migration matrix, e.g. MC
        weight_train : array of event-weights for train sample
        df_test : test dataframe used for unfolding target, e.g. data
        weight_test : array of event-weights for test sample
        name_var_true : string of index in the dataframe for a true varibale, e.g. 'true_q2'
        name_var_reco : string of index in the dataframe for a reconstructed varibale, e.g. 'reco_q2'
        show_var : string of variable to be shown in plot, e.g. '$q^{2}$'
        bins: an array of binning
        reco_bin_error (optional) : measured bin-wiese uncertainty, default is statistical error based on bin content
        reco_cov (optional) : measured covariance matrix, default is statistical covariance
        kcovtoy (optional) : flag provided by ROOUNFOLD. Default is False and the full covariance matrix 'reco_cov' propagated through unfolding. If True, the error propagation is based on toys generated internally by RooUnfold.
        mc_stat_err (optional) : relate to ROOUNFOLD::IncludeSystematics(). "3"[kAll, default] is to include the effect of MC statistical uncertainties on the migration matrix which is is evaluated with internal toys. "0"[kNoSystematics] is to exclude this effect and then only the propagated measurement error will be obtained. "2" [kAlphas] is for only counting the MC statistical uncertainties of measured distribtuon and migration matrix.
        
        
        Result:
        result_df : dataframe with columns=['bin_index', 'truth_central', 'truth_stat_error', 'measured_central', 'measured_error', 'unfolded_central', 'unfolded_error', 'coverage_perbin']
        result_cov : post-unfolding covariance matrix
        chi2: provided by RooUnfold Chi2 (const Hist* hTrue,RooUnfolding::ErrorTreatment DoChi2=RooUnfolding::kCovariance)
        """

        
        if(kcovtoy):
            self.witherror = ROOT.RooUnfold.kCovToys  #  error propagation based on toys generated internally by RooUnfold
        else:
            self.witherror = ROOT.RooUnfold.kCovariance  # [default] error propagation based on full covariance matrix
        
        
        self.mc_stat_err = mc_stat_err
        
        self.hist_train_true, self.hist_train_measure, self.hist_respon, self.hist_train_Adet, self.hist_test_true,  self.hist_test_measure = df_to_roounf(
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
        self.kcovtoy = kcovtoy
        
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
        
        self.unf_chi2 = -1
        self.unf_coverage_perbin = -1  #coverage probability of the unfolded result within one sigma for each bin
        
        
        if len(self.reco_bin_error)==self.nbins:
            for x in range(self.nbins):
                self.hist_test_measure.SetBinError(x+1, self.reco_bin_error[x])
        
    def do_Ids(self, para):
        if para is None or para <0 : print('ERROR: Ids method requires a iteration number (>=0).')
        elif para>=0 :
            unf = ROOT.RooUnfoldIds(self.hist_respon, self.hist_test_measure, para)
            unf.IncludeSystematics(self.mc_stat_err)
            if self.reco_cov != 'False':
                unf.SetMeasuredCov(ndarr_to_tmatrix(self.reco_cov))
            self.unfres = unf.Hunfold(self.witherror)
            self.unfres_cov = unf.Eunfold(self.witherror)
            self.unf_result_cov()
            self.unf_chi2 = unf.Chi2(self.hist_test_true, self.witherror)
            self.unf_coverage_perbin = tvec_to_arr(unf.CoverageProbV(1))
            self.unf_result_df()
                    
    
    def do_Svd(self, para):
        if para is None or para <0 : print('ERROR: Svd method require a regularisation number (<= nbins, 0 using default nbins/2).')
        elif para > self.hist_test_measure.GetNbinsX(): print('ERROR: Svd method do not work when regularisation number > nbins.')
        elif para>= 0 & para <= self.hist_test_measure.GetNbinsX():
            unf = ROOT.RooUnfoldSvd(self.hist_respon, self.hist_test_measure, para)
            unf.IncludeSystematics(self.mc_stat_err)
            if self.reco_cov != 'False':
                unf.SetMeasuredCov(ndarr_to_tmatrix(self.reco_cov))
            self.unfres = unf.Hunfold(self.witherror)
            self.unfres_cov = unf.Eunfold(self.witherror)
            self.unf_result_cov()
            self.unf_chi2 = unf.Chi2(self.hist_test_true, self.witherror)
            self.unf_coverage_perbin = tvec_to_arr(unf.CoverageProbV(1))
            self.unf_result_df()
            
    def do_Bayes(self, para):
        if para is None :
            print('ERROR: Bayes method requires a iteration number.')
        else:
            unf = ROOT.RooUnfoldBayes(self.hist_respon, self.hist_test_measure, para)
            unf.IncludeSystematics(self.mc_stat_err)
            if self.reco_cov != 'False':
                unf.SetMeasuredCov(ndarr_to_tmatrix(self.reco_cov))
            self.unfres = unf.Hunfold(self.witherror)
            self.unfres_cov = unf.Eunfold(self.witherror)
            self.unf_result_cov()
            self.unf_chi2 = unf.Chi2(self.hist_test_true, self.witherror)
            self.unf_coverage_perbin = tvec_to_arr(unf.CoverageProbV(1))
            self.unf_result_df()
            
            
    def do_Invert(self):
        unf = ROOT.RooUnfoldInvert(self.hist_respon, self.hist_test_measure)
        unf.IncludeSystematics(self.mc_stat_err)
        if self.reco_cov != 'False':
            unf.SetMeasuredCov(ndarr_to_tmatrix(self.reco_cov))
        self.unfres = unf.Hunfold(self.witherror)
        self.unfres_cov = unf.Eunfold(self.witherror)
        self.unf_result_cov()
        self.unf_chi2 = unf.Chi2(self.hist_test_true, self.witherror)
        self.unf_coverage_perbin = tvec_to_arr(unf.CoverageProbV(1))
        self.unf_result_df()
        
        
    def do_TUnfold(self, para=None):
        if para is None :
            print('Info: Regularisation parameter tau is not indicated and then will be optimised internally.')
            unf = ROOT.RooUnfoldTUnfold(self.hist_respon, self.hist_test_measure)
        else:
            unf = ROOT.RooUnfoldTUnfold(self.hist_respon, self.hist_test_measure, para)
        unf.IncludeSystematics(self.mc_stat_err)
        if self.reco_cov != 'False':
            unf.SetMeasuredCov(ndarr_to_tmatrix(self.reco_cov))
        self.unfres = unf.Hunfold(self.witherror)
        self.unfres_cov = unf.Eunfold(self.witherror)
        self.unf_result_cov()
        self.unf_chi2 = unf.Chi2(self.hist_test_true, self.witherror)
        self.unf_coverage_perbin = tvec_to_arr(unf.CoverageProbV(1))
        self.unf_result_df()
        
    
    def do_BinByBin(self):
        unf = ROOT.RooUnfoldBinByBin(self.hist_respon, self.hist_test_measure)
        unf.IncludeSystematics(self.mc_stat_err)
        if self.reco_cov != 'False':
            unf.SetMeasuredCov(ndarr_to_tmatrix(self.reco_cov))
        self.unfres = unf.Hunfold(self.witherror)
        self.unfres_cov = unf.Eunfold(self.witherror)
        self.unf_result_cov()
        self.unf_chi2 = unf.Chi2(self.hist_test_true, self.witherror)
        self.unf_coverage_perbin = tvec_to_arr(unf.CoverageProbV(1))
        self.unf_result_df()
        
        
    def do_GP(self, para=None):
        unf = ROOT.RooUnfoldGP(self.hist_respon, self.hist_test_measure, kernel = 1) ## only 1 is available in RooUnfold for now
        if para is not None :
            unf.SetRegParm(para)
        unf.IncludeSystematics(self.mc_stat_err)
        if self.reco_cov != 'False':
            unf.SetMeasuredCov(ndarr_to_tmatrix(self.reco_cov))
        self.unfres = unf.Hunfold(self.witherror)
        self.unfres_cov = unf.Eunfold(self.witherror)
        self.unf_result_cov()
        self.unf_chi2 = unf.Chi2(self.hist_test_true, self.witherror)
        self.unf_coverage_perbin = tvec_to_arr(unf.CoverageProbV(1))
        self.unf_result_df()
        
        
        
    def unf_result_df(self):
        self.result_df['truth_central'], self.result_df['truth_stat_error'] = th1_to_arr(self.hist_test_true)
        self.result_df['measured_central'], self.result_df['measured_error'] = th1_to_arr(self.hist_test_measure)
        self.result_df['unfolded_central'], self.result_df['unfolded_error'] = th1_to_arr(self.unfres)
        self.result_df['coverage_perbin'] = self.unf_coverage_perbin
        
    def unf_result_cov(self):
        """  Convert TMatrixT to ndarray
        """
        for x in range(self.nbins):
            for y in range(self.nbins):
                self.result_cov[x][y] = self.unfres_cov[x][y]
                
    
    def set_hist_measure_cov(self, reco_cov):
        self.reco_cov = reco_cov
        
    def set_hist_measure(self, central_array, error_array='False'):
        self.hist_test_measure = arr_to_th1_withErr(self.bins, central_array, error_array)
        
                   
    
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
        h : averaged coverage probability
        
        """
        self.bias_a, self.bias_b, self.bias_c, self.bias_d, self.bias_e, self.bias_f, self.bias_g = study_complex_errors(th1_to_arr(self.unfres)[0], th1_to_arr(self.hist_test_true)[0], self.result_cov)
        self.bias_h = sum(self.unf_coverage_perbin)/len(self.unf_coverage_perbin)
        
        
            

    # Unfolding methods.
    def algorithm(method):
      alg = None
      if method == "bayes":
        alg= ROOT.RooUnfolding.kBayes;
      elif method == "bbb":
        alg= ROOT.RooUnfolding.kBinByBin;
      elif method == "inv":
        alg= ROOT.RooUnfolding.kInvert;
      elif method == "svd":
        alg= ROOT.RooUnfolding.kSVD;
      elif method == "root":
        alg= ROOT.RooUnfolding.kTUnfold;
      elif method == "ids":
        alg= ROOT.RooUnfolding.kIDS;
      elif method == "gp":
        alg= ROOT.RooUnfolding.kGP;
      else:
        print("The passed unfolding method does not match any of the supported methods. Please pass one of the following methods:")
        print("bayes")
        print("bbb")
        print("inv")
        print("svd")
        print("root")
        print("ids")
        print("gp")
        exit(0)
      return alg

  
    def try_RooFitUnfold(self, method, para=None):
        """
        A simple example of using RooFitUnfold based on ROOT.RooUnfoldSpec()
  
        """
        spec = ROOT.RooUnfoldSpec("unfold","unfold", self.hist_train_true,"train_truth", self.hist_train_measure, "train_reco",
                                  self.hist_train_Adet, 0, self.hist_test_measure, False, 0.0005, False)

        if not para == None:
            unfolding = spec.makeFunc(algorithm(method), para)
        else:
            unfolding = spec.makeFunc(algorithm(method))
            
        ROOT.gDirectory.Clear()
        
        test_truth = spec.makeHistogram(self.hist_test_true)
        test_truth_func = test_truth.func()
        test_truth_func.SetName("test_truth")
        ws = ROOT.RooWorkspace("workspace","workspace")
        getattr(ws,"import")(unfolding)
        getattr(ws,"import")(test_truth_func)
        src = str(unfolding.getStringAttribute("source"))
        func = ws.function(src)
        func.unfolding()
        
        # plot
        unfoldfunc = ws.obj("unfold")
        test_truth_ws = ws.obj("test_truth")
        train_truth = ws.obj("train_truth")
        plot_truth = train_truth.frame()
        test_truth_ws.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name("test_truth_graph"))
        unfoldfunc.plotOn(plot_truth,ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.Name("unfolded_data_graph"))
        canvas_truth = ROOT.TCanvas("unfolded","unfolded")
        plot_truth.Draw()
        leg_truth = ROOT.TLegend(0.1, 0.8, 0.3, 0.9)
        leg_truth.AddEntry( plot_truth.findObject("test_truth_graph"), "Truth", "l" )
        leg_truth.AddEntry( plot_truth.findObject("unfolded_data_graph"), "Unfolded", "l" )
        leg_truth.Draw()
        canvas_truth.SaveAs("unfolded.pdf")




    
    
       
