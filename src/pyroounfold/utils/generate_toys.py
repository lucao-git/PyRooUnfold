#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Plotting
    
This provides some tools to easily plot data and other outputs from unfolding.

"""

import numpy as np
import pandas as pd


def get_toys_BinnedData(n_arr, size=1000, bin_err='False', cov='False',  poisson=False):

    if cov!='False':
        print('Generate toys abased on multivariate Gaussian of input covariance.')
        toys = np.random.multivariate_normal(n_arr, cov, size)
    else:
        toys=np.array([])
        for i in range(0,size):
            if poisson==False:
                if bin_err!='False':
                    print('Generate toys abased on Gaussian of input bin-error.')
                    toy = np.random.normal(n_arr, bin_err)
                if bin_err=='False': # stat. error
                    print('No input error found. Generate toys abased on Gaussian of statistical bin-error.')
                    toy = np.random.normal(n_arr, np.sqrt(n_arr))
            else:
                    print('Generate toys abased on Poisson of statistical error.')
                    toy = np.random.poisson(n_arr)
            toys=np.append(toys, toy, axis=0)
    toys=toys.reshape(size, len(n_arr))
    df = pd.DataFrame(data=toys, columns=['bin_index_'+str(i) for i in range(0,len(n_arr))])

    return df
    
def get_toys_fromDataFrame_toBinnedData(df, variable, bins_var, mean_w, size=1000, cov='False', bin_err='False', poisson=False):
    n,_ = np.histogram(df[variable], bins=bins_var, density=False, weights=mean_w)
    
    return get_toys_BinnedData(n, size, bin_err, cov, poisson)



