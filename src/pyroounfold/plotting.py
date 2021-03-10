#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Plotting
    
This provides some tools to easily plot data and other outputs from unfolding.

"""

import numpy as np
import matplotlib.pyplot as plt
from pyroounfold.utils.roo_convertor import th1_to_arr
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyroounfold.plotstyle as ps
from pyroounfold.plotstyle.colors import PaperColors as p


def get_bin_centers(bins_array):
    """ Extract bin center from a defined bins array
        
        Args:
        bins_array (array): defined bins (length = number of bins + 1)
        
        Returns:
        centers : an array of center postion of each bin (length = number of bins)
        """
    
    centers = [a + (b-a)/2 for a, b in zip(bins_array[0:-1], bins_array[1:])]
    return np.asarray(centers)



def get_bin_widths(bins_array):
    """ Extract bin width from a defined bins array
        
        Args:
        bins_array (array): defined bins (length = number of bins + 1)
        
        Returns:
        widths : an array of width for each bin (length = number of bins)
        """
    
    widths = [b-a for a, b in zip(bins_array[0:-1], bins_array[1:])]
    return np.asarray(widths)


@ps.WG1_decorator
def plot_compare_single_run(df, bins, x_title, y_title):
    """ Plot for comparing unfolded results, truth and measured distributions
        
        Args:
        df : dataframe with 'truth_central', 'measured_central', 'unfolded_central' and 'unfolded_error'
        bins : defined bins array (length = number of bins + 1)
        x_title : string for x-axis title
        y_title : string for y-axis title
        
        Returns:
        fig : one plot for all compared distributions
        """
    fig = plt.figure(figsize=(9.0, 5.5))
    plt.plot(get_bin_centers(bins), df['truth_central'], marker='o', color='greenyellow', ls='', label='true')
    fig = plt.plot(get_bin_centers(bins), df['measured_central'], marker='o', color='red', ls='', label='reco')
    #plt.bar(get_bin_centers(bins), height=reco_err, width=get_bin_widths(bins),bottom=reco_cen - reco_err/2, fill=False, label='stat. error')
    
    
    plt.errorbar(x=get_bin_centers(bins), y=df['unfolded_central'], yerr=df['unfolded_error'], xerr=0, marker='.', color='black',ls='', elinewidth=0.5, ecolor='black', label='unfolded')
    
    plt.xlabel(x_title)
    plt.ylabel(y_title)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.show()
    plt.close()
    return fig

@ps.WG1_decorator
def plot_compare_para(nerr, para, name, err_list, str_list, c_list, label, last_para_index=None):
    """ Plot for comparing errors with various parameters in SVD or IDS method
        
        Args:
        nerr : number of error types to be compared, used to define subplots
        para : array of tested parameters
        name : name of parameter, e.g. 'Regularisation k' for SVD, 'Iteration n' for IDS
        err_list : list of errors
        str_list : list of error lables or expressions in the same order of errors stored in err_list
        c_list : color list for each type of error
        label : string indicating plotted variable
        last_para_index : (int) index of last parameter for comparasion, default (=None) is comparing all in the 'para' list
        
        Returns:
        fig : a group of plots comparing input errors
        
        """
    if last_para_index is not None:
        para = para[:last_para_index]
        err_list = [err_list[x][:last_para_index] for x in range(nerr)]
        str_list = str_list[:last_para_index]
        c_list = c_list[:last_para_index]
    
    fig, axs = plt.subplots(1, nerr, figsize=(nerr*4,4))
    for i in range(nerr):
        axs[i].scatter(x= para, y = err_list[i], marker = 'o', c = c_list[i], label=label)
        #axs[i].legend(loc='center')
        axs[i].set_xlabel(name)
        axs[i].set_ylabel(str_list[i])
        axs[i].xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    return fig


@ps.WG1_decorator
def get_migration(true_data, reco_data, weight, bin_var, name_var, txt_offset=0, txt_fontsize=12, fname_mig=None, fname_n=None):
    """ Get migration matrix (also called response matrix) and entries-number matrix
        
        Args:
        true_data : array of variable in MC truth
        reco_data : array of reconstructed variable
        weight : event weight, will apply on ture_data and reco_data
        bin_var : binning array of variable
        name_var : string for variable definition, e.g."$p$"
        txt_offset : modify the text position on 2D plot
        
        Returns:
        mig_matrix (2D-array): migration matrix
        n_matrix (2D-array): number of entries matrix
        mig_fig : plot of migration matrix in percentage
        n_fig : plot of number of entries matrix
        
        To Do:
        ---- could use get_bin_centers() to automatically locate text instead of manully seting by txt_offset
        """
    rc('text', usetex=True)
    plt.rc('text', usetex=True)
    bin_layout = range(len(bin_var)-1)
    n_matrix, _, _ = np.histogram2d(x=reco_data, y=true_data, bins=bin_var, weights=weight)
    x = np.linalg.norm(n_matrix, ord=1, axis=0)
    mig_matrix = n_matrix/np.expand_dims(x, axis=0)
    mig_fig, ax = plt.subplots(figsize=(8,8))
    ax.set_aspect("equal")
    im = ax.imshow(mig_matrix*100, origin='low',cmap=plt.cm.GnBu)
    ax.set_xlabel("Bin index of true " + name_var)
    ax.set_ylabel("Bin index of reconstracted " + name_var)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    ax.minorticks_off()
    plt.colorbar(im, cax=cax)
    
    for i in range(len(bin_var)-1):
        for j in range(len(bin_var)-1):
            ax.text(bin_layout[j]+txt_offset, bin_layout[i]+txt_offset, round(mig_matrix[i,j]*100, 1), fontsize = txt_fontsize,
                    color="black", ha="center", va="center", fontweight="bold")
    if fname_mig:
            plt.savefig(fname_mig)
    plt.show()
    plt.close()

    n_fig, ax2 = plt.subplots(figsize=(8,8))
    ax2.set_aspect("equal")
    im2 = ax2.imshow(n_matrix, origin='low',cmap=plt.cm.RdPu)
    ax2.set_xlabel("Bin index of true " + name_var )
    ax2.set_ylabel("Bin index of reconstructed " + name_var )
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="5%", pad=0.1)
    ax2.minorticks_off()
    plt.colorbar(im2, cax=cax2)
    for i in range(len(bin_var)-1):
        for j in range(len(bin_var)-1):
            ax2.text(bin_layout[j]+txt_offset, bin_layout[i]+txt_offset, int(n_matrix[i,j]),fontsize = txt_fontsize,
                 color="gray", ha="center", va="center", fontweight="bold")
    if fname_n:
            plt.savefig(fname_n)
    plt.show()
    plt.close()

    return mig_matrix, n_matrix, mig_fig, n_fig


# compare unfolded reslut obtained in SVD and inv
@ps.WG1_decorator
def plot_unf_scan(hist_test_reco, hist_test_true,  respon_matrix, unf_cen_svd,
                     unf_err_svd, unf_cen_inv,unf_err_inv,unfold_bin, xlab, leg, case, yup=None, fname=None):
    
    mea_cen, mea_err = th1_to_arr(hist_test_reco)
    mc_cen, mc_err = th1_to_arr(hist_test_true)


    fig = plt.figure(figsize=(9.0, 5.5))
    
    plt.hist(get_bin_centers(unfold_bin), bins=unfold_bin, weights=mea_cen, color=p.p_blue,
              density=False, histtype='step', lw=1.5, label=r'Reco MC '+case+' with $\sigma^{tot}$')
    plt.bar(get_bin_centers(unfold_bin),
           height= 2*mea_err,
           width=get_bin_widths(unfold_bin),
           bottom=mea_cen - mea_err, #label='stat. uncertainty',
           alpha=0.5, color=p.p_light_blue
           )
    
    plt.bar(get_bin_centers(unfold_bin),
           height= 2*unf_err_inv,
           width= get_bin_widths(unfold_bin),
           bottom= unf_cen_inv - unf_err_inv,
           edgecolor=p.p_orange, fill=False, lw=1.5, label='Matrix-inversion unfolded'
           )
    
    plt.scatter(get_bin_centers(unfold_bin),unf_cen_inv, color=p.p_orange, marker='s')

    plt.hist(get_bin_centers(unfold_bin), bins=unfold_bin, weights=mc_cen, color=p.p_red, density=False,
             histtype='step', lw=1.5, label=r'True MC '+case)

    
    plt.errorbar(x=[], y=[], xerr=[],yerr=[], marker='.', markersize=0.1, color='black',ls='',
                 elinewidth=0.5, ecolor='black',label='SVD unfolded, k=(2..' + str(unfold_bin.size-1)+')')

    if xlab == r'$E^{B}_{\ell}$ [GeV]':
        ms = 1.2
    else:
        ms = 3
    for x in range (unfold_bin.size -2):  # loop all k
        plt.errorbar(x=unfold_bin[:-1]+get_bin_widths(unfold_bin)/(unfold_bin.size -1)*(x+1),
                     y=unf_cen_svd[x], yerr=unf_err_svd[x],
                     xerr=0, marker='.', markersize=ms, color='black', capsize=0.0, ls='', lw=0.8, elinewidth=0.5, ecolor='black')

    plt.xlabel(xlab)
    plt.ylabel('Events')
    plt.xlim(unfold_bin[0], unfold_bin[-1])

    plt.legend(loc=leg)
    plt.tight_layout()
    if yup:
        plt.ylim(0,yup)
    if fname:
            plt.savefig(fname)
            
            
            
@ps.WG1_decorator
def plot_unf_scan_perBinWidth(hist_test_reco, hist_test_true,  respon_matrix, unf_cen_svd,
                     unf_err_svd, unf_cen_inv,unf_err_inv,unfold_bin, xlab, leg, case, yup=None, fname=None):
    
    mea_cen, mea_err = th1_to_arr(hist_test_reco)
    mc_cen, mc_err = th1_to_arr(hist_test_true)


    fig = plt.figure(figsize=(9.0, 5.5))
    
    plt.hist(get_bin_centers(unfold_bin), bins=unfold_bin, weights=mea_cen/get_bin_widths(unfold_bin), color=p.p_blue,
              density=False, histtype='step', lw=1.5, label=r'Reco MC '+case+' with $\sigma^{exp}_{stat}$')
    plt.bar(get_bin_centers(unfold_bin),
           height= 2*mea_err/get_bin_widths(unfold_bin),
           width=get_bin_widths(unfold_bin),
           bottom=mea_cen/get_bin_widths(unfold_bin) - mea_err/get_bin_widths(unfold_bin), #label='stat. uncertainty',
           alpha=0.5, color=p.p_light_blue
           )
    
    plt.bar(get_bin_centers(unfold_bin),
           height= 2*unf_err_inv/get_bin_widths(unfold_bin),
           width=get_bin_widths(unfold_bin),
           bottom=unf_cen_inv/get_bin_widths(unfold_bin) - unf_err_inv/get_bin_widths(unfold_bin),
           edgecolor=p.p_orange, fill=False, lw=1.5, label='Matrix-inversion unfolded'
           )

    plt.hist(get_bin_centers(unfold_bin), bins=unfold_bin, weights=mc_cen/get_bin_widths(unfold_bin), color=p.p_red, density=False,
             histtype='step', lw=1.5, label=r'True MC '+ case)

    
    plt.errorbar(x=[], y=[], xerr=[],yerr=[], marker='.', markersize=0.1, color='black',ls='',
                 elinewidth=0.5, ecolor='black',label='SVD unfolded, k=(2..' + str(unfold_bin.size-1)+')')

    if xlab == r'$E^{B}_{\ell}$ [GeV]':
        ms = 1.2
    else:
        ms = 3
    for x in range (unfold_bin.size -2):  # loop all k
        plt.errorbar(x=unfold_bin[:-1]+get_bin_widths(unfold_bin)/(unfold_bin.size -1)*(x+1),
                     y=unf_cen_svd[x]/get_bin_widths(unfold_bin), yerr=unf_err_svd[x]/get_bin_widths(unfold_bin),
                     xerr=0, marker='.', markersize=ms, color='black', capsize=0.0, ls='', lw=0.8, elinewidth=0.5, ecolor='black')

    plt.xlabel(xlab)
    plt.ylabel('Events/ bin width')
    plt.xlim(unfold_bin[0], unfold_bin[-1])

    plt.legend(loc=leg)
    plt.tight_layout()
    if yup:
        plt.ylim(0,yup)
    if fname:
            plt.savefig(fname)
            

# compare fom obtained in SVD and inv
@ps.WG1_decorator
def plot_unf_fom(k_arr, svd_fom_dict, inv_fom_dict, variable, leg='center', fname=None):
    
    err_list = [svd_fom_dict['fom_'+ i] for i in ['a', 'b', 'c', 'd', 'e', 'f']]
    errors_inv = [inv_fom_dict['fom_'+ i] for i in ['a', 'b', 'c', 'd', 'e', 'f']]
    
    err_names = ['$\sum{|b_{i}|}$',
             '$\sum{|b_{i}| / N_{i}^{true}}$',
             '$\sum{b_{i}}$',
             '$\sqrt{\sum Cov_{i,j}}$',
             '$\sum{|b_{i}|}/\sqrt{\sum Cov_{i,j}}$',
             '$\sqrt{(\sum{|b_{i}|})^{2} + \sum Cov_{i,j}}$' ]

    
    fig, axs = plt.subplots(2,3, figsize=(30, 15))
    i=0
    
    for loc_x in [0,1]:
        for loc_y in [0,1,2]:
            axs[loc_x, loc_y].scatter(x= k_arr, y = err_list[i], marker = 'o', s=200, c = p.p_deep_blue, label='SVD')
            axs[loc_x, loc_y].hlines(errors_inv[i], k_arr[0]-0.5, k_arr[-1]+0.5, lw=2, color='#737373', label='Matrix-inv.')
            axs[loc_x, loc_y].set_xlabel('k')
            axs[loc_x, loc_y].set_ylabel(err_names[i])
            axs[loc_x, loc_y].xaxis.set_major_locator(MaxNLocator(integer=True))
            axs[loc_x, loc_y].set_xlim(k_arr[0]-0.5, k_arr[-1]+0.5)
            i=i+1
            
    axs[0, 0].plot([], [], ' ', label=variable)

    axs[0, 0].legend(loc=leg, fontsize=26)
    
    plt.tight_layout()
    
    if fname:
        plt.savefig(fname)
            
            
# compare unfolded result obtained in SVD and inv
@ps.WG1_decorator
def plot_unf_toys(toy, inv_cen,inv_err, xlab, leg, case, violin_widths=0.2, yup=None, fname=None):
    
    inv_err = np.asarray(inv_err)
    inv_cen = np.asarray(inv_cen)
    mea_cen, _ = th1_to_arr(toy.hist_test_measure)
    mea_err = toy.reco_bin_error
    mc_cen, _ = th1_to_arr(toy.hist_test_true)
    
    unf_cen = [toy.result_df.loc[toy.result_df.bin_index==i,'unfolded_central'] for i in range(0, len(toy.bins)-1)]
    
    unf_cen_mean = toy.result_cen_mean
    unf_err = toy.result_cen_err


    fig = plt.figure(figsize=(9.0, 5.5))
    
    violin_parts = plt.violinplot(unf_cen, get_bin_centers(toy.bins), points=60, widths=violin_widths, showmeans=False,
                      showextrema=False, showmedians=False, bw_method=0.5)
    
    for pc in violin_parts['bodies']:
        pc.set_facecolor('gray')
        pc.set_edgecolor('gray')
    
    
    plt.hist(get_bin_centers(toy.bins), bins=toy.bins, weights=mea_cen, color=p.p_blue,
              density=False, histtype='step', lw=1.5, label=r'Reco MC '+case+' with $\sigma^{tot}$')
   
    plt.hist(get_bin_centers(toy.bins), bins=toy.bins, weights=mc_cen, color=p.p_red, density=False,
             histtype='step', lw=1.5, label=r'True MC '+case)
 
    plt.bar(get_bin_centers(toy.bins),
           height= 2*mea_err,
           width=get_bin_widths(toy.bins),
           bottom=mea_cen - mea_err,
           alpha=0.5, color=p.p_light_blue
           )
    
    plt.bar(get_bin_centers(toy.bins),
           height= 2*inv_err,
           width= get_bin_widths(toy.bins),
           bottom= inv_cen - inv_err,
           edgecolor=p.p_orange, fill=False, lw=1.5, label='Matrix-inversion unfolded'
           )
    
    plt.scatter(get_bin_centers(toy.bins),inv_cen, color=p.p_orange, marker='s')

    
    plt.errorbar(x=get_bin_centers(toy.bins), y=unf_cen_mean, xerr=np.zeros(len(toy.bins)-1), yerr=unf_err,
                 marker='o', markersize=2, color='black',ls='',
                 elinewidth=0.6, ecolor='black',label=r'Mean and $\sigma_{RMS}$ of toys')
    
    if xlab == r'$E^{B}_{\ell}$ [GeV]':
        ms = 1.2
    else:
        ms = 3
    
    plt.xlabel(xlab)
    plt.ylabel('Events')
    plt.xlim(toy.bins[0], toy.bins[-1])

    plt.legend(loc=leg)
    plt.tight_layout()
    if yup:
        plt.ylim(0,yup)
    if fname:
            plt.savefig(fname)
            
            

@ps.WG1_decorator
def plot_unf_fom_my(k_arr, svd_fom_dict, inv_fom_dict, variable, leg='center', fname=None):
    
    err_list = [svd_fom_dict['fom_'+ i] for i in ['a',  'c', 'e', 'g']]
    errors_inv = [inv_fom_dict['fom_'+ i] for i in ['a',  'c', 'e', 'g']]
    
    err_names = ['$\sum{|b_{i}|}$',
             '$\sum{b_{i}}$',
             '$\sum{|b_{i}|}/\sqrt{\sum Cov_{ij}}$',
             '$\sum{|b_{i}|/\sqrt{Cov_{ii}}}$' ]

    
    fig, axs = plt.subplots(2,2, figsize=(12, 8))
    i=0
    
    for loc_x in [0,1]:
        for loc_y in [0,1]:
            axs[loc_x, loc_y].hlines(errors_inv[i], k_arr[0]-0.5, k_arr[-1]+0.5, lw=2, color=p.p_orange, label='Matrix-inv.')
            axs[loc_x, loc_y].scatter(x= k_arr, y = err_list[i], marker = 'o', s=50, c = p.p_deep_blue, label='SVD')
            axs[loc_x, loc_y].set_xlabel('k')
            axs[loc_x, loc_y].set_ylabel(err_names[i])
            axs[loc_x, loc_y].xaxis.set_major_locator(MaxNLocator(integer=True))
            axs[loc_x, loc_y].set_xlim(k_arr[0]-0.5, k_arr[-1]+0.5)
            i=i+1
            
    axs[0, 0].plot([], [], ' ', label=variable)

    axs[0, 0].legend(loc=leg)
    
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.4)
    
    if fname:
        plt.savefig(fname)
            
