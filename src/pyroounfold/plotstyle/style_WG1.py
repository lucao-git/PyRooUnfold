# -*- coding: utf-8 -*-
"""
style_WG1.py
======

Context decorator for producing figures for bonn thesis.

EPS format for image file are used, because it is high quality and working
in actual physical size rahter than pixel unit.

cadidate colormaps:
    # cmap = cm.viridis
    # cmap = cm.RdYlGn
    # cmap = cm.RdYlBu
    # cmap = cm.binary
    # cmap = cm.Greys
    # cmap = cm.coolwarm
    # cmap = cm.hot
    # cmap = cm.gnuplot
    # cmap = cm.gnuplot2
    # cmap = cm.inferno
    # cmap = cm.magma

rcparams:
https://matplotlib.org/tutorials/introductory/customizing.html

"""

from .decorator import MPLdecorator
from .colors import *
from .layout import GOLDEN_RATIO
import matplotlib.font_manager

__all__ = ['WG1_decorator', ]

# Constants from APS Authour Guidelines. https://journals.aps.org/prl/authors
# Figures should have a width of a 8.6 cm or 3 3/8 in, the width of a single manuscript column...
width_single_column = 3.375
width_double_column = 6.75

# Default ratio for a single plot figure
# I prefer a little higher than goden ratio, from 0.618 to about 0.68
height_width_ratio = GOLDEN_RATIO * 1.1  # = height / width

_width = width_single_column
_height = width_single_column * height_width_ratio

_params = {'font.family': 'DejaVu Sans',
           'font.serif': ['Times New Roman'],
           # ['Helvetica', 'Computer Modern Sans serif'] are not default installed,
           #'font.sans-serif': ['Helvetica'],
           #'font.size': 8,
           #'text.usetex': True,
           # To force LaTeX use Helvetica fonts.
#           'text.latex.preamble': [r'\usepackage{siunitx}',
#                                   r'\sisetup{detect-all}',
#                                   r'\usepackage{helvet}',
#                                   r'\usepackage[eulergreek,EULERGREEK]{sansmath}',
#                                   r'\sansmath'],
           #    'axes.prop_cycle': monochrome,
           #    'axes.prop_cycle': bonn_cycler_1,
           #'axes.prop_cycle': bonn_cycler_2,
           #    'axes.prop_cycle': default_color_cycler,

           #'image.cmap': 'RdYlBu',
           'axes.labelsize': 22,
           'axes.prop_cycle': paper_color_cycler,
           'axes.formatter.limits': (-4, 4),
           'axes.formatter.use_mathtext': True,
           'axes.titlesize': 'large',
           'axes.labelpad': 6.0, #6.0,
           'axes.linewidth': 1.2,
           
           'errorbar.capsize': 2,
           

           'figure.figsize': (_width, _height),
           # 'figure.subplot.left' : 0.125,
           # 'figure.subplot.right' : 0.95,
           # 'figure.subplot.bottom' : 0.1,
           # 'figure.subplot.top' : 0.95,

           'savefig.dpi': 600,
           'savefig.format': 'pdf',  # 'eps'#'png'
           # 'savefig.bbox': 'tight',
           # this will crop white spaces around images that will make
           # width/height no longer the same as the specified one.

           
           'legend.frameon': False,
           'legend.fontsize': 15,
#           'legend.numpoints': 1,
#           'legend.handlelength': 2,
#           'legend.scatterpoints': 1,
#           'legend.labelspacing': 0.5,
#           'legend.markerscale': 0.9,
#           'legend.handletextpad': 0.5,  # pad between handle and text
#           'legend.borderaxespad': 0.5,  # pad between legend and axes
#           'legend.borderpad': 0.5,  # pad between legend and legend content
#           'legend.columnspacing': 1,  # pad between each legend column

           # 'text.fontsize' : 8,
           'xtick.labelsize': 15,
           'ytick.labelsize': 15,
           'xtick.minor.visible':True,
           'ytick.minor.visible':True,
           'xtick.top':True,
           'ytick.right':True,
           'xtick.major.size': 6,
           'ytick.major.size': 6,
           'xtick.minor.size': 3,
           'ytick.minor.size': 3,
           'xtick.direction':'in',
           'ytick.direction':'in',
           'xtick.major.pad': 8,
           'ytick.major.pad': 8,
           'lines.linewidth': 1.5,
           # 'lines.markeredgewidth' : 0,
           # 0 will make line-type markers, such as '+', 'x', invisible
           }

WG1_decorator = MPLdecorator(_params)
