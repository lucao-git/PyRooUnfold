# -*- coding: utf-8 -*-
"""
aps.py
======

Context decorator for producing figures which is ready to publish
in American Chemical Society.

EPS format for image file are used, because it is high quality and working
in actual physical size rahter than pixel unit.

"""

from .decorator import MPLdecorator
from .colors import default_color_cycler,custom_cycler
from .layout import GOLDEN_RATIO

__all__ = ['aps_decorator', ]

# Constants from APS Authour Guidelines. https://journals.aps.org/prl/authors
# Figures should have a width of a 8.6 cm or 3 3/8 in, the width of a single manuscript column...
width_single_column = 3.375
width_double_column = 6.75

# Default ratio for a single plot figure
# I prefer a little higher than goden ratio, from 0.618 to about 0.68
height_width_ratio = GOLDEN_RATIO * 1.1  # = height / width

_width = width_single_column
_height = width_single_column * height_width_ratio

_params = {'font.family': 'sans-serif',
           'font.serif': ['Times', 'Computer Modern Roman'],
           'font.sans-serif': ['Helvetica'], # ['Helvetica', 'Computer Modern Sans serif'] are not default installed,
           'font.size': 8,
           'text.usetex': True,
           # To force LaTeX use Helvetica fonts.
           'text.latex.preamble': [r'\usepackage{siunitx}',
                                   r'\sisetup{detect-all}',
                                   r'\usepackage{helvet}',
                                   r'\usepackage[eulergreek,EULERGREEK]{sansmath}',
                                   r'\sansmath'],
           'axes.prop_cycle': default_color_cycler,
           'axes.labelsize': 8,
           'axes.linewidth': 1,

           'figure.figsize': (_width, _height),
           # 'figure.subplot.left' : 0.125,
           # 'figure.subplot.right' : 0.95,
           # 'figure.subplot.bottom' : 0.1,
           # 'figure.subplot.top' : 0.95,

           'savefig.dpi': 600,
           'savefig.format': 'eps',
           # 'savefig.bbox': 'tight',
           # this will crop white spaces around images that will make
           # width/height no longer the same as the specified one.

           'legend.fontsize': 7.5,
           'legend.frameon': False,
           'legend.numpoints': 1,
           'legend.handlelength': 2,
           'legend.scatterpoints': 1,
           'legend.labelspacing': 0.5,
           'legend.markerscale': 0.9,
           'legend.handletextpad': 0.5,  # pad between handle and text
           'legend.borderaxespad': 0.5,  # pad between legend and axes
           'legend.borderpad': 0.5,  # pad between legend and legend content
           'legend.columnspacing': 1,  # pad between each legend column

           # 'text.fontsize' : 8,
           'xtick.labelsize': 8,
           'ytick.labelsize': 8,

           'lines.linewidth': 0.5,
           'lines.markersize': 3,
           # 'lines.markeredgewidth' : 0,
           # 0 will make line-type markers, such as '+', 'x', invisible
           }

aps_decorator = MPLdecorator(_params)
