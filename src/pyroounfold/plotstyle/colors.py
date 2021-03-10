# -*- coding: utf-8 -*-
"""
colors.py
=========

A set of colors and color maps as provided by python package `palettable` (https://jiffyclub.github.io/palettable).

color wiki:
https://encycolorpedia.com/
"""

import numpy as np

from palettable.colorbrewer.qualitative import Set1_9
from palettable.tableau import Tableau_10, Tableau_20
from cycler import cycler
__all__ = ['paper_color_cycler', 'tango_color_cycler','default_color_cycler', 'bonn_cycler_1','bonn_cycler_2','monobar_cycle','monochrome']
# Set some commonly used colors
almost_black = '#262626'
light_grey = np.array([float(248) / float(255)] * 3)

# ColorBrewer
# by Drs. Cynthia Brewer and Mark Harrower of Pennsylvania
# State University. For more information on ColorBrewer, see:
# - Flash-based interactive map:
# http://colorbrewer2.org/
# - A quick visual reference to every ColorBrewer scale:
# http://bl.ocks.org/mbostock/5577023
# https://jiffyclub.github.io/palettable/colorbrewer
#
# ColorBrewer scale 'Qualitative.Set1'.
# This one has nice "traditional" colors like reds and blues
# It is suitable for multi-line plots and category items.
brewer_set1 = Set1_9.mpl_colors
# Remove the sixth color (yellow) which is too bright
brewer_set1.pop(5)
# Swap the red and blue to let blue come first
brewer_set1[0], brewer_set1[1] = brewer_set1[1], brewer_set1[0]
# Add a decent black color to this list
brewer_set1.append(almost_black)

# Tableau 10 & 20 Classic
# See: https://jiffyclub.github.io/palettable/tableau/
# Another great qualitative set of colors suitable for multi-line plots.
# 10 and 20 are number of colors in the list.
tableau_10 = Tableau_10.mpl_colors
# Add a decent black color
tableau_10.append(almost_black)
# Swap orange and red
tableau_10[1], tableau_10[3] = tableau_10[3], tableau_10[1]
# swap orange and purple
# now table_au has similar sequence to brewer_set1
tableau_10[3], tableau_10[4] = tableau_10[4], tableau_10[3]
# This is 20-color Tableau which contains light version of tableau_10
tableau_20 = Tableau_20.mpl_colors


"""
Use `cycler` to customize styles in multiline plot.

Example:
# https://matplotlib.org/3.1.1/tutorials/intermediate/color_cycle.html
custom_cycler = (cycler(color=['c', 'm', 'y', 'k']) +
                 cycler(lw=[1, 2, 3, 4]))

fig, (ax0, ax1) = plt.subplots(nrows=2)
ax0.plot(yy)
ax0.set_title('Set default color cycle to rgby')
ax1.set_prop_cycle(custom_cycler)
ax1.plot(yy)
ax1.set_title('Set axes color cycle to cmyk')

# Add a bit more space between the two plots.
fig.subplots_adjust(hspace=0.3)
plt.show()

ref:
linestyle = ['-', '--', ':', '-.']
marker:
https://matplotlib.org/3.1.1/api/markers_api.html
'marker', ['^',',', '.','s']
"""

default_color_cycler = cycler('color', tableau_10)
custom_cycler = (cycler(color=['c', 'm', 'y', 'k'])
                 + cycler(linestyle=['-', '--', ':', '-.'])
                 + cycler(lw=[0.5, 1.0, 1.5, 2.0]))
bonn_blue = '#07529a'
bonn_yellow = '#eab90c'
bonn_grey ='#909085'
almost_black = '#262626'
bonn_red = 'r'
bonn_tblue = '#81D8D0'
bonn_green = '#15f571'
bonn_magenta= '#ea0a8e'
bonn_cycler_1 = (cycler(color = [almost_black,bonn_blue,bonn_yellow,bonn_red,bonn_green,bonn_magenta,bonn_tblue])
# + cycler(linestyle=['-', '-', '-', '-','-'])
)
bonn_cycler_2 = (cycler(color = [almost_black,bonn_blue,bonn_yellow,bonn_red,bonn_green,bonn_magenta,bonn_tblue])
+ cycler(linestyle=['-','-','-', '--', '-', '-.','-']))
# monochrome
# http://olsgaard.dk/monochrome-black-white-plots-in-matplotlib.html
monobar_cycle = (cycler('hatch', ['///', '--', '...','\///', 'xxx', '\\\\']) * cycler('color', 'w')*cycler('zorder', [10]))
monochrome = (cycler('color', ['k']) * cycler('linestyle', ['-', '--', ':', '=.']) * cycler('marker', ['^',',', '.']))



class TangoColors(object):
#"""
#Provides the Tango colors.
#"""
    scarlet_red_light = '#ef2929'
    scarlet_red = '#cc0000'
    scarlet_red_dark = '#a40000'

    aluminium_light = '#eeeeec'
    aluminium = '#d3d7cf'
    aluminium_dark = '#babdb6'

    butter_light = '#fce94f'
    butter = '#edd400'
    butter_dark = '#c4a000'

    chameleon_light = '#8ae234'
    chameleon = '#73d216'
    chameleon_dark = '#4e9a06'

    orange_light = '#fcaf3e'
    orange = '#f57900'
    orange_dark = '#ce5c00'

    chocolate_light = '#e9b96e'
    chocolate = '#c17d11'
    chocolate_dark = '#8f5902'

    sky_blue_light = '#729fcf'
    sky_blue = '#3465a4'
    sky_blue_dark = '#204a87'

    plum_light = '#ad7fa8'
    plum = '#75507b'
    plum_dark = '#5c3566'

    slate_light = '#888a85'
    slate = '#555753'
    slate_dark = '#2e3436'

    default_colors = [
        sky_blue,
        orange,
        chameleon,
        scarlet_red,
        plum,
        chocolate,
        butter,
        slate,
        aluminium,
    ]


tango_color_cycler = cycler("color", TangoColors.default_colors)


class PaperColors(object):
#"""
#Provides the Tango colors.
#"""

#    p_red = '#d53e4f'
#    p_orange = '#fc8d59'
#    p_yellow = '#fee08b'
#    p_butter = '#fee08b'
#    p_spring = '#e6f598'
#    p_green = '#99d594'
#    p_blue = '#3288bd'
#
#
#    default_colors = [
#    p_red,
#    p_orange,
#    p_yellow,
#    p_butter,
#    p_spring,
#    p_green,
#    p_blue
#    ]

    p_red = '#d53e4f'
    p_light_red = '#fcbba1'
    p_orange = '#fdae61'
    p_yellow = '#ffffbf'
    p_gray_blue = '#c6dbef'
    p_light_blue = '#6baed6'
    p_blue = '#3182bd'
    p_deep_blue = '#2171b5'


    default_colors = [
    p_red,
    p_orange,
    p_yellow,
    p_gray_blue,
    p_light_blue,
    p_blue,
    p_deep_blue
    ]

paper_color_cycler = cycler("color", PaperColors.default_colors)
