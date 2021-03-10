# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = 'PyrooUnfold'
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound

from .decorator import *
from .colors import *
from .utilities import *
from .layout import *
from .style_aps import *
from .style_thesis import *
from .style_presentation import *
from .style_WG1 import *
