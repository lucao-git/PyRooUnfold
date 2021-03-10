"""
from
https://github.com/MarkusPrim/WG1Template.git
"""

import os
from typing import Union, Optional, Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
__all__ = ['export',]

def export(
        fig: plt.Figure,
        filename: Union[str, os.fspath],
        target_dir: str = "plots/",
        file_formats: Tuple[str] = (".pdf", ".png")
) -> None:
    """
    Convenience function for saving a matplotlib figure.

    :param fig: A matplotlib figure.
    :param filename: Filename of the plot without .pdf suffix.
    :param file_formats: Tuple of file formats specifying the format
    figure will be saved as.
    :param target_dir: Directory where the plot will be saved in.
    Default is './plots/'.
    :return: None

    example:
    export(fig, 'simple', 'examples')
    """
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)

    for file_format in file_formats:
        fig.savefig(os.path.join(target_dir, f'{filename}{file_format}'), bbox_inches="tight")