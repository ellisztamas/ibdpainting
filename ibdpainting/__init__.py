"""Top-level package for ibdpainting."""

__author__ = """Tom Ellis"""
__email__ = 'thomas.ellis@gmi.oeaw.ac.at'
__version__ = '0.4.0'

import sys

# Only import everything if we're not running from command line
if not sys.argv[0].endswith('ibdpainting'):
    from ibdpainting.load_genotype_data import load_genotype_data
    from ibdpainting.geneticDistance import *
    from ibdpainting.ibd_table import ibd_table
    from ibdpainting.ibd_scores import ibd_scores
    from ibdpainting.plot_ibd_table import plot_ibd_table

    __all__ = [
        'load_genotype_data',
        'ibd_table',
        'ibd_scores',
        'plot_ibd_table'
    ]