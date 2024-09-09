"""Top-level package for ibdpainting."""

__author__ = """Tom Ellis"""
__email__ = 'thomas.ellis@gmi.oeaw.ac.at'
__version__ = '0.2.0'

from ibdpainting.load_genotype_data import load_genotype_data
from ibdpainting.geneticDistance import *
from ibdpainting.ibd_table import ibd_table
from ibdpainting.ibd_scores import ibd_scores
from ibdpainting.plot_ibd_table import plot_ibd_table