# -*- coding = utf-8 -*-
# @Time: 2024/7/31 14:11
# @File: __init__.PY
"""Backbone Rigid Invariant (BRI) - A Python Package for Protein Structure Analysis

The BRI package provides tools for analyzing protein structures using geometric
descriptors based on backbone rigid invariants. It enables:

- Extraction and cleaning of protein structure data from PDB/mmCIF files
- Computation of backbone rigid invariants (BRI)
- Statistical analysis of BRI distributions
- Structural comparison using BRI descriptors
- Visualization of BRI patterns through diagrams and barcodes

Main Components
--------------
===================    ===================================
Module                 Description
===================    ===================================
:mod:`.base`           Core utilities and mathematical functions
:mod:`.filter`         Structure validation and cleaning tools
:mod:`.invariant`      BRI computation and analysis
:mod:`.pdbx2df`        Structure data handling and representation
:mod:`.invariant_compare` Structural comparison using BRI
===================    ===================================

Example Usage
------------
>>> from bri import Entry, Chain
>>> # Load a structure
>>> entry = Entry("1abc")
>>> # Compute BRI
>>> invariants = entry.get_entry_invariant()
>>> # Visualize BRI
>>> chain = entry[0]
>>> chain.generate_BID()

For more information and documentation, visit:
https://backbone-rigid-invariant.readthedocs.io
"""

__version__ = "1.2.2"
__author__ = "Ziqiu Jiang"
__maintainer__ = "Ziqiu Jiang"
__email__ = "jiangziqiu02@gmail.com"
__license__ = "CC-BY-NC-SA-4.0"
__copyright__ = "Copyright 2024, Ziqiu Jiang"


from .pdbx2df import Chain, Entry, MiniChain, MiniEntry
from .invariant_compare import (
    group_invariant_compare,
    get_theo_lipschitz_constant,
    get_e_lipschitz_constant,
)
from .invariant import *
from .base import base_util
from .filter import integrated_chainwise_filter, entry_integrated_cleaning, minientry_integrated_cleaning

# https://average-minimum-distance.readthedocs.io.
