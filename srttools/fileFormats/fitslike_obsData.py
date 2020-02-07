# -*- coding: utf-8 -*-
"""Observation data container"""

import fitslike_keywords
import fitslike_component

class Fitslike_obsdata(fitslike_component.Fitslike_component):
    """Container for observation genereic data"""
    
    def __init__(self):
        """
        Obs data generation starting from input file dict and a dedicated
        parser.
        Observation data resides on m_obsdata dictionary.

        Returns
        -------
        None.
        
        """
        fitslike_component.Fitslike_component.\
        __init__(self, fitslike_keywords.fitslike_obsdata_keywords)

