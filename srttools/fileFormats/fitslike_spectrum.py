#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spectrum data container"""

import fitslike_keywords
import fitslike_component


class Fitslike_spectrum(fitslike_component.Fitslike_component):
    """Container for spectrum fitslike data"""

    def __init__(self, p_componentName, p_keywords):
        """
        Spectrum data generation starting from input file dict and a dedicated
        parser.        

        Parameters
        p_componentName: string fitslike component name
        p_keywords: Keyword list

        Returns
        -------
        None.

        """
        fitslike_component.Fitslike_component.\
        __init__(self, p_componentName, p_keywords)