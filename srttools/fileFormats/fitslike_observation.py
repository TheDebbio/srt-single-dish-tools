# -*- coding: utf-8 -*-
"""Observation data container"""

import fitslike_keywords
import fitslike_component


class Fitslike_observation(fitslike_component.Fitslike_component):
    """Container for observation fitslike data"""

    def __init__(self, p_componentName, p_keywords):
        """
        Obs data generation starting from input file dict and a dedicated
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
