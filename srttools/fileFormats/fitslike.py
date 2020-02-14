# -*- coding: utf-8 -*-
"""
Root representation of fits like data set

It's the contatiner for a generic scan data set representation
"""
import fitslike_commons
import logging


class Fitslike():
    """Base container to represent generic scan data

    It's a semantic representation of scan data fitted in fits like files
    It means important data are grouped by semantic field:
        - Generic observation data
        - Rf inputs
        - Data and their context
        - Map
    """

    def __init__(self, p_representation):
        """
        Istantiate Fitslike object
        Takes an already processed representation from *fits input file.
        This high level rep might be processed again here.
        
        
            Parameters:
                p_representation: dict
                    Generic fits like representation
                    
        """
        self.m_commons = fitslike_commons.Fitslike_commons()
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        self.m_inputRepr = p_representation.copy()
        

    def dump(self):
        """Object contents dumping"""        
        print(self.m_inputRepr)
        
    def dump_keys(self):
        """Nested keys dumping"""
        
        def _recursive(p_dict):
            """Iteration """
            for l_key, l_value in p_dict.items():
                if type(l_value) is dict:
                    yield from _recursive(l_value)
                else:
                    yield (l_key, l_value)
        
        # iterations
        for l_key, l_value in _recursive(self.m_inputRepr):
            print(l_key)
            
            