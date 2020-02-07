#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:26:18 2020

@author: debbio
"""
import json
import os
import logging

class Keyword_json():
    """Json keyword reader
    
    It has to instantiated for every input files
    eg fitszilla, mbfits..
    """
    
    def __init__(self, p_fitsType):
        """
        Open and decodes json file for the selected type

        Parameters
        ----------
        p_fitsType : fits type
            'fitszilla',
            ...

        Returns
        -------
        None

        """
        self.m_jsonDefs = None
        self.m_components = None
        l_jsonFileName = p_fitsType + '.json'
        if os.path.exists(l_jsonFileName):
            with open(l_jsonFileName, 'r') as l_jsonFile:
                try:
                    #Json definition loaded
                    self.m_jsonDefs = json.load(l_jsonFile)
                except IOError:
                    logging.error('Json file definition not found for '
                                  + p_fitsType)
                    return
        # Decoding json definitions
        self.m_components = self.m_jsonDefs['components']
        
    def input_tables_list(self):
        """
        Returns input file table name list

        Returns
        -------
        Input file table name list

        """        
        if 'tables' in self.m_jsonDefs.keys():
            return self.m_jsonDefs['tables']
        return []
        
    def fitslike_keywords(self, p_component):
        """
        Returns fitslike representation keyword for selected component

        Parameters
        ----------
        p_component : p_component
            Choosen component 
            eg. 'observation', 'map', ..

        Returns
        -------
        Keyword list.

        """        
        if p_component in self.m_jsonDefs['components'].keys():
            return self.m_jsonDefs['components'][p_component]['fitslike']
        return []

    def parser(self, p_component):
        """
        Returns appropriate parser for p_component
        Parser binds input file keywords to fitslike keywords

        Parameters
        ----------
         p_component : p_component
            Choosen component
            eg. 'observation', 'map', ..

        Returns
        -------
        Parser dictionary

        """        
        if p_component in self.m_jsonDefs['components'].keys():
            return self.m_jsonDefs['components'][p_component]['parser']
        return []
        

class Keyword_helper():
    """Helper keyword handling"""
        
    @staticmethod
    def get_table_fitszilla(p_parseDict, p_key):
        """Return table index for given parse dictionary and given key"""
        if p_key not in p_parseDict.keys():
            return None
        l_tableName = p_parseDict[p_key][0]        
        if l_tableName in awarness_fitszilla_tables.keys():
            return awarness_fitszilla_tables[l_tableName]
        else:
            return None
        
    @staticmethod
    def is_header_or_data(p_parseDict, p_key):
        """Decode awarness parsign dictionary locating where the keyword is:
            in header or in data table
            
            Returns True for 'HEADER', false for 'DATA'
        """
        if p_key not in p_parseDict.keys():
            return None
        if p_parseDict[p_key][1] == 'HEADER':
            return True
        else:
            return False
  
    @staticmethod
    def get_zilla_key(p_parseDict, p_key):
        """Get fitszilla keyword for selected dictionary keyword"""
        if p_key not in p_parseDict.keys():
            return None
        return p_parseDict[p_key][2]