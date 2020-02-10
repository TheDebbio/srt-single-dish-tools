# -*- coding: utf-8 -*-
"""Fits like parser

    Fits format awarness, it handles fits data toward fits like representation
"""
import astropy.io.fits as fits
import fitslike_keywords
import pdb


class Awarness_fitszilla():
    """fitszilla data parser"""

    def __init__(self, l_fitszilla):
        """
        Store fitzilla file

        Parameters
        ----------
        l_fitszilla : astropy fits
        Fitszilla to be handled before it gets represented as fitslike

        Returns
        -------
        None.

        """
        self.m_jsonRepr = fitslike_keywords.Keyword_json('fitszilla')
        self.m_components = self.m_jsonRepr.fitslike_components()        
        self.m_parsingDicts = {}
        for l_component in self.m_components:
            self.m_parsingDicts[l_component] = \
                self.m_jsonRepr.parser(l_component)
        self.m_fitszilla = l_fitszilla

    def parse(self):
        """Fitszilla parsing

        Creates an intermediate data representation matching given rule
        dictionaries.
        For every coponent dictionary looks for corresponding matches inside
        fitszilla input file

        Returns
        -------
        A dictionary = {
            fitslike_key : fitszilla_value
            ... : ...
        }
        """
        l_intermediate = {}        
        for l_component in self.m_parsingDicts.keys():            
            for l_key in self.m_parsingDicts[l_component].keys():
                l_intermediate[l_key] = None
                l_keyParserEntry = self.m_parsingDicts[l_component][l_key]
                # Decoding entry
                l_tableName, l_headerOrData, l_inputKeyword  =\
                    fitslike_keywords.Keyword_helper.\
                    decode_entries(l_keyParserEntry)
                if l_tableName is None:
                    return
                # todo semplificare la gestione tabella
                l_table = fitslike_keywords.Keyword_helper.\
                    get_table_index(self.m_jsonRepr.input_tables_dict(),
                                     l_tableName)
                if l_table is None:
                    return
                l_fitsTable = self.m_fitszilla[l_table]
                # is in header or in data ?
                l_headerOrData = fitslike_keywords.Keyword_helper.\
                    is_header_or_data(l_headerOrData)
                # looking for keyword or data
                l_fitsTableContent = None
                if l_headerOrData is True:
                    l_fitsTableContent = l_fitsTable.header
                else:
                    l_fitsTableContent = l_fitsTable.data                
                l_intermediate[l_key] = l_fitsTableContent[l_inputKeyword]
        return l_intermediate
