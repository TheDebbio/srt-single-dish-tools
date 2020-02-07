# -*- coding: utf-8 -*-
"""Fits like parser

    Fits format awarness, it handles fits data toward fits like representation
"""

import astropy.io.fits as fits
import fitslike_keywords

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
        self.m_parsingDicts = [fitslike_keywords.awarness_fitszilla_obsdata]
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
        for l_parsingDict in self.m_parsingDicts:
            for l_key in l_parsingDict.keys():
                l_intermediate[l_key] = None
                # table number
                l_table = \
                    fitslike_keywords.Keyword_helper.\
                        get_table_fitszilla(l_parsingDict, l_key)
                if l_table is None:
                    return
                l_fitsTable = self.m_fitszilla[l_table]
                # is in header or in data ?
                l_headerOrData = fitslike_keywords.Keyword_helper.\
                    is_header_or_data(l_parsingDict, l_key)
                # looking for keyword or data
                l_fitsTableContent = None
                if l_headerOrData is True:
                    l_fitsTableContent = l_fitsTable.header
                else:
                    l_fitsTableContent = l_fitsTable.data
                l_fitsKey = fitslike_keywords.Keyword_helper.\
                    get_zilla_key(l_parsingDict, l_key)
                l_intermediate[l_key] = l_fitsTableContent[l_fitsKey]
        return l_intermediate
