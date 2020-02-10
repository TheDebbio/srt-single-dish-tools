# -*- coding: utf-8 -*-
"""
Root representation of fits like data set

It's the contatiner for a generic scan data set representation
"""
import fitslike_commons
import fitslike_obsData
import os
import json
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

    def __init__(self):
        """Composing Fitslike object"""
        self.m_fitslikeJson = None
        self._load_json()
        self.m_componentData = None
        self._build_components()

    def _load_json(self):
        """
        Private, handles loading json fitslike defs

        Returns
        -------
        None

        """
        if os.path.exists('fitslike.json'):
            with open('fitslike.json', 'r') as l_jsonFitslike:
                try:
                    self.m_fitslikeJson = json.load(l_jsonFitslike)
                except IOError:
                    logging.error('Json file definition not' +
                                  'found for fitslike')

    def _build_components(self):
        """
        Private, fitslike component builder

        Returns
        -------
        None.

        """
        l_obsKeywords = self.m_fitslikeJson['observation']
        self.m_obsData = fitslike_obsData.Fitslike_obsdata('observation',
                                                           l_obsKeywords)
        self.m_componentData = {'observation': self.m_obsData}

    def update(self, p_key, p_value):
        """
        Ask to every component to update the key, if the key is owned by the
        component itself.
        Every key added is unique hence it's usefull to parse additional
        tables from multi-table fits input files

        Parameters
        ----------
        p_dict : dictionary
            { column : data - type }

        p_parser : parser object
            Parser per estrarre la keyword e trasformarla in rappresentazione
            generica

        Returns
        -------
        None.

        """
        for l_component in self.m_componentData.values():
            l_component.update(p_key, p_value)

    def dump(self):
        """Object contents dumping"""
        for l_component in self.m_componentData.values():
            l_component.dump()
