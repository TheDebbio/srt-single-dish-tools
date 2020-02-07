# -*- coding: utf-8 -*-
"""
Root representation of fits like data set

It's the contatiner for a generic scan data set representation
"""
import fitslike_commons
import fitslike_obsData


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
        self.m_obsdata = fitslike_obsData.Fitslike_obsdata()
        self.m_componentdata = {'observation' : self.m_obsdata}
        
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
        for l_component in self.m_componentdata.values():
            l_component.update(p_key, p_value)

    def dump(self):
        """Object contents dumping"""
        for l_component in self.m_componentdata.values():
            l_component.dump()
                    