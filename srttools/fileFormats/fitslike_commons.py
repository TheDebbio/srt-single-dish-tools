# -*- coding: utf-8 -*-
"""Fitslike commons helper"""


from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import astropy.units as unit


class Fitslike_commons():
    """Helper class"""

    @staticmethod
    def logger_name():
        """
        Getter logger name

        Returns
        -------
        String with logger name.

        """
        return 'fitlike_logger'
    
    @staticmethod
    def build_attribute_dict(p_lstAttributes):
        """Build fitlslike component attribute dictionary
        Parameters
        ----------
        p_lstAttributes : string list
            List with dictionary entries

        Returns
        -------
        Attribute dictionary with 'value' and  'type' entries
        """
        l_outDict = {}
        for l_attribute in p_lstAttributes:
            l_outDict[l_attribute] = {'value': '', 'type': ''}
        return l_outDict.copy()
    
    @staticmethod
    def dump_attribute(p_key, p_attributeValue):
        """fitslike component value dumping
        component's dict :
            {'key' : {'value' :'', 'type':''}}
        It prints key and value type part
        """
        print(p_key + " " + str(p_attributeValue))
        
    @staticmethod
    def get_site_location(p_site):
        """
        Getter site earth location
        """
        locations = {
            'srt': EarthLocation(4865182.7660, 791922.6890, 4035137.1740,
                                  unit=unit.m),
             'medicina': EarthLocation(Angle("11:38:49", unit.deg),
                                       Angle("44:31:15", unit.deg),
                                       25 * unit.meter),
             'greenwich': EarthLocation(lat=51.477*unit.deg, lon=0*unit.deg)}
        return locations[p_site]

    @staticmethod
    def filter_input_files(p_type):
        """
        Input file extension filter based         
        
        Parameters
        ----------
        p_type : string
            input data file type

        Returns
        -------
        None.

        """    
        l_inputDict={
            'fitszilla': r"[\s\S]*\.fits"
            }
        try:
            return l_inputDict[p_type]
        except KeyError:
            return ''
        return ''