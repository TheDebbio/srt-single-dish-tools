# -*- coding: utf-8 -*-
"""Fitslike commons helper"""


from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import astropy.units as unit
import numpy as np

keywords= {
    "key_on":"SIGNAL",
    "key_off": "REFERENCE",
    "key_off_cal": "REFCAL",
    "key_on_cal": "REFSIG",
    "key_sig_cal": "SIGCAL",
    "keys_on": ["SIGNAL", "SIGCAL", "REFSIG"],
    "keys_off": ["REFCAL", "REFERENCE"],
    "keys_cal_on": ["REFCAL", "SIGCAL", "REFSIG"],
    "keys_cal_off": ["REFERENCE", "SIGNAL"]
    
}

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
        return 'root'
    
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
        return locations[p_site.lower()]

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
    
    @staticmethod
    def calculate_weather(p_tmp, p_u):
        """Get the meters of H2O, using the formula from old converter.
        Unsure if this is correct in all cases"""
        RS = 8.314472
        mw = 0.018015
        md = 0.0289644
        eps0 = mw / md
        k1 = 77.60
        k2 = 70.4
        k3 = 3.739E5    
        p_tmp = p_tmp - 273.15
        H = (np.log10(p_u) - 2.0) / 0.4343 + (17.62 * p_tmp) / (243.12 + p_tmp)    
        DPT = 243.12 * H / (17.62 - H)    
        p_tmp = p_tmp + 273.15    
        Tm = 0.673 * p_tmp + 83.0
        C = 1E6 * mw / (k2 - k1 * eps0 + k3 / Tm) / RS
        e0 = np.exp(1.81 + 17.27 * DPT / (DPT + 237.5))
        ZWDS = 0.002277 * (0.005 + 1255 / p_tmp) * e0    
        return ZWDS * C * 100.