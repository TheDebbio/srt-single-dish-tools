# -*- coding: utf-8 -*-
"""Fitslike commons helper"""


from astropy.coordinates import EarthLocation, AltAz, Angle, ICRS
import astropy.units as unit
import numpy as np
import re

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
            'srt': EarthLocation(4865182.7660, 791922.6890, 4035137.1740, unit=unit.m),
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
    def telescope_label_from_channel_name(p_chName):
        """
        Generation telescope name for feed
        """            
        
        def interpret_chan_name(chan_name):
            """Get feed, polarization and baseband info from chan name.
        
            Examples
            >>> feed, polar, baseband = interpret_chan_name('blablabal')
            >>> feed  # None
            >>> polar  # None
            >>> baseband  # None
            >>> feed, polar, baseband = interpret_chan_name('Ch0')
            >>> feed
            0
            >>> polar  # None
            >>> baseband  # None
            >>> feed, polar, baseband = interpret_chan_name('Feed1_LCP')
            >>> feed
            1
            >>> polar
            'LCP'
            >>> baseband  # None
            >>> feed, polar, baseband = interpret_chan_name('Feed2_LCP_3')
            >>> feed
            2
            >>> polar
            'LCP'
            >>> baseband
            3
            """
            chan_re = re.compile(r'^Ch([0-9]+)$'
                         r'|^Feed([0-9]+)_([a-zA-Z]+)$'
                         r'|^Feed([0-9]+)_([a-zA-Z]+)_([0-9]+)$')
            
            matchobj = chan_re.match(chan_name)
            if not matchobj:
                return None, None, None
        
            matches = [matchobj.group(i) for i in range(7)]
            polar, baseband = None, None
            if matches[6] is not None:
                baseband = int(matchobj.group(6))
                polar = matchobj.group(5)
                feed = int(matchobj.group(4))
            elif matches[3] is not None:
                polar = matchobj.group(3)
                feed = int(matchobj.group(2))
            else:
                feed = int(matchobj.group(1))
    
            return feed, polar, baseband
        
        _, polar, _ = interpret_chan_name(p_chName)
    
        if polar.startswith('L'):
            return 'LL'
        elif polar.startswith('R'):
            return 'RR'
        elif polar.startswith('Q'):
            return 'LR'
        elif polar.startswith('U'):
            return 'RL'
        else:
            raise ValueError('Unrecognized polarization')
    
    @staticmethod 
    def class_telescope_name(p_scan):
        """
        Generazione strigna telescope classfits
        @todo da definire..tradurre con le funziona sopra
        """
        return '{}-{}-{}'.format(
                        p_scan['scheduled']['antenna'],
                        p_scan['frontend']['feed'],                        
                        p_scan['frontend']['polarizations']
                        )

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