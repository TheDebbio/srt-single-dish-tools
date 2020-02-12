# -*- coding: utf-8 -*-
"""Fitslike objects testing
"""

import fitslike_commons
import fitslike
import awarness_fitszilla
import pytest
from astropy.io import fits
import pdb
import logging

nodding_dir = '/home/debbio/discos/Scan/SARDARA_Nodding/20181210-232201-24-18-W3OH/'
nodding_zilla_1_8 = '20181210-232307-24-18-W3OH_001_008.fits'


class TestFitsLike_arch():
    """Fitslike architecture test class"""

    @staticmethod
    def test_instance():
        """Instancing fitslike object"""
        l_fits = fitslike.Fitslike()
        #l_fits.dump()

    @staticmethod
    def test_Awareness_fitszilla_parse(p_file):
        """Parse a fitszilla through Awarness_fitszilla"""
        l_fits = fits.open(p_file)        
        l_aware = awarness_fitszilla.Awarness_fitszilla(l_fits)
        #pdb.set_trace()
        l_intermediateDict = l_aware.parse()
        #print(l_intermediateDict)
        #pdb.set_trace()
        l_processedDict = l_aware.process() 
        pdb.set_trace()
        #print(l_processedDict)

if __name__ == "__main__":
    l_commons = fitslike_commons.Fitslike_commons()    
    l_logger = logging.getLogger(l_commons.logger_name())
    l_formatter =logging.Formatter('[%(filename)s:%(lineno)s - %(funcName)s() ] %(message)s')
    l_logCh = logging.StreamHandler()            
    l_logCh.setFormatter(l_formatter)
    l_logger.setLevel(logging.INFO)
    l_logger.info("Testing fitslike unit")
    if not len(l_logger.handlers):
        l_logger.addHandler(l_logCh)
    # Fitslike    
    l_unitArch = TestFitsLike_arch()
    l_unitArch.test_instance()
    # Awareness
    l_file = nodding_dir + nodding_zilla_1_8
    l_unitArch.test_Awareness_fitszilla_parse(l_file)
    
