#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 18:55:45 2020

@author: debbio

Fitslike Handler objects testing
"""

import pdb
import logging
import fitslike_commons
import fitslike_handler

#nodding_dir = '/home/debbio/discos/Scan/SARDARA_Nodding/20181210-232201-24-18-W3OH_small/'
#nodding_dir = '/home/debbio/discos/Scan/SARDARA_Nodding/very_small/'
nodding_dir = '/home/debbio/discos/Scan/SARDARA_Nodding/20181210-232201-24-18-W3OH/'
output_path= '/home/debbio/discos/Scan/SARDARA_Nodding/debbio-20181210-232201-24-18-W3OH/'
nodding_zilla_1_8 = '20181210-232307-24-18-W3OH_001_008.fits'

class TestFitslike_handler():
    """Fitslike handler unit test"""
    @staticmethod
    def test_scan():
        """Scan input file test"""
        l_fh= fitslike_handler.Fitslike_handler('fitszilla', 'nod')        
        l_fh.setOutputPath(output_path)
        l_fh.scan_data(nodding_dir)
        l_fh.group_on_off_cal()
        l_fh.normalize()
        l_fh.ClassFitsAdaptations()
        l_fh.classfitsWrite()
        #l_fh._on_off_match()
        #pdb.set_trace()
        
if __name__ == "__main__":
    l_commons = fitslike_commons.Fitslike_commons()    
    l_logger = logging.getLogger(l_commons.logger_name())
    l_formatter =logging.Formatter('[%(filename)s:%(lineno)s - %(funcName)s() ] %(message)s')
    l_logCh = logging.StreamHandler()            
    l_logCh.setFormatter(l_formatter)
    l_logger.setLevel(logging.INFO)
    l_logger.info("Testing fitslike handler unit")
    if not len(l_logger.handlers):
        l_logger.addHandler(l_logCh)
    # Fitslike    
    l_fhu= TestFitslike_handler()
    l_fhu.test_scan()
    
    


