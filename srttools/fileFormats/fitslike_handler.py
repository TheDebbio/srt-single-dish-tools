# -*- coding: utf-8 -*-

import os
import re
import logging
import pdb
from astropy.io import fits
import fitslike_commons
import fitslike
import awarness_fitszilla
from multiprocessing import Pool


g_subscans= [] # Subscans module data collector

def _envelope_subscan(p_logger, p_ftype, p_path):
      """
      Embeds Fitslike subscan parsing 

      Parameters
      ----------
      
      p_path : string
          fits file type.
      p_path : string
          path fits file.

      Returns
      -------
      ....

      """
      "todo comporre il parsing attraverso il fitslike ?"      
      if  p_ftype== 'fitszilla':            
          p_logger.info("fitszilla parsing: " + p_path)
          l_fits = fits.open(p_path)
          l_aware= awarness_fitszilla.Awarness_fitszilla(l_fits) 
          l_aware.parse()
          l_repr= l_aware.process()
          #l_fitslike= fitslike.Fitslike(l_repr)
          #return l_fitslike.get_inputRepr()
          return l_repr
      return {}
  
def _subscan_callback(l_dict):
    """
    Call back at the end of subscan parsing process
    Adds subscan dict to subscan data

    Parameters
    ----------
    l_dict : dict
        Subscan fitslike dict.

    Returns
    -------
    None.

    """    
    g_subscans.append(l_dict)                    

class Fitslike_handler():
    """
    Subscan handler
    It takes fitslike subscan data carrying:
        - data
        - coordinates
        - observation data ( son off cal, site etc..)
    It applies normalization and calibration 
    Providing data ready to be studied or to be written to disc (through
    appropriate output conversion module)    
    """
    
    def __init__(self, p_inputType):
        """
        Internal struct definition

        Parameters
        ----------
        p_dataDir : string
            Input scan data type
            Fitszilla etc..

        Returns
        -------
        None.

        """
        self.m_commons = fitslike_commons.Fitslike_commons()
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        self.m_inputType = p_inputType
        self.m_fitslikes= g_subscans
        self.m_dataDir =""
        self.m_pool= Pool()
     
        
    def scan_data(self, p_dataDir):
        """
        Takes data input directory and launchs subscan data conversion
        to fitslike

        Parameters
        ----------
        p_dataDir : string
            Input data directory path

        Returns
        -------
        None.

        """        
        # Parsing dir
        self.m_dataDir = p_dataDir        
        if not os.path.isdir(p_dataDir):
            self.m_logger.error("Input data dir is not valid")
            return
        # Processing
        l_filt= self.m_commons.filter_input_files(self.m_inputType)
        l_inputFiles= []
        for dirpath, dirnames, filenames in os.walk(self.m_dataDir):
            l_inputFiles.extend(filenames)                    
        l_inputFiles= [f for f in l_inputFiles if re.match(l_filt,f)]
        # Split parsing multiprocessing        
        for l_fPath in l_inputFiles:
            self.m_logger.info("Subscan name: " + l_fPath)
            self.m_pool.apply_async(
                _envelope_subscan,
                args= (self.m_logger,
                        'fitszilla',
                        p_dataDir + l_fPath),
                 callback= _subscan_callback
                )
        self.m_pool.close()
        self.m_pool.join()
        print(len(self.m_fitslikes))
        pdb.set_trace()
        
        
        