# -*- coding: utf-8 -*-

import os
import re
import logging
import pdb
from astropy.io import fits
import numpy as np
import astropy.units as unit
import fitslike_commons
import sys
import fitslike
import awarness_fitszilla
from multiprocessing import Pool
from fitslike_commons import keywords as kws


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
      p_logger.info("Scanning " + p_path)            
          
      if  p_ftype== 'fitszilla':            
          if "summary" in p_path:
              "todo summary.fits"
              return {}
          else:
              p_logger.info("fitszilla parsing: " + p_path)
              l_fits = fits.open(p_path)
              l_aware= awarness_fitszilla.Awarness_fitszilla(l_fits) 
              l_aware.parse()
              l_repr= l_aware.process()
              l_fits.close()
              l_fitslike= fitslike.Fitslike(l_repr)
              l_fitslike.data_channel_integration()               
              l_aware = None
              l_fits =None                  
              return l_fitslike.get_inputRepr()          
      return {}
                  
class Fitslike_handler():
    """    
    Subscan master operation handler.
    Based on input data type starts appropiate parsing/processing
    
    Reading
        Takes a path as input directory with source'subscans
        Spawns multiprocessing parsing of input subscan files
        In particular add a process to Pool for every fits input file through
        global(at module level) functions:
            _envelope_subscan
            _subscan_callback
        This approch comes in handy to overcome  Pickle serialization 
        limitations/complexity when dealing with class methods/isntances
        _subscan_callback imposes parsing data are returned as dictionary
        i.e. as python built-in type
        
    
    Processing
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
        self.m_subscans=[]
        self.m_dataDir =""        
        self.m_group_on_off_cal={
            'on':[],
            'off':[],
            'cal':[],
            }
     
        
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
        global g_subscans
        
        def chunker(p_seq, p_size):
            return (p_seq[pos:pos + p_size] for pos in range(0,len(p_seq), p_size))
        
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
        "todo Memory Leak ad ogni pool ! se prelevi il risutlato"
        " Questione non risolta"
        self.m_results=[]
        l_poolSize= 3
        for l_group in chunker(l_inputFiles, l_poolSize ):
            l_results=[]
            self.m_pool=Pool(l_poolSize)            
            for l_fPath in l_group:                            
                l_results.append( self.m_pool.apply_async(
                        _envelope_subscan,
                        [self.m_logger,
                        'fitszilla',
                        p_dataDir + l_fPath
                        ])      
                    )                          
            self.m_pool.close()
            self.m_pool.join()                        
            self.m_subscans= self.m_subscans + [x.get() for x in l_results]
        self.m_logger.info("subscan numbers " + str(len(self.m_subscans)))        
        
            
        
    def group_on_off_cal(self):
        """
        Creates a new dict for :
            - On
            - Off
            - Cal
        g_on_off_cal_dict will store data groups per type
        Subscan collection parsing  per feed
        """        
        
        for l_subscan in self.m_subscans:
            for l_chx in l_subscan:
                # regroup by ref                 
                if l_subscan[l_chx]['scheduled']['signal'] == kws["key_on"]:
                   self.m_group_on_off_cal['on'].append(l_subscan[l_chx])
                
                if l_subscan[l_chx]['scheduled']['signal'] == kws["key_off"]:
                    self.m_group_on_off_cal['off'].append(l_subscan[l_chx])
                        
                if l_subscan[l_chx]['scheduled']['signal'] == kws["key_cal"]:
                    self.m_group_on_off_cal['cal'].append(l_subscan[l_chx])      
  
    def _on_off_match(self):                  
        """
        Look the best off for every on.
        It works per feed.
        Best match is intended where shortest distance occours
        Adds a key reference to best off ch_xoff at ch_xon  
        Looks for 
            ['ch_x']['integrated_data']
        Adds to ch_x new key ['best_off'] as reference to best off channel
        """
        for l_on in self.m_group_on_off_cal['on']:
            l_best_off={}
            l_dist_min = 10E6 * unit.deg
            l_el1= l_on['integrated_data']['data_el']
            l_az1= l_on['integrated_data']['data_az']            
            for l_off in self.m_group_on_off_cal['off']:            
                l_el2= l_off['integrated_data']['data_el']
                l_az2= l_off['integrated_data']['data_az']               
                l_dist_loc= np.sqrt(np.square(l_el1-l_el2)+\
                                    np.square(l_az1-l_az2))
                if ( l_dist_loc < l_dist_min):
                    l_dist_min= l_dist_loc
                    l_best_off= l_off
            l_on['best_off'] = l_best_off
            "todo eseguire il calcolo on meno off"
            
                
  