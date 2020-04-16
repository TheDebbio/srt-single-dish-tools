# -*- coding: utf-8 -*-

import os
import re
import logging
from astropy.io import fits
from collections.abc import Iterable
import numpy as np
import astropy.units as unit
from astropy.time import Time
import astropy.constants as const
import fitslike_commons
import sys
import fitslike
import awarness_fitszilla
from multiprocessing import Pool
from fitslike_commons import keywords as kws
import os
import pdb

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
        l_path, l_filename= os.path.split(p_path)
        p_logger.info("fitszilla parsing: " + p_path)
        l_fits = fits.open(p_path)
        l_aware= awarness_fitszilla.Awarness_fitszilla(l_fits, p_path) 
        l_aware.parse()
        l_repr= l_aware.process()
        l_fits.close()
        l_fitslike= fitslike.Fitslike(l_repr)       
        if not l_fitslike.is_summary():
            l_fitslike.data_channel_integration()               
        l_aware = None
        l_fits =None                  
        l_repr= l_fitslike.get_inputRepr()          
        l_repr['file_name']= l_filename
        return l_repr
                  
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
    
    def __init__(self, p_inputType, p_scanType):
        """
        Internal struct definition

        Parameters
        ----------
        p_dataDir : string
            Input scan data type
            Fitszilla etc..
        p_dataDir : string
            Scan type, on-off, nodding, map

        Returns
        -------
        None.

        """
        self.m_commons = fitslike_commons.Fitslike_commons()
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        self.m_inputType = p_inputType
        self.m_scanType= p_scanType
        self.m_subscans=[]
        self.m_dataDir =""        
        self.m_group_on_off_cal={}
        self.m_summary= {}
     
    def get_summary(self):
        """Getter summary"""
        return self.m_summary
    
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
            - Cal_on
            - Cal_off
            For every ch_x
            
        it means:
            ch_x{
                on:[]
                off:[]
                cal_on:[]
                cal_off:[]
            }
            
        g_on_off_cal_dict will store data groups per type
        Subscan collection parsing per feed
        Output groups are separeted basing on scan type,
        but generally per feed
        """
        
        def _is_cal(p_subScan):
            """
            Genera il flag calibrazione attiva

            Parameters
            ----------
            p_sub : dict
                subscan ch_x

            Returns
            -------
            true calibrazione attiva

            """            
            l_signal= p_subScan['scheduled']['signal']
            l_flagCal= p_subScan['spectrum']['flag_cal']            
            if l_signal in kws['keys_cal_on']:
                return True
                " @todo keyword subtype ? "
                " Controllar anche flag_cal any"
            elif np.any(l_flagCal):
                return True
            return False
        
        def _is_on(p_subScan):
            """
            Verifica se la subscan è un on o un off

            Parameters
            ----------
            p_subscan : dict
                subscan ch_x

            Returns
            -------
            None.

            """
            l_signal= p_subScan['scheduled']['signal']            
            l_feed= p_subScan['frontend']['feed']
            if l_signal in kws['keys_on']:
                if l_feed == 0:
                    return True
                else:
                    return False
            else:
                if l_feed == 0:
                    return False
                else:
                    return True
            if l_signal == None :
                return p_subScan['ccordinates']['az_offset'] > 1e-4 * unit.rad
                                
        " ordino le subscan in base al file name "        
        self.m_subscans= sorted(self.m_subscans,\
                                key= lambda item:('file_name' not in item, item.get('file_name', None)))        
            
        for l_subscan in self.m_subscans:     
            if 'file_name' in l_subscan:
                self.m_logger.info(l_subscan['file_name'])
            " separo il summary "
            if 'summary' in l_subscan.keys():
                self.m_summary= l_subscan
                self.m_logger.info("summary.fits excluded from on off cal")
                continue
            " raggruppo le scan non summary.fits "                
            for l_chx in l_subscan:
                if 'ch_' not in l_chx:                                     
                    continue
                if self.m_scanType == 'on_off' or self.m_scanType == 'nod':
                    """                    
                    on_data = table[~cal_on & onsource]
                    off_data = table[~cal_on & ~onsource]
                    calon_data = table[cal_on & onsource]
                    caloff_data = table[cal_on & ~onsource]
                    """ 
                    l_feed= l_subscan[l_chx]['frontend']['feed']                                                                      
                    l_isCal= _is_cal(l_subscan[l_chx])
                    l_isOn= _is_on(l_subscan[l_chx])                
                    " list per feed  ch_x setup "
                    l_onOffdict= {
                        'on':[], 'off':[],'cal_on':[], 'cal_off':[]
                        }
                    try:
                        if l_chx not in self.m_group_on_off_cal.keys():
                            self.m_group_on_off_cal[l_chx]=l_onOffdict                            
                    except KeyError as e :
                        pdb.set_trace()
                        self.m_logger.error("key error: " + str(e))
                        continue
                    
                    " Grouping feed sub scans on of cal on off"
                    if l_isOn and not l_isCal:                            
                        self.m_group_on_off_cal[l_chx]['on'].append(l_subscan[l_chx])           
                        self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is on ')
                    elif not l_isOn and not l_isCal:                        
                        self.m_group_on_off_cal[l_chx]['off'].append(l_subscan[l_chx])
                        self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx +  ' is off ')
                    elif l_isOn and  l_isCal:
                        self.m_group_on_off_cal[l_chx]['cal_on'].append(l_subscan[l_chx])
                        self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is cal_on ')
                    elif not l_isOn and  l_isCal:
                        self.m_group_on_off_cal[l_chx]['cal_off'].append(l_subscan[l_chx])                                            
                        self.m_logger.info('feed ' + str(l_feed)+ ' ' + l_chx + ' is cal off ')
                if self.m_scanType == 'map':
                    " @todo on off in caso di mappe "
                    pass
    
    def normalize(self):
        """                
        On - Off - Cal normalize 
        
        cal on off, off , on are grouped:
            nod, on off: left right polarization
                group_on_off_cal['on']['0']['ch_0']=on ch_0
                group_on_off_cal['on']['0']['ch_1']=on ch_1
                and so on
            map : ?
            
        In general terms calculations resemble the following:
            
            off, average
            cal_on, average
            cal_off, average
            
            (on - offAvg) / offAvg
            
            (calOn - on[0])/ on[0] [it should be onRef]
            (calOff - off[0])/ off[0] [it should be onRef]
            
            Check calOn calOff Normalized , nan, inf 
            
            if [calOnNorm, calOffNorm] is a long array. average / median the array 
            
            calibrationFactor= 1 / [calOnNorm, calOffNorm] * cal_temp
    
            normalized= on * calibrationFactor
            
        """         
        for ch in self.m_group_on_off_cal.keys():
            l_group= self.m_group_on_off_cal[ch]        
            l_calMarkTemp= l_group['on'][0]['frontend']['cal_mark_temp']                         
            l_offAvg= sum(v['integrated_data']['spectrum'] for v in l_group['off']) / len(l_group['off'])            
            if len(l_group['cal_on']):
                l_calOnAvg= sum(v['integrated_data']['spectrum'] for v in l_group['cal_on']) / len(l_group['cal_on'])
            if len(l_group['cal_off']):
                l_calOffAvg= sum(v['integrated_data']['spectrum'] for v in l_group['cal_off']) / len(l_group['cal_off'])
            
            "calcolo"
                        
            l_group['on_off']=[]            
            l_group['on_off_cal']=[]
            for elOn in l_group['on']:
                on= elOn['integrated_data']['spectrum']
                on_off= (on - l_offAvg)/l_offAvg
                l_group['on_off'].append(on_off)
                        
            #pdb.set_trace()
            cal = np.array([l_calOnAvg, l_calOffAvg])
            good = (cal != 0) & ~np.isnan(cal) & ~np.isinf(cal)
            cal = cal[good]
            if len(cal) > 0:
                meancal = np.median(cal) if len(cal) > 30 else np.mean(cal)    
                calibration_factor = 1 / meancal * l_calMarkTemp
            else:   
                return None, ""
            " Calibrated spectrum added to chx"            
            for elOn in l_group['on']:                
                elOn['integrated_data']['calibrated'] = elOn['integrated_data']['spectrum'] * \
                                    calibration_factor
                
                
    def _on_off_match(self):                  
        """
        Look the best off for every on.
        It works per feed.
        Best match is intended closest scan in time
        Adds a key reference to best off ch_xoff at ch_xon  
        Looks for 
            ['ch_x']['integrated_data']
        Adds to ch_x new key ['best_off'] as reference to best off channel
        """
        for l_on in self.m_group_on_off_cal['on']:            
            l_best_off={}
            l_dist_min = 10E6 
            l_onTime= l_on['integrated_data']['data_time']            
            for l_off in self.m_group_on_off_cal['off']:            
                l_offTime= l_off['integrated_data']['data_time']               
                l_dist= l_onTime - l_offTime
                l_dist= abs(l_dist)
                if ( l_dist < l_dist_min ):
                    l_dist_min= l_dist
                    l_best_off= l_off
            l_on['best_off'] = l_best_off
            " calcolo on meno off "
            

    def ClassFitsAdaptations(self):
        """
        Generazione struttura dati secondo la definizione del classfist
        
        info base
        
        coordinate comandate in az, el o ra, dec
        coordinate osservate in crdelt2,3
        spettri separati per polarizzazione e per feed
        un file per ogni uno
        """
        " data in group on off cal sono divisi per  "
        " chx "
        "   on off cal "
        "       [chx...]"
        " @todo inserire il campo cal is on ? quindi diversificare on ed off ?"      
        " @todo gestire i dati in caso di campo singolo spettro mediato o serie di spettri"
        for l_ch in self.m_group_on_off_cal:              
            " single entry on classfits table"
            for l_ch in self.m_group_on_off_cal[l_ch]['on']:
                l_ch['classfits']={}
                " ch by ch "                                
                try:
                    " Lavoro con  i dati integrati "
                    " ut "                                                                           
                    l_tMjd= l_ch['integrated_data']['data_mjd'].mjd
                    l_ch['classfits']['UT']= ( l_tMjd - np.floor(l_tMjd)) * 86400
                    " date "                    
                    
                    l_ch['classfits']['DATE']= l_ch['integrated_data']['data_mjd'].strftime('%d/%m/%y') 
                    " lsts "
                    l_lsts= l_ch['integrated_data']['data_mjd'].sidereal_time('apparent', \
                              fitslike_commons.Fitslike_commons.\
                                  get_site_location(l_ch['scheduled']['antenna']).lon)
                    l_lsts= l_lsts.value * unit.hr                                     
                    " infos "                    
                    l_ch['classfits']['OBJECT']= l_ch['scheduled']['source']
                    l_ch['classfits']['LINE']= "F{}-{:3.3f}-MHz"\
                        .format(l_ch['frontend']['feed'], l_ch['backend']['bandwith'])
                    l_ch['classfits']['TELESCOP']= self.m_commons.class_telescope_name(l_ch)
                    " temp "
                    l_ch['classfits']['TSYS']= 1.0
                    l_ch['classfits']['CALTEMP']= l_ch['frontend']['cal_mark_temp']
                    " time "
                    l_ch['classfits']['LSTS'] = l_lsts.to('s').value     
                    l_ch['classfits']['OBSTIME']= l_ch['integrated_data']['data_integration']
                    "  "                    
                    l_ch['classfits']['CDELT1']= (l_ch['frontend']['bandwidth'] / 
                                                l_ch['backend']['bins']).to('Hz')                    
                    " freq and velocity "                    
                    l_ch['classfits']['RESTFREQ']= self.m_summary['summary']['restfreq']
                    l_ch['classfits']['VELOCITY']= l_ch['scheduled']['vlsr']                    
                    l_df= (l_ch['frontend']['bandwidth'] / l_ch['backend']['bins']).to('Hz')
                    l_ch['classfits']['CDELT1']= l_df
                    l_deltav= - l_df/ l_ch['classfits']['RESTFREQ'] * const.c
                    l_ch['classfits']['DELTAV']= l_deltav.to('m/s').value
                    " Objects Coordinates "                    
                    l_ch['classfits']['CDELT2'] = l_ch['scheduled']['ra_offset'].to(unit.deg).value
                    l_ch['classfits']['CDELT3'] = l_ch['scheduled']['dec_offset'].to(unit.deg).value
                    l_ch['classfits']['AZIMUTH']= l_ch['integrated_data']['data_az']
                    l_ch['classfits']['ELEVATION']= l_ch['integrated_data']['data_el']
                    l_ch['classfits']['CRVAL2']= l_ch['integrated_data']['data_ra']
                    l_ch['classfits']['CRVAL3']= l_ch['integrated_data']['data_dec']
                    " data "                                        
                    l_ch['classfits']['MAXIS1'] = l_ch['integrated_data']['data_integration']    
                    l_ch['classfits']['MAXIS1'] = l_ch['backend']['bins']
                    l_ch['classfits']['spectrum']= l_ch['integrated_data']['spectrum']
                except:
                    pdb.set_trace()                    
            
        def ClassfitsWrite(self):
            """
            Scrittura file con calcolo header
            header prende i dati dai dati generici ricavati dalla scansione scansione 
            """
            " creazione tabella "
            l_hdData= self.m_obs_general_data
            " header "
            l_hdu= fits.PrimaryHDU() 
            l_hdu.header['CTYPE1']= "FREQ"
            l_hdu.header['CRVAL1']= 0
            l_hdu.header['CRVAL2']= l_hdData['ra']
            l_hdu.header['CRVAL3']= l_hdData['dec']                        
            l_hdu.header['OBJECT'] = l_hdData['SOURCE']
            l_hdu.header['SOURCE'] = l_hdData['SOURCE']
            l_hdu.header['DATE-RED'] = Time.now().to_datetime().strftime('%d/%m/%y')
            l_hdu.header['LINE'] = l_hdData['LINE']
            l_hdu.header['CDELT1'] = l_hdData['CDELT1']
            l_hdu.header['RESTFREQ'] = l_hdData['RESTFREQ']
            l_hdu.header['MAXIS1'] = l_hdData['MAXIS1']
            " data "
            " @todo verificare "                    
                
                
  