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
import shutil
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
        self.m_outputPath= ''
     
    def setOutputPath(self, p_path):
        """ output path setter"""
        self.m_outputPath= p_path
        
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
            feed{
             ch_x{
                on:[]
                off:[]
                cal_on:[]
                cal_off:[]
                }
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
            for l_feed in l_subscan:      
                #pdb.set_trace()
                for l_chx in l_subscan[l_feed]:
                    if 'ch_' not in l_chx:                                     
                        continue
                    if self.m_scanType == 'on_off' or self.m_scanType == 'nod':
                        """                    
                        on_data = table[~cal_on & onsource]
                        off_data = table[~cal_on & ~onsource]
                        calon_data = table[cal_on & onsource]
                        caloff_data = table[cal_on & ~onsource]
                        """ 
                        l_feed= l_subscan[l_feed][l_chx]['frontend']['feed']                                                                      
                        l_isCal= _is_cal(l_subscan[l_feed][l_chx])
                        l_isOn= _is_on(l_subscan[l_feed][l_chx])                
                        " list per feed ch_x setup "
                        l_onOffdict= {
                            'on':[], 'off':[],'cal_on':[], 'cal_off':[]
                            }
                        try:
                            if l_feed not in self.m_group_on_off_cal.keys():
                                self.m_group_on_off_cal[l_feed]={}
                            if l_chx not in self.m_group_on_off_cal[l_feed].keys():
                                self.m_group_on_off_cal[l_feed][l_chx]=l_onOffdict                            
                        except KeyError as e:
                            pdb.set_trace()
                            self.m_logger.error("key error: " + str(e))
                            continue
                                                
                        " Grouping feed sub scans on of cal on off"
                        if l_isOn and not l_isCal:                            
                            self.m_group_on_off_cal[l_feed][l_chx]['on'].append(l_subscan[l_feed][l_chx])           
                            self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is on ')
                        elif not l_isOn and not l_isCal:                        
                            self.m_group_on_off_cal[l_feed][l_chx]['off'].append(l_subscan[l_feed][l_chx])
                            self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx +  ' is off ')
                        elif l_isOn and  l_isCal:
                            self.m_group_on_off_cal[l_feed][l_chx]['cal_on'].append(l_subscan[l_feed][l_chx])
                            self.m_logger.info('feed ' + str(l_feed) + ' ' + l_chx + ' is cal_on ')
                        elif not l_isOn and  l_isCal:
                            self.m_group_on_off_cal[l_feed][l_chx]['cal_off'].append(l_subscan[l_feed][l_chx])                                            
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
        for l_feed in self.m_group_on_off_cal.keys():
            for ch in self.m_group_on_off_cal[l_feed].keys():
                l_group= self.m_group_on_off_cal[l_feed][ch]        
                l_calMarkTemp= l_group['on'][0]['frontend']['cal_mark_temp']                         
                #l_offAvg= sum(v['integrated_data']['spectrum'] for v in l_group['off']) / len(l_group['off'])
                #l_offAvg= np.mean( l_group['off'], axis= 0)
                " Avg off "
                l_offAvgData= []
                for el in l_group['off']:
                    l_offAvgData.append(el['integrated_data']['spectrum'])
                l_offAvg= np.mean(l_offAvgData, axis= 0)
                " Avg call on"
                l_CalOnAvg= None
                if len(l_group['cal_on']):
                    #l_calOnAvg= sum(v['integrated_data']['spectrum'] for v in l_group['cal_on']) / len(l_group['cal_on'])
                    l_CalOnAvgData= []
                    for el in l_group['cal_on']:
                        l_CalOnAvgData.append(el['integrated_data']['spectrum'])                    
                    l_calOnAvg= np.mean(l_CalOnAvgData, axis= 0)
                    l_calOnAvg = (l_calOnAvg - l_offAvg) / l_offAvg
                " Avg cal off "
                l_calOffAvg= None
                if len(l_group['cal_off']):
                    #l_calOffAvg= sum(v['integrated_data']['spectrum'] for v in l_group['cal_off']) / len(l_group['cal_off'])
                    l_CalOffAvgData= []
                    for el in l_group['cal_off']:
                        l_CalOffAvgData.append(el['integrated_data']['spectrum'])
                    l_calOffAvg= np.mean(l_CalOffAvgData, axis= 0)
                    l_calOffAvg = (l_calOffAvg - l_offAvg) / l_offAvg    
                                        
                " On - Off "                                
                for elOn in l_group['on']:
                    on= elOn['integrated_data']['spectrum']
                    on_off= (on - l_offAvg)/l_offAvg
                    elOn['integrated_data']['spectrum_on_off']= on_off
                " Rescale with cal mark temp "
                cal = np.array([l_calOnAvg, l_calOffAvg])                
                good = (cal != 0) & ~np.isnan(cal) & ~np.isinf(cal)
                cal = cal[good]                
                if len(cal) > 0:
                    meancal = np.median(cal) if len(cal) > 30 else np.mean(cal)    
                    calibration_factor = 1 / meancal * l_calMarkTemp
                    pdb.set_trace()
                    print("calibration_factor: " + str(calibration_factor))
                    print("meancal: " + str(meancal))
                    print("cal mark temp: " + str(l_calMarkTemp))
                else:   
                    return None, ""
                " Calibrated spectrum added to chx"                            
                for elOn in l_group['on']:                
                    elOn['integrated_data']['calibrated'] = elOn['integrated_data']['spectrum_on_off'] * \
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
            

    def ClassFitsAdaptations(self, p_wichGroup):
        """
        Generazione struttura dati secondo la definizione del classfist
                
        data in group on off cal sono divisi per
        self.m_group_on_off_cal['ch_0']['on'][0].keys()
        
        info base
        
        coordinate comandate in az, el o ra, dec
        coordinate osservate in crdelt2,3
        spettri separati per polarizzazione e per feed
        un file per ogni uno
        
        Paramters
        --------
        
        p_wichGroup: string
            'on', 'off', 'cal' ... scelta del gruppo con cui lavorare
            
        """

        " @todo inserire il campo cal is on ? quindi diversificare on ed off ?"      
        " @todo gestire i dati in caso di campo singolo spettro mediato o serie di spettri"        
        for l_feed in self.m_group_on_off_cal:        
            for l_chx in self.m_group_on_off_cal[l_feed]:              
                " single entry on classfits table "
                for l_ch in self.m_group_on_off_cal[l_feed][l_chx][p_wichGroup]:
                    #pdb.set_trace()
                    " Generic observation data copy to dedicated dict, more copies "
                    " below during calculations "
                    self.m_obs_general_data={}
                    self.m_obs_general_data['ra']= l_ch['scheduled']['ra'].to(unit.deg).value
                    self.m_obs_general_data['dec']= l_ch['scheduled']['dec'].to(unit.deg).value
                    self.m_obs_general_data['source']= l_ch['scheduled']['source']
                    self.m_obs_general_data['date-red']= Time.now().to_datetime().strftime('%d/%m/%y')
                    " classfits new dict filling process "
                    l_ch['classfits']={}
                    " ch by ch "
                    try:
                        " Lavoro con  i dati integrati "
                        " ut "                                                                           
                        l_tMjd= l_ch['integrated_data']['data_mjd'].mjd
                        l_ch['classfits']['UT']= ( l_tMjd - np.floor(l_tMjd)) * 86400
                        " date "                                            
                        l_ch['classfits']['DATE-OBS']= l_ch['integrated_data']['data_mjd'].strftime('%d/%m/%y') 
                        " lsts "
                        l_lsts= l_ch['integrated_data']['data_mjd'].sidereal_time('apparent', \
                                  fitslike_commons.Fitslike_commons.\
                                      get_site_location(l_ch['scheduled']['antenna']).lon)
                        l_lsts= l_lsts.value * unit.hr                              
                        " infos "                    
                        l_ch['classfits']['OBJECT']= l_ch['scheduled']['source']
                        l_ch['classfits']['LINE']= "F{}-{:3.3f}-MHz"\
                            .format(l_ch['frontend']['feed'], l_ch['backend']['bandwidth'])
                        self.m_obs_general_data['line']= l_ch['classfits']['LINE']
                        l_ch['classfits']['TELESCOP']= self.m_commons.class_telescope_name(l_ch)
                        l_mH2O= l_ch['integrated_data']['weather']
                        l_ch['classfits']['MH2O']= l_mH2O
                        " temp "
                        l_ch['classfits']['TSYS']= 1.0
                        l_ch['classfits']['CALTEMP']= l_ch['frontend']['cal_mark_temp'].value
                        " time "
                        l_ch['classfits']['LST'] = l_lsts.to('s').value     
                        l_ch['classfits']['OBSTIME']= l_ch['integrated_data']['data_integration']
                        "  "                    
                        l_ch['classfits']['CDELT1']= (l_ch['frontend']['bandwidth'] / 
                                                    l_ch['backend']['bins']).to('Hz')         
                        " freq and velocity "                    
                        l_ch['classfits']['RESTFREQ']= self.m_summary['summary']['restfreq'].to(unit.Hz).value                                                
                        self.m_obs_general_data['restfreq']= l_ch['classfits']['RESTFREQ']                              
                        l_ch['classfits']['VELOCITY']= l_ch['scheduled']['vlsr'].to("m/s").value                              
                        l_df= (l_ch['backend']['bandwidth'] / l_ch['backend']['bins']).to('Hz')
                        l_ch['classfits']['CDELT1']= l_df.value
                        self.m_obs_general_data['cdelt1']= l_ch['classfits']['CDELT1']
                        l_deltav= - l_df/ l_ch['classfits']['RESTFREQ'] * const.c
                        l_ch['classfits']['DELTAV']= l_deltav.value
                        " LOG test "
                        #self.m_logger.warn("RESTFREQ {}".format(l_ch['classfits']['RESTFREQ']))                                                
                        " Objects Coordinates "
                        l_ch['classfits']['CDELT2'] = l_ch['scheduled']['ra_offset'].to(unit.deg).value
                        l_ch['classfits']['CDELT3'] = l_ch['scheduled']['dec_offset'].to(unit.deg).value
                        l_ch['classfits']['AZIMUTH']= l_ch['integrated_data']['data_az'].to(unit.deg).value
                        l_ch['classfits']['ELEVATIO']= l_ch['integrated_data']['data_el'].to(unit.deg).value
                        l_ch['classfits']['CRVAL2']= l_ch['integrated_data']['data_ra'].to(unit.deg).value
                        l_ch['classfits']['CRVAL3']= l_ch['integrated_data']['data_dec'].to(unit.deg).value
                        " data "
                        l_ch['classfits']['OBSTIME'] = l_ch['integrated_data']['data_integration']    
                        l_ch['classfits']['MAXIS1'] = l_ch['backend']['bins']
                        self.m_obs_general_data['maxis1']= l_ch['classfits']['MAXIS1']
                        l_ch['classfits']['SPECTRUM_CAL']= l_ch['integrated_data']['calibrated']
                        l_ch['classfits']['SPECTRUM_ON']= l_ch['integrated_data']['spectrum']
                        l_ch['classfits']['SPECTRUM_ON_OFF']= l_ch['integrated_data']['spectrum_on_off']
                        l_ch['classfits']['CRPIX1']=  l_ch['backend']['bins'] // 2 + 1           
                    except Exception as e:
                        self.m_logger.error("Error preparing class data: " +str(e))
                    
            
    def classfitsWrite(self, p_group, p_on_what):
        """
        Scrittura file con calcolo header
        header prende i dati dai dati generici ricavati dalla scansione scansione      
        astropy fits works per column, i have to transpose all data while 
        generating new fits cols ( input data are per row basis
                                  
        Parameters
        ---------
        
        p_group: string
            which group, on ,off , cal_on, cal_off
        p_on_what: string
            wich kind of normalized data to use in case of 'on' group
            it assumes 'on', 'on_off', 'cal'        
                                          
        """        
        " clear - create destination folder "
        if not os.path.exists(self.m_outputPath):
            os.makedirs(self.m_outputPath)            
        " for every feed "
        for l_feed in self.m_group_on_off_cal:            
            l_outFileName= self.m_outputPath+ "feed_{}_{}_{}.fits".format(l_feed, p_group, p_on_what)
            self.m_logger.info("Preparing classfits file : " + l_outFileName)
            l_newCols=[]
            " for every column expressed in classfits definition.."
            for classCol in self.m_commons.getClassfitsColumnsZip():
                " [ name, form, unit ] column by column data building "                
                " fill one column looking into every feed[on], and builds column data "                     
                l_colData=[]                                
                l_columnFound = False                 
                " conditionals, some fields needs dedicated approach"
                for l_chx in self.m_group_on_off_cal[l_feed]:                    
                    for l_ch in self.m_group_on_off_cal[l_feed][l_chx][p_group]:
                        " converted fits data matches with classfits columns? "
                        " some columns need special care"
                        l_inferredCol= classCol[0]
                        "  we can choose between on on-off and calibrated spectrum for 'on' group "
                        if classCol[0] == "SPECTRUM":
                            l_inferredCol += '_'+ p_on_what.upper()                                
                        " "                                 
                        if l_inferredCol in l_ch['classfits'].keys():
                            " found match, add data to column data "
                            l_colData.append(l_ch['classfits'][l_inferredCol])
                            l_columnFound= True
                try:
                    " adding column to classfits if fitszilla representation matches it"                    
                    if l_columnFound:                        
                        " some fields needs dedicated approach"
                        if classCol[0] == "SPECTRUM":
                            l_newCols.append(fits.Column(name= classCol[0], format= "{}D".format(len(l_colData[0])),\
                                                         unit= classCol[2], array= l_colData))
                        else:                    
                            l_newCols.append(fits.Column(name= classCol[0], format= classCol[1],\
                                                unit= classCol[2], array= l_colData) )
                except Exception as e:
                    self.m_logger.error("classfits column creation exception: "+ str(e))
                    self.m_logger.error("column: " +str(classCol))
                    self.m_logger.error("column data: " + str(l_colData))
                    pdb.set_trace()
                                                                                
            l_hdData= self.m_obs_general_data
            " header "
            l_hdu= fits.PrimaryHDU() 
            try:
                l_hdu.header['CTYPE1']= "FREQ"
                l_hdu.header['CRVAL1']= 0
                l_hdu.header['CRVAL2']= l_hdData['ra']
                l_hdu.header['CRVAL3']= l_hdData['dec']
                l_hdu.header['OBJECT'] = l_hdData['source']
                l_hdu.header['SOURCE'] = l_hdData['source']
                l_hdu.header['DATE-RED'] = l_hdData['date-red']
                l_hdu.header['LINE'] = l_hdData['line']
                l_hdu.header['CDELT1'] = l_hdData['cdelt1']
                l_hdu.header['RESTFREQ'] = l_hdData['restfreq']
                l_hdu.header['MAXIS1'] = l_hdData['maxis1']
            except KeyError as e:
                self.m_logger.error("Exception filling " + l_outFileName + " header data: "+ str(e))                    
                " data "
            try:
                l_cdefs= fits.ColDefs(l_newCols)
                l_hdu= fits.BinTableHDU().from_columns(l_cdefs)                                        
            except Exception as e:                
                self.m_logger.error("Exception creating classfits model file")
                self.m_logger.error("classfits file: " + l_outFileName)
                self.m_logger.error(str(e))
                return
                #pdb.set_trace()
            try:
                if os.path.exists(l_outFileName):
                    os.remove(l_outFileName)
                l_hdu.writeto(l_outFileName)
            except Exception as e:                
                self.m_logger.error("Exception writings file")
                self.m_logger.error("classfits file: " + l_outFileName)
                self.m_logger.error(str(e))
