# -*- coding: utf-8 -*-
"""Fits like parser

    Fits format awarness, it handles fits data toward fits like representation
"""
import astropy.io.fits as fits
import fitslike_keywords
import fitslike_commons
import pdb
import logging


class Awarness_fitszilla():
    """fitszilla data parser"""

    def __init__(self, l_fitszilla):
        """
        Store fitzilla file

        Parameters
        ----------
        l_fitszilla : astropy fits
        Fitszilla to be handled before it gets represented as fitslike

        Returns
        -------
        None.

        """
        self.m_commons = fitslike_commons.Fitslike_commons() 
        self.m_jsonRepr = fitslike_keywords.Keyword_json('fitszilla')
        self.m_components = self.m_jsonRepr.fitslike_components()        
        self.m_parsingDicts = {}
        self.m_logger = logging.getLogger(self.m_commons.logger_name())
        self.m_logger.info("Fitszilla input file dedocding")
        for l_component in self.m_components:
            self.m_logger.info("fitszilla registering component :  %s ",
                               l_component)
            self.m_parsingDicts[l_component] = \
                self.m_jsonRepr.parser(l_component)        
        self.m_fitszilla = l_fitszilla

    def parse(self):
        """Fitszilla parsing

        Parse fitzilla input file following parsing dictionaried
        It only extract keywords from fitzilla without processing them.

        Returns
        -------
        A dictionary = {
            fitslike_key : fitszilla_value
            ... : ...
        }
        """
        self.m_intermediate = {}
        for l_component in self.m_parsingDicts.keys():            
            for l_key in self.m_parsingDicts[l_component].keys():
                self.m_intermediate[l_key] = None
                l_keyParserEntry = self.m_parsingDicts[l_component][l_key]
                # Decoding entry
                l_tableName, l_headerOrData, l_inputKeyword  =\
                    fitslike_keywords.Keyword_helper.\
                    decode_entries(l_keyParserEntry)
                if l_tableName is None:
                    self.m_logger.error("empty table name %s, %s, %s",\
                                 l_tableName, l_headerOrData, l_inputKeyword)
                    return
                # todo semplificare la gestione tabella
                l_table = fitslike_keywords.Keyword_helper.\
                    get_table_index(self.m_jsonRepr.input_tables_dict(),
                                     l_tableName)
                if l_table is None:
                    self.m_logger.error("table not found %s", l_tableName)
                    return
                l_fitsTable = self.m_fitszilla[l_table]
                # is in header or in data ?
                l_headerOrData = fitslike_keywords.Keyword_helper.\
                    is_header_or_data(l_headerOrData)
                # looking for keyword or data
                l_fitsTableContent = None
                if l_headerOrData is True:
                    l_fitsTableContent = l_fitsTable.header
                else:
                    l_fitsTableContent = l_fitsTable.data     
                # fitszilla keyword look up
                try:
                    self.m_intermediate[l_key] = l_fitsTableContent[l_inputKeyword]
                except:
                    self.m_logger.error("Missing [table, keyword] %s: %s",
                                        l_tableName, l_inputKeyword)                                    
        return self.m_intermediate
    
    def process(self):
        """
        Parsed keyword processing
        
        This class knows how to prepare data to fit fitslike.
        Those data are store into a processed representation.
        Processing order is mandatory, coordinates should be processed after 
        spectrum.
        

        Returns
        -------
        Processed data representation.

        """
        self.m_processedRepr = {}
        self._process_spectrum()
        self._process_coordinates()
        return self.m_processedRepr
        
        
    def _process_spectrum(self):
        """
        Spectrum keyword processing.
        Fills spectrum final representation
        
        frontend[backend_id]{
            "fe_feeds"
			"fe_if"
			"fe_polarizations",
			"fe_be_id"
			"fe_frequency"
			"fe_bandwith"
			"fe_local_oscillator"
			"fe_cal_mark_temp"
        }
        
        backend[frontend_id]{            
            "be_id"
			"be_bins"
			"be_sample_rate"
			"be_bandwith"
			"be_frequency"
			"be_integration"
			"be_data_type"
        }    
        
        ch_X: {
            frontend: frontend[be_id]
            backend: backend[fe_id]
            spectrum
        }
        
        Returns
        -------
        None.

        """
        # zip front end
        l_frontEnds= {}
        l_zipFrontEnds = zip(self.m_intermediate['fe_be_id'],
                    self.m_intermediate['fe_feeds'],
                    self.m_intermediate['fe_if'],
                    self.m_intermediate['fe_polarizations'],
                    self.m_intermediate['fe_frequency'],
                    self.m_intermediate['fe_bandwidth'],
                    self.m_intermediate['fe_local_oscillator'],
                    self.m_intermediate['fe_cal_mark_temp'])
        # create dict[backend_id]= front end
        for l_zipFe in l_zipFrontEnds:
            l_frontEnds[l_zipFe[0]]= l_zipFe      
        #  zip backend
        l_backEnds= {}
        l_zipBackend= zip(self.m_intermediate['be_id'],
                    self.m_intermediate['be_bins'],
                    self.m_intermediate['be_sample_rate'],
                    self.m_intermediate['be_bandwidth'],
                    self.m_intermediate['be_frequency'],                    
                    self.m_intermediate['be_data_type'])        
        # create dict[backend_id]= back end
        for l_zipBe in l_zipBackend:
            l_backEnds[l_zipBe[0]]= l_zipBe
        #pdb.set_trace()
        # Creates chX_feed_pol: frontend, backend, spectrum            
        for l_elBe in l_backEnds.keys():            
            l_innerDict= {}
            l_innerDict['backend']= l_backEnds[l_elBe]
            l_innerDict['frontend']= l_frontEnds[l_elBe]
            l_innerDict['spectrum']= self.m_intermediate['ch'+str(l_elBe)]
            self.m_processedRepr['ch_'+str(l_elBe)] = l_innerDict.copy()

    def _process_coordinates(self):
        """
        Coordinate data processing
        
        Every data table entry fitszilla coordinates refers to central feed.

        feed offsets:
            "fe_x_offset"
			"fe_y_offset"
            
        coordinates : [commons to every data table entry]           	
			"data_ra"
			"data_dec"
			"data_az"
			"data_el"
			"data_par_angle"
			"data_derot_angle"
        
            - Converts radians to arcsec ?
            - Apply offset and rotation angle to every feed in az, el
            - Convert final az, el in ra, dec
            
        Replicate coordinates for every ch_x (feed only)
        Apply feed offset (az, el)
        Apply derot_angle to every feed (az, el)
        Infer ra, dec for every feed
            
        Returns
        -------
        None.

        """
        l_coordinatesDict= {
            'data_az': self.m_intermediate['data_az'],
            'data_el': self.m_intermediate['data_el'],
            'data_derot_angle': self.m_intermediate['data_derot_angle']
            }
        pdb.set_trace()
        # reduce feeds removing duplicate (left/right)
        l_feedCoordinatesDict= dict.fromkeys(self.m_intermediate['fe_feeds'])
        # copy usefull coord. for every feed
        for l_feeds in l_feedCoordinatesDict.keys():
            l_feedCoordinatesDict[l_feeds]= l_coordinatesDict.copy()
        # Apply offset + rotation
        "todo"
        # Calculate ra, dec
        "todo"
            
        