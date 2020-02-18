#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 09:31:06 2020

@author: debbio
"""


"""Calibration details

Normalizzazione :
    (ON-OFF)/OFF

Calibrazione :
    CAL + OFF
    CAL + ON
    dipende dai dati acquisiti
    La marca conferisce la scala fisica al dato normalizzato    

Dal codice di Bachetti:
    - Prende la scan già orientata al classfits
    - ordina i dati in base 'mjd', 'telescope', 'line'
    - fa un sottogruppo 'telescope', 'line'
    - ricerca la ripetizione di un pattern ( find_cycle ) 
      sulle keyword 'SIGNAL', 'CAL_IS_ON' (?)
    - crea un gruppo con le ntries del pattern trovato
    - normalizza on - off e cal
        - on-source - off-source (off più vicina in posizione)
        - signal/on-source o reference/off-source ?
        - se abbiamo cal su off (OFFCAL) => off/(offcal-off)
        - se abbiamo cal su on (ONCAL) => on/(oncal-on)
    

Steps:

Keywords:
    'SIGNAL': preso dall'header di fitszilla con valori:
        - SIGNAL (ON)
        - REFERENCE (OFF)
        - REFCAL (OFF + CAL)
        - REFSIG (ON + CAL)
        - NONE (OFF)
        
    Nel multifeed, fitszilla dovrebbe in futuro splittare file per ogni feed
    quindi ogni table ha la sua info a riguardo!?
        

"""