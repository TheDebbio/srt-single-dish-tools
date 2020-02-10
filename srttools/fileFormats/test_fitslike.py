# -*- coding: utf-8 -*-
"""Fitslike objects testing
"""

import fitslike_commons
import fitslike_obsData
import fitslike
import awarness_fitszilla
import pytest
from astropy.io import fits
import pdb

nodding_dir = '/home/debbio/discos/Scan/SARDARA_Nodding/20181210-232201-24-18-W3OH/'
nodding_zilla_1_8 = '20181210-232307-24-18-W3OH_001_008.fits'


class TestFitsLike_arch():
    """Fitslike architecture test class"""

    @staticmethod
    def test_instance():
        """Instancing fitslike object"""
        l_fits = fitslike.Fitslike()
        l_fits.dump()

    @staticmethod
    def test_Awareness_fitszilla_parse(p_file):
        """Parse a fitszilla through Awarness_fitszilla"""
        l_fits = fits.open(p_file)        
        l_aware = awarness_fitszilla.Awarness_fitszilla(l_fits)
        #pdb.set_trace()
        l_intermediateDict = l_aware.parse()
        print(l_intermediateDict)


if __name__ == "__main__":
    l_unitArch = TestFitsLike_arch()
    l_unitArch.test_instance()
    #
    l_file = nodding_dir + nodding_zilla_1_8
    l_unitArch.test_Awareness_fitszilla_parse(l_file)

