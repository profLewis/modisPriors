#!/usr/bin/env python

# Python version of albedo2.pro
# NSS 23/11/12
# Lewis 11 June 2013: major modifications

'''Python code to scan over a set of MODIS albedo product files and calculate a weighted mean and measure of variation. Stage 1.'''

__author__      = "P.Lewis, NSS"
__copyright__   = "Copyright 2013, UCL"
__license__ = "GPL"
__version__ = "1.0.7"
__maintainer__ = "P. Lewis"
__email__ = "p.lewis@ucl.ac.uk"
__status__ = "Production"

import numpy as np
import sys,os,datetime,math,ast,glob,resource
from pyhdf import SD
#import gdal

from optparse import OptionParser
#the two imports below are needed if you want to save the output in ENVI format, but they conflict with the netcdf code, so commented out.
#from osgeo import gdal
#import osgeo.gdalconst as gdc
from subprocess import check_call
import netCDF4 as nc
import logging
from scipy import stats


class modis_climatology_utils(object):
  '''
  Modis climatology utilities
  '''
  def __init__(self):
    # make the op directory
    if os.path.exists(self.ip.opdir) == 0:
            os.makedirs(self.ip.opdir)

  def insensitive_glob(self,pattern):
    """ From: http://stackoverflow.com/\
	questions/8151300/ignore-case-in-glob-on-linux

    form a case insensitive glob version of a filename (or other) pattern

    Parameters
    ----------
    pattern : a string (filename)
          String to process

    Returns
    -------
    npattern : a string (filename) that is case insensitive

    Examples
    --------

    insenitive = insensitive_glob('/data/fbloggs/*/*.hdf')

    >>> insensitive_glob('/data/*/plewis/h25v06/DATA/*.HDF')[0]
    '/data/geospatial_11/plewis/h25v06/data/\
	MCD43A1.A2001001.h25v06.005.2006360202416.hdf'

    """
    def either(c):
        return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(self.either,pattern)))



  def processArgs(self,args=None,parser=None):

    usage = "usage: %prog [options]"

    #process input arguements
    parser = parser or OptionParser(usage=usage)
    prog = %prog

    # options
    parser.add_option('--logfile',dest='logfile',type='string',\
			default='%s.log'%(prog),\
                      	help="set log file name")
    parser.add_option('--logdir',dest='logdir',type='string',\
			default='logs',\
                      	help="set log directory name")
    parser.add_option('--srcdir',dest='srcdir',type='string',\
			default='modis',\
                      	help="Source (MODIS MCD43) data directory")
    parser.add_option('--tile',dest='tile',type='string',\
			default='h18v03',\
                      	help="MODIS tile ID")
    parser.add_option('--backupscale',dest='backupscale',type='string',\
			default='[1.,0.7,0.49,0.343]',\
                      	help="Array defining the scale \
			to map MODIS QA flags to, e.g. [1.,0.7,0.49,0.343]")
    parser.add_option('--opdir',dest='opdir',type='string',\
			default='results',\
                      	help="Output directory")
    parser.add_option('--compression',dest='compression',\
			action='store_true',default=True,\
                      	help='Compress output file')
    parser.add_option('--type',dest='type',default=["NoSnow"],\
                      	help='Specify data types to process: \
			list of e.g. ["NoSnow", "Snow", "All"]')
    parser.add_option('--shrink',dest='shrink',type='int',default=1,\
                      help="Spatial shrink factor (integer: default 1)")
    parser.add_option('--sdims',dest='sdims',type='string',\
			default='[-1,-1,-1,-1]',\
                      	help='image subsection: \
			default [-1,-1,-1,-1]i or [l0,nl,s0,ns]')
    parser.add_option('--bands',dest='bands',type='string',\
			default='[7,8,9]',\
			help='list of bands to process. Default [7,8,9]')

  
