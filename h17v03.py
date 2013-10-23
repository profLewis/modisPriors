# put the local directory files/python in path
import sys
sys.path.append('files/python')
import pylab as plt
import numpy as np

# try to import the modules we will need
# if any of these fail, you need to install them
# e.g. 
# easy_install --user pyhdf

import pyhdf
import netCDF4 
import cloud_sptheme
import optparse 
import subprocess 
import genPnw
import pweave
import ast,os,glob,shutil,sys
import logging
import tempfile
import os

# local utilities that should be in files/python
import globAlbedo
import plotter
import albedo_pix
import prior2Fast

# set up for a single tile:

srcdir = './files'
sdims  = None
years  = None
focus  = None
shrink = 2
scale  = None
bands  = range(0,7)
tile   = 'h17v03'
stage  = '1'
atype  = ['NoSnow']

from globAlbedo import globAlbedo
self = globAlbedo(srcdir,sdims=sdims,years=years,\
                  tile=tile,focus=focus,\
                  stage=stage,type=atype,\
                  shrink=shrink,scale=scale,bands=bands)

self.prep()

self.stage1()

