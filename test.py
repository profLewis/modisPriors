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
import optparse 
import subprocess 
import ast,os,glob,shutil,sys
import logging
import tempfile
import os
# local utilities that should be in files/python
import globAlbedo
import albedo_pix

# set up for a single tile:

srcdir = './files'
sdims  = None
years  = None
focus  = None
shrink = 2
scale  = None
bands  = [0]
tile   = 'h17v03'
stage  = '1'
atype  = ['NoSnow']

import pdb;pdb.set_trace()

self = globAlbedo.globAlbedo(srcdir,sdims=sdims,years=years,\
                  tile=tile,focus=focus,\
                  stage=stage,type=atype,\
                  shrink=shrink,scale=scale,bands=bands)


