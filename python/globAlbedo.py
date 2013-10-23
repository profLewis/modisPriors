#!/bin/env python

__author__      = "P.Lewis"
__copyright__   = "Copyright 2013, UCL"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "P. Lewis"
__email__ = "p.lewis@ucl.ac.uk"
__status__ = "Production"

from optparse import OptionParser
import numpy as np
import netCDF4 as nc
import ast,os,glob,shutil,sys
from subprocess import check_call
import logging


'''
Script to run globAlbedo approach for processing MODIS climatology
'''

genPnw = '''
Generate data pages
---------------------

<<results = "rst">>=
from genPnw import *
genIndex('results',tile='TILE',stage=STAGE)
@

'''



class globAlbedo():
  '''
  Class for globAlbedo processing utilities
  '''
  def __init__(self,srcdir,stage=None,tile=None,sdims=None,\
               type=None,doyList=None,years=None,focus=None,shrink=1,scale=None,bands=None):
    '''
    Class initialisation.

    Parameters
    ===========

    srcdir : string
             directory containing MCD43 files (either in this directory or
             one level below this). E.g. /data/geospatial_11/plewis/$tile/data
             then the data could be in /data/geospatial_11/plewis/TILE/data 
             or /data/geospatial_11/plewis/TILE/data/*
             The string can contain the patterns TILE and/or STAGE, which
             are then replaced by the actual tile and stage strings

    Optional:

    stage : string 
            '1' or '2' depending on which stage you want to run
            If not specified (default None) then this is pulled from 
            the current directory name, assuming it is
            of the form e.g. stage1_h18v04   

    tile : string 
            Tile identifier, e.g. 'h18v04'
            If not specified (default None) then this is pulled from 
            the current directory name, assuming it is
            of the form e.g. stage1_h18v04   

    sdims : array e.g. [-1,-1,-1,-1] (which is the default, meaning everything)
            which specifies [l0,nl,s0,ns]

    years : array e.g. [2000,2001,2002,2003,2004] defining
            which years data to process

    focus : array e.g. [0.1,0.5,1.0,0.5,0.1] of weights to apply to the data in stage 2 processing

    shrink : spatial shrink factor (default 1)

    scale : array of weighting factors to assign to the MODIS MCD43 QA values
            default 0.618**np.arange(4)

    bands : array which bands to extract from MCD43
            default [7,8,9] which are the broadband values

    doyList : list of doys to process
    '''
    from datetime import date

    self.here = os.path.abspath(os.curdir)
    # try to find the tile and stage from this
    self.thisdir = self.here.split(os.sep)[-1]
    # assuming its of the form stage1_h18v04
    self.stage = stage or self.thisdir.split('_')[0].replace('stage','')
    self.tile  = tile  or self.thisdir.split('_')[1]
    print 'processing in',srcdir,'for tile',self.tile,'stage',self.stage

    self.srcdir = srcdir.replace('TILE',self.tile).replace('STAGE',self.stage)
    self.sdims = sdims or [-1,-1,-1,-1]
    self.years = years or np.arange(2000,date.today().year+1)
    self.focus = focus or np.ones_like(self.years).astype(float)
    self.shrink = shrink
    self.scale = scale or 0.618**np.arange(4)
    self.bands = bands or [7,8,9]
    self.doyList = doyList or ['%03s'%i for i in xrange(1,366,8)]
    self.type = type or ['Snow', 'NoSnow', 'SnowAndNoSnow']

    print 'processing',type
    try:
      # list the srcdir to ensure automounting
      cmd = 'ls %s > /dev/null'%self.srcdir
      os.system(cmd)
    except:
      sys.stderr('Failed to list contents of source directory %s\n'%self.srcdir)
      return

    self.whereFrom = os.path.abspath(os.path.dirname(__file__))
    # make sure whereever this file comes from is in pythonpath
    os.sys.path.insert(0,self.whereFrom)

  def stage1(self):
    '''
    Run stage 1 processing
    '''
    from albedo_pix import processArgs,insensitive_glob,albedo_pix

    # -clean 
    args = [sys.argv[0],'--logfile=%s.log'%self.tile,'--dontclean',\
            '--shrink=%d'%int(self.shrink),\
            '--tile=%s'%self.tile,'--srcdir=%s'%self.srcdir] 

    # pass through any from the cmd line, which can override these
    [args.append(i) for i in sys.argv[1:]]
    # use the cmd line parser from albedo_pix 
    try:
      args.remove('--pylab')
    except: pass
    try:
      args.remove('--pylab=inline')
    except: pass
    try:
      args.remove('notebook')
    except: pass
    opts, args = processArgs(args=args)
    opts.backupscale = self.scale
    opts.sdims       = self.sdims
    opts.bands       = self.bands
    opts.years       = self.years
    opts.focus       = self.focus
    opts.type        = np.atleast_1d(self.type)
    opts.dontwithsnow = opts.dontnosnow = opts.dontsnow = True
    if 'SnowAndNoSnow' in opts.type:
      opts.dontwithsnow = False
    if 'Snow' in opts.type:
      opts.dontsnow = False
    if 'NoSnow' in opts.type:
      opts.dontnosnow = False

    print opts.type

    opts.doyList = np.array(self.doyList)
    #np.array(ast.literal_eval(opts.doyList))

    #class call
    self.albedo = albedo_pix(opts)

    self.albedo.runAll()

  def snowtype(self,this):
    '''
    Translation for snow type representation
    '''
    if this == '--withsnow':
      return '--type=SnowAndNoSnow'
    elif this == '--snow':
      return '--type=Snow'
    elif this == '--nosnow':
      return '--type=NoSnow'
    else:
      return this

  def stage2(self):
    '''
    Run stage 2 processing
    '''
    from prior2Fast import prior2,processArgs

    #sys.argv = [self.snowtype(i) for i in sys.argv]

    args = [sys.argv[0],'--tmp=tmp','--opdir=results','--srcdir=results','--clean','--tile=%s'%self.tile,\
            '--logfile=%s.log'%self.tile]
    [args.append(i) for i in sys.argv[1:]]
    # pass through any from the cmd line, which can override these
    [args.append(i) for i in sys.argv[1:]]
    # use the cmd line parser from albedo_pix 
    try:
      args.remove('--pylab')
    except: pass
    try:
      args.remove('--pylab=inline')
    except: pass
    try:
      args.remove('notebook')
    except: pass
    opts, args = processArgs(args=args)
    opts.backupscale = self.scale
    opts.sdims       = self.sdims
    opts.bands       = self.bands
    opts.years       = self.years
    opts.focus       = self.focus
    opts.type        = np.atleast_1d(self.type)
    opts.dontwithsnow = opts.dontnosnow = opts.dontsnow = True
    if 'SnowAndNoSnow' in opts.type:
      opts.dontwithsnow = False
    if 'Snow' in opts.type:
      opts.dontsnow = False
    if 'NoSnow' in opts.type:
      opts.dontnosnow = False

    print opts.type

    self.prior2 = prior2(opts)
    self.prior2.processAll()


  def prep(self):
    '''
    Create correct links and directory structure
    '''
    # cheat -- do this in unix for the moment
    cmd = '(mkdir -p source;mkdir -p build/html;pushd build/html;rm -f python logs results;ln -s ../../logs logs;ln -s ../../python python;ln -s ../../results results;popd) >& /dev/null'
    try:
      os.system(cmd)
    except:
      sys.stderr.write('error creating sane directory structure in %s'%self.here)
      return
    # sort path in conf.py
    this = open(self.whereFrom + os.sep + 'conf.py').readlines()
    this = [i.replace('theCode',self.whereFrom) for i in this]
    open('source' + os.sep + 'conf.py','w').writelines(this)
    # reporting using Pweave
    cmd = 'sphinx-apidoc -f -o source %s'%self.whereFrom
    try:
      os.system(cmd)
    except:
      sys.stderr.write('error creating sane API files in %s'%self.here)

    # create source/genPnw.Pnw
    try:
      f = open('source/genPnw.Pnw','w')
    except:
      sys.stderr.write('error opening genPnw.Pnw for writing in %s'%self.here)
      return
    # copy over genPnw.py
    cmd = 'cp files/python/genPnw.py source'
    try:
      os.system(cmd)
    except:
      sys.stderr.write('error creating sane API files in %s: cant find genPnw.py in bin as expected'%self.here)

    try:
      f.write(genPnw.replace('TILE',self.tile).replace('STAGE',self.stage))
      f.close()
    except:
      sys.stderr.write('error writing to genPnw.Pnw in %s'%self.here)
      return

  def report(self):
    '''
    run reporting modules
    '''
    import pweave
    
    os.chdir('source')
    try:
      pweave.pweave('genPnw.Pnw')
    except:
      sys.stderr.write('error creating data files in %s'%self.here + '/source')
      return
    import glob
    files = glob.glob('K*Pnw')
    for file in files:
      print file
      try:
        pweave.pweave(file)
      except:
        sys.stderr.write('error Pweaving data file %s in %s'%(file,self.here + '/source'))
    os.chdir(self.here)
    cmd = 'make clean html'
    try:
      os.system(cmd)
    except:
      sys.stderr.write('error making html')
    return


from albedo_pix import processArgs

if __name__ == '__main__':
  srcdir = './files'

  sdims=None
  years=None
  focus=None
  shrink=1
  scale=None
  bands=range(0,7)
  type=['NoSnow']
  stage = '1'
  tile = 'h18v03'


  isRun = False
  isReport = False

  if '--runStage' in sys.argv:
    isRun = True
    sys.argv.remove('--runStage')

  if '--report' in sys.argv:
    isReport = True
    sys.argv.remove('--report')

  
  opts, args = processArgs()

  tile = opts.tile
  opts.stage = opts.stage

  # set up this utility passing through any obvious options
  self = globAlbedo(srcdir,sdims=sdims,years=years,\
             stage=stage,type=type,tile=tile,\
             focus=focus,shrink=shrink,scale=scale,bands=bands)

  # prep
  self.prep()

  # run the stage or just report?
  if isRun:
    if self.stage == '1':
      self.stage1()
    elif self.stage == '2':
      self.stage2()

  if isReport:
    self.report()  
  

