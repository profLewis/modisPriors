#!/usr/bin/env python

#Python version of prior2_adam_2.pro
#NSS 15/01/13

'''Python code to scan over a set of MODIS albedo product files and calculate a weighted mean and measure of variation. Stage 2.'''

__author__      = "P.Lewis, NSS"
__copyright__   = "Copyright 2013, UCL"
__license__ = "GPL"
__version__ = "1.0.3"
__maintainer__ = "P. Lewis"
__email__ = "p.lewis@ucl.ac.uk"
__status__ = "Production"

from optparse import OptionParser
import numpy as np
import netCDF4 as nc
import ast,os,glob,shutil
from subprocess import check_call
import logging

class prior2():
    """
    Stage 2 of MODIS climatology processing.

    In this stage, we gather all of the stage 1 processed files for all dates for a 
    given tile. These might, for instance, be in a directory ``results``.

    
    """
    def __init__(self,inputs):
        #constructor
        self.inputs = inputs

        try:
          from socket import gethostname, gethostbyname
          self.clientip = gethostbyname(gethostname()) 
        except:
          self.clientip = ''
        try:
          import getpass
          self.user = getpass.getuser()
        except:
          self.user = ''
        if os.path.exists(inputs.logdir) == 0:
          os.makedirs(inputs.logdir)

        self.log = inputs.logdir + '/' + inputs.logfile
        logging.basicConfig(filename=self.log,\
                     filemode='w+',level=logging.DEBUG,\
                     format='%(asctime)-15s %(clientip)s %(user)-8s %(message)s')
        self.d = {'clientip': self.clientip, 'user': self.user}
        logging.info("__version__ %s"%__version__,extra=self.d)
        print "logging to %s"%self.log
        # set an id for files
        try:
          import uuid
          self.id = unicode(uuid.uuid4())
        except:
          self.id = ""
        
        #inputs.ipdir = inputs.srcdir + '/' + inputs.opdir_
        #inputs.opdir = inputs.srcdir + '/' + inputs.opdir + '/' + inputs.opdir_

        inputs.ipdir = inputs.srcdir
        inputs.opdir = inputs.opdir
        try:
          if os.path.exists(inputs.srcdir) == 0:
            os.makedirs(inputs.srcdir)

          #create a symlink to the priors 1 directory
          #if os.path.exists(inputs.ipdir) == 0:
          #  os.symlink(inputs.prior1dir,inputs.ipdir)

          if os.path.exists(inputs.opdir) == 0:
            os.makedirs(inputs.opdir)
          if inputs.tmp == 'None':
            import tempfile
            f = tempfile.NamedTemporaryFile()
            inputs.tmp = f.name
            f.close()
            try:
              os.unlink(f.name)
            except:
              pass

          if os.path.exists(inputs.tmp) == 0:
            os.makedirs(inputs.tmp)
        except:
          logging.error("Error configuring directories: %s %s %s or %s "%(inputs.srcdir,inputs.ipdir,\
              inputs.opdir, inputs.tmp),extra=self.d)
        
       
    def getValidFiles(self,opdir,tile,version,type,tmp):
        """
        Obtain a list of valid files by looking in the directory ``opdir``
        for files of the pattern opdir + ``'/Kernels.*' + version + '.' + tile + '.' + type + '*.nc'``
        or possibly compressed versions (.gz suffix expected)

        """
        comps = glob.glob(opdir + '/Kernels.*' + version + '.' + tile + '.' + type + '.*.nc.gz')
        flats = glob.glob(opdir + '/Kernels.*' + version + '.' + tile + '.' + type + '.*.nc')
        if len(comps) == 0 and len(flats) == 0:
          comps = glob.glob(opdir + '/*/Kernels.*' + version + '.' + tile + '.' + type + '.*.nc.gz')
          flats = glob.glob(opdir + '/*/Kernels.*' + version + '.' + tile + '.' + type + '.*.nc')

        doy = []
        dataFiles = []

        for flat in flats:
          doyy = flat.split('/')[-1].split('.')[1]
          doy.append(doyy)
          dataFiles.append(flat)

        # scan over the compressed files
        for comp in comps:
          doyy = comp.split('/')[-1].split('.')[1]
          if doyy not in doy:
            # dont process if there is a flat file
            try:
              import os
              dataFile = '%s/%s'%(tmp,os.path.basename(comp.replace('.gz','')))
              if glob.glob(dataFile) != []:
                logging.info('uncompressed data '+dataFile+' exists',extra=self.d)
              else:
                import os
                logging.info('uncompressing %s'%dataFile,extra=self.d)
                cmd = "zcat %s > %s"%(comp,dataFile)
                os.system(cmd)
              doy.append(doyy)
              dataFiles.append(dataFile)
            except:
              try:
                logging.warning("file compression error %s"%dataFile,extra=self.d)
              except:
                logging.warning("file compression error: couldn't form filename",extra=self.d)

        # order the entries
        order = np.argsort(doy)
        doy = np.array(doy)[order]
        dataFiles = np.array(dataFiles)[order]
  
        return dataFiles,doy


    def getFileInfo(self,inimage):
        logging.info('testing image %s'%inimage,extra=self.d)
        ncfile = nc.Dataset(inimage,'r')
        ns = len(ncfile.dimensions['ns'])
        nl = len(ncfile.dimensions['nl'])
        nb = len(ncfile.variables.values())

        descrip = ncfile.description
        bnames = ncfile.variables.keys()
        dataType = ncfile.variables[bnames[0]].dtype
        
        logging.info('ns: %d nl: %d nb: %d'%(ns,nl,nb),extra=self.d)
        ncfile.close()
        
        return ns,nl,nb,descrip,dataType,bnames

    def readAll(self,b,ns,nl,nBands,bnames,dataFiles,weight):
        nDoys = len(dataFiles)
        # loop over bands
        n = np.zeros((nDoys,ns,nl),dtype=float)
        mean = np.zeros((nDoys,ns,nl),dtype=float)
        sd = np.zeros((nDoys,ns,nl),dtype=float)
        mask = np.zeros((nDoys,ns,nl),dtype=float)
        nSum = np.zeros((ns,nl),dtype=float)
        for i in range(nDoys):
          inimage = dataFiles[i]
          try:
            logging.info('opening %s'%inimage,extra=self.d)
            ncfile = nc.Dataset(inimage,'r')
            n_ = n[i,:,:] = (ncfile.variables['Weighted number of samples'])[:]
            land = (ncfile.variables['land mask'])[:]
            ok = True
          except:
            logging.error('error reading n and land mask ... %s'%inimage,extra=self.d)   
            ok = False
          try:
            mean[i,:,:] = (ncfile.variables[bnames[2*b]])[:]
            sd_ = (ncfile.variables[bnames[2*b+1]])[:]
            nSum += n[i,:,:] * weight[i]
          except:
            ok = False
            logging.error('error reading mean/sd ... %s'%inimage,extra=self.d)   
          try:
            ncfile.close()
          except:
            logging.error('error closing... %s'%inimage,extra=self.d)
     
          if ok: 
            # correct any apparently bum std values 
            w = np.where((land>0)&(sd_<=0)&(n_>0))
            sd_[w] = np.sqrt(self.minvar)

            # convert std to std err (see atbd)
            if self.inputs.stderr:
              w = np.where(n_>0)
              sd_[w] /= np.sqrt(n_[w])

            w1 = np.where((n_>0) & (sd_ <= np.sqrt(self.minvar)*1000))
            w = np.where(n_>0) 
            nChange = np.array(w1).shape[1]

            if self.inputs.normalise:
              # we shouldnt normally be here
              nsd = np.zeros((ns,nl),dtype=float)
              nsd[w] = np.sqrt(float(self.nYears)/n_[w])*sd_[w]
            else:
              nsd = sd_

            # set nominally low values to the median
            if len(w[0]) == 0:
              median = np.sqrt(self.maxvar)
            else:
              median = np.median(nsd[w])
              #median = np.median(nsd[w][nsd[w]>median])
            if np.isnan(median):
              median = 1.0
            if median == 0:
              median = np.sqrt(self.maxvar)
            if nChange > 0:
              logging.info("  %d %d Modifying %d low threshold std err to upper quartile (%.2f)"%(i,b,nChange,median),extra=self.d)
            nsd[w1] = median
            # deal with low sample number terms
            if self.inputs.nYears:
              w1 = np.where((n_<0.01*self.nYears)&(n_>0)&(nsd<median))
            else:
              w1 = np.where((n_<np.max(n_)*0.15)&(n_>0)&(nsd<median))
            nChange = np.array(w1).shape[1]
            if nChange > 0:
              logging.info("  %d %d Modifying %d low sample number entries to upper quartile (%.2f)"%(i,b,nChange,median),extra=self.d)
            nsd[w1] = median
            nsd[nsd>np.sqrt(self.maxvar)] = np.sqrt(self.maxvar)
            sd_ = nsd
            w = np.where(sd_>0)
            sd_[w] = weight[i]/(sd_[w]*sd_[w])
            sd[i,:,:] = sd_
            mask[i,:,:] = (sd_>0)*weight[i]
        return n,mean,sd,mask,nSum

    def weightIt(self,data,wt):
        [n,mean,sd,mask,nSum] = data
        xn = (n.T * wt).T
        xsd = (sd.T * wt).T
        xmask = (mask.T * wt).T
        xnSum = np.sum(xmask,axis=0)
        return xn,mean,xsd,xmask,xnSum

    def processAll(self):
        inputs = self.inputs
        #import pdb;pdb.set_trace()
        try:
          inputs.type = inputs.type[0]
        except:
          pass
        #import pdb;pdb.set_trace()
        # look in opdir for all files
        dataFiles, doy = self.getValidFiles(inputs.ipdir,inputs.tile,inputs.version,inputs.type,inputs.tmp)
        #import pdb;pdb.set_trace()    
        #get file info
        ns,nl,nb,descrip,dataType,bnames = self.getFileInfo(dataFiles[0])
        nBands = (nb - 2) / 2
        try:
          nYears = np.array(inputs.focus).sum()
          self.nYears = nYears
          logging.info("source data derived from %f years of samples"%nYears,extra=self.d)
          logging.info("  %s"%str(inputs.focus),extra=self.d)
        except:
          try:
            years = eval(descrip.split('years')[1].split('version')[0].replace(' [','[').replace('] ',']').replace(' ',','))
            nYears = len(years)
            logging.info("source data derived from %d years of samples"%nYears,extra=self.d)
            logging.info("  %s"%str(years),extra=self.d)
          except:
            nYears = 10
            logging.error("couldn't find years from data description",extra=self.d)
            logging.error("%s"%descrip,extra=self.d)
            logging.error("  setting to %d",extra=self.d)

          self.nYears = self.inputs.nYears or nYears
        self.minvar = 1e-20
        self.maxvar = 1.0

        doyflt = doy.astype(float)
        self.hasDoy = doy.copy()
        self.allDoys = np.array(['%03d'%i for i in np.arange(365/8+1)*8+1])
        doyflt = self.allDoys.astype(float)
        hasDoyflt = self.hasDoy.astype(float)
        # we have to process data for *all* doys 
        allMask = np.zeros((ns,nl)).astype(bool)
        nSum  = np.zeros((len(self.allDoys),ns,nl)).astype(float)

        ncfiles = []
        timer = 0
        ones = np.ones(len(self.allDoys))
        # loop over bands and store the info for all days for this band
        #import pdb;pdb.set_trace()
        for b in xrange(nBands):
          logging.info("... storing band %d\n"%(b),extra=self.d)
          # store data for all doys for this band
          n,meani,invVar,mask,nSumi = self.readAll(b,ns,nl,nBands,bnames,dataFiles,ones)

          # loop over days for processing
          for count,d in enumerate(range(len(self.allDoys))):
            logging.info("\n\n-------- processing band %d doy %d / %d %s -------\n"%(b,d,len(self.allDoys),self.allDoys[d]),extra=self.d)
            # sort the temporal filter, assuring correct wrap around
            dd = np.abs(doyflt[d] - hasDoyflt)
            ww = np.where(dd>=365/2)
            dd[ww] = 365-dd[ww]
            weight = np.exp(-dd/inputs.halflife)
        
            outimage = inputs.opdir + '/Kernels.'+ doy[d] + '.' + inputs.version + '.' + inputs.tile + '.background.' + inputs.type + '.nc'
            #names = bnames[0:nBands]
            #names2 = [n.replace('MEAN','SD') for n in names]
            #names[nBands:2*nBands] = names2
            #names.append('N samples')
            #names.append('Mask')

            # create the output files if necessary
            ok = True
            if b == 0:
              try:
                logging.info("opening output image %s"%outimage,extra=self.d)
                ncfile = nc.Dataset(outimage,'w', format ='NETCDF4')
                ncfile.createDimension('ns',ns)
                ncfile.createDimension('nl',nl)
                setattr(ncfile,'description',descrip)
                logging.info(" ... done",extra=self.d)
                ncfiles.append(ncfile)
              except:
                logging.error("... failed",extra=self.d)
                ok = False
            else:
              ncfile = ncfiles[count]

              # check for file errors
              if len(bnames) != len(np.unique(bnames)):
                logging.error("... failed: inconsistent band names",extra=self.d)
                logging.error("    %s"%bnames,extra=self.d)
                logging.error("    %s"%np.unique(bnames),extra=self.d)
                ok = False

            # process this band of data for this doy
            if ok:
              xn,xmeani,xinvVar,xmask,xnSumi = self.weightIt([n,meani,invVar,mask,nSumi],weight)
              # weight it
              
              invPostVar = np.sum(xinvVar,axis=0)
              sumN = np.sum(xmask,axis=0)
              w = np.where(invPostVar>0)
              # this is the posterior variance (when multiplied by the sum of the weights)
              num = np.sum(xinvVar*xmeani,axis=0)
              denom = invPostVar #np.sum(invVar,axis=0)
              mean = np.zeros_like(num)
              mean[w] = num[w]/denom[w]
              postVar = np.zeros_like(num)
              postVar[w] = 1./invPostVar[w]
              sd = np.sqrt(postVar)
              sd[sd>1.0] = 1.0
              # now store these
              allMask = (sumN>0) | allMask
              # write data
              try:
                logging.info("saving data %s %.4f"%(bnames[2*b],np.median(mean[w])),extra=self.d)
                logging.info("saving data %s %.4f"%(bnames[2*b+1],np.median(sd[w])),extra=self.d)
              except:
                pass
              data = ncfile.createVariable(bnames[2*b],'f8',('ns','nl'))
              data[:] = mean
              data = ncfile.createVariable(bnames[2*b+1],'f8',('ns','nl'))
              data[:] = sd
              nSum[count] += nSumi
              ncfile.sync()

            if b == nBands-1:
              logging.info("saving data %s %s"%('N samples','Mask'),extra=self.d)
              data = ncfile.createVariable('N samples','f8',('ns','nl'))
              data[:] = nSum[count]/float(nBands)
              data = ncfile.createVariable('Mask','f8',('ns','nl'))
              data[:] = (nSum[count]>0).astype(float)
              try:
                ncfile.close()
              except:
                pass
        if inputs.clean == True:
          try:
            shutil.rmtree('%s'%inputs.tmp)
          except:
            logging.warning("error attempting to clean tmp directory %s"%inputs.tmp,extra=self.d)

    
def processArgs(args=None,parser=None):

    usage = "usage: %prog [options]"
   
    args = args or sys.argv 
    #process input arguements
    parser = parser or OptionParser(usage=usage)
    prog = 'logger'

    parser.add_option('--srcdir',dest='srcdir',type='string',default='.',\
      help="the source directory, i.e. where the stage 1 files are (default: '.'")
    parser.add_option('--tile',dest='tile',type='string',default='h18v03',\
      help="the MODIS tile name (default: h18v03)")
    #parser.add_option('--prior1dir',dest='prior1dir',type='string',default='results',\
    #  help="directory containing stage 1 priors")
    parser.add_option('--clean',dest='clean',action='store_true',default=True,\
      help="clean any tmp files (default True)")
    parser.add_option('--noclean',dest='clean',action='store_false',\
      help="don't clean any tmp files")
    parser.add_option('--opdir',dest='opdir',type='string',default='.',\
      help="The output directory (default: processed)")
    parser.add_option('--halflife',dest='halflife',type='float',default=11.54,\
      help="half life decay of the double exponential function used in temporal weighting")
    parser.add_option('--version',dest='version',type='string',default='005',\
      help="MODIS data collection number as a string (default 005)")
    parser.add_option('--opdir_',dest='opdir_',type='string',default='processed',\
      help="The output directory (default: processed)")
    parser.add_option('--type',dest='type',type='string',default='SnowAndNoSnow',\
      help="The data type: SnowAndNoSnow, Snow or NoSnow (default: SnowAndNoSnow)")
    parser.add_option('--tmp',dest='tmp',type='string',default='None',\
      help="temporary directory (default: None -> unique directorty in /tmp)")
    parser.add_option('--logfile',dest='logfile',type='string',default='%s.log'%(prog),\
                      help="set log file name") 
    parser.add_option('--logdir',dest='logdir',type='string',default='logs',\
                      help="set log directory name")
    parser.add_option('--normalise','--normalize',dest='normalise',action='store_true',default=False,\
      help="Normalise the uncertainty information (default: False)")
    parser.add_option('--nonormalise','--nonormalize',dest='normalise',action='store_false',\
      help="Don't normalise the uncertainty information")
    parser.add_option('--sdscale',dest='nYears',type='float',default=None,\
                      help="set the scale for sd normalisation. If this is not set"+\
                           " the number of years used in processing stage 1"+\
                           " is used as an expectation of the number of samples" +\
                           "NB only relevant if --normalise used")
    parser.add_option('--stderr',dest='stderr',action='store_true',default=True,\
      help="Convert std dev to std err (default: True)")
    parser.add_option('--nostderr',dest='stderr',action='store_false',\
      help="Don't convert std dev to std err")


    #parser.add_option('--hack',dest='hack',action='store_true',default=False,\
    #  help="Hack for incorrect file info (default: False)")
    #parser.add_option('--nohack',dest='hack',action='store_false',\
    #  help="Don't use hack for incorrect file info")


    return parser.parse_args(args or sys.argv)

if __name__ == '__main__':
    
    opts, args = processArgs()

    prior2 = prior2(opts)
    prior2.processAll()
    
