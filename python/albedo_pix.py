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
from optparse import OptionParser
#the two imports below are needed if you want to save the output in ENVI format, but they conflict with the netcdf code, so commented out.
#from osgeo import gdal
#import osgeo.gdalconst as gdc
from subprocess import check_call
import netCDF4 as nc
import logging
from scipy import stats

def insensitive_glob(pattern):
    """ From: http://stackoverflow.com/questions/8151300/ignore-case-in-glob-on-linux

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
    '/data/geospatial_11/plewis/h25v06/data/MCD43A1.A2001001.h25v06.005.2006360202416.hdf'

    """
    def either(c):
        return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
    return glob.glob(''.join(map(either,pattern)))

class albedo_pix():
    """
    Calculate mean and variation in albedo for albedo climatology.

    Initialise as:

    ::

        al = albedo_pix(inputs)


    where ``inputs`` should contain (e.g. ``inputs.clean = True``):

        ``doyList`` : list of strings identifying which days to process
                      e.g. ['001','009'] or set to None to select all
                      available
        ``clean`` : boolean flag to determine whether previously calculated netCDF files
                       should be read rather than re-calculating from the MODIS hdf files.
        ``rst`` : ReST directory name (e.g. `source`). 
                     Set to the string 'None' for no ReST output.
        ``logfile`` : log file name (string)

        ``logdir`` : directory for log file (string)

        ``srcdir`` : Source (MODIS MCD43) data directory

        ``tile``   : MODIS tile ID

        ``backupscale`` : Array defining the scale to map MODIS QA flags to, e.g. [1.,0.7,0.49,0.343]

        ``opdir`` : Output directory

        ``compression`` : Compress output file (boolean flag)

        ``shrink`` : spatial shrink factor (integer) 

        ``sdims`` : image subsection: default [-1,-1,-1,-1] or [l0,nl,s0,ns]

        ``bands`` : list of bands to process. e.g. [7,8,9]

        ``years`` : list of years to process e.g. 
                [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013]

        ``focus`` : list or weights for years to process (defaults to 1.0)

        ``version`` : MODIS collection number (as string). e.g. 005

        ``product`` : product name e.g. mcd43a (case insensitive)

        ``part`` : How many parts to split the tile into e.g. 1

        ``dontsnow`` : don't output snow data

        ``dontnosnow`` : don't output non-snow data

        ``dontwithsnow`` : don't output snow and non-snow (combined) data

        ``rstbase``  : base name for ReST files

    """    
    def __init__(self,inputs):
        """
        Class constructor


        Parameters
        ----------
        inputs : structure (e.g. from parser) containing settings

        """

        #constructor
        self.ip = inputs
        
        self.a1Files = None
        self.a2Files = None
        self.weighting = None

        # threshold sd to 1.0 
        self.maxvar = 1.0
        self.minvar = 1e-20
       
        # logging
        if os.path.exists(self.ip.logdir) == 0:
          os.makedirs(self.ip.logdir)
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
        self.log = self.ip.logdir + '/' + self.ip.logfile
        logging.basicConfig(filename=self.log,\
                     filemode='w+',level=logging.DEBUG,\
                     format='%(asctime)-15s %(clientip)s %(user)-8s %(message)s')
        self.d = {'clientip': self.clientip, 'user': self.user}
        logging.info("__version__ %s"%__version__,extra=self.d)
        if len(inputs.focus) != len(inputs.years):
          logging.warning('warning: inconsistent number of years defined for --years and --focus',extra=self.d)
          logging.warning('setting focus to 1.0 for all years',extra=self.d)
          inputs.focus = np.ones_like(inputs.years.astype(float))
        print "logging to %s"%self.log
        # set an id for files
        try:
          import uuid
          self.id = unicode(uuid.uuid4())
        except:
          self.id = ""
        # create error file directory
        self.errFiles = './logs'
        d = os.path.dirname(self.errFiles)
        if not os.path.exists(d):
          os.makedirs(d)

    def getDates(self,list1):
        """
        Process a list of date/year entries into a list of days of year (doyList) 
        and years (yearlist)
        
        Parameters
        ----------

        list1 : list of MODIS hdf file names from which we will extract the date information
        
        Returns
        -------

        doyList :  unique, sorted list of DOY 
        yearlist : unique, sorted list of years

        """
        doyList=[]
        yearList=[]
        
        for file in list1:
            fileBaseName = os.path.basename(file)
            year=fileBaseName[9:13]
            doy=fileBaseName[13:16]
            if doyList.count(doy) == 0:
                doyList.append(doy)
            if yearList.count(year) == 0:
                yearList.append(year)

        yearList = np.sort(np.unique(yearList))
        doyList = np.sort(np.unique(doyList))
        logging.info(str(yearList),extra=self.d)
        logging.info(str(doyList),extra=self.d)

        return doyList, yearList

    def idSamples(self,nb):
        """
        Return a string array describing the data storage for the variance/covariance structures

        Parameters
        ----------

        nb  : long
              number of bands

        Returns
        -------

        retval : string array 
                 containing output band text descriptions


        """
        params = ['F0','F1','F2']
        this= []
        for i in range(nb):
            for j in range(3):
              this.append('MEAN: BAND ' + str(i) + ' PARAMETER ' + params[j])
              this.append('SD: BAND ' + str(i) + ' PARAMETER ' + params[j])

        return this



    def getValidFiles(self,yearList,doy):

        """       
        Return a list of files that can be processed 


        Parameters
        ----------
        
        yearlist  : string array  
                    years to process

        doy       : string       
                    which day to process 


        Returns
        -------

        dictionary : {A1FILES:string array,A2FILES:string array}

        where:

        A1FILES : data files
        A2FILES : QA files

        """
        if(self.a1Files==None):
            product = self.ip.product
            tile = self.ip.tile
            version = self.ip.version
            srcdir = self.ip.srcdir
            doy = '%03d'%int(doy)
            self.a1Files = []
            self.a2Files = []
            self.weighting = []
            for i,year in enumerate(yearList):
                year = str(yearList[i])
                # Lewis modified for different structure
                a2 = insensitive_glob(srcdir+'/%s2/%s/%s/*/'%(product,year,tile)+\
                         '*'+product+'2.A'+year+str(doy)+'.'+tile+'.'+version+'.*.hdf')
                a1 = insensitive_glob(srcdir+'/%s1/%s/%s/*/'%(product,year,tile)+\
                         '*'+product+'1.A'+year+str(doy)+'.'+tile+'.'+version+'.*.hdf')

                if len(a1) == 0 or len(a2) == 0: 
                  a2 = insensitive_glob(srcdir+'/'+'*'+product+'2.A'+year+str(doy)+'.'+tile+'.'+version+'.*.hdf')
                  a1 = insensitive_glob(srcdir+'/'+'*'+product+'1.A'+year+str(doy)+'.'+tile+'.'+version+'.*.hdf')

                if len(a2) == 0 or len(a1) == 0:
                  a2.extend(insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'2/%s/*.A%s%s.%s.*hdf'%(year,year,str(doy),self.ip.tile)))
                  a1.extend(insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'1/%s/*.A%s%s.%s.*hdf'%(year,year,str(doy),self.ip.tile)))


                if len(a1) == 0 or len(a2) == 0 :
                    logging.info( '========',extra=self.d)
                    logging.info( 'inconsistent or missing data for doy \
                              %s year %s tile %s version %s'%(str(doy),str(year),tile,version),extra=self.d)
                    logging.info( 'A1: %s'%a1,extra=self.d)
                    logging.info( 'A2: %s'%a2,extra=self.d)
                    logging.info( '========',extra=self.d)
                else:
                    self.a1Files.append(a1[0])
                    self.a2Files.append(a2[0])
                    self.weighting.append(self.ip.focus[i])
                    logging.info( a1[0],extra=self.d)
       
        return self.a1Files,self.a2Files,self.weighting


    def sortDims(self,ip,op):
        """
        Change the dimensions of the dataset if you want to process a sub-image
        

        Parameters
        ----------
        
        ip  :  Long array[4] 
               containing [s0, ns,l0,nl]

        op  :  Long array[4]
               original data dimensions [s0,send,l0,lend]


        where send, lend are the end sample and line numbers required
        
        Returns
        -------

        op  :  Long array[4]
               DIMS format containing [s0,ns,l0,nl]

        where send, lend are the end sample and line numbers required

        """
        s0 = op[0]
        ns = op[1] - op[0]
        l0 = op[2]
        nl = op[3] - op[2]
        s0 = max(ip[0],s0)

        if ip[1] != -1:
            ns = ip[1]
        if ip[3] != -1:
            nl = ip[3]

        l0 = max(ip[2],l0)
        s0 = min(max(0,s0),op[1])
        l0 = min(max(0,l0),op[3])

        ns = max(min(ns,op[1] - op[0] - s0),1)
        nl = max(min(nl,op[3] - op[2] - l0),1)

        op = [int(s0), int(ns), int(l0), int(nl)]

        return op



    def getModisAlbedo(self,fileName,QaFile,sdim,backupscale):
        """
        Extract data from MODIS data file (C5)

        Parameters
        ----------

        fileName    : MCD43 data file (A1)
                      The HDF filename containing the required data (e.g. MCD43A1.A2008001.h19v08.005.2008020042141.hdf)
        QaFile      : MCD43 QA file (A2)  
                      The HDF filename for the associated QA datafile(e.g. MCD43A2.A2008001.h19v08.005.2008020042141.hdf)
        sdim        : array [-1,-1,-1,-1] used to ectract subset / subsample. 
                      The format is [s0,ns,l0,nl]
        backupscale : translation coefficients for QA flags
                      e.g. [1.,0.7,0.49,0.3]
 
        Returns
        -------

        dictionary : containing

        err      : long              
                   error code (set to 1 if file error on qa data and ERRCODE also set as string
        nSnow    : long              
                   number of snow samples found
        isSnow   : long[ns,nl]       
                   spatial flag for snow (1), no snow (0), no data (-1)
        n        : float[ns,nl]     
                   spatial flag for weighted number fo samples (according to backupscale)
        data : float[nb,3,ns,nl] 
                   extracted dataset of model parameters
        nb       : long              
                   number of bands
        ns       : long             
                   number of samples
        nl       : long              
                   number of lines
        land     : long[ns,nl]       
                   land category (15 = do not process)

        Notes
        -----
        0  :    Shallow ocean
        1  :   Land (Nothing else but land)
        2  :   Ocean and lake shorelines
        3  :  Shallow inland water
        4  :  Ephemeral water
        5  :  Deep inland water
        6  :  Moderate or continental ocean
        7 :   Deep ocean 
         
        see https://lpdaac.usgs.gov/products/modis_products_table/mcd43a2
                                    
        """
        bands = self.ip.bands

        duff = long(32767)
        oneScale = 1000
        scale = 1.0/oneScale
        nBands = len(self.ip.bands)
        mask = 3
        #get the data size
        logging.info( '...reading qa... %s'%QaFile,extra=self.d)

        err=0
        openQA=0;
        try:
            hdf = SD.SD(QaFile)
            openQA=1
            sds_1 = hdf.select(0)
            nl,ns = sds_1.dimensions().values()
            fileDims = [0,ns,0,nl]
            #take a subset if input sdims values require it
            s0,ns,l0,nl = self.sortDims(sdim,fileDims)
            #get the data
            QA = sds_1.get(start=[s0,l0],count=[ns,nl])
    
            sds_2 = hdf.select(1)
            QA2 = sds_2.get(start=[s0,l0],count=[ns,nl])
    
            sds_3 = hdf.select(2)
            QA3 = sds_3.get(start=[s0,l0],count=[ns,nl])
    
            sds_4 = hdf.select(3)
            QA4 = sds_4.get(start=[s0,l0],count=[ns,nl])
            
            hdf.end()
    
            logging.info( 'done ...reading data...',extra=self.d)
            data = []
            goodData = np.ones((ns,nl),dtype=np.int8)
            typeOfInversion = np.zeros((ns,nl),dtype=np.dtype('f4')) - 1.0
            typeOfSnow = np.zeros((ns,nl),dtype=np.int8) - 1
       
            hdf = SD.SD(fileName)
            
            for i in range(nBands):
                logging.info( '  ... band %d'%i,extra=self.d)
                # Lewis: ensure this is an int
                sds = hdf.select(int(bands[i]))  
                ithis = sds.get(start=[s0,l0,0],count=[ns,nl,3])
                #filter out duff values
                for j in range(3):
                    badValues = np.where(ithis[:,:,j] == duff)
                    goodData[badValues] = 0
                    
                data.append(ithis.astype(np.float32))

            for i in range(nBands):
                data[i][goodData==0] = 0.
            hdf.end()
            
            logging.info( 'done ...sorting mask...',extra=self.d)
            #get the land mask (taking only 'Land (Nothing else but land)')
            land = (QA3 >> 4 & 15) & (goodData == 1)
            w = np.where(land != 1)
            land[w] = 0
            bandqa = QA4 & mask
    
            for k in range(1,7):
                QA4 = QA4/16
                bandqb = QA4 & mask
                w=np.where(bandqb > bandqa)
                bandqa[w] = bandqb[w]
    
            """     
            set wanted data to 1
            QA == 0 --> 'full inversion'
            QA == 1 --> 'mag inversion'
            QA2 ==0 --> 'no snow'
            QA2 ==1 --> 'snow'
            sort the QA info to assign a weight (specified by backupscale)
            for 0 : best quality, full inversion (WoDs, RMSE majority good)
            1 : good quality, full inversion 
            2 : Magnitude inversion (numobs >=7) 
            3 : Magnitude inversion (numobs >=3&<7) 
            where the QA is determined as the maximum (ie poorest quality) over the wavebands  
            """  
            if len(backupscale) == 2:
                bandqa = QA
            for i in range(len(backupscale)):
                if backupscale[i] != 0.0:
                    # VERSION 2: only take land 1
                    w = np.where((bandqa == i) & (land == 1) & (goodData == 1))                
                    typeOfInversion[w] = backupscale[i]
                    logging.info( 'weight: %.3f %d'%(backupscale[i],len(w[0])),extra=self.d)
                else:
                    logging.info( 'weight: %.3f %d'%(backupscale[i],0),extra=self.d)
    
    
            w = np.where((QA2 == 1) & (typeOfInversion > 0))
            typeOfSnow[w] = 1
            nSnow = len(w[0])
        
    
            w = np.where((QA2 == 0) & (typeOfInversion > 0))
            typeOfSnow[w] = 0
    
            w = np.where(typeOfInversion < 0)
            typeOfInversion[w] = 0.0
            N = typeOfInversion
    
            logging.info('done',extra=self.d)
        
        except:
            fileErr = None
            
            if(openQA==1):
                filenameTab=fileName.replace("\\", "/").split("/");
                filename=filenameTab[len(filenameTab)-1]
                fileErr = open(self.errFiles + '/' + 'errFile.'+filename, 'w')
            else:
                filenameTab=QaFile.replace("\\", "/").split("/");
                filename=filenameTab[len(filenameTab)-1]
                fileErr = open(self.errFiles + '/' + 'errFile.'+filename, 'w')
            
            fileErr.close();
            
            logging.warning( 'error in %s'%fileName,extra=self.d)
            err=1
             
        if err==0:
            return dict(err=err, nSnow=nSnow, isSnow=typeOfSnow, N=N, data=np.asarray(data)*scale, nb=nBands, ns=ns, nl=nl, land=land)
        else:
            return dict(err=err, nSnow=None, isSnow=None, N=None, data=None, nb=None, ns=None, nl=None, land=None)
        

    def allocateData(self,nb,ns,nl,nsets,flag):
        """

        Allocate data for calculating image statistics

        Parameters
        ----------
        
        ns  : long 
              number of samples
        nl  : long 
              number of lines
        nb  : long 
              number of bands
        nsets : long
              number of data sets      
 
        Returns
        ------- 

        Dictionary containing:
        
        sum  : sum 
               (weighted by w) for band of nb, kernel of 3, sample of ns, line of nl
        sum2 : sum^2 
               (weighted by w^2) for each combination of band/kernel, sample of ns, line of nl
        n    : sum of weight 
               for sample of ns, line of nl
        nSamples : number of samples 
                   per pixel

        """
        # allocate the arrays for mean and variance
        if (flag == 0):
            
            data = np.zeros((nsets,nb,ns,nl,3),dtype=np.float32)           
            n = np.zeros((nsets,nb,ns,nl,3),dtype=np.float32)
            nSamples = np.zeros((nsets,nb,ns,nl),dtype=np.int)
            logging.info( 'numb of bands %d, numb of samples %d, numb of lines %d'%(nb, ns, nl),extra=self.d)
 
        else:
            data = -1
            n = -1
            nSamples = -1

        return dict(n=n, data=data, nSamples=nSamples)


    def incrementSamples(self,sumData,samples,isSnow,index):
        """
        Increment information in the sumdata dictionary with data from a MODIS image in samples dictionary

        The data here are in samples (a dictionary).

        The snow mask is in data['isSnow'] != -1
        i.e. this is set to -1 for 'no data'
        It is set to 1 if a pixel is 'snow' and 0 if snow free

        Parameters
        ----------

        sumData : dictionary
                  Containing ``sum``, ``sum2``, ``n`` and ``nSamples`` that we wish to accumulate into
        samples : dictionary
                  Containing ``sum``, ``sum2``, ``n`` and ``nSamples`` that are the values to be added to sumData
        isSnow : integer
                 Code for snow processing type. The flag isSnow is used to determine the type of coverage:
          0 : no snow only
          1 : snow only
          2 : snow and no snow together

         index : integer
                 dataset index

        Returns
        --------

        sumData : dictionary
                  With ``data``, ``n`` and ``nSamples`` in the arrays

        """
        ns = samples['ns']
        nl = samples['nl']
        nb = samples['nb']

       
        # generate some negative masks 
        if (isSnow == 0):
            # no snow only
            w = np.where(samples['isSnow'] != 0)
        elif (isSnow == 1):
            # snow only
            w = np.where(samples['isSnow'] != 1)
        else:
            w = np.where(samples['isSnow'] == -1)

        # so w is a mask of where we *dont* have data (that we want)
        samplesN = samples['N'].copy()
        samplesN[w] = 0

        # sum f0 over all bands
        f0sum = samples['data'][:,:,:,0].sum(axis=0)
        # find where == 0  as f0 == 0 is likely dodgy     
        w = np.where((f0sum<=0) & (samplesN>0))
        if len(w)>0 and len(w[0])>0:
          logging.info('N %.2f'%samplesN.sum(),extra=self.d)
          logging.info("deleting %d samples that are zero"%len(w[0]),extra=self.d)
          samplesN[w] = 0
          logging.info('N %.2f'%samplesN.sum(),extra=self.d)
        else:
          logging.info('N %.2f'%samplesN.sum(),extra=self.d)
        # save some time on duffers
        if samplesN.sum() == 0:
          logging.info('No samples here ...',extra=self.d)
          return sumData

        weight = np.zeros((nb,ns,nl,3),dtype=np.float32)

        for i in range(nb):
            for j in range(3):
                weight[i,:,:,j] = samplesN

        # shrink the data
        sweight = np.zeros((nb,ns/self.ip.shrink,nl/self.ip.shrink,3),dtype=float)
        sdata = np.zeros((nb,ns/self.ip.shrink,nl/self.ip.shrink,3),dtype=float)

        # so sweightdata is the observations multiplied by the weight
        sweightdata = samples['data']*weight
        for i in xrange(nb):
          for j in xrange(3):
            sweight[i,:,:,j] = self.shrunk(weight[i,:,:,j],ns,nl,self.ip.shrink)
            sdata[i,:,:,j] = self.shrunk(sweightdata[i,:,:,j],ns,nl,self.ip.shrink)
            ww = np.where(sweight[i,:,:,j]>0)
            sdata[i,:,:,j][ww] /= sweight[i,:,:,j][ww]
        # now sdata is re-normalised so its just the data again

        # store the data
        sumData['n'][index,...] = sweight
        sumData['nSamples'][index,...] = (sweight[:,:,:,0] > 0).astype(int)
        sumData['data'][index,...] = sdata

        return sumData



    def processAlbedo(self,yearList,doy,sdmins):
        """
         For some given doy, read and process MODIS kernel data 

        Parameters
        ----------

         yearlist  :  string array
                      candidate years to process
         doy       :  long 
                      day of year to process
         SDIMS     :  long array[4]
                      [s0,ns,l0,ns] or [-1,-1,-1,-1] for full dataset
         

         Returns
         -------

         processed              : boolean
                                 True if data read ok otherwise False
         totalSnow              :  long
                                   total number of snow samples
         sumdataNoSnow          :  sumdata-type dict 
                                   for no snow information
         sumdataSnow            :  sumdata-type dict
                                   for snow information
         sumdataWithSnow        : sumdata-type dict 
                                  for snow and no-snow information
         nb                     : long  
                                   number of bands
         ns                     : long
                                   number of samples
         nl                     : long
                                   number of lines
         tile                   : string
                                  tile name
         version                : string
                                   version name (005)
         doy                    : string
                                   doy string to process e.g. 001
         land                   : float[ns,nl]
                                  land codes 
        
         Where:
 
          land category (15 = do not process)
                                0          Shallow ocean
                                1         Land (Nothing else but land)
                                2         Ocean and lake shorelines
                                3         Shallow inland water
                                4         Ephemeral water
                                5         Deep inland water
                                6         Moderate or continental ocean
                                7         Deep ocean
         see https://lpdaac.usgs.gov/lpdaac/products/modis_products_table/brdf_albedo_quality/16_day_l3_global_500m/v5/combined

         See Also
         ---------
         self.increment_samples()
                       

        """
        self.a1Files = None
        self.a2Files = None
        a1Files, a2Files, weighting = self.getValidFiles(yearList,doy[0])

        logging.info('doy %s'%str(doy[0]),extra=self.d)
        for i in xrange(len(a1Files)):
           logging.info('  %d %s %s'%(i,str(a1Files[i]),str(a2Files[i])),extra=self.d)

        foundOne = False
        thisOne = 0
        # try to file at least one file that works ...
        while not foundOne: 
          thisData = self.getModisAlbedo(a1Files[thisOne],a2Files[thisOne],\
                          sdmins,np.asarray(self.ip.backupscale)*weighting[thisOne])
          if thisData['err'] != 0:
            thisOne += 1
            logging.warning('error in getModisAlbedo for %s %s'%(a1Files[thisOne],a2Files[thisOne]),extra=self.d)
            # try another one?
            if thisOne == len(a1Files):
              logging.error('error in getModisAlbedo: No valid data files found')
              thisData['err'] = 1
              return False,0,0,0,0,nb, ns, nl, 0
          else:
            foundOne = True

        ns = thisData['ns']
        nl = thisData['nl']
        nb = thisData['nb']
        totalSnow = 0
        #set up arrays for sum, n and sum2
        logging.info('data allocation',extra=self.d)
        dontsnow = self.ip.dontsnow
        dontnosnow = self.ip.dontnosnow 
        dontwithsnow = self.ip.dontwithsnow 
 
        nsets = len(a1Files)
        if nsets == 0:
          logging.error('error in file specification: zero length list of files a1Files',extra=self.d)
          thisData['err'] = 1
          return False,0,0,0,0,nb, ns, nl, 0

        try:
          sumDataSnow = self.allocateData(nb,ns/self.ip.shrink,nl/self.ip.shrink,nsets,dontsnow)
          sumDataNoSnow = self.allocateData(nb,ns/self.ip.shrink,nl/self.ip.shrink,nsets,dontnosnow)
          sumDataWithSnow = self.allocateData(nb,ns/self.ip.shrink,nl/self.ip.shrink,nsets,dontwithsnow)
        except:
          logging.error('error in memory allocation: nb %d ns %d nl %d nsets %s'%(nb,ns/self.ip.shrink,nl/self.ip.shrink,nsets),extra=self.d)
          return False,0,0,0,0,nb, ns, nl, 0 

        land = np.zeros((ns,nl),dtype='bool') 
 
        if dontsnow == 0 or dontwithsnow == 0:
            totalSnow = thisData['nSnow']
            logging.info('n snow %d'%totalSnow,extra=self.d)

        for i in range(len(a1Files)):
            logging.info( 'file %d/%d'%(i,len(a1Files)),extra=self.d)
            logging.info( 'doy %s %s'%(str(doy[0]),str(a1Files[i])),extra=self.d)
            #only read if i > 1 as we have read first file above
            if (i != thisOne):
                thisData = self.getModisAlbedo(a1Files[i],a2Files[i],sdmins,\
                                     np.asarray(self.ip.backupscale)*weighting[i])
            if (thisData['err'] != 0):
                logging.warning( 'warning opening file: %s'%str(a1Files[i]),extra=self.d)
            else:
                # Lewis: sort the land info
                land = (land | (thisData['land'] == 1))
                
                logging.info( '... incrementing samples',extra=self.d)
                if (dontsnow == 0):
                    sumDataSnow = self.incrementSamples(sumDataSnow,thisData,1,i)
                if (dontnosnow == 0):
                    sumDataNoSnow = self.incrementSamples(sumDataNoSnow,thisData,0,i)
                if (dontwithsnow == 0):
                    sumDataWithSnow = self.incrementSamples(sumDataWithSnow,thisData,2,i)
                if (dontsnow == 0 or dontwithsnow == 0):
                    totalSnow += thisData['nSnow']
                logging.info( 'done',extra=self.d)
        
    
        return True,totalSnow, sumDataNoSnow, sumDataSnow, sumDataWithSnow, nb, ns, nl, land


    def calculateStats(self,sumData):
        """
        Given data from sumdata dict, calculate mean and var/covar information

        Parameters
        -----------
        sumData : sumdata dict

        Returns
        -------

        n : float array (ns,nl)
            containing weight for sample
        meanData : float array (nb,ns,nl,3)
                   weighted mean 
        sdData : float array (nb,ns,nl,3)
                   weighted std dev


        See Also
        --------

        self.increment_samples()

        """
        focus = np.array(self.weighting)
        logging.info("sorting data arrays ...",extra=self.d)
        n = sumData['n']
        data = sumData['data']
        nSamples = sumData['nSamples']
        logging.info("...done",extra=self.d)

        sumN = np.sum(nSamples,axis=0)
        ww = np.where(sumN>0)

        #if not (focus == 1).all():
        #  # weight the n terms
        #  for i in xrange(len(focus)):
        #    n[i,...] *= focus[i]

        total = np.sum(data*n,axis=0)
        ntot = np.sum(n,axis=0)
        ntot2 = np.sum(n**2,axis=0)
        meanData = np.zeros_like(total)
        meanData[ww] = total[ww]/ntot[ww]
        diff = data - meanData
        d2 = n*(diff**2)
        var = np.sum(d2,axis=0)
        #var[ww] /= ntot[ww]
        # small number correction: see ATBD
        num = (ntot**2 - ntot2)
        # set all valid to min err
        store = np.zeros_like(var)
        store[ww] = np.sqrt(self.minvar)
        store[var>0] = var[var>0]
        var = store
        # now fill in
        num = ntot**2 - ntot2
        ww = np.where(num>0)
        var[ww] = ntot[ww] * var[ww]/num[ww]
        # sqrt
        sdData = np.sqrt(var)
        ww = np.where(sumN>0)
        samp = sdData[ww]
        sqmax = np.sqrt(self.maxvar)
        sdData[ww][samp==0.] = np.sqrt(self.minvar)
        sdData[ww][samp>sqmax] = sqmax
          
        return ntot, meanData, sdData


    def rebin(self,a,shape):
        """
        python equivalent of IDL rebin, copied from
        stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
        """
        sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
        return a.reshape(sh).mean(-1).mean(1)


    def shrunk(self,fdata,ns,nl,shrink,sdata=None):
        """
        shrink but account for zeros

        if we specify sdata (sd) we use this in the scaling
        
        var = sdata^2
        vscale = 1./var
        fs = fdata * vscale
        
        Then shrink (mean) and normalise by vscale

        Otherwise, we normalise the mean to not count zero values in the inoput data.

        Parameters
        -----------
        fdata : float array
                data to be shrunk
        ns  : integer
               number of samples
        nl : integer
               number of lines
        shrink : integer
               shrink factor (1 for do nothing)

        sdata : float array
                sd of data


        Returns
        --------
        odata : float array
                derived from fdata but shrunk by `shrink` factor
                where shrinking produces a local mean, ignoring zero values

        Or, if ``sdata != None``, then in addition

        osdata : float array
                 Shrunk std dev array

        """
        if shrink == 1:
          if sdata != None:
            return fdata,sdata
          else:
            return fdata
        if sdata != None:
          var = sdata**2
          var[fdata==0] = 0.
          vscale = np.zeros_like(var)
          w = np.where(var>0)
          vscale[w] = 1./var[w]
          idata = fdata * vscale
          shrunkData = self.rebin(idata,(ns/shrink,nl/shrink))
          shrunkVar = self.rebin(vscale,(ns/shrink,nl/shrink))
          w = np.where(shrunkVar>0)
          odata = np.zeros_like(shrunkData)
          odata[w] = shrunkData[w]/shrunkVar[w]
          mdata = np.zeros_like(fdata)
          w1 = np.where(var>0)
          mdata[w1] = 1
          n = self.rebin(mdata,(ns/shrink,nl/shrink))
          osdata = np.zeros_like(shrunkData)
          osdata[w] = np.sqrt(n[w]/(shrunkVar[w]))
          return odata,osdata
        # else, rescale and account for zeros
        odata = self.rebin(fdata,(ns/shrink,nl/shrink))
        w = np.where(fdata == 0)
        mdata = np.ones_like(fdata)
        mdata[w] = 0
        n = self.rebin(mdata,(ns/shrink,nl/shrink))
        w2 = np.where((n > 0) & (n < 1))
        odata[w2] /= n[w2]

        return odata

    def writeRST(self,ns,nl,nb,mean,sd,n,land,snowType,p,doy,shrink=None,order=0):
        """
        Write relevant images and stats to an rst file (for further processing with sphinx)

        Parameters
        -----------
        ns : integer
             number of samples
        nl : integer
             number of lines 
        nb : integer
             number of bands
        mean : float array 
               containing mean data
        sd : float array 
               containing std dev data
        n  : float array 
               containing sample weight data
        land : float array 
               containing land / sea mask
        snowType : string 
                e.g. SnowAndNoSnow, used in the output filenames
        p : integer 
            index for part (0 default)
        doy : string
            DOY string e.g. 001

        Options
        --------

        shrink : integer
                 shrink factor applied
        order : integer
                 Specify a variation in the array ordering of `mean` and `sd`.
                 If `order` is 0, we have ``mean[band,sample,line,kernel_number]``
                 , otherwise it is ``mean[band*3+kernel_number,sample,line]``      

        Notes
        -----

        The filenames for graphs are derived from:

        filename =  'Kernels_' + doy + '_' + self.ip.version + '_' + self.ip.tile + '_' 
                                                        + snowType +'_'+ p + 'XXXX.png'

        And the ReST file:

        rstFile = self.ip.rstdir + '/Kernels_' + doy + '_' + self.ip.version + '_' + 
                                  self.ip.tile + '_' + snowType +'_'+ p + 'XXXX.rst'

        """
        import pylab as plt

        p = '%02d'%p
        shrink = shrink or self.ip.shrink
        doy = '%03d'%int(doy)


        kernels = self.ip.basename or Kernels
        filename = '%s_'%kernels + doy + '_' + self.ip.version + '_' +\
                                 self.ip.tile + '_' + snowType +'_'+ p + 'XXXX.png'
        rstFile = self.ip.rstdir + '/%s_'%kernels + doy + '_' + self.ip.version + '_' +\
                                 self.ip.tile + '_' + snowType +'_'+ p + 'XXXX.rst'

        try:
          f0 = open(rstFile.replace('XXXX','F0'),'w')
          f1 = open(rstFile.replace('XXXX','F1'),'w')
          f2 = open(rstFile.replace('XXXX','F2'),'w')
        except:
          logging.info( 'rst write failure: %s'%rstFile,extra=self.d)
          return


        header = '''
Stage 1 processing
====================

Stage 1 prior processing for DOY %s\n
        
MODIS version: %s
        
Tile: %s

Snow processing: %s

Process part: %s

Generated with:

.. code-block:: csh

  %s

        '''%(doy,self.ip.version,self.ip.tile,snowType,p,''.join( [i.replace('[',"'[").replace(']',"]'") + ' ' for i in sys.argv]))
        bNames = self.idSamples(nb)
        bNames.append('Weighted number of samples')
        bNames.append('land mask')

        logging.info( 'writing %s'%filename,extra=self.d)
        if nb == 2:
            defBands = [4,1,1]
        else:
            defBands = [1,10,7]
        descrip = snowType + ' MODIS Mean/SD ' + doy + ' over the years ' + str(self.yearList)  + \
            ' version ' + self.ip.version + ' tile '+ self.ip.tile + \
            ' using input MODIS bands '
        for band in self.ip.bands:
            descrip = descrip + str(band) + ' '

        for f in [f0,f1,f2]:
          f.write(header + '\n')
          f.write('`Log file <%s>`_\n\n\n'%self.log)
          f.write(descrip + '\n')

        fs = [f0,f1,f2]
        try:
          count = self.count
        except:
          count = 0

        #shrunkLand = self.shrunk(land,ns,nl,shrink)
        shrunkLand = land
        fig = plt.figure(count);plt.clf()
        ax = fig.add_subplot(111)
        im=ax.imshow(shrunkLand,interpolation='nearest');
        ax.set_aspect('equal')
        fn = filename.replace('XXXX',"LandMask_%s"%self.id)
        fig.savefig( self.ip.rstdir + '/' + fn)
        count += 1

        for f in [f0,f1,f2]:
          f.write("\n\n.. image:: %s\n\n\n\nLand mask\n"%fn)
          f.write("\n\nNumber of land pixels: %d\nProportion of land pixels:%.3f\n\n"%(np.sum(shrunkLand),np.mean(shrunkLand)))

          f.write("\n\n\nMean and SD results\n--------------------------\n")

        if order == 0:
          shrunkN = n[0,:,:,0]
          #shrunkN = self.shrunk(n[0,:,:,0],ns,nl,shrink)
        else:
          shrunkN = n[:,:]
          #shrunkN = self.shrunk(n[:,:],ns,nl,shrink)

        for i in range(nb):
            for j in range(3):
                if order == 0:
                  shrunkMean = mean[i,:,:,j]
                  shrunkSd = sd[i,:,:,j]
                  #shrunkMean,shrunkSd = self.shrunk(mean[i,:,:,j],ns,nl,shrink,sdata=sd[i,:,:,j])
                else:
                  shrunkMean = mean[i*3+j,:,:]
                  shrunkSd = sd[i*3+j,:,:]
                  #shrunkMean,shrunkSd = self.shrunk(mean[i*3+j,:,:],ns,nl,shrink,sdata=sd[i*3+j,:,:])

                fig = plt.figure(count);plt.clf()
                ax = fig.add_subplot(111)
                im=ax.imshow(self.fix(shrunkMean),interpolation='nearest');fig.colorbar(im)
                ax.set_aspect('equal')
                fn = filename.replace('XXXX',"Mean_band_%03d_F%d_%s"%(i,j,self.id))
                fig.savefig( self.ip.rstdir + '/' + fn)
                fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'mean: band %03d F%d'%(i,j)))
                count += 1
                # histogram stackoverflow.com/questions/5328556/histogram-matplotlib
                data = shrunkMean[shrunkMean>0]
                try:
                  hist,bins = np.histogram(self.fix(data),bins=50)
                  width = 0.7*(bins[1]-bins[0])
                  center = (bins[:-1]+bins[1:])/2
                  plt.figure(count);plt.clf()
                  plt.bar(center,hist,align='center',width=width)
                  fn = filename.replace('XXXX',"hist_Mean_band_%03d_F%d"%(i,j))
                  plt.savefig( self.ip.rstdir + '/' + fn)
                  fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'mean: band %03d F%d'%(i,j)))
                  fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                except:
                  if len(data) == 0:
                    data = [0.]
                  fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                  logging.warning('failed to create mean histogram',extra=self.d)
                count += 1

                # STD DEV
                fig = plt.figure(count);plt.clf()
                ax = fig.add_subplot(111)
                im=ax.imshow(self.fix(shrunkSd),interpolation='nearest');fig.colorbar(im)
                ax.set_aspect('equal')
                fn = filename.replace('XXXX',"SD_band_%03d_F%d_%s"%(i,j,self.id))
                fig.savefig( self.ip.rstdir + '/' + fn)
                fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'sd: band %03d F%d'%(i,j)))
                count += 1
                # histogram stackoverflow.com/questions/5328556/histogram-matplotlib
                data = shrunkSd[shrunkSd>0]
                try:
                  hist,bins = np.histogram(self.fix(data),bins=50)
                  width = 0.7*(bins[1]-bins[0])
                  center = (bins[:-1]+bins[1:])/2
                  plt.figure(count);plt.clf()
                  plt.bar(center,hist,align='center',width=width)
                  fn = filename.replace('XXXX',"hist_SD_band_%03d_F%d"%(i,j))
                  plt.savefig( self.ip.rstdir + '/' + fn)
                  fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'sd: band %03d F%d'%(i,j)))
                  fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                except:
                  if len(data) == 0:
                    data = [0.]
                    fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                  logging.warning('failed to create sd histogram',extra=self.d)
                count += 1

                # coefficient of variation
                fig = plt.figure(count);plt.clf()
                ax = fig.add_subplot(111)
                cv = np.zeros_like(shrunkSd)
                w = np.where(shrunkMean>0)
                cv[w] = shrunkSd[w]/shrunkMean[w]
                im=ax.imshow(self.fix(cv),interpolation='nearest');fig.colorbar(im)
                ax.set_aspect('equal')
                fn = filename.replace('XXXX',"CV_band_%03d_F%d_%s"%(i,j,self.id))
                fig.savefig( self.ip.rstdir + '/' + fn)
                fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'cv: band %03d F%d'%(i,j)))
                count += 1
                # histogram stackoverflow.com/questions/5328556/histogram-matplotlib
                data = cv[cv>0]
                try:
                  hist,bins = np.histogram(self.fix(data),bins=50)
                  width = 0.7*(bins[1]-bins[0])
                  center = (bins[:-1]+bins[1:])/2
                  plt.figure(count);plt.clf()
                  plt.bar(center,hist,align='center',width=width)
                  fn = filename.replace('XXXX',"hist_CV_band_%03d_F%d"%(i,j))
                  plt.savefig( self.ip.rstdir + '/' + fn)
                  fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'cv: band %03d F%d'%(i,j)))
                  fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                except:
                  if len(data) == 0:
                    data = [0.]
                    fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                  logging.warning('failed to create cv histogram',extra=self.d)
                count += 1

                '''
                # normalise uncertainty
                fig = plt.figure(count);plt.clf()
                ax = fig.add_subplot(111)
                nsd = np.zeros_like(shrunkSd)
                nYears = len(self.yearList)
                # The idea here is that if we had nYears samples
                # each with std dev S, then we would expect std dev to be S/sqrt(nYears)
                # in fact, we have shrunkN samples, so we might consider
                # raising the std dev by a factor sqrt(nYears/shrunkN)
                w1 = np.where((shrunkN>0) & (shrunkSd <= np.sqrt(self.minvar)*1000))
                w = np.where(shrunkN>0)
                nChange = np.array(w1).shape[1]
                nsd[w] = np.sqrt(float(nYears)/shrunkN[w])*shrunkSd[w]
                # set nominally low values to the median
                if len(w[0]) == 0:
                  median = np.sqrt(self.maxvar)
                else:
                  median = np.median(nsd[w])
                  median = np.median(nsd[w][nsd[w]>median])
                if np.isnan(median):
                  median = 1.0
                if median == 0:
                  median = np.sqrt(self.maxvar)
                if nChange > 0:
                  logging.info("  %d %d Modifying %d low threshold std err to upper quartile (%.2f)"%(i,j,nChange,median),extra=self.d) 
                nsd[w1] = median
                # deal with low sample number terms
                w1 = np.where((shrunkN<2)&(shrunkN>0)&(nsd<median))
                nChange = np.array(w1).shape[1]
                if nChange > 0:
                  logging.info("  %d %d Modifying %d low sample number entries to upper quartile (%.2f)"%(i,j,nChange,median),extra=self.d)
                nsd[w1] = median
                nsd[nsd>np.sqrt(self.maxvar)] = np.sqrt(self.maxvar)
                im=ax.imshow(self.fix(nsd),interpolation='nearest');fig.colorbar(im)
                ax.set_aspect('equal')
                fn = filename.replace('XXXX',"NSD_band_%03d_F%d_%s"%(i,j,self.id))
                fig.savefig( self.ip.rstdir + '/' + fn)
                fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'normalised std dev: band %03d F%d'%(i,j)))
                fs[j].write("\n\nThe normalisation is a scaling of std dev by sqrt(Nyears/Nsamples) where \
                            Nyears = %d here (threshold at %.2f)\n\n\n"%(nYears,self.maxvar))
                count += 1
                # histogram stackoverflow.com/questions/5328556/histogram-matplotlib
                data = nsd[nsd>0]
                try:
                  hist,bins = np.histogram(self.fix(data),bins=50)
                  width = 0.7*(bins[1]-bins[0])
                  center = (bins[:-1]+bins[1:])/2
                  plt.figure(count);plt.clf()
                  plt.bar(center,hist,align='center',width=width)
                  fn = filename.replace('XXXX',"hist_NSD_band_%03d_F%d"%(i,j))
                  plt.savefig( self.ip.rstdir + '/' + fn)
                  fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'nsd: band %03d F%d'%(i,j)))
                  fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                except:
                  if len(data) == 0:
                    data = [0.]
                    fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                  logging.warning('failed to create nsd histogram',extra=self.d)
                count += 1

                # normalised coefficient of variation
                fig = plt.figure(count);plt.clf()
                ax = fig.add_subplot(111)
                cv = np.zeros_like(shrunkSd)
                w = np.where(shrunkMean>0)
                cv[w] = nsd[w]/shrunkMean[w]
                im=ax.imshow(self.fix(cv),interpolation='nearest');fig.colorbar(im)
                ax.set_aspect('equal')
                fn = filename.replace('XXXX',"NSDCV_band_%03d_F%d_%s"%(i,j,self.id))
                fig.savefig( self.ip.rstdir + '/' + fn)
                fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'cv from normalised std dev: band %03d F%d'%(i,j)))
                count += 1
                # histogram stackoverflow.com/questions/5328556/histogram-matplotlib
                data = cv[cv>0]
                try:
                  hist,bins = np.histogram(self.fix(data),bins=50)
                  width = 0.7*(bins[1]-bins[0])
                  center = (bins[:-1]+bins[1:])/2
                  plt.figure(count);plt.clf()
                  plt.bar(center,hist,align='center',width=width)
                  fn = filename.replace('XXXX',"hist_NSDCV_band_%03d_F%d"%(i,j))
                  plt.savefig( self.ip.rstdir + '/' + fn)
                  fs[j].write("\n\n.. image:: %s\n\n\n\n%s\n"%(fn,'cv from normalised std dev: band %03d F%d'%(i,j)))
                  fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                except:
                  if len(data) == 0:
                    data = [0.]
                    fs[j].write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                             %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
                  logging.warning('failed to create nsd cv histogram',extra=self.d)
                count += 1
                '''
        for f in fs:
          f.write("\n\n\nSamples\n-----------------\n")

        fig = plt.figure(count);plt.clf()
        ax = fig.add_subplot(111)
        im=ax.imshow(shrunkN,interpolation='nearest');fig.colorbar(im)
        ax.set_aspect('equal')
        fn = filename.replace('XXXX',"Nsamp_%s"%self.id)
        fig.savefig( self.ip.rstdir + '/' + fn)
        for f in fs:
          f.write("\n\n.. image:: %s\n    :width: 80%%\n\n\nNumber of samples (weighted)\n"%fn)
        count += 1
        # histogram stackoverflow.com/questions/5328556/histogram-matplotlib
        try:
          data = shrunkN[shrunkN>0]
          hist,bins = np.histogram(data,bins=50)
          width = 0.7*(bins[1]-bins[0])
          center = (bins[:-1]+bins[1:])/2
          plt.figure(count);plt.clf()
          plt.bar(center,hist,align='center',width=width)
          fn = filename.replace('XXXX',"hist")
          plt.savefig( self.ip.rstdir + '/' + fn)
          for f in fs:
            f.write("\n\n.. image:: %s\n\n\n\n"%(fn))
            f.write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                            %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
        except:
          if len(data) == 0:
            data = [0.]
            f.write("\n\nMean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f\n"\
                            %(np.mean(data),np.median(data),np.min(data),np.max(data),stats.sem(data,axis=None,ddof=1)))
          logging.warning('failed to create N histogram',extra=self.d)

        # count incrementer
        self.count = count + 1
        

    def fix(self,data):
        '''
        Threshold the data at 99th centile for plotting
        quantising data to 0.001 bins for histogram calculation

        If this fails for any reason, we return the original data

        Parameters
        -----------

        data : float array
               Data to be manipulated
              

        Returns
        -------
 
        data : float array
               Data having been manipulated

 
        '''
        try:
          (n,b) = np.histogram(data.flatten(),np.arange(0,1,0.001)) 
          nc = n.cumsum()
          N = nc[-1]*0.99 
          w = np.where(nc >= N)[0][0]
          max = b[w]
          d = data.copy()
          d[d>=max]=max
          return d
        except:
          logging.warning('failed to threshold histogram',extra=self.d)
          return data
        
    def readNetCdf(self,nb,snowType,p,doy,filename=None):
        """

        read mean and var/covar of MODIS albedo datasets to file (NetCDF format)
        filenames are of form
        OPDIR + '/' + 'Kernels.' + doy + '.' + version + '.' + tile  + '.' + Snowtype

        (or override name structure with filename option)

        Parameters
        ----------

        nb : integer
             number of bands
        snowType : string
               e.g. SnowAndNoSnow, used in filename
        p        : int -- part of image
        filename : string
               Override for defining the filename to read
        doy   : string
               dot string e.g. 001   
 
        Returns
        -------- 
 
        mean : float array
               mean array
        sd : float array
               sd  array
        n : float array
               sample weight
        l : float array
               land mask

        """
        p = '%02d'%p

        filename = filename or self.ip.opdir + '/Kernels.' + '%03d'%int(doy) + '.' + self.ip.version + '.' +\
                                 self.ip.tile + '.' + snowType +'.'+ p + '.nc'
        logging.info( 'reading %s'%filename,extra=self.d)
        try:
          ncfile = nc.Dataset(filename,'r')
        except:
          # maybe there is a compressed version?
          try:
            from subprocess import call
            call(['gzip','-df',filename+'.gz'])
          except:
             logging.info( "Failed to uncompress output file %s"%filename,extra=self.d)
        try:
          ncfile = nc.Dataset(filename,'r')
        except:
          logging.warning("Failed to read Netcdf file %s"%filename,extra=self.d)
          
        bNames = self.idSamples(nb)
        meanNames = np.array(bNames)[np.where(['MEAN: ' in i for i in bNames])].tolist()
        sdNames = np.array(bNames)[np.where(['SD: ' in i for i in bNames])].tolist()
        nNames = ['Weighted number of samples']
        lNames = ['land mask']
        mean = []
        for b in meanNames:
          mean.append(ncfile.variables[b][:])   
        sd = []
        for b in sdNames:
          sd.append(ncfile.variables[b][:])
        n = ncfile.variables[nNames[0][:]] 
        l = ncfile.variables[lNames[0][:]]
        return np.array(mean),np.array(sd),np.array(n),np.array(l)


    def writeNetCdf(self,ns,nl,nb,mean,sd,n,land,snowType,p,doy):
        """

        write mean and var/covar of MODIS albedo datasets to file (NetCDF format)
        filenames are of form
        OPDIR + '/' + 'Kernels.' + doy + '.' + version + '.' + tile  + '.' + Snowtype

        If self.ip.compression is set, an attenpt is made to compress the netCDF file.

        Parameters
        ----------

        ns : integer
             number of samples
        nl : integer
             number of lines 
        nb : integer
             number of bands
        mean : float array
               mean array
        sd : float array
               sd  array
        n : float array
               sample weight
        land : float array
               land mask
        snowType : string
               e.g. SnowAndNoSnow, used in filename
        p : integer
             image part number
        doy  : string
             DOY string e.g. 009

        Returns
        -------

        None 

        """

        p = '%02d'%p
        shrink = self.ip.shrink
        doy = '%03d'%int(doy)

        filename = self.ip.opdir + '/Kernels.' + doy + '.' + self.ip.version + '.' +\
                                 self.ip.tile + '.' + snowType +'.'+ p + '.nc'

        bNames = self.idSamples(nb)
        bNames.append('Weighted number of samples')
        bNames.append('land mask')
     
        logging.info( 'writing %s'%filename,extra=self.d)
        if nb == 2:
            defBands = [4,1,1]
        elif nb == 7:
            defBands = [1,14,19]
        else:
            defBands = [1,10,7]
        descrip = snowType + ' MODIS Mean/SD ' + doy + ' over the years ' + str(self.yearList)  + \
            ' version ' + self.ip.version + ' tile '+ self.ip.tile + \
            ' using input MODIS bands '
        for band in self.ip.bands:
            descrip = descrip + str(band) + ' '

        ncfile = nc.Dataset(filename,'w',format = 'NETCDF4')
        ncfile.createDimension('ns',ns)
        ncfile.createDimension('nl',nl)
        
        count = 0
        for i in range(nb):
            for j in range(3):
                shrunkMean = mean[i,:,:,j]
                shrunkSd = sd[i,:,:,j]
                #shrunkMean,shrunkSd = self.shrunk(mean[i,:,:,j],ns,nl,shrink,sdata=sd[i,:,:,j])
                data = ncfile.createVariable(bNames[count],'f4',('ns','nl'),zlib=True,least_significant_digit=10)
                data[:] = shrunkMean
                count = count +1
                data = ncfile.createVariable(bNames[count],'f4',('ns','nl'),zlib=True,least_significant_digit=10)
                data[:] = shrunkSd
                count = count + 1
      
        shrunkN = n[0,:,:,0] 
        #shrunkN = self.shrunk(n[0,:,:,0],ns,nl,shrink)
        data = ncfile.createVariable(bNames[count],'f4',('ns','nl'),zlib=True,least_significant_digit=10)
        data[:] = shrunkN
        count = count + 1
        
        shrunkLand = land
        #shrunkLand = self.shrunk(land,ns,nl,shrink)
        data = ncfile.createVariable(bNames[count],'f4',('ns','nl'),zlib=True,least_significant_digit=2)
        data[:] = shrunkLand

        setattr(ncfile,'description',descrip)
        setattr(ncfile,'data ignore value',-1.0)
        setattr(ncfile,'default bands',defBands)
 
        try: 
          ncfile.close()
        except:
          pass
        if(False and (self.ip.compression == 1)):
          try:
            from subprocess import call
            call(['gzip','-f',filename])
          except:
            logging.info( "Failed to compress output file %s"%filename,extra=self.d)
    

    def runAll(self):
        if os.path.exists(self.ip.opdir) == 0:
            os.makedirs(self.ip.opdir)

        self.years = self.ip.years 
        if len(self.years) == 0:
          self.years = map(int,self.yearList)

        # Lewis: try some options on data location
        self.a2list = insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'2/*/%s/*/*hdf'%self.ip.tile)
        self.a1list = insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'1/*/%s/*/*hdf'%self.ip.tile)

        if len(self.a2list) == 0 or len(self.a1list) == 0:
          self.a2list = insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'2*hdf')
          self.a1list = insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'1*hdf')
        if len(self.a2list) == 0 or len(self.a1list) == 0:
          for year in self.years:
            self.a2list.extend(insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'2/%s/*.%s.*hdf'%(year,self.ip.tile)))
            self.a1list.extend(insensitive_glob(self.ip.srcdir+'/'+'*'+self.ip.product+'1/%s/*.%s.*hdf'%(year,self.ip.tile)))

        self.doyList, self.yearList = self.getDates(self.a1list)

        self.nDays = len(self.doyList)
        self.nYears = len(self.years)
         
        logging.info( 'N Years = %d; N Days = %d'%(self.nYears,self.nDays),extra=self.d)

        logging.info(str(self.doyList),extra=self.d)
        logging.info('file compression: %s'%str(self.ip.compression),extra=self.d)

        isRST = False
        if self.ip.rstdir != 'None':
          try:
            import pylab as plt
            isRST = True
            # try to make directory
            if os.path.exists(self.ip.rstdir) == 0:
              os.makedirs(self.ip.rstdir)
          except:
            logging.info("failed to load pylab or create required directory: required for rst output",extra=self.d)
            self.ip.rstdir = 'None'
            isRST = False

        #loop over each day of interest
        if (self.ip.doyList == [None]).all():
          doyList = self.doyList
        else:
          doyList = self.ip.doyList

        ok = True
        for doy in doyList:
          logging.info("DOY %s"%str(doy),extra=self.d)
          dimy=2400
          biny=dimy/self.ip.part

          # try reading the data from CDF files
          if self.ip.clean == False:
            for p in range(self.ip.part):
              nb = len(self.ip.bands)
              # try reading the data unless --clean
              try:
                if (self.ip.dontwithsnow == False):
                  mean, sdData, n, land = self.readNetCdf(nb,'SnowAndNoSnow',p,doy)
                  (nb,ns,nl) = mean.shape;nb/=3
                  # write RST file, but dont shrink the data again
                  if isRST:
                    self.writeRST(ns,nl,nb,mean,sdData,n,land,'SnowAndNoSnow',p,doy,shrink=1,order=1)
              except:
                ok = False
              # try reading the data unless --clean
              try:
                if (self.ip.dontsnow == False):
                  mean, sdData, n, land = self.readNetCdf(nb,'Snow',p,doy)
                  (nb,ns,nl) = mean.shape;nb/=3
                  if isRST:
                    self.writeRST(ns,nl,nb,mean,sdData,n,land,'Snow',p,doy,shrink=1,order=1)
              except:
                ok = False
              # try reading the data unless --clean
              try:
                if (self.ip.dontnosnow == False):
                  mean, sdData, n, land = self.readNetCdf(nb,'NoSnow',p,doy)
                  (nb,ns,nl) = mean.shape;nb/=3
                  if isRST:
                    self.writeRST(ns,nl,nb,mean,sdData,n,land,'NoSnow',p,doy,shrink=1,order=1)
              except:
                ok = False

          # don't do multiple parts as yet .. use --sdmims instead
          if not ok: 
           for p in range(self.ip.part):
            logging.info( "\n\n-------- partition %d/%d --------"%(p,self.ip.part),extra=self.d)
            y0=p*biny
            #sdimsL=[-1, -1, -1, -1]
            sdimsL = self.ip.sdims
            if self.ip.part>1:
              sdimsL[2]=y0
              sdimsL[3]=biny
            logging.info("sub region: %s"%str(sdimsL),extra=self.d)

            processed,totalSnow,sumDataNoSnow,sumDataSnow,sumDataWithSnow,nb,ns,nl,land = \
                                self.processAlbedo(self.years,[doy],sdimsL)
          
            if processed: 
              land = self.shrunk(land,ns,nl,self.ip.shrink)
              ns /= self.ip.shrink
              nl /= self.ip.shrink

              logging.info('... calculating stats',extra=self.d)
              n = np.zeros((nb,ns,nl,3),dtype=np.float32)
              mean = np.zeros((nb,ns,nl,3),dtype=np.float32)
              sdData = np.zeros((nb,ns,nl,3),dtype=np.float32)

              if (self.ip.dontwithsnow == False):
                n, mean, sdData = self.calculateStats(sumDataWithSnow)
                self.writeNetCdf(ns,nl,nb,mean,sdData,n,land,'SnowAndNoSnow',p,doy)
                if isRST:
                  self.writeRST(ns,nl,nb,mean,sdData,n,land,'SnowAndNoSnow',p,doy)
 
              if self.ip.dontnosnow == False:
                n, mean, sdData = self.calculateStats(sumDataNoSnow)  
                self.writeNetCdf(ns,nl,nb,mean,sdData,n,land,'NoSnow',p,doy)
                if isRST:       
                  self.writeRST(ns,nl,nb,mean,sdData,n,land,'NoSnow',p,doy)
 
              if self.ip.dontsnow == False:
                n, mean, sdData = self.calculateStats(sumDataSnow)
                self.writeNetCdf(ns,nl,nb,mean,sdData,n,land,'Snow', p,doy)
                if isRST:       
                  self.writeRST(ns,nl,nb,mean,sdData,n,land,'Snow',p,doy)


def processArgs(args=None,parser=None):

    usage = "usage: %prog [options]"

    #process input arguements
    parser = parser or OptionParser(usage=usage)
    prog = 'logger'

    # Lewis: add defaults here

    parser.add_option('--clean',dest='clean',action='store_true',default=False,\
                      help="ignore previous cdf files")
    parser.add_option('--dontclean',dest='clean',action='store_false',\
                      help="read previous cdf files if available (default action)")
    parser.add_option('--logfile',dest='logfile',type='string',default='%s.log'%(prog),\
                      help="set log file name") 
    parser.add_option('--logdir',dest='logdir',type='string',default='logs',\
                      help="set log directory name")
    parser.add_option('--rst',dest='rstdir',type='string',default='None',\
                      help="Write out an rst file with the processed data")
    parser.add_option('--srcdir',dest='srcdir',type='string',default='modis',\
                      help="Source (MODIS MCD43) data directory")
    parser.add_option('--tile',dest='tile',type='string',default='h18v03',\
                      help="MODIS tile ID")
    parser.add_option('--backupscale',dest='backupscale',type='string',default='[1.,0.7,0.49,0.343]',\
                      help="Array defining the scale to map MODIS QA flags to, e.g. [1.,0.7,0.49,0.343]")
    parser.add_option('--opdir',dest='opdir',type='string',default='results',\
                      help="Output directory")
    parser.add_option('--compression',dest='compression',action='store_true',default=True,\
                      help='Compress output file')
    parser.add_option('--dontsnow',dest='dontsnow',action='store_true',default=True,\
                      help="Don't output snow prior")
    parser.add_option('--dontwithsnow',dest='dontwithsnow',action='store_true',default=False,\
                      help="Dont output prior with snow and no snow data. Default False (i.e. do process this)") 
    parser.add_option('--dontnosnow',dest='dontnosnow',action='store_true',default=True,\
                      help="Don't output prior with no snow (i.e. snow free) data")
    parser.add_option('--nocompression',dest='compression',action='store_false',\
                      help="Don't compress output file")
    parser.add_option('--snow',dest='dontsnow',action='store_false',\
                      help="Do output snow prior")
    parser.add_option('--withsnow',dest='dontwithsnow',action='store_false',\
                      help="Do output prior with snow and no snow data")
    parser.add_option('--nosnow',dest='dontnosnow',action='store_false',\
                      help="Do output prior with no snow (i.e. snow free) data")
    parser.add_option('--shrink',dest='shrink',type='int',default=1,\
                      help="Spatial shrink factor (integer: default 1)")
    parser.add_option('--sdims',dest='sdims',type='string',default='[-1,-1,-1,-1]',\
                      help='image subsection: default [-1,-1,-1,-1]i or [l0,nl,s0,ns]')
    parser.add_option('--bands',dest='bands',type='string',default='[7,8,9]',help='list of bands to process. Default [7,8,9]')
    parser.add_option('--years',dest='years',type='string',default=\
                              '[2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013]',help='list of years to process')
    parser.add_option('--focus',dest='focus',type='string',default='[1.,1,1,1,1,1,1,1,1,1,1,1,1]',help='weighting for each year used')
    parser.add_option('--version',dest='version',type='string',default='005',help='MODIS collection number (as string). Default 005')
    parser.add_option('--product',dest='product',type='string',default='mcd43a',help='product name (default mcd43a)')
    parser.add_option('--part',dest='part',type='int',default=1,help='How many parts to split the tile into (default 1)')
    parser.add_option('--rstbase',dest='basename',type='string',default='Kernels',help='root name for ReST files')
    parser.add_option('--doy',dest='doyList',type='string',default='None',help='list of doys to process e.g. "[1,9]". N.B. do not put leading zeros on the dates (e.g. do not use e.g. [001,009])')

    
    return parser.parse_args(args or sys.argv)

if __name__ == '__main__':
    
    opts, args = processArgs()

    # Lewis: need str at times on the eval 
    opts.backupscale = np.array(ast.literal_eval(opts.backupscale))
    opts.sdims       = np.array(ast.literal_eval(opts.sdims))
    opts.bands       = np.array(ast.literal_eval(opts.bands))
    opts.years       = np.array(ast.literal_eval(opts.years))
    opts.focus       = np.array(ast.literal_eval(opts.focus))
    opts.doyList = np.array(ast.literal_eval(opts.doyList)) 

    #class call
    albedo = albedo_pix(opts)

    albedo.runAll()

        
