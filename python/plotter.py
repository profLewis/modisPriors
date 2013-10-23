import netCDF4 as nc
import numpy as np
import pylab as plt
import os

'''
Data plotting utility: class plotter
'''
 
inimage = 'results/Kernels.153.005.h18v03.background.Snow.nc'
opdir = 'source'

class plotter():
  '''
  Data plotter utility
  '''
  def __init__(self,inimage,n=0,plotData=False,plotStats=False,plotHist=False):
    '''
    Initialise and run plotting routines
    '''
    self.plotter(inimage,n=n,plotData=plotData,plotStats=plotStats,plotHist=plotHist)


  def fix(self,data,iter=0):
    '''
    Threshold dataset at 99% cumulative sum for easier display
    '''
    if iter > 10 or (data>0).size == 0 or np.std(data[data>0]) == 0:
      return data
    try:
      fdata = data.flatten()
      (n,b) = np.histogram(fdata[fdata>0])
      nc = n.cumsum()
      N = nc[-1]*0.99
      w = np.where(nc >= N)[0][0]
      if w == 0:
        w = 1
        max = b[w]
        d = data.copy()
        d[d>=max]=max
        return self.fix(d,iter=iter+1)
      max = b[w]
      d = data.copy()
      d[d>=max]=max
      if np.std(d) == 0:
        return data
      return d
    except:
      return data
  
  def plotter(self,inimage,opdir=None,n=None,count=0,pause=False,plotStats=False,plotData=False,plotHist=False):
    '''
    Read datafile inimage and plot data
  
    If the data file is compressed, we uncompress it.
  
    Parameters
    ----------
  
    inimage : name of input image
    
    opdir : set to write plots to png files (don't use this with Pweave wrapper)
    n     : dataset number
    count : plot number
    pause : set True to step through a big loop (don't normally use)
    plotStats : plot stats if True
    plotData : plot data if True
    plotHist : plot histogram if True 
   
    '''
    tmpFile = None
    if opdir:
      if os.path.exists(opdir) == 0:
        os.makedirs(opdir)
    try:
      ncfile = nc.Dataset(inimage,'r')
    except:
      import tempfile,os
      tmpFile = tempfile.NamedTemporaryFile(delete=False).name
      cmd = "zcat %s > %s"%(inimage+'.gz',tmpFile)
      try:
        os.system(cmd)
        ncfile = nc.Dataset(tmpFile,'r')
      except:
        cmd = "cat %s > %s"%(inimage+'.gz',tmpFile)
        try:
          os.system(cmd)
          ncfile = nc.Dataset(tmpFile,'r')
        except:
          from sys import stderr,exit
          stderr.write("Error reading or uncompressing file %s\n"%inimage)
          stderr.write("cmd: file %s\n"%cmd)
          exit(0)
  
    descrip = ncfile.description
    bnames = ncfile.variables.keys()
    if n != None:
      bnames = [bnames[n]]
    for i in xrange(len(bnames)):
      try:
        data = (ncfile.variables[bnames[i]])[:]
        if(plotData):  self.dataPlot(data,bnames[i],count,opdir)
        if(plotHist):  self.histPlot(data,bnames[i],count,opdir)
        if(plotStats): self.statsPrint(data,bnames[i],count,opdir)
        #if pause:
          #import pdb;pdb.set_trace()
      except:
        pass
    try:
      ncfile.close()
    except:
      pass
    if tmpFile:
      import os
      os.remove(tmpFile)
  
  def histPlot(self,data,name,count,opdir,type=None):
    '''
    Plot a histogram of data
  
    Parameters
    -----------
  
    data : numpy dataset
    name : plot title
    count : plot number 
    opdir : None (no save) or filename to save png image to
    type: unused
  
    '''
    try:
      hist,bins = np.histogram(self.fix(data[data>0],iter=0),bins=50)
      width = 0.7*(bins[1]-bins[0])
      center = (bins[:-1]+bins[1:])/2
      plt.figure(count);plt.clf();count+=1
      plt.bar(center,hist,align='center',width=width)
      plt.title(name)
      if opdir:
        fig.savefig('%s/hist_%s_%d.png'%(opdir,name.replace(' ','_'),count))
      else:
        fig.show()
    except:
      pass
  
  def statsPrint(self,data,name,count,opdir,noZero=True,type=None):
    '''
    Print summary statistics
  
    Parameters
    -----------
  
    data : numpy dataset
    noZero : set True for discounting zero (default True)
    name : plot title (unused)
    count : plot number (unused)
    opdir : None (no save) or filename to save png image to (unused)
    type: unused
  
    '''
    print name + "\n"
    if noZero:
      data = data[data>0]
    print 'Mean: %.3f Median: %.3f Min: %.3f Max: %.3f SD: %.3f'%\
        (np.mean(data),np.median(data),np.min(data),np.max(data),np.std(data))
  
  def dataPlot(self,data,name,count,opdir,type=None):
    '''
    Plot data as image
  
    Parameters
    -----------
  
    data : numpy dataset
    name : plot title
    count : plot number 
    opdir : None (no save) or filename to save png image to
    type: unused
  
    '''
  
    plt.clf()
    plt.imshow(self.fix(data,iter=0),interpolation="nearest")
    plt.colorbar()
    plt.title(name)
    if opdir:
      plt.savefig('%s/%s_%d.png'%(opdir,name.replace(' ','_'),count))
    else:
      plt.show()
  
