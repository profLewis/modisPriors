import netCDF4 as nc
import numpy as np
import pylab as plt
import os
from sys import stderr

baseStr0 = '''
Stage ZZZ results for YYY
------------------------------------------------------------------------------------------------------------------------------------------------


<<fig=True>>=
import pylab as plt
from plotter import plotter
inimage = "XXX"
plotter(inimage,n=NN,plotData=True)
@

<<results = 'rst', echo = False>>=
plotter(inimage,n=NN,plotStats=True)
@

<<fig=True>>=
plotter(inimage,n=NN,plotHist=True)
@

'''


baseStr = '''
<<fig=True>>=
plotter(inimage,n=NN,plotData=True)
@


<<results = 'rst', echo = False>>=
plotter(inimage,n=NN,plotStats=True)
@

<<fig=True>>=
plotter(inimage,n=NN,plotHist=True)
@

'''

baseStr1 = '''
<<fig=True>>=
plotter(inimage,n=NN,plotData=True)
@


<<results = 'rst', echo = False>>=
plotter(inimage,n=NN,plotStats=True)
@

'''




def nImages(inimage,srcdir='results'):
  inimage = srcdir + '/' + inimage
  try:
    tmpFile = None
    ncfile = nc.Dataset(inimage,'r')
    stderr.write("Read %s\n"%inimage)
  except:
    # try uncompressing it
    import os
    import tempfile
    tmpFile = tempfile.NamedTemporaryFile(delete=False).name
    cmd = "zcat %s > %s"%(inimage+'.gz',tmpFile)
    try:
      os.system(cmd)
      ncfile = nc.Dataset(tmpFile,'r')
      stderr.write("Read %s uncompressed to %s\n"%(inimage,tmpFile))
    except:
      #from sys import exit
      stderr.write("Error reading or uncompressing file %s.gz\n"%inimage)
      #stderr.write("cmd: file %s\n"%cmd)
      #exit(0)
      try:
        cmd = "cat %s > %s"%(inimage+'.gz',tmpFile)
        os.system(cmd)
        ncfile = nc.Dataset(tmpFile,'r')
        stderr.write("Read %s copied to %s\n"%(inimage,tmpFile))
      except:
        from sys import exit
        stderr.write("Error reading file %s\n"%inimage)
        stderr.write("cmd: file %s\n"%cmd)
        exit(0)


  bnames = ncfile.variables.keys()
  ncfile.close()
  return bnames,tmpFile

def genPage(inimage,srcdir='results',stage=1):
  import os
  ofile = inimage.replace(srcdir,'source').replace('.','_').replace('_nc','.Pnw')
  if os.path.exists('source') == 0:
    os.makedirs('source')
  bnames,tmp = nImages(inimage,srcdir=srcdir)
  doy = int(inimage.split('/')[-1].split('.')[1]) 
  f = open(ofile,'w')

  for i in xrange(len(bnames)):
    if i == 0:
      f.write(baseStr0.replace('ZZZ','%d'%stage).replace('NN',str(i))\
           .replace('XXX','%s/%s'%(srcdir,inimage))\
           .replace('YYY','DOY %03i'%int(doy)))
    elif i == len(bnames)-1:
      f.write(baseStr1.replace('NN',str(i)))
    else:
      f.write(baseStr.replace('NN',str(i)))
  if tmp:
    import os
    os.remove(tmp)

# generate index page
indexTxt1 = '''
Stage ZZZ processing TILE
====================================================================================

Contents:

.. toctree::
   :maxdepth: 4

   Introduction<intro>
   globAlbedo
   prior2
   albedo_pix
   genPnw
   plotter

'''

indexTxt2 = '''

Indices and tables
==================



* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


'''

def genIndex(srcdir,tile='h18v04',stage=1):
   import glob
   if stage == 1:
     where = srcdir+'/*.%s.*.??.nc*'%tile
   else:
     where = srcdir+'/*.%s.background.*.nc*'%tile

   stderr.write("searching in %s\n"%where)
   files = glob.glob(where)
   files = [i.split('/')[-1].replace(".gz","") for i in files]
   # write index
   stderr.write("files: %s\n"%str(files))
   file = 'index.rst'
   f = open(file,'w')
   f.write(indexTxt1.replace('TILE',tile).replace('ZZZ','%d'%stage))
   # 
   for file in np.sort(files):
     genPage(file,srcdir=srcdir,stage=stage) 
     f.write('   %s<%s>\n'%(file,file.replace(".nc","").replace(".","_")))
     stderr.write(" ... %s done.\n"%file)
   f.write(indexTxt2)
   f.close() 


#inimage = 'results/Kernels.001.005.h18v03.background.Snow.nc'
#genPage(inimage)    
