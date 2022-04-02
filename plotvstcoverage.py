import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import numpy.ma as ma
import sys
import astropy.io.fits as pf
from astropy.table import vstack
from astropy.table.table import Table
import pyvo
from pyvo.dal import tap
from datetime import datetime,timedelta


class VSTVISTAcoverage:
    
    '''
    Plot 9-band coverage by VST and VISTA in given region of sky
    Methods:
    On initialisation load the archive of VST/VISTA observations.
      with argument newdb=True, always get new info from ESO archive,
          otherwise look for local files first
    plotmaps(ra,dec,radius) then performs the search 
    It returns a 9-panel plot of total exposure.
    '''

    def __init__(self,newdb=False):
        '''
           If newdb is True, query ESO archive of observations
           Otherwise first check for existing files before querying.
           Query results are written out as fits tables.
           The resulting vstobs and vistaobs attributes are Tables
        '''

        def getobsfromeso(instrument='OMEGACAM'):

            print('Querying ESO archive for '+instrument+' observations...')
            ESO_TAP_OBS = "http://archive.eso.org/tap_obs"
            tapobs = tap.TAPService(ESO_TAP_OBS)
            
            query="""
            select dp_id,ob_name,object,prog_id,ra,dec,date_obs,exposure,
            filter_path,access_url
            from dbo.raw
            where dp_cat='SCIENCE' and instrument='"""+instrument+"'"
            
            obs = tapobs.search(query=query,maxrec=3000000).to_table()
            
            print('Obtained '+instrument+' observations from TAP query to ESO')
            # turn 'object' data types into strings to fitsio can deal with them
            for col in obs.columns:
                if obs[col].dtype=='O':
                    obs.replace_column(col,obs[col].astype('str'))
            obs.sort(keys='date_obs')
            fitsname=instrument+'.fits'
            if instrument=='OMEGACAM':   fitsname='vst_obs.fits'
            if instrument=='VIRCAM':     fitsname='vista_obs.fits'
            obs.write(fitsname,overwrite=True)
            print('Observations saved to fits table',fitsname)
            return obs

        if newdb:
            self.vstobs=getobsfromeso(instrument='OMEGACAM')
            self.vistaobs=getobsfromeso(instrument='VIRCAM')
        else:
            try:
            # first check if have a table downloaded
                self.vstobs=Table(pf.open('vst_obs.fits')[1].data)
                self.vistaobs=Table(pf.open('vista_obs.fits')[1].data)
                print('Read fits tables with VST and VISTA observations.')
            except:
                print('No fits tables of VST and VISTA observations found.')
                self.vstobs=getobsfromeso(instrument='OMEGACAM')
                self.vistaobs=getobsfromeso(instrument='VIRCAM')
    
        # Now sort both tables by date_obs
#        self.vstobs.sort(keys='date_obs')
#        self.vistaobs.sort(keys='date_obs')
        # turn dates into datetime objects, add as extra DATE column
        # note that VISTA date_obs has 4 sig figs in seconds, should be 3 or 6
        # this breaks fromisoformat()
        # so here only use first 3 decimals
        for v in [self.vstobs,self.vistaobs]:
            v.add_column([datetime.fromisoformat(x[:23])
                for x in v['date_obs']],name='DATE')
        self.firstdate=min(self.vstobs['DATE'][0],
                               self.vistaobs['DATE'][0]).date()
        self.lastdate=max(self.vstobs['DATE'][-1],
                              self.vistaobs['DATE'][-1]).date()
        print('Observations span',self.firstdate,'to',self.lastdate)
        # set reference level for exposure times to KiDS/VIKING values
        self.vstexposure={'u':1000.,'g':900.,'r':1800.,'i':1200.}
        self.vistaexposure={'Z':480.,'Y':400.,'J':400.,'H':300.,'K':480.}
        print('Reference exposure times set to KiDS/VIKING')
        self.cmap='Spectral'

    def plotmaps(self,ra0,dec0,radius=3.,label='',plotcircle=True):
        '''
        Make a 9-panel plot of exposure maps of all VST and VISTA observations
        in an area centered on ra0,dec0, with given radius (all in deg)
        The label is used for the title of the plot, and can be the field name
        if plotcircle is True, a circle of the specified radius is overplotted
        '''
        binperdeg=100.
        ra=self.vstobs['ra']
        ra=(ra-ra0+540) % 360 -180 +ra0     # wrap ra to as close to ra0 as possible
        de=self.vstobs['dec']
        cosfac=1/np.cos(dec0*np.pi/180)
        dra=0.5*cosfac  # 0.5 deg margins
        dde=0.5
        ralo,rahi=ra0-radius*cosfac-dra, ra0+radius*cosfac+dra
        delo,dehi=dec0-radius-dde,dec0+radius+dde
        nra=int(2*radius*binperdeg)+1
        nde=int(2*radius*binperdeg)+1
        ra1d=np.linspace(ralo,rahi,nra)
        de1d=np.linspace(delo,dehi,nde)
        x,y=np.meshgrid(ra1d,de1d)
        # use all observations with ctr in field or within 0.5 deg of it
        use=(ra>ralo-dra) & (ra<rahi+dra) & (de>delo-dde) & (de<dehi+dde)
        ra=ra[use]
        de=de[use]
        texp=self.vstobs['exposure'][use]
        filt=self.vstobs['filter_path'][use]
        self.exposuremap={}
        for filter in ['u','g','r','i']:
            z=np.zeros((nde,nra))
            for i in range(len(ra)):
                if filt[i][0]==filter.upper():
                    inpix=(x>ra[i]-dra) & (x<ra[i]+dra) & (y>de[i]-dde) & (y<de[i]+dde) 
                    z[inpix]+=texp[i]
                    plt.plot(ra[i],de[i],'ko',markersize=100./nde)
            plt.title(label+' '+filter)
            plt.xlabel('RA')
            plt.ylabel('DEC')
            mexp=self.vstexposure[filter]
            plt.imshow(ma.masked_equal(z,0.),origin='lower',
                           extent=(ralo,rahi,delo,dehi),aspect=cosfac,
                           norm=colors.LogNorm(vmin=mexp/30,vmax=mexp*30),
                           cmap=self.cmap)
            if plotcircle:
                th=np.linspace(0,2*np.pi,100)
                xcirc=np.cos(th)*radius*cosfac+ra0
                ycirc=np.sin(th)*radius+dec0
                plt.plot(xcirc,ycirc,'b-')
            plt.colorbar(label='Exposure time [s]')
            plt.xlim(rahi,ralo)
            plt.ylim(delo,dehi)
            plt.savefig(label+'_'+filter+'.png')
            plt.clf()
            self.exposuremap[filter]=z
        ra=self.vistaobs['ra']
        ra=(ra-ra0+540) % 360 -180 +ra0     # wrap ra to as close to ra0 as possible
        de=self.vistaobs['dec']
        cosfac=1/np.cos(dec0*np.pi/180)
        dra=0.5*cosfac  # 0.5 deg margins
        dde=0.5
        ralo,rahi=ra0-radius*cosfac-dra, ra0+radius*cosfac+dra
        delo,dehi=dec0-radius-dde,dec0+radius+dde
        nra=int(2*radius*binperdeg)+1
        nde=int(2*radius*binperdeg)+1
        ra1d=np.linspace(ralo,rahi,nra)
        de1d=np.linspace(delo,dehi,nde)
        x,y=np.meshgrid(ra1d,de1d)
        # use all observations with ctr in field or within 0.5 deg of it
        use=(ra>ralo-dra) & (ra<rahi+dra) & (de>delo-dde) & (de<dehi+dde)
        ra=ra[use]
        de=de[use]
        texp=self.vistaobs['exposure'][use]
        filt=self.vistaobs['filter_path'][use]
        for filter in ['Z','Y','J','H','K']:
            z=np.zeros((nde,nra))
            for i in range(len(ra)):
                if filt[i][0]==filter.upper():    # assume 1.5x1 deg VISTA tile is aligned East-West
                    inpix=(x>ra[i]-dra*1.5) & (x<ra[i]+dra*1.5) & (y>de[i]-dde) & (y<de[i]+dde) 
                    z[inpix]+=texp[i] * (2/3)   # Allow for dither losses in standard VISTA tiling
                    plt.plot(ra[i],de[i],'ko',markersize=100./nde)
            plt.title(label+' '+filter)
            plt.xlabel('RA')
            plt.ylabel('DEC')
            mexp=self.vistaexposure[filter]
            plt.imshow(ma.masked_equal(z,0.),origin='lower',
                           extent=(ralo,rahi,delo,dehi),aspect=cosfac,
                           norm=colors.LogNorm(vmin=mexp/30,vmax=mexp*30),
                           cmap=self.cmap)
            if plotcircle:
                th=np.linspace(0,2*np.pi,100)
                xcirc=np.cos(th)*radius*cosfac+ra0
                ycirc=np.sin(th)*radius+dec0
                plt.plot(xcirc,ycirc,'b-')
            plt.colorbar(label='Exposure time [s]')
            plt.xlim(rahi,ralo)
            plt.ylim(delo,dehi)
            plt.savefig(label+'_'+filter+'.png')
            plt.clf()
            self.exposuremap[filter]=z
