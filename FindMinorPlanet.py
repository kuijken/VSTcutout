from astroquery.mpc import MPC
from astroquery.jplhorizons import Horizons
from sys import argv
import astropy.io.fits as pf
from astropy.table import vstack
from astropy.table.table import Table
import numpy as np
from datetime import datetime,timedelta

class MinorPlanetSearch:
    
    '''
    Look for minor planets that intersected VST and VISTA fields of view
    Methods:
    On initialisation load the archive of VST/VISTA observations.
      with argument newdb=True, always get new info from ESO archive,
          otherwise look for local files first
    search(designation) then performs the search by first getting an ephemeris
      and then comparing this with the observations archive.
      It returns a table of observations that might contain the target
    '''

    def __init__(self,newdb=False):
        '''
           If newdb is True, query ESO archive of observations
           Otherwise first check for existing files before querying.
           Query results are written out as fits tables.
           The resulting vstobs and vistaobs attributes are Tables
        '''

        def getobsfromeso(instrument='OMEGACAM'):
            import pyvo
            from pyvo.dal import tap
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
                print('Read fits tables with VST observations.')
            except:
                print('No fits tables of VST observations found.')
                self.vstobs=getobsfromeso(instrument='OMEGACAM')
            try:
            # first check if have a table downloaded
                self.vistaobs=Table(pf.open('vista_obs.fits')[1].data)
                print('Read fits tables with VISTA observations.')
            except:
                print('No fits tables of VISTA observations found.')
                self.vistaobs=getobsfromeso(instrument='VIRCAM')
    
        # Now sort both tables by date_obs
        self.vstobs.sort(keys='date_obs')
        self.vistaobs.sort(keys='date_obs')
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

        self.paranal,self.vistapeak='309','W91'     # IAU obs location codes

    def search(self,asteroid:str,ephemservice='MPC',id_type='smallbody'):
        '''
            Search for object with MPC ID and return possible observations
            ephemservice can be 'MPC' or 'JPL'
        '''

        if (ephemservice=='MPC') & (id_type == 'smallbody'):
            d=self.firstdate
            Neph=1000  # can request a maximum of 1441 dates in MPEC
            eph=[]
            while d < self.lastdate:
                eph.append(
                    MPC.get_ephemeris(asteroid, location=self.paranal, start=str(d),
                                    step='1d', number=Neph)
                    )
                # fix sporadic formatting error in MPC
                if eph[-1]['Direction'].dtype != 'float64':  
                    eph[-1].replace_column('Direction',[float(x[-5:]) for x in eph[-1]['Direction']])
                d+=timedelta(Neph)
            eph=vstack(eph)
            eph.write(asteroid+'.eph.fits',overwrite=True)
            print('MPC empemerides written to',asteroid+'.eph.fits')
            eph.add_column([x.to_datetime() for x in eph['Date']],name='DATE')
        elif ephemservice=='JPL':
            eph=Horizons(id=asteroid, location=self.paranal, epochs={
                'start':self.firstdate.isoformat(),
                'stop': self.lastdate.isoformat(),
                'step': '1d'},id_type=id_type).ephemerides()
            eph.write(asteroid+'.eph.fits',overwrite=True)
            print('JPL empemerides written to',asteroid+'.eph.fits')
            eph.add_column([datetime.strptime(x,'%Y-%b-%d %H:%M')
                        for x in eph['datetime_str']],name='DATE')
            eph.rename_column('DEC','Dec')
            if 'V' in eph.columns:
                pass
            elif 'mag' in eph.columns:
                eph.rename_column('mag','V')
            elif 'Tmag' in eph.columns:
                eph.rename_column('Tmag','V')
            else:
                print('Found no mag in JPL ephemeris columns:',eph.columns)
                eph.add_column(0*eph['RA'],name='V')
        else:
            print('No valid ephem service given: MPC [default] or JPL')
            return [99,99]

        # Interpolate where the asteroid was at the time of every observation
        refdate=datetime(2000,1,1)
        ephtime=[(x-refdate).total_seconds() for x in eph['DATE']]
        #first check VST
        obstime=[(x-refdate).total_seconds() for x in self.vstobs['DATE']]
        ephra= np.interp(obstime,ephtime,eph['RA'])
        ephdec=np.interp(obstime,ephtime,eph['Dec'])
        obsra= self.vstobs['ra']
        obsdec=self.vstobs['dec']
        dra= (obsra-ephra) * np.cos(ephdec*np.pi/180)
        ddec=obsdec-ephdec
        infield=(dra**2 < 0.55**2) & (ddec**2<0.55**2)
        invst=self.vstobs[infield]
        tobs=np.array(obstime)[infield]
        invst.add_column(ephra[infield],name='RAtarg')
        invst.add_column(ephdec[infield],name='DECtarg')
        invst.add_column(dra[infield]*np.cos(obsdec[infield]*np.pi/180),name='dRA')
        invst.add_column(ddec[infield],name='dDEC')
        invst.add_column(np.round(np.interp(tobs,ephtime,eph['V']),1),name='Vmag')
        print()
        print(len(invst),'VST matches:')
        print(invst)

        #then check VISTA
        obstime=[(x-refdate).total_seconds() for x in self.vistaobs['DATE']]
        ephra= np.interp(obstime,ephtime,eph['RA'])
        ephdec=np.interp(obstime,ephtime,eph['Dec'])
        obsra= self.vistaobs['ra']
        obsdec=self.vistaobs['dec']
        dra= (obsra-ephra) * np.cos(ephdec*np.pi/180)
        ddec=obsdec-ephdec
        infield=(dra**2 < 0.8**2) & (ddec**2<0.8**2)
        invista=self.vistaobs[infield]
        tobs=np.array(obstime)[infield]
        invista.add_column(ephra[infield],name='RAtarg')
        invista.add_column(ephdec[infield],name='DECtarg')
        invista.add_column(dra[infield]*np.cos(obsdec[infield]*np.pi/180),name='dRA')
        invista.add_column(ddec[infield],name='dDEC')
        invista.add_column(np.round(np.interp(tobs,ephtime,eph['V']),1),name='Vmag')
        print()
        print(len(invista),'VISTA matches:')
        print(invista)

        invst.write(asteroid+'.vst.csv',format='csv',overwrite=True)
        invista.write(asteroid+'.vista.csv',format='csv',overwrite=True)
        print('Result written out to csv files',asteroid+'.vst.csv and .vista.csv')
        return invst,invista
