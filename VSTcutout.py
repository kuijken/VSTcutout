import astropy.wcs as wcs
import astropy.io.fits as pf
import numpy as np
import numpy.random as rnd
from scipy.ndimage.filters import uniform_filter1d
from pyextract import pysex
from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import os
import requests
import cgi

class VSTcutout:
    '''
    Extract pixel array cutout around a WCS coordinates. 
    Can work with raw OMEGACAM data, includes crude astrometric correction
       to work around pointing errors - accuracy ~ 5 arcsec
    '''
    def __init__(self,image,fixastrom='MEF',refcatname='GSC2.3'):
        '''
        Load a fits file <image> (can be URL) and its WCS.
        fixastrom determines whether to tweak the WCS.
         if fixastrom=True, then a source catalogue is extracted and compared
           to GSC2.3, and the mean astrometric offset is applied to world
           coordinates
         if fixastrom='MEF' then only do this for MEF files, not for
           single-extension files.
         any other value means NO astrometry correction is applied

        attributes are :

        XOFF,YOFF: (=0 if fixastrom=False) average offset image-ref in RA,DEC
        xoff,yoff: arrays of astrom offsets, one value per extension
        xofffit,yofffit: linear fit of xoff,yoff values vs RA,DEC
        images    list with all image extensions, (or 1 if not a MEF)
        wcs       list with the corresponding WCS
        over/underscanx/y to mark start of physical pixels
        X/Ylo/hi  are extreme coordinates of the image
        srccat    the source catalogues generated for astrometry
        srcXref   a 2D image of the crosscorrelation with the reference cat
        srcXrefbins[0],[1] corresponding bin edges in RA,Dec
        '''

        if (image[:4]=='http'): # if URL then first download it
            headers={}
            response=requests.get(image, stream=True, headers=headers)
            contentdisposition = response.headers.get('Content-Disposition')
            if contentdisposition == None:    # if file cannot be returned
                self.fitsfile=None
                return
            value, params = cgi.parse_header(contentdisposition)
            filename = params["filename"]
            if not os.path.exists(filename):
                if response.status_code == 200:
                    print('Downloading image',filename)
                    with open(filename, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=50000):
                            f.write(chunk)
            else:
                print('File',filename,'was already present')
        else:
            filename=image
        # if image is Z compressed, decompress it first (fitsio can't handle it)
        if (filename[-2:]=='.Z'):
            os.system('uncompress '+filename)
            filename=filename[:-2]
        f=pf.open(filename)
        self.fitsfile=filename
        self.MEF=len(f)>1
        if self.MEF:
            self.images=f[1:]
            # assume under/overscans are same for all extensions !
            # define start and end indices of physical pixel area
            self.underscanx=f[1].header['HIERARCH ESO DET OUT1 PRSCX']
            self.overscanx=f[1].header['HIERARCH ESO DET OUT1 NX']+self.underscanx
            self.underscany=f[1].header['HIERARCH ESO DET OUT1 PRSCY']
            self.overscany=f[1].header['HIERARCH ESO DET OUT1 NY']+self.underscany
        else:
            # if single extension fits file, assume there is no overscan
            self.images=f
            self.underscanx=0
            self.overscanx=f[0].data.shape[1]
            self.underscany=0
            self.overscany=f[0].data.shape[0]
        self.wcs=[wcs.WCS(ext.header) for ext in self.images]
        #get corners of the image
        xlo,xhi,ylo,yhi=10000,-10000,10000,-10000
        for w in self.wcs:
            if not w.is_celestial:
                print('*****WCS is not a celestial system - abort')
                self.fitsfile=None
                f.close()
                return
            xcorner,ycorner=w.calc_footprint().T
            xlo=min(xlo,xcorner.min())
            xhi=max(xhi,xcorner.max())
            ylo=min(ylo,ycorner.min())
            yhi=max(yhi,ycorner.max())
        # check for crossing RA=0, assume this happened if width>180deg
        if xhi-xlo > 180: # redo RA edges after moving RA cut to 180 deg
            xlo,xhi=10000,-10000
            for w in self.wcs:
                xcorner=w.calc_footprint()[:,0]
                xcorner=(xcorner+540) % 360  - 180   # cut RA coords at 180
                xlo=min(xlo,xcorner.min())
                xhi=max(xhi,xcorner.max())
        self.Xlo=xlo
        self.Xhi=xhi
        self.Ylo=ylo
        self.Yhi=yhi
        ra=(xhi+xlo)/2
        dec=(yhi+ylo)/2
        #search GSC catalogue for overplotting
        Vizier.ROW_LIMIT=50000
        r=Vizier.query_region(coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg)),
                    width=(xhi-xlo)*u.deg,height=(yhi-ylo)*u.deg,
                    catalog=refcatname)[0]
        print(len(r),'GSC stars found.')
        self.refstars=r
        if (fixastrom==True) | ((fixastrom=='MEF') & self.MEF):
            # make world coordinate catalogue using pyextractor,
            # and crosscorrelate with a reference catalogue.
            # record pixel offsets as self.XOFF, self.YOFF
            # Note pysex does not run on compressed fits files, so decompress first.
            #    (make separate copy and delete it after SExtractor has run)
            if filename[-3:]=='.fz':
                print('Funpacking the image to run SExtractor')
                os.system('funpack '+filename)
                fitsfile=filename[:-3]
            elif filename[-2:]=='.Z':
                print('Uncompressing the image to run SExtractor')
                os.system('gunzip -k '+filename)
                fitsfile=filename[:-2]
            elif filename[-3:]=='.gz':
                print('Gunzipping the image to run SExtractor')
                os.system('gunzip -k '+filename)
                fitsfile=filename[:-3]
            else:
                fitsfile=filename
            self.srccat=pysex.run(imageref=fitsfile,
                                 params=['X_WORLD','Y_WORLD'],
                                 conf_args={'DETECT_THRESH':10})
            if filename != fitsfile:
                os.system('rm -f '+fitsfile)
            Noff=401   # binning for course astrometric offset search
            box=0.1  # ± range of offsets to probe
            h=np.zeros((Noff,Noff))
            for ext in range(len(self.images)):
                src=self.srccat[2*ext+2].data
                xw=src['X_WORLD']  # source positions in extension ext
                yw=src['Y_WORLD']
                dx=np.array([x-self.refstars['RAJ2000'] for x in xw]).flatten()
                dy=np.array([y-self.refstars['DEJ2000'] for y in yw]).flatten()
                rx=rnd.randn(len(dx)) / Noff * box *2
                ry=rnd.randn(len(dx)) / Noff * box *2
                hh,hbins=np.histogramdd((dx+rx,dy+ry),bins=Noff,
                                         range=[[-box,box],[-box,box]])
                h+=hh
            #plt.imshow(h.T,extent=[-box,box,-box,box],origin='lower')
            #plt.xlabel('Delta(RA) [deg]')
            #plt.ylabel('Delta(DEC) [deg]')
            #plt.title('Astrometry Xcorr with '+refcatname)
            #plt.savefig(fitsfile+'_astrom_wide.png')
            #plt.clf()
            self.srcXref=h
            self.srcXrefbins=hbins
            pk=np.unravel_index(h.argmax(),h.shape)
            xoff=hbins[0][0]+(hbins[0][1]-hbins[0][0])*(pk[0]+0.5)
            yoff=hbins[1][0]+(hbins[1][1]-hbins[1][0])*(pk[1]+0.5)
            print('Offsets found for RA, DEC:',xoff,yoff,'deg [rough]')
            # Refined crosscorrelation with smaller bins in smaller range
            # measure peak as mean of histogram above median after adding
            # Gaussian noise to smooth distribution
            Noff=25   # binning for fine astrometric offset search
            box=0.002  # refined plus/min range of offsets to probe (deg)
            Nrand=20
            h=np.zeros((Noff,Noff))
            self.xoff=np.zeros(len(self.images))
            self.yoff=np.zeros(len(self.images))
            badoff=np.zeros(len(self.images),dtype=bool)
            for ext in range(len(self.images)):
                src=self.srccat[2*ext+2].data
                xw=src['X_WORLD']  # source positions in extension ext
                yw=src['Y_WORLD']
                dx=np.array([x-self.refstars['RAJ2000'] for x in xw]).flatten()
                dy=np.array([y-self.refstars['DEJ2000'] for y in yw]).flatten()
                keep=(np.abs(dx-xoff)<box*1.5) & (np.abs(dy-yoff)<box*1.5)
                dx=dx[keep]
                dy=dy[keep]
                for irand in range(Nrand):  # smooth histogram by adding noise
                    rx=rnd.randn(len(dx)) / Noff * box *2
                    ry=rnd.randn(len(dx)) / Noff * box *2
                    hh,hbins=np.histogramdd((dx+rx-xoff,dy+ry-yoff),bins=Noff,
                                         range=[[-box,box],[-box,box]])
                    h+=hh/Nrand
                xoffext,yoffext=iterweightedctr(dx,dy,1./3600)
                self.xoff[ext]=xoffext # offset image-ref for this ext [deg]
                self.yoff[ext]=yoffext
                badoff[ext]=False   # flag that can be used to indicate poor offsets
            plt.imshow(h.T,extent=[xoff-box,xoff+box,yoff-box,yoff+box],
                           origin='lower')
            plt.xlabel('Delta(RA) [deg]')
            plt.ylabel('Delta(DEC) [deg]')
            plt.title('Astrometry Xcorr with '+refcatname+' (zoom)')
            plt.colorbar(label='Ref stars / offset bin')
            plt.savefig(fitsfile+'_astrom_zoom.png')
            plt.clf()
            self.srcXrefzoom=h
            self.srcXrefzoombins=hbins
            pk=np.unravel_index(h.argmax(),h.shape)
            xoff += hbins[0][0]+(hbins[0][1]-hbins[0][0])*(pk[0]+0.5)
            yoff += hbins[1][0]+(hbins[1][1]-hbins[1][0])*(pk[1]+0.5)
            self.XOFF=xoff   # best-fit offset image-reference, in deg.
            self.YOFF=yoff
            self.xoff[badoff]=self.XOFF  # fix poorly determined offsets to global value
            self.yoff[badoff]=self.YOFF
            self.refcatname=refcatname
            print('Offsets found for RA, DEC:',xoff,yoff,'deg')
            # Make a fit of the offsets vs field position -
            #   eg in OMEGACAM there is an indication for a small rotation,
            #   and a linear fit of xoff, yoff vs wcs[ext].wcs.crpix[0]
            #   and ..[1] can take this out.
            # also provides a way to interpolate over poorly determined chips.
            # Note differential refraction can be corrected with a linear fit.
            # provide this fit as self.xoff_fit, self.yoff_fit
            xy00=np.array([w.wcs_pix2world(0,0,0) for w in self.wcs])  # lower left corners of each ext
            A=np.vstack([np.ones(len(xy00)), xy00[:,0],xy00[:,1]]).T
            cxoff=np.linalg.lstsq(A,self.xoff,rcond=None)[0]
            cyoff=np.linalg.lstsq(A,self.yoff,rcond=None)[0]
            modxoff=np.dot(A,cxoff)
            modyoff=np.dot(A,cyoff)
            self.xofffit=modxoff
            self.yofffit=modyoff
        else:
            self.XOFF=0
            self.YOFF=0
            self.xoff=np.zeros(len(self.images))
            self.yoff=np.zeros(len(self.images))
            self.xofffit=np.zeros(len(self.images))
            self.yofffit=np.zeros(len(self.images))
        print('RA range %8.4f %8.4f ; DEC range %8.4f %8.4f' % (self.Xlo,self.Xhi,self.Ylo,self.Yhi))


    def cutout(self,ra,dec,width=100,label='',savefile='',offset='avg'):
        '''
        cut out width x width pixels around (ra,dec) coordinate, 
            scanning all FITS extensions
        preserve orientation of the pixels (no interpolations).
        apply offsets to all world coordinates (only if VSTcutout was called with fixastrom=True):
            offset='avg' : XOFF,YOFF   (average over all extensions)
            offset='ext' : xoff,yoff   (measured value per extensions)
            offset='fit' : xofffit,yofffit   (linear fit to measured value per extensions)
            offset= anything else: just use header WCS directly without correction
        No coaddition, simple projection of any pixels that cover cutout area
        Makes a figure:
            label is used in the title of the plot (eg source name)
            if savefile is given, use that as filename, otherwise display 
        The cutout pixels are returned by the method as a numpy array.
        '''

        if self.fitsfile == None:
            print('No fits file to make cutout from!')
            return np.zeros((0,0))
        
        out=np.zeros((width,width))
        if (ra+self.XOFF<self.Xlo) | (ra+self.XOFF>self.Xhi) | \
                      (dec+self.YOFF<self.Ylo) | (dec+self.YOFF>self.Yhi):
            print('coordinate %8.4f %8.4f falls outside the image - no cutout!'
                      % (ra,dec))
            return out
        xc,yc=width//2,width//2
        xref=np.zeros((0))
        yref=np.zeros((0))
        if offset=='ext':
            print('offset parameter is ext: so using WCS with per-extension shifts')
        elif offset=='fit':
            print('offset parameter is fit: so using linear fit to per-extension WCS offsets')
        elif offset=='avg':
            print('offset parameter is avg: so using WCS with best-fit global shift')
        else:
            print('offset parameter is not ext,fit,avg, so using WCS without correction')
        for im in range(len(self.images)):
            w=self.wcs[im]
            ny,nx=w.array_shape
            if offset=='ext':
                x,y=w.wcs_world2pix(ra+self.xoff[im],dec+self.yoff[im],1)
            elif offset=='fit':
                x,y=w.wcs_world2pix(ra+self.xofffit[im],dec+self.yofffit[im],1)
            elif offset=='avg':
                x,y=w.wcs_world2pix(ra+self.XOFF,dec+self.YOFF,1)
            else:
                x,y=w.wcs_world2pix(ra,dec,1)
            ix,iy=int(np.round(x)),int(np.round(y))
            #min, max pixels that need to be copied to cutout
            xlo=min(max(self.underscanx,ix-xc),self.overscanx)
            ylo=min(max(self.underscany,iy-yc),self.overscany)
            xhi=max(min(self.overscanx,ix-xc+width),self.underscanx)
            yhi=max(min(self.overscany,iy-yc+width),self.underscany)
            if self.MEF:  # assume MEF files are raw and need debiasing
                bias=np.median(self.images[im].data[ylo:yhi,:self.underscanx-5],axis=1).reshape(yhi-ylo,1)
                bias=uniform_filter1d(bias,size=21,mode='nearest')
            else:
                bias=0
            # copy raw pixels; could also include overscan correction, FF, BG subtraction here.
            out[ylo-iy+yc:yhi-iy+yc,xlo-ix+xc:xhi-ix+xc]=self.images[im].data[ylo:yhi,xlo:xhi] - bias
            if (xlo<xhi) & (ylo<yhi):   # if thumbnail on this fits extension
                # then get ref.stars in pixel coordinates using this WCS
                if offset=='ext':
                    refra=self.refstars['RAJ2000']+self.xoff[im]
                    refdec=self.refstars['DEJ2000']+self.yoff[im]
                elif offset=='fit':
                    refra=self.refstars['RAJ2000']+self.xofffit[im]
                    refdec=self.refstars['DEJ2000']+self.yofffit[im]
                elif offset=='avg':
                    refra=self.refstars['RAJ2000']+self.XOFF
                    refdec=self.refstars['DEJ2000']+self.YOFF
                else:
                    refra=self.refstars['RAJ2000']
                    refdec=self.refstars['DEJ2000']
                xref,yref=w.wcs_world2pix(refra,refdec,1)
                xref-= ix-xc
                yref-= iy-yc
        # set scale limits from image statistics - 2x interquartile range
        zm=np.median(out[out!=0])
        z1=np.median(out[(out!=0) & (out<zm)])
        z2=np.median(out[(out!=0) & (out>zm)])
        plt.imshow(out,origin='lower',vmin=2*z1-z2,vmax=2*z2-z1,cmap='Greys_r')
        plt.title(label+(' RA %9.5f,  DEC %9.5f' % (ra,dec)))
        plt.colorbar()
        infield=(xref>0) & (xref<width) & (yref>0) & (yref<width)
        plt.plot(xref[infield],yref[infield],'r+')
        plt.plot([xc*1.1,xc*1.2],[yc*1.1,yc*1.2],'w')
        plt.plot([xc*1.1,xc*1.2],[yc*0.9,yc*0.8],'w')
        plt.plot([xc*0.9,xc*0.8],[yc*1.1,yc*1.2],'w')
        plt.plot([xc*0.9,xc*0.8],[yc*0.9,yc*0.8],'w')
        plt.plot([xc*1.1,xc*1.2],[yc,yc],'k')
        plt.plot([xc,xc],[yc*0.9,yc*0.8],'k')
        plt.plot([xc,xc],[yc*1.1,yc*1.2],'k')
        plt.plot([xc*0.9,xc*0.8],[yc,yc],'k')
        plt.xlim(0,width)
        plt.ylim(0,width)
        if savefile == '':
            plt.show()
        else:
            plt.savefig(savefile)
            plt.clf()
            print('Plot saved in',savefile)

        return out

def iterweightedctr(x,y,sig,cap=4,verbose=False):
        '''
        x,y are two lists of N coordinates
        Find the centre of the distrubution by adaptively fitting Gaussian
        sig is fixed rms width of Gaussian weight, truncated at cap sigma
        Each iteration adjusts weights based on new ctr and then
           recomputes the ctr
        '''
        xc,yc=np.mean(x),np.mean(y)
        if verbose:
            print('Iterated centres:')
            print('-', xc,yc)
        for iter in range(10):
            rsig2=np.minimum(cap**2,((x-xc)**2+(y-yc)**2)/sig**2)
                 # square rad in units of filter sigma, capped at 4
            w=np.exp(-0.5*rsig2)-np.exp(-0.5*cap**2)
            wtot=w.sum()
            xc,yc=(w*x).sum()/wtot, (w*y).sum()/wtot
            if verbose:
                print(iter, xc,yc)
        return xc,yc
    
