import VSTcutout
import FindMinorPlanet

s=FindMinorPlanet.MinorPlanetSearch()   # initialise the observations db

vst,vista=s.search('22191')  # get empeherides for '22191 achucarro'
                             # return matches in VST and VISTA observations
                             #    as tables (and write out csv files)

# make a cutout with the first VST match:
ra=vst['RAtarg'][0]
dec=vst['DECtarg'][0]
url=vst['access_url'][0]

# download the image from the archive, calculate the astrometry corrections
v=VSTcutout.VSTcutout(url)

# display a 300x300 pixel cutout around the RA,DEC of the target
cutout=v.cutout(ra,dec,300,offset='fit')

# to make more cutouts from the same image, keep calling v.cutout
# to load a new image from the archive, rerun VSTcutout.VSTcutout(url)
