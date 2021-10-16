'''
Example use of FindMinorPlanet.
Read in a table of Minor planet ID's, and use the MPC or JPL empeherides to find a match with the ESO VST and VISTA observations.
'''


import FindMinorPlanet as mp
import astropy.io.fits as pf
import os

DESTNO=pf.open('DESY4TNO_apjsab6bd8_table2.fits')[1].data['MPC']
#DESTNO=pf.open('DESY4TNO_apjsab6bd8_table3.fits')[1].data['MPC']   # table of 'others'

# Initialise the VST and VISTA observations databases
s=mp.MinorPlanetSearch(newdb=False)

# go through all objects in the DESY4 TNO table and look for matches
# the outputs are csv files listing observations that may show the objects
# these can be further processed with the VSTcutout package

for id in DESTNO:
    if not os.path.exists(id+'.eph.fits'):
        try:
            vst,vista=s.search(id)
        except:
            try:
                vst,vista=s.search(id,'JPL')
            except:
                print('No ephemerides for',id)
