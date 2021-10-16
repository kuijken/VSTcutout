# VSTcutout
Scripts to make cutouts from the ESO VST raw data archive

These scripts were developed to look for chance observations of named asteroids in VST and VISTA data.
There is a tool to get the ESO observations archive as a fits file, and then match that to MPC or JPL 
ephemerides for any object. See FindMinorPlanet.py for more details.

A separate tool provides pixel cutouts from raw VST images, given its coordinates. It includes optional 
tweaking of the astrometry in the raw fits headers to give good alignments. See VSTcutout.py for mode details.

An example script is provided.
<img width="668" alt="Screenshot 2021-10-16 at 16 05 08" src="https://user-images.githubusercontent.com/6078683/137590490-b0e7c78c-a81b-441f-a8b5-c78cb167623d.png">
