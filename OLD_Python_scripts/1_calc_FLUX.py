import numpy as np 
from astropy.io import fits 
import matplotlib.pyplot as plt

name_fibA    = np.genfromtxt('fitsort_fibA.dat',dtype ='string', skip_header=1)
name_fibB    = np.genfromtxt('fitsort_fibB.dat',dtype ='string', skip_header=1)

#print np.size(name_fibA, axis=0)
for ii in xrange(0,np.size(name_fibA,axis=0)):
	nameA = name_fibA[ii,0]
	nameB = name_fibB[ii,0]
	hdulist = fits.open(nameA)
	scidata = hdulist[0].data
	fluxA = np.sum(scidata[69,:]/np.size(scidata[69,:]))
	hdulist.close()

	hdulist = fits.open(nameB)
	scidata = hdulist[0].data
	fluxB = np.sum(scidata[69,:]/np.size(scidata[69,:]))
	hdulist.close()
	print name_fibA[ii,0], fluxA, np.sqrt(np.abs(fluxA)), fluxB, np.sqrt(np.abs(fluxB))

