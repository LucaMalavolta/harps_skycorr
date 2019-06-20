import numpy as np
import matplotlib.pyplot as plt

input_sky = np.genfromtxt('line_sky.txt',dtype ='string', skip_header=2)
name_sky = input_sky[:,0]
data_sky = np.asarray(input_sky[:,1:],dtype=np.double)

input_pro = np.genfromtxt('line_profile.txt',dtype ='string', skip_header=2)
name_pro = input_pro[:,0]
data_pro = np.asarray(input_pro[:,1:],dtype=np.double)

# data_pro[:,0]  JD
# data_pro[:,1]  RV [km/s]
# data_pro[:,2]  Bis
# data_pro[:,3]  FWHM

input_fibA = np.genfromtxt('fitsort_fibA.dat',dtype ='string', skip_header=1)
name_fibA = input_fibA[:,0]
data_fibA = np.asarray(input_fibA[:,1:],dtype=np.double)

# data_fibA[:,0]  RV
# data_fibA[:,1]  FWHM
# data_fibA[:,2]  Contrast
# data_fibA[:,3]  CCF Noise [km/s]
# data_fibA[:,4]  V mag
# data_fibA[:,5]  Seeing
# data_fibA[:,6]  BERV
# data_fibA[:,7]  SN
# data_fibA[:,8]  SN
# data_fibA[:,9]  SN

input_fibB = np.genfromtxt('fitsort_fibB.dat',dtype ='string', skip_header=1)
name_fibB = input_fibB[:,0]
data_fibB = np.asarray(input_fibB[:,1:],dtype=np.double)

input_flux = np.genfromtxt('flux_fits.dat',dtype ='string')
name_flux = input_flux[:,0]
data_flux = np.asarray(input_flux[:,1:],dtype=np.double)


data_x1out = np.ones(np.size(name_fibA))
data_x2out = np.ones(np.size(name_fibA))
data_x3out = np.ones(np.size(name_fibA))

data_y1out = np.ones(np.size(name_fibA))
data_y2out = np.ones(np.size(name_fibA))
data_y3out = np.ones(np.size(name_fibA))

for ii in xrange(np.size(name_fibA)):
   cond_fib=np.in1d(name_fibB,name_fibA[ii])
   cond_pro=np.in1d(name_pro,name_fibA[ii])
   cond_sky=np.in1d(name_sky,name_fibA[ii])
   cond_flx=np.in1d(name_flux,name_fibA[ii])

   data_y1out[ii]=  data_flux[cond_flx,2]/data_flux[cond_flx,0]
   data_y2out[ii]=((data_sky[cond_sky,1]-data_pro[cond_pro,1])*1000.)
   # RV[sky]-RV[DRS], m/s
   data_y3out[ii] = (data_fibA[ii,3]*1000.)

cond_flag2 = (np.abs(data_y2out)>2.*data_y3out)

plt.scatter(data_fibA[:,3]*1000, data_y1out)
plt.scatter(data_fibA[cond_flag2,3]*1000, data_y1out[cond_flag2], c='r')
plt.axvline(5)
plt.xlabel('CCF_noise fibA [m/s]')
plt.ylabel('<flux B> / <flux A> ')
plt.show()
