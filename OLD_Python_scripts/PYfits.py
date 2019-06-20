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

print np.size(name_fibA)
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
   #data_y1out[ii]=((data_sky[cond_pro,1]-data_pro[cond_sky,1])*1000.)
   data_y1out[ii]=((data_sky[cond_sky,1]-data_pro[cond_pro,1])*1000.)
   # RV[sky]-RV[DRS], m/s



   data_y2out[ii]= np.abs(data_y1out[ii]/(data_fibA[ii,3]*1000.))
   # (RV[sky]-RV[DRS])/CCF_noise

   if data_y2out[ii] > 2:
       print 'COND 2', data_pro[cond_pro],data_fibA[ii,3], data_y1out[ii],data_fibA[ii,6]-data_fibA[ii,0]
   data_y3out[ii]= np.abs(data_fibB[cond_fib,9])/data_fibA[ii,9]


   data_x1out[ii]=(data_fibA[ii,3]*1000.)
   data_x2out[ii]=(data_fibA[ii,6]-data_fibA[ii,0])
   data_x3out[ii]=(data_flux[cond_flx,2]/data_flux[cond_flx,1])

cond_flag1 = (np.abs(data_y1out)>1.*data_x1out)
cond_flag2 = (np.abs(data_y1out)>2.*data_x1out)
cond_flag3 = (np.abs(data_y1out)>3.*data_x1out)
cond_flag4 = (data_x3out > 15.) & cond_flag2
cond_flag5 = (data_y3out > 0.07)
cond_10ms = (data_x1out<100.)
cond_10ms_sel = (data_x1out<100.) & (cond_flag2)


for ii in xrange(-45,35,10):
	cond_rvber = (data_x2out>ii) & (data_x2out<=ii+10) & cond_10ms
	conf_flag_sel = cond_flag2 & cond_rvber
	print ii+5, np.sum(conf_flag_sel), np.sum(cond_rvber)


xx=np.arange(0,100,1)


print np.sum(cond_10ms), np.sum(cond_10ms_sel)
plt.xlim(0.,15.)
plt.ylim(-25.,25.)
plt.xlabel('CCF noise [m/s]')
plt.ylabel('RV[sky_cor]-RV[DRS] [m/s]')
plt.plot(data_x1out,data_y1out,'ro')
plt.plot(data_x1out[cond_flag2],data_y1out[cond_flag2],'bo')
plt.plot(data_x1out[cond_flag3],data_y1out[cond_flag3],'go')
plt.plot(data_x1out[cond_flag4],data_y1out[cond_flag4],'co')
plt.plot(xx,xx,'c')
plt.plot(xx,-xx,'c')
plt.plot(xx,2*xx,'m')
plt.plot(xx,-2*xx,'m')
plt.show()


plt.ylim(-25.,25.)
plt.plot(data_x2out,data_y2out,'ro')
plt.plot(data_x2out[cond_flag2],data_y2out[cond_flag2],'bo')
plt.plot(data_x2out[cond_flag3],data_y2out[cond_flag3],'go')
plt.plot(data_x2out[cond_flag4],data_y2out[cond_flag4],'co')
plt.plot(data_x2out[cond_flag5],data_y2out[cond_flag5],'yo')
plt.show()

def left_out():
    plt.ylim(-25.,25.)
    plt.plot(data_x2out,data_y2out,'ro')
    plt.plot(data_x2out[cond_flag2],data_y2out[cond_flag2],'bo')
    plt.plot(data_x2out[cond_flag3],data_y2out[cond_flag3],'go')
    plt.plot(data_x2out[cond_flag4],data_y2out[cond_flag4],'co')
    plt.plot(data_x2out[cond_flag5],data_y2out[cond_flag5],'yo')
    plt.show()

    #plt.ylim(-0.1,1.)
    plt.xlim(0.,15.)
    plt.plot(data_y2out,data_y3out,'ro')
    plt.plot(data_y2out[cond_flag2],data_y3out[cond_flag2],'bo')
    plt.plot(data_y2out[cond_flag3],data_y3out[cond_flag3],'go')
    plt.plot(data_y2out[cond_flag4],data_y3out[cond_flag4],'co')
    plt.plot(data_y2out[cond_flag5],data_y3out[cond_flag5],'yo')

    plt.show()

    plt.xlim(-2,10)
    plt.ylim(-25.,25.)
    plt.plot(data_x3out,data_y2out,'ro')
    plt.plot(data_x3out[cond_flag2],data_y2out[cond_flag2],'bo')
    plt.plot(data_x3out[cond_flag3],data_y2out[cond_flag3],'go')
    plt.plot(data_x3out[cond_flag5],data_y2out[cond_flag5],'yo')
    plt.show()
