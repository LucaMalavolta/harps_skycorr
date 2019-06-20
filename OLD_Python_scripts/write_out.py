import numpy as np

input_fibA = np.genfromtxt('fitsort_fibA.dat',dtype ='string', skip_header=1)
name_fibA = input_fibA[:,0]
data_fibA = np.asarray(input_fibA[:,1:],dtype=np.double)

# data_fibA[:,0]  RV
# data_fibA[:,1]  FWHM
# data_fibA[:,2]  Contrast
# data_fibA[:,3]  CCF Noise [km/s]

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

rv_diff = data_pro[:,1]-data_sky[:,1]
rv_ratio = np.abs(rv_diff)/data_fibA[:,3]

f1_out = open('EP212521166_skycorrected.dat', 'w')
f1_out.write('file_root\tjdb\trv_DRS\tnoise\tfwhm_DRS\tbis_DRS\trv_skycor\tnoise\tfwhm_skycor\tbis_skycor\trv_DRS-rv_skycor\trvdiff_noise_ratio\n')
f1_out.write('---------\t---\t------\t-----\t--------\t-------\t---------\t-----\t-----------\t----------\t----------------\t------------------\n')
for ii in xrange(np.size(data_fibA[:,0])):
    f1_out.write('{0:s}\t{1:f}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:f}\t{9:f}\t{10:f}\t{11:f}\n'.format(
                input_fibA[ii,0],data_pro[ii,0],
                data_pro[ii,1],data_fibA[ii,3],data_pro[ii,2],data_pro[ii,3],
                data_sky[ii,1],data_fibA[ii,3],data_sky[ii,2],data_sky[ii,3],
                rv_diff[ii], rv_ratio[ii]))
f1_out.close()
