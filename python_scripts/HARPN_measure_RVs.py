import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit
from astropy.io import fits

def gauss_p0(x, y):
    return np.asarray([-1*((max(y)-min(y))/max(y)),\
                x[y.argmin()],8*(x[1]-x[0]),1.])

def gauss_f(x, depth, x0, sigma, continuum):
    return continuum + depth * np.exp(- (x-x0)**2 / (2*sigma**2))

#   noise=ccdsigdet*sqrt(ic_extmeanzone)
ic_extmeanzone=6
ic_plate_scale=0.82
ccf_step = 0.25
ccf_oversampling=ic_plate_scale/ccf_step


parser = argparse.ArgumentParser(prog='HARPN_measure_RVs.py', description='Analyze chromatic RVs')
parser.add_argument('file_list', type=str, nargs=1, help='input file list')
parser.add_argument('output_dir', type=str, nargs=1, help='directory with skycorr output ')
parser.add_argument('date_output', type=int, nargs=1, help='use the date to search for output files (0=False)')

args = parser.parse_args()
file_list_name = args.file_list[0]
output_dir = args.output_dir[0]
date_output = args.date_output[0]

file_list = np.genfromtxt(file_list_name, dtype=str)

# open the first file to create the proper arrays
if date_output==0:
    CCF_hdu = fits.open(output_dir + '/' + file_list[0,1] + '_ccf_' + file_list[0,2] + '_A_chromatic.fits')
else:
    CCF_hdu = fits.open(output_dir + '/' + file_list[0,0] + '/' + file_list[0,1] + '_ccf_' + file_list[0,2] + '_A_chromatic.fits')

n_chromatic = CCF_hdu[0].header['NAXIS2'] - 1
n_files = len(file_list)

results_drs = {
    'bjd': np.zeros(n_files, dtype=np.double),
    'rv': np.zeros([n_files,n_chromatic+1], dtype=np.double),
    'noise': np.zeros([n_files,n_chromatic+1], dtype=np.double),
    'contrast': np.zeros([n_files,n_chromatic+1], dtype=np.double),
    'fwhm': np.zeros([n_files,n_chromatic+1], dtype=np.double)
}

results_sky = {
    'bjd': np.zeros(n_files, dtype=np.double),
    'rv': np.zeros([n_files,n_chromatic+1], dtype=np.double),
    'noise': np.zeros([n_files,n_chromatic+1], dtype=np.double),
    'contrast': np.zeros([n_files,n_chromatic+1], dtype=np.double),
    'fwhm': np.zeros([n_files,n_chromatic+1], dtype=np.double)
}

for i_file in range(0,n_files):
    file_date, file_rad_obs, file_mask = file_list[i_file,:3]
    print(' Analyzing:', file_date, file_rad_obs, file_mask)

    if date_output==0:
        file_rad = output_dir + '/' + file_rad_obs + '_ccf_' + file_mask + '_A'
    else:
        file_rad = output_dir + '/' + file_date + '/' + file_rad_obs + '_ccf_' + file_mask + '_A'

    if file_rad_obs[:5] == 'HARPS':
        instr_key = 'ESO'
    elif file_rad_obs[:5] == 'HARPN':
        instr_key = 'TNG'
    else:
        print('Instrument not supported')
        quit()

    CCF_hdu = fits.open(file_rad + ".fits")
    RV_CCF = CCF_hdu[0].header['CRVAL1'] + np.arange(0, CCF_hdu[0].header['NAXIS1'], 1) * CCF_hdu[0].header['CDELT1']

    results_drs['bjd'][i_file] = CCF_hdu[0].header['HIERARCH '+instr_key+' DRS BJD']
    results_sky['bjd'][i_file] = CCF_hdu[0].header['HIERARCH '+instr_key+' DRS BJD']
    noise = CCF_hdu[0].header['HIERARCH '+instr_key+' DRS CCD SIGDET']*np.sqrt(ic_extmeanzone)

    CCF_hdu.close()

    CCF_table = np.genfromtxt(file_rad+".tbl", skip_header=2)

    for ccf_type in ["", "_skyR"]:
    # Chromatic RVs

        if not os.path.exists(file_rad + ccf_type + ".fits"):
            continue

        CCF_hdu = fits.open(file_rad + ccf_type + ".fits")
        CCF_data = CCF_hdu[0].data
        CCF_hdu.close()

        CCF_hdu = fits.open(file_rad + "_chromatic"+ccf_type+".fits")


        for i_ccf in range(0,n_chromatic+1):
            average_CCF = CCF_hdu[0].data[i_ccf,:]
            normalized_CCF=average_CCF/max(average_CCF)
            popt, pcov = curve_fit(gauss_f, RV_CCF, normalized_CCF, p0=gauss_p0(RV_CCF, normalized_CCF))
            # plt.scatter(RV_CCF, normalized_CCF, s=2, c='C0')
            # plt.plot(RV_CCF, gauss_f(RV_CCF, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f , c=%5.3f' % tuple(popt))

            if i_ccf == 0:
                order_list = np.asarray(CCF_table[:,0], dtype=np.int16)
            else:
                order_list = np.genfromtxt(file_rad+"_chromatic_"+repr(i_ccf)+".dat", skip_header=1, usecols=0, dtype=np.int16)

            ## CCF correction
            CCF_noise_temp = average_CCF*0.0000
            for ii in order_list:

                if CCF_table[ii,2] < 1: continue
                CCF_noise_temp += CCF_data[ii,:] +  noise**2 * np.max(CCF_data[ii,:])/CCF_table[ii,1] *1.225

            CCF_noise_tot=np.sqrt(CCF_noise_temp) # ic_plate_scale
            CCF_slope=(average_CCF[2:]-average_CCF[:-2])/(RV_CCF[2:]-RV_CCF[:-2])

            indexlist=map(int,np.arange(np.round(len(CCF_slope)/ccf_oversampling))*ccf_oversampling)

            qq=np.zeros(len(CCF_slope))
            for i in indexlist:
                qq[i]=1
            factor = 1. / 1.404
            ccf_noise=factor * (np.sum(np.compress(qq,CCF_slope**2)/np.compress(qq,CCF_noise_tot[1:-1]**2)))**(-0.5)


            if ccf_type == "" :
                results_drs['rv'][i_file, i_ccf] =  popt[1]
                results_drs['contrast'][i_file, i_ccf] = (-popt[0])/popt[3]*100.
                results_drs['fwhm'][i_file, i_ccf] = popt[2] * np.sqrt(8 * np.log(2))
                results_drs['noise'][i_file, i_ccf] =  ccf_noise
            else:
                results_sky['rv'][i_file, i_ccf] =  popt[1]
                results_sky['contrast'][i_file, i_ccf] = (-popt[0])/popt[3]*100.
                results_sky['fwhm'][i_file, i_ccf] = popt[2] * np.sqrt(8 * np.log(2))
                results_sky['noise'][i_file, i_ccf] =  ccf_noise

        CCF_hdu.close()

file_rad = file_list_name[:-5] + "_chromatic"

fileout_drs = open(file_rad + "_drs.dat", "w")
fileout_sky = open(file_rad + "_sky.dat", "w")

fileout_drs.write('# date rad bjd rv_drs noise_drs contrast_drs fwhm_drs ')
for i_ccf in range(1,n_chromatic+1):
    fileout_drs.write('rv_{0:1.0f} noise_{0:1.0f} contrast_{0:1.0f} fwhm_{0:1.0f} '.format(i_ccf))
fileout_drs.write('\n')

fileout_sky.write('# date rad bjd rv_sky noise_sky contrast_sky fwhm_sky ')
for i_ccf in range(1,n_chromatic+1):
    fileout_sky.write('rv_{0:1.0f} noise_{0:1.0f} contrast_{0:1.0f} fwhm_{0:1.0f} '.format(i_ccf))
fileout_sky.write('\n')

for i_file in range(0,n_files):
    file_date, file_rad_obs, file_mask = file_list[i_file,:3]

    fileout_drs.write('{0:s} {1:s} {2:15.8f} '.format(file_date, file_rad_obs, results_drs['bjd'][i_file]))
    for i_ccf in range(0,n_chromatic+1):
        fileout_drs.write('{0:12.5f} {1:12.5f} {2:12.5f} {3:12.5f} '.format(results_drs['rv'][i_file, i_ccf],
            results_drs['noise'][i_file, i_ccf],
            results_drs['contrast'][i_file, i_ccf],
            results_drs['fwhm'][i_file, i_ccf]))

    fileout_drs.write('\n')


    fileout_sky.write('{0:s} {1:s} {2:15.8f} '.format(file_date, file_rad_obs, results_sky['bjd'][i_file]))
    for i_ccf in range(0,n_chromatic+1):
        fileout_sky.write('{0:12.5f} {1:12.5f} {2:12.5f} {3:12.5f} '.format(results_sky['rv'][i_file, i_ccf],
            results_sky['noise'][i_file, i_ccf],
            results_sky['contrast'][i_file, i_ccf],
            results_sky['fwhm'][i_file, i_ccf]))
    fileout_sky.write('\n')
