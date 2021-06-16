import matplotlib.pyplot as plt
#import seaborn
import numpy as np
import argparse
import os
import sys

def single_plot(target_name, plot_dirs):
    data = np.genfromtxt(target_name + '_sky_contamination.dat', skip_header=1, delimiter='\t',
        names=['TARG_NAME','RA-DEG','DEC-DEG','file_rad','file_mask','BJD','ccf_DRS_RV',
        'ccf_DRS_FWHM','ccf_DRS_contrast','ccf_DRS_continuum','ccf_COR_RV','ccf_COR_FWHM',
        'ccf_COR_contrast','ccf_COR_continuum','CCF_NOISE','SCI_SN20','SCI_SN40','SCI_SN60',
        'SKY_SN20','SKY_SN40','SKY_SN60','fluxA','fluxB','BERV','MOON_dist','MOON_El_app',
        'MOON_airmass','MOON_APmag','MOON_surf_bright','MOON_illum','MOON_RV_Obs','MOON_flag'])

    data_str = np.genfromtxt(target_name+'_sky_contamination.dat', skip_header=1, delimiter='\t',dtype=str,
        usecols=(0,3))

    NOISE_ratio = (data['ccf_COR_RV']-data['ccf_DRS_RV'])/data['CCF_NOISE']
    Delta_RV = data['BERV']-data['ccf_DRS_RV']

    sel_nofibB =  (data['MOON_flag']<0.5)
    sel_wtfibB =  (data['MOON_flag']>0.5)

    flag_fiberB = np.zeros(len(sel_wtfibB), dtype=np.int16)
    flag_fiberB[sel_wtfibB] = 1

    NOISE_ratio[sel_nofibB] = 0.000

    sel_contaminated = (np.abs(NOISE_ratio)> 2.) & (sel_wtfibB)
    sel_uncontaminated = (np.abs(NOISE_ratio) < 2.) & (sel_wtfibB)

    sel_likely_contaminated = (np.abs(Delta_RV)<15.) & (data['MOON_illum']>95) & (data['MOON_dist']<60) & (sel_nofibB)

    sel_2outof3_contaminated = (~sel_contaminated) & (~sel_likely_contaminated) & (((np.abs(Delta_RV)<15.) & (data['MOON_illum']>95)) | \
     ((data['MOON_illum']>95) & (data['MOON_dist']<60)) | (((np.abs(Delta_RV)<15.) & (data['MOON_dist']<60))))

    flag_contamination = np.zeros(len(sel_contaminated), dtype=np.int16)
    flag_contamination[sel_2outof3_contaminated] = 1
    flag_contamination[sel_likely_contaminated] = 2
    flag_contamination[sel_contaminated] = 3

    fileout = open(target_name + '_skycorr_summary_selected.dat', 'w')
    fileout.write('# CONTAMINATED ONLY (flag 3): \n')
    fileout.write('# Target\tFILE\tBJD\tDelta_V_vs_Noise\tBERV-RV\tMOON_dist\tMOON_illum\tMOON_El\tfiberB\tflag\n')
    fileout.write('# ------\t----\t---\t----------------\t-------\t---------\t----------\t-------\t------\t----\n')
    for i_sel, f_sel in enumerate(sel_contaminated):
        if f_sel:
            fileout.write('{0:s}\t{1:s}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:d}\t{9:d}\n'.format(
            data_str[i_sel,0], data_str[i_sel,1], data['BJD'][i_sel], \
            NOISE_ratio[i_sel],Delta_RV[i_sel], \
            data['MOON_dist'][i_sel], data['MOON_illum'][i_sel], data['MOON_El_app'][i_sel], flag_fiberB[i_sel], 3))
    fileout.write('#\n')
    fileout.write('#\n')

    fileout.write('# NO FIBER B BUT LIKELY CONTAMINATED (flag 2):\n')
    fileout.write('# Target\tFILE\tBJD\tDelta_V_vs_Noise\tBERV-RV\tMOON_dist\tMOON_illum\tMOON_El\tfiberB\tflag\n')
    fileout.write('# ------\t----\t---\t----------------\t-------\t---------\t----------\t-------\t------\t----\n')
    for i_sel, f_sel in enumerate(sel_likely_contaminated):
        if f_sel:
            fileout.write('{0:s}\t{1:s}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:d}\t{9:d}\n'.format(
            data_str[i_sel,0], data_str[i_sel,1], data['BJD'][i_sel], \
            0.0000,Delta_RV[i_sel], \
            data['MOON_dist'][i_sel], data['MOON_illum'][i_sel], data['MOON_El_app'][i_sel], flag_fiberB[i_sel], 2))
    fileout.write('#\n')
    fileout.write('#\n')

    fileout.write('# TWO CONTAMINATION CRITERIA OUT OF THREE (flag 1):\n')
    fileout.write('# Target\tFILE\tBJD\tDelta_V_vs_Noise\tBERV-RV\tMOON_dist\tMOON_illum\tMOON_El\tfiberB\tflag\n')
    fileout.write('# ------\t----\t---\t----------------\t-------\t---------\t----------\t-------\t------\t----\n')
    for i_sel, f_sel in enumerate(sel_2outof3_contaminated):
        if f_sel:
            fileout.write('{0:s}\t{1:s}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:d}\t{9:d}\n'.format(
            data_str[i_sel,0], data_str[i_sel,1], data['BJD'][i_sel], \
            NOISE_ratio[i_sel],Delta_RV[i_sel], \
            data['MOON_dist'][i_sel], data['MOON_illum'][i_sel], data['MOON_El_app'][i_sel], flag_fiberB[i_sel], 1))

    fileout.close()

    fileout = open(target_name + '_skycorr_summary_all.dat', 'w')
    fileout.write('# Target\tFILE\tBJD\tDelta_V_vs_Noise\tBERV-RV\tMOON_dist\tMOON_illum\tMOON_El\tfiberB\tflag\n')
    fileout.write('# ------\t----\t---\t----------------\t-------\t---------\t----------\t-------\t------\t----\n')
    for i_sel, f_sel in enumerate(flag_contamination):
            #if f_sel:
            fileout.write('{0:s}\t{1:s}\t{2:f}\t{3:f}\t{4:f}\t{5:f}\t{6:f}\t{7:f}\t{8:d}\t{9:d}\n'.format(
            data_str[i_sel,0], data_str[i_sel,1], data['BJD'][i_sel], \
            NOISE_ratio[i_sel],Delta_RV[i_sel], \
            data['MOON_dist'][i_sel], data['MOON_illum'][i_sel], data['MOON_El_app'][i_sel], flag_fiberB[i_sel], f_sel))

    fileout.close()

    plt.xlabel('MOON_illum')
    plt.ylabel('BERV-RV')
    plt.axvline(95.)
    plt.xlim(0,100)
    plt.scatter(data['MOON_illum'][sel_uncontaminated],Delta_RV[sel_uncontaminated], c='C0', s=2)
    plt.scatter(data['MOON_illum'][sel_contaminated],Delta_RV[sel_contaminated], c='C1', s=2)
    plt.scatter(data['MOON_illum'][sel_nofibB],Delta_RV[sel_nofibB], c='C2', s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOONillum_DeltaRV.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('SCI_SN60')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['SCI_SN60'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_SCI_SN60.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('SKY_SN60')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['SKY_SN60'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_SKY_SN60.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('flux_B')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['fluxB'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_fluxB.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('MOON_airmass')
    plt.ylim(3.,1.)
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['MOON_airmass'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOON_airmass.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('MOON_dist')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['MOON_dist'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOON_dist.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('MOON_El_app')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['MOON_El_app'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOON_El_app.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('MOON_APmag')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['MOON_APmag'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOON_APmag.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('MOON_surf_bright')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['MOON_surf_bright'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOON_surf_bright.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('MOON_illum')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['MOON_illum'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOON_illum.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('FWHM [km/s]')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['ccf_DRS_FWHM'],s=5,c='r')
    plt.scatter(NOISE_ratio,data['ccf_COR_FWHM'],s=5,c='b')
    plt.savefig(plot_dirs+'/'+target_name+'_ccf_DRS_FWHM.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('BERV-RV')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,Delta_RV, s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_BERV_RV.pdf')
    plt.close()

    plt.xlabel('CCF Noise ratio')
    plt.ylabel('MOON_RV_Obs')
    plt.axvline(-2.)
    plt.axvline(2.)
    plt.scatter(NOISE_ratio,data['MOON_RV_Obs'], s=2)
    plt.savefig(plot_dirs+'/'+target_name+'_MOON_RV_Obs.pdf')
    plt.close()

parser = argparse.ArgumentParser(prog='GET_vals', description='Analyze SkyCorr output')
parser.add_argument('file_list', type=str, nargs=1, help='input file list')


args = parser.parse_args()
file_list_name = args.file_list[0]
target_name = file_list_name[:-5]
plot_dirs = target_name + '_skycorr_plots'

os.system("mkdir -p " + plot_dirs)
single_plot(target_name, plot_dirs)
