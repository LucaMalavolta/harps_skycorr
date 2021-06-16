import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from astropy.coordinates import SkyCoord
import os

def residuals_gaussf(p,x,y,err=1.000):

	"""
	gaussian(x, A, mean, FWHM, cont)

	p[0] = RV;
	p[1] = FWHM;
	p[2] = contrast;
	p[3] = continuum

	"""

	sigma = np.abs(p[1]) /( 2 * np.sqrt(2 * np.log(2)) )
	tau = -((x - p[0])**2) / (2*(sigma**2))
	model = (1 - p[2]/100. * np.exp(tau))*p[3]
	return (y-model)/err

# Retrival of Moon tabulated data


parser = argparse.ArgumentParser(prog='GET_vals', description='Analyze SkyCorr output')
parser.add_argument('file_list', type=str, nargs=1, help='input file list')
parser.add_argument('output_dir', type=str, nargs=1, help='directory with skycorr output ')
parser.add_argument('date_output', type=int, nargs=1, help='use the date to search for output files (0=False)')

args = parser.parse_args()
file_list_name = args.file_list[0]
output_dir = args.output_dir[0]
date_output = args.date_output[0]

target = file_list_name[:-5]

print(' ** Reading MOON data tables - it may be slow!')
try:
	MOON_TNG = np.genfromtxt('/home/malavolta/CODE/harps_skycorr/tables/JPL_moon_TNG_2012_2020.csv',
		names = "BJD, flag1, flag2, RA, DEC, RA_app, DEC_app, Az_app, El_app, airmass, mag_ex, APmag, surf_bright, illum, RV_Sun, RV_Obs, none",
					 delimiter=',', skip_header=167)
except:
	MOON_TNG = np.genfromtxt('/Users/malavolta/Astro/CODE/harps_skycorr/tables/JPL_moon_TNG_2012_2020.csv',
		names = "BJD, flag1, flag2, RA, DEC, RA_app, DEC_app, Az_app, El_app, airmass, mag_ex, APmag, surf_bright, illum, RV_Sun, RV_Obs, none",
					 delimiter=',', skip_header=167)

try:
	MOON_ESO = np.genfromtxt('/home/malavolta/CODE/harps_skycorr/tables/JPL_moon_LaSilla_2012_2017.csv',
		names = "BJD, flag1, flag2, RA, DEC, RA_app, DEC_app, Az_app, El_app, airmass, mag_ex, APmag, surf_bright, illum, RV_Sun, RV_Obs, none",
					 delimiter=',', skip_header=167)
except:
	MOON_ESO = np.genfromtxt('/Users/malavolta/Astro/CODE/harps_skycorr/tables/JPL_moon_LaSilla_2012_2017.csv',
		names = "BJD, flag1, flag2, RA, DEC, RA_app, DEC_app, Az_app, El_app, airmass, mag_ex, APmag, surf_bright, illum, RV_Sun, RV_Obs, none",
					 delimiter=',', skip_header=167)

MOON={
	'TNG': MOON_TNG,
	'ESO': MOON_ESO
}

print()


name_fibA = np.genfromtxt(file_list_name,dtype =str, skip_header=0)
fileout = open(target + '_sky_contamination.dat','w')

fileout.write('#TARG_NAME\tRA-DEG\tDEC-DEG\tfile_rad\tfile_mask\tBJD\t'\
		'ccf_DRS_RV\tccf_DRS_FWHM\tccf_DRS_contrast\tccf_DRS_continuum\t'\
		'ccf_COR_RV\tccf_COR_FWHM\tccf_COR_contrast\tccf_COR_continuum\t'\
		'CCF_NOISE\tSCI_SN20\tSCI_SN40\tSCI_SN60\tSKY_SN20\tSKY_SN40\tSKY_SN60\t'\
		'fluxA\tfluxB\tBERV\tMOON_dits\tMOON_El_app\tMOON_airmass\tMOON_APmag\tMOON_surf_bright\tMOON_illum\tMOON_RV_Obs\n')

for file_date, file_rad,file_mask,file_fib in zip(name_fibA[:,0], name_fibA[:,1],name_fibA[:,2],name_fibA[:,3]):

	if date_output==0:
		name_DRS = output_dir + '/' + file_rad + '_ccf_' + file_mask + '_' + file_fib + '.fits'
		name_COR = output_dir + '/' + file_rad + '_ccf_' + file_mask + '_' + file_fib + '_skyR.fits'
		name_SKY = output_dir + '/' + file_rad + '_ccf_' + file_mask + '_B.fits'
	else:
		name_DRS = output_dir + '/' + file_date + '/' + file_rad + '_ccf_' + file_mask + '_' + file_fib + '.fits'
		name_COR = output_dir + '/' + file_date + '/' + file_rad + '_ccf_' + file_mask + '_' + file_fib + '_skyR.fits'
		name_SKY = output_dir + '/' + file_date + '/' + file_rad + '_ccf_' + file_mask + '_B.fits'

	if file_rad[:5] == 'HARPS':
		instr_key = 'ESO'
	elif file_rad[:5] == 'HARPN':
		instr_key = 'TNG'

	hdulist = fits.open(name_DRS)

	#header saved here for future reference
	header_DRS =  hdulist[0].header

	# Get coordinates - we assume that proper motions are irrilevant for our analysis
	try:
		RA_deg = header_DRS['RA-DEG']
		DEC_deg = header_DRS['DEC-DEG']
	except:
		RA_deg = header_DRS['RA-RAD'] * 180.00 / np.pi
		DEC_deg = header_DRS['DEC-RAD'] * 180.00 / np.pi # weird choice of using DEC in hours


	BJD = header_DRS['HIERARCH '+instr_key+' DRS BJD']
	BERV = header_DRS['HIERARCH '+instr_key+' DRS BERV']
	CCF_noise = header_DRS['HIERARCH '+instr_key+' DRS CCF NOISE']

	""" Other keywords that mey be useful
	header_DRS['HIERARCH TNG DRS SPE EXT SN20']
	header_DRS['HIERARCH TNG DRS SPE EXT SN40']
	header_DRS['HIERARCH TNG DRS SPE EXT SN60']
	header_DRS['HIERARCH TNG OBS TARG NAME']
	"""

	""" Total amount of flux for fib A"""
	fluxA = np.sum(hdulist[0].data[69,:]/np.size(hdulist[0].data[69,:]))


	""" extracting Moon info
		1) identification of the closest (in time) tabulated data
	"""
	idx = np.argmin(np.abs(MOON[instr_key]['BJD']-BJD))

	""" 2) Picking the two closest tabulated data (one before and one after the observed BJD)
	"""
	if BJD>MOON[instr_key]['BJD'][idx]:
		idx1 = idx
		idx2 = idx+1
	else:
		idx1 = idx-1
		idx2 = idx

	""" 3) Moon distance computed by linear interpolation of the two closest points
	"""
	target = SkyCoord(ra=RA_deg, dec=DEC_deg, frame='icrs', unit='deg')
	dist_idx1 = target.separation(SkyCoord(ra=MOON[instr_key]['RA_app'][idx1], dec=MOON[instr_key]['DEC_app'][idx1],frame='icrs', unit='deg')).deg
	dist_idx2 = target.separation(SkyCoord(ra=MOON[instr_key]['RA_app'][idx2], dec=MOON[instr_key]['DEC_app'][idx2],frame='icrs', unit='deg')).deg

	target_dist = (dist_idx2-dist_idx1)/(MOON[instr_key]['BJD'][idx2]-MOON[instr_key]['BJD'][idx1])*(BJD-MOON[instr_key]['BJD'][idx1]) + dist_idx1

	""" 4) Other MOON values that could be useful
		MOON[instr_key]['El_app']
		MOON[instr_key]['airmass']
		MOON[instr_key]['APmag']
		MOON[instr_key]['surf_bright']
		MOON[instr_key]['illum']
		MOON[instr_key]['RV_Obs']
	"""

	#extracting the CCF
	ccf_x = np.arange(0., header_DRS['CDELT1']*header_DRS['NAXIS1'],
						header_DRS['CDELT1']) + header_DRS['CRVAL1']
	ccf_y = hdulist[0].data[-1,:]
	ccf_e = np.ones(np.size(ccf_x))

	#Getting the starting points for RV measurement
	p0 = np.zeros(4)
	p0[0] = header_DRS['HIERARCH '+instr_key+' DRS CCF RV']
	p0[1] = header_DRS['HIERARCH '+instr_key+' DRS CCF FWHM']
	p0[2] = header_DRS['HIERARCH '+instr_key+' DRS CCF CONTRAST']
	p0[3] = np.amax(ccf_y)

	plsq = leastsq(residuals_gaussf, p0, args=(ccf_x, ccf_y, ccf_e))

	ccf_pams_DRS = plsq[0][:]

	hdulist.close()
	#Repeating for the sky-corrected CCF

	try:

		hdulist = fits.open(name_COR)
		#extracting the CCF
		ccf_y = hdulist[0].data[-1,:]
		ccf_e = np.ones(np.size(ccf_x))

		plsq = leastsq(residuals_gaussf, p0, args=(ccf_x, ccf_y, ccf_e))
		ccf_pams_COR = plsq[0][:]

		hdulist.close()

		""" Retrival of a few interesting keywords and data from fib B """
		hdulist = fits.open(name_SKY)
		header_SKY =  hdulist[0].header
		fluxB = np.sum(hdulist[0].data[69,:]/np.size(hdulist[0].data[69,:]))
		hdulist.close()

		header_sky_SN20 = header_SKY['HIERARCH '+instr_key+' DRS SPE EXT SN20']
		header_sky_SN40 = header_SKY['HIERARCH '+instr_key+' DRS SPE EXT SN40']
		header_sky_SN60 = header_SKY['HIERARCH '+instr_key+' DRS SPE EXT SN60']
		MOON_flag = 1

	except:
		ccf_pams_COR = np.zeros(4)
		fluxB = 0.0
		header_sky_SN20 = 0.0
		header_sky_SN40 = 0.0
		header_sky_SN60 = 0.0
		MOON_flag = 0


	#fileout.write('{ 0:s}\t{ 1:f}\t{ 2:f}\t{ 3:s}\t{ 4:s}\t{ 5:f}\t{ 6:f}\t{ 7:f}\t{ 8:f}\t{ 9:f}\t{10:f}\t{11:f}\t{12:f}\t{13:f}\t{14:f}\t{15:f}\t{16:f}\t{17:f}\t{18:f}\t{19:f}\t{20:f}\t{21:f}\t{22:f}\t{23:f}\t{24:f}\t{25:f}\t{26:f}\t{27:f}\t{28:f}\t{29:f}\t{30:f}\n'.format(
	fileout.write('{0:s}\t{1:f}\t{2:f}\t{3:s}\t{4:s}\t{5:f}\t{6:f}\t{7:f}\t{8:f}\t{9:f}\t'\
				  '{10:f}\t{11:f}\t{12:f}\t{13:f}\t{14:f}\t{15:f}\t{16:f}\t{17:f}\t{18:f}\t{19:f}\t'\
				  '{20:f}\t{21:f}\t{22:f}\t{23:f}\t{24:f}\t{25:f}\t{26:f}\t{27:f}\t{28:f}\t{29:f}\t{30:f}\t{31:d}\n'.format(
		header_DRS['HIERARCH '+instr_key+' OBS TARG NAME'], RA_deg, DEC_deg,
		file_rad, file_mask, header_DRS['HIERARCH TNG DRS BJD'],
		ccf_pams_DRS[0], ccf_pams_DRS[1], ccf_pams_DRS[2], ccf_pams_DRS[3],
		ccf_pams_COR[0], ccf_pams_COR[1], ccf_pams_COR[2], ccf_pams_COR[3],
		header_DRS['HIERARCH '+instr_key+' DRS CCF NOISE'],
		header_DRS['HIERARCH '+instr_key+' DRS SPE EXT SN20'],
		header_DRS['HIERARCH '+instr_key+' DRS SPE EXT SN40'],
		header_DRS['HIERARCH '+instr_key+' DRS SPE EXT SN60'],
		header_sky_SN20,
		header_sky_SN40,
		header_sky_SN60,
		fluxA, fluxB,
		header_DRS['HIERARCH '+instr_key+' DRS BERV'],target_dist,
		MOON[instr_key]['El_app'][idx],MOON[instr_key]['airmass'][idx],MOON[instr_key]['APmag'][idx],
		MOON[instr_key]['surf_bright'][idx],MOON[instr_key]['illum'][idx],MOON[instr_key]['RV_Obs'][idx], MOON_flag))
	#print target_dist, (ccf_pams_COR-ccf_pams_DRS)*1000., CCF_noise*1000. , BERV-ccf_pams_DRS[0]

fileout.close()
# hdulist[0].header['HIERARCH TNG DRS CCF NOISE']
