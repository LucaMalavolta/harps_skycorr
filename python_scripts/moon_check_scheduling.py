import numpy as np
import matplotlib.pyplot as plt

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u


print(' ** Reading MOON data tables - it may be slow!')


MOON_TNG = np.genfromtxt('JPL_moon_TNG_2012_2022.csv', names = "BJD, flag1, flag2, RA, DEC, RA_app, DEC_app, Az_app, El_app, airmass, mag_ex, APmag, surf_bright, illum, RV_Sun, RV_Obs, none", delimiter=',', skip_header=167)


night_date = '2021-01-23' # night of observation
night_start = '20:00:00.000' # beginning of the observations
night_duration = 4.00 # in hours

epoch_start = Time([night_date+'T'+night_start], format='isot', scale='utc')
epoch_interval = [epoch_start.jd[0], epoch_start.jd[0]+night_duration/24.]




#TNG_location = EarthLocation.from_geodetic(lat=0.0000*u.deg, lon=0.0000*u.deg, height=0*u.m)
TNG_location = EarthLocation.of_site('Roque de los Muchachos')

target = SkyCoord(ra=37.765299*u.deg, dec=8.381625*u.deg)
#target = SkyCoord("23:13:58.76", "+08:45:40.57", unit=(u.hourangle, u.deg), frame='icrs')

target_rv = 12.63074




""" extracting Moon info
  1) identification of the closest (in time) tabulated data
  Mixing BJD with UTC, don't tell anyone I did this, I'll deny!!!
"""
idx_start = np.argmin(np.abs(MOON_TNG['BJD']-epoch_interval[0]))
idx_end = np.argmin(np.abs(MOON_TNG['BJD']-epoch_interval[1]))

""" 2) Picking the two closest tabulated data (one before and one after the observed BJD)
"""
if epoch_interval[0]>MOON_TNG['BJD'][idx_start]:
    idx_start -= 1

if epoch_interval[1]<MOON_TNG['BJD'][idx_end]:
    idx_end += 1


idx_length = idx_end - idx_start + 1


print('            BJD                      iso    dist   elev   illu              berv           Delta_RV    ')

for ii in range(0, idx_length, 1):
    idx = idx_start + ii

    bjd_time =  Time(MOON_TNG['BJD'][idx], format='jd', scale='utc')
    berv = target.radial_velocity_correction(obstime=bjd_time, location=TNG_location)
    deltaRV = target.radial_velocity_correction(obstime=bjd_time, location=TNG_location) - target_rv*(u.km/u.s)
    moon_distance = target.separation(SkyCoord(ra=MOON_TNG['RA_app'][idx], dec=MOON_TNG['DEC_app'][idx],frame='icrs', unit='deg')).deg

    print('{0:15.4f}  {1}  {2:6.2f} {3:6.2f} {4:6.2f}  {5} {6} '.format(MOON_TNG['BJD'][idx],
                                                               bjd_time.isot,
                                                               moon_distance, MOON_TNG['El_app'][idx],
                                                               MOON_TNG['illum'][idx],
                                                                berv.to_string(u.km/u.s), deltaRV.to_string(u.km/u.s), ))
