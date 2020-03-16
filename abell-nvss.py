#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 12:02:02 2020

@author: Kabelo Kesebonye
"""

import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
from pandas import DataFrame
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import binned_statistic, chisquare
from scipy import stats


cross_cat = pd.read_csv('./1.44GHz/Abell194_NVSS_1.44GHz.csv')
#cross_cat = pd.read_csv('./1.38GHz/Abell194_NVSS_1.38GHz.csv')

meerKAT_ra = np.asarray(cross_cat['RA'])
meerKAT_dec = np.asarray(cross_cat['DEC'])
meerKAT_flux = np.asarray(cross_cat['Total_flux'])
meerKAT_flux = meerKAT_flux*1e3

nvss_ra = np.asarray(cross_cat['_RAJ2000'])
nvss_dec = np.asarray(cross_cat['_DEJ2000'])
nvss_flux = np.asarray(cross_cat['S1.4'])

flux_ratio = meerKAT_flux/nvss_flux

rcParams['font.family'] = 'monospace'
rcParams['font.size'] = '12'

'''
FLUX VS FLUX PLOT
'''
#fig = plt.figure(figsize=(10,10))
#
#plt.scatter(meerKAT_flux, nvss_flux, s=15,color='royalblue')
#plt.xscale('log')
#plt.yscale('log')
##plt.ylim(10**0,10**3)
#plt.xlabel('Abell 194 MeerKAT fluxes (mJy)')
#plt.ylabel('Abell 194 NVSS fluxes (mJy)')


meerKAT_coord = SkyCoord(ra=meerKAT_ra*u.deg, dec=meerKAT_dec*u.deg, frame='fk5')
centre = SkyCoord(ra=21.45*u.deg, dec=-1.373056*u.deg, frame='fk5')
distance = meerKAT_coord.separation(centre)
distance = distance.to(u.arcminute)

mean_stat = binned_statistic(distance.value, flux_ratio, statistic='mean',
                             bins=7, range=(0, 70))
mean_count = binned_statistic(distance.value, flux_ratio, statistic='count',
                             bins=7, range=(0, 70))
bin_count = mean_count.statistic
bin_dist = mean_stat.statistic
bin_dist = np.insert(bin_dist,0,[1])

fratio_dist = DataFrame({'ratios' : flux_ratio,
                         'distances' : distance.value})
fratio_sorted = fratio_dist.sort_values(by=['distances'])
fratios_sorted = np.asarray(fratio_sorted['ratios'])


#stdev = np.std(flux_ratio, ddof = 1)
#stderr = stdev/np.sqrt(len(flux_ratio))

array = [fratios_sorted[0:2],fratios_sorted[2:8],fratios_sorted[8:20],
         fratios_sorted[20:53],fratios_sorted[53:87],fratios_sorted[87:127],
         fratios_sorted[127:152]]

std_err = [0]
for i in array:
    err = stats.sem(i)
    std_err.append(err)

'''
PRIMARY BEAM MODEL EQUATION
'''
rho = np.linspace(0,90,170)
theta_b = (1.28/1.5)**-1*57.5
angle = (rho/theta_b)

x = 1.189*np.pi*angle
y = 1-4*(1.189*angle)**2
z = x/y

a_b = (np.cos(x)/y)**2

PB_mean_stat = binned_statistic(rho, a_b, statistic='mean',
                             bins=7, range=(0, 70))

'''
CHI-SQUARE CALCULATIONS
'''
fit_test = chisquare(mean_stat.statistic, PB_mean_stat.statistic)
PB_mean = np.insert(PB_mean_stat.statistic,0,[1])
'''
BINNED FLUXES PLOT
'''
fig = plt.figure(figsize=(10,10))

#plt.scatter(mean_stat.bin_edges, bin_dist, s=35,color='blue')
plt.errorbar(mean_stat.bin_edges-5,bin_dist,yerr=std_err,fmt='o',ms=6,color='blue')
plt.plot(rho,a_b,color='black')
#plt.scatter(PB_mean_stat.bin_edges-5,PB_mean,color='black')
plt.xlim(0,70)
plt.ylim(0,1.6)
#plt.yscale('log')
plt.xlabel('Distance from field centre (arcmin)')
plt.ylabel('Mean flux ratios (MeerKAT/NVSS)')
#plt.text(7, 0.8, '({:.0f})'.format(bin_count[0]))
#plt.text(17, 1.26, '({:.0f})'.format(bin_count[1]))
#plt.text(26, 0.86, '({:.0f})'.format(bin_count[2]))
#plt.text(36, 0.42, '({:.0f})'.format(bin_count[3]))
#plt.text(46, 0.19, '({:.0f})'.format(bin_count[4]))
#plt.text(56, 0.07, '({:.0f})'.format(bin_count[5]))
#plt.text(66, 0.027, '({:.0f})'.format(bin_count[6]))
plt.grid()
plt.show()

'''
DISTANCE VS FLUX RATIOS PLOT
'''
#fig = plt.figure(figsize=(10,10))
#
#plt.scatter(distance, flux_ratio, s=15,color='teal')
#plt.errorbar(distance.value, flux_ratio, yerr=std_err, ms=4,fmt='o',color='teal')
#plt.xlabel('Distance from field centre (arcmin)')
#plt.yscale('log')
#plt.ylim(10**-3,10**1)
#plt.ylabel('flux ratio (MeerKAT/NVSS)')
#plt.show()

#plt.scatter(meerKAT_ra,meerKAT_dec,s=45,color ='red', label='MeerKAT')
#plt.scatter(nvss_ra,nvss_dec,s=10,color ='cyan',label='NVSS')
#plt.legend(loc='upper right', fontsize = 13.5)
#plt.xlabel('RA')
#plt.ylabel('DEC')


meerKAT_fluxerr = np.asarray(cross_cat['E_Total_flux'])
meerKAT_fluxerr = meerKAT_fluxerr*1e3
nvss_fluxerr = np.asarray(cross_cat['e_S1.4'])

'''
DISTANCE VS FLUXES PLOT
'''
#fig = plt.figure(figsize=(10,10))
#
#plt.errorbar(distance.value, meerKAT_flux, yerr=meerKAT_fluxerr, ms=4, fmt='o',color='red',label='MeerKAT')
#plt.errorbar(distance.value, nvss_flux, yerr=nvss_fluxerr, ms=4,fmt='o',color='teal',label='NVSS')
#plt.xlabel('Distance from field centre (arcmin)')
#plt.yscale('log')
#plt.ylabel('Source fluxes (mJy)')
#plt.legend(loc='best',fontsize=12.5)
#plt.show()








