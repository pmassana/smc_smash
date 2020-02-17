import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from astropy.io import fits
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec


def read_fits(name):
    hdulist = fits.open(name)
    field_data = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()
    return field_data


SC_data = np.loadtxt('/user/HS128/pm00518/Documents/PhD/CMD_SMASH/SC_16Feb_2020_coordinates.dat').T
#SCs, SCs_error = SC_data[], SC_data[]

belokurov_data = read_fits('/scratch/CloudsInArms_data/lmc_den.fits')

#Plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':18})

#fig = plt.figure(1, figsize=(7,10.5))
#gs = GridSpec(1, 2)
