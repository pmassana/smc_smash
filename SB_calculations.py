import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from astropy.io import fits
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec
from colorspacious import cspace_convert
from astropy.coordinates import SkyCoord
import matplotlib.patches as patches


vertex_coord = [(0.1,23.), (0.1,22.1), (0.45,22.1), (0.45,23.)]



#Isochrone data
isochrones_data = np.loadtxt('/scratch/Isochrones/PARSEC_table.dat')

old_metallicity = -1.
log_age_old = 9.9
distance_modulus = 18.9

met_isochrones = isochrones_data[np.around(isochrones_data.T[1], 2)==old_metallicity]
old_isochrone = met_isochrones[np.around(met_isochrones.T[2], 2)==log_age_old]

g_i_old = old_isochrone[:,12] - old_isochrone[:,14]
g_old = old_isochrone[:,12]
old_points = np.stack([g_i_old, g_old+distance_modulus])

msto_box = Path(vertex_coord)
msto_mask = msto_box.contains_points(old_points.T)

int_IMF = old_isochrone[:,4]
Mass = old_isochrone[:,5]
logL = old_isochrone[:,6]

total_IMF = np.max(int_IMF) - np.min(int_IMF)
cumulative_percentage_IMF_msto = (int_IMF[msto_mask] - np.min(int_IMF))/total_IMF
int_IMF_msto = np.max(int_IMF[msto_mask]) - np.min(int_IMF[msto_mask])
percentage_IMF_msto = cumulative_percentage_IMF_msto[1:]-cumulative_percentage_IMF_msto[:-1]
mean_L_msto = np.power(10, (logL[msto_mask][1:]+logL[msto_mask][:-1])/2.)


cumulative_percentage_IMF = (int_IMF - np.min(int_IMF))/total_IMF
percentage_IMF = cumulative_percentage_IMF[1:]-cumulative_percentage_IMF[:-1]
mean_L = np.power(10, (logL[1:]+logL[:-1])/2.)

print(np.sum(percentage_IMF_msto*mean_L_msto)/np.sum(percentage_IMF*mean_L))
