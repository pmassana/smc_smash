import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import astropy
import sys
from matplotlib.path import Path
import emcee
import corner

def read_fits (name):
    hdulist = fits.open(name)
    field_data = hdulist[1].data
    header = hdulist[1].header
    hdulist.close()
    return field_data, header

def dereddener(mag, filter_index, EBV, field_number, filterset = 'decam'):

    if field_number == '5' or field_number == '6' or field_number == '9' or field_number == '10' or field_number == '11':
        EBV[np.where(EBV > 0.3)] = 0.037 #This line is to ensure we don't have unrealistic extinction values in the SMC centre. The limit is debatable.
    if filterset == 'decam':
        return mag - extinction_coeff_decam[filter_index]*EBV
    else:
        return mag - extinction_coeff_sdss[filter_index]*EBV

def isochrone_plotter(isochrone_points, mag_g, mag_i, distance_modulus, field_number, radius, PA):
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams.update({'font.size':10})
    fig = plt.figure(1, figsize=(4,5))
    ax = fig.add_subplot(111)

    ax.scatter(isochrone_points[0], isochrone_points[1]+distance_modulus, zorder = 10, c='red', s=3.,
        label='PARSEC isochrone')
#    ax.plot(isochrone_points[0], isochrone_points[1]+distance_modulus, zorder = 10, c='red', lw=4.,
#        label='PARSEC isochrone')

    #ax.hist2d(mag_g[mask_bool] - mag_i[mask_bool], mag_g[mask_bool], bins = [xbins, ybins], cmin = 1, zorder = 1)
    counts, xedges, yedges, imag = ax.hist2d(mag_g - mag_i, mag_g, bins = [xbins, ybins], cmin = 1, zorder = 1, cmap='Greys',
        norm=colors.LogNorm(), label = 'SMASH')

    ax.annotate(f'Field {field_number} \n $r = {radius:.1f}$ \n PA$\,={PA:.1f}$', (1.2,16.5))
    #vertex_list = [list(x) for x in vertex_coord]
    #vertex_array = np.array(vertex_list)
    #polygon = matplotlib.patches.Polygon(vertex_array, alpha = 0.3, edgecolor = 'black', facecolor = 'orange', linestyle = '-', linewidth = '2.5', zorder = 2)
    #ax.add_patch(polygon)
    ax.set_xlabel('$g-i$')
    ax.set_ylabel('$g$')
    plt.gca().invert_yaxis()
    #cbax = fig.add_axes([0.05,0.9,0.8,0.05])
    plt.colorbar(imag, orientation = 'horizontal', panchor=(0.5, 0.0), aspect=90)
    plt.minorticks_on()
    ax.legend(loc='upper left')
    #plt.grid(which='both')
    plt.savefig('/user/HS128/pm00518/Documents/PhD/CMD_SMASH/Plots/CMD_%s_isochrone.pdf' % (field_number), bbox_inches='tight', dpi=100)
    plt.close()


#field_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','139','141','142','143','144','145','147','149','150','176','177','178','181']
#field_list = ['13']
field_list_file = np.loadtxt('/user/HS128/pm00518/Documents/PhD/CMD_SMASH/SB_data_errors_dartmouth_decam.dat')
#print(field_list.T)
#sys.exit()
field_list = field_list_file.T[0].astype(int)
field_radius = field_list_file.T[1]
field_PA = field_list_file.T[2]

xbinsize = 0.02
ybinsize = 0.06
xbins = np.arange( start = -1, stop = 2, step = xbinsize)
ybins = np.arange( start = 15, stop = 24, step = ybinsize)

# Respectively u,g,r,i,z. Ex: u0 = u - 4.239*E(B-V)
extinction_coeff_sdss = np.array([4.239, 3.303, 2.285, 1.698, 1.263])
# Respectively u,g,r,i,z,Y. Ex: A_u = R_u*EBV_SFD98 ; u0 = u - A_u
extinction_coeff_decam = np.array([3.9631, 3.1863, 2.1401, 1.5690, 1.1957, 1.0476])

isochrones_data = np.loadtxt('/scratch/Isochrones/PARSEC_table.dat')

metallicity = -1.
log_age = 9.5

met_isochrones = isochrones_data[np.around(isochrones_data.T[1], 2)==metallicity]
final_isochrone = met_isochrones[np.around(met_isochrones.T[2], 2)==log_age]
g_i_isochrone = final_isochrone[:,12] - final_isochrone[:,14]
g_isochrone = final_isochrone[:,12]
isochrone_points = np.stack([g_i_isochrone, g_isochrone])

for field, radius, PA in zip(field_list, field_radius, field_PA):
    if field != 13:
        continue
    if field == 4 or field == 9 or field == 1:
        name_field_data = '/scratch/SMASH_DATA/Field%s_allobj_stars.fits' % (field)
    else:
        name_field_data = '/scratch/SMASH_DATA/Field%s_allobj_stars.fits.gz' % (field)
    smash_data, header = read_fits(name_field_data)
    #smash_data, smash_header = read_fits('/scratch/SMASH_DATA/Field%s_allobj_stars.fits.gz' % (field))
    mag_g = dereddener(smash_data['G'], 1, smash_data['EBV'], field)
    mag_i = dereddener(smash_data['I'], 3, smash_data['EBV'], field)
    distance_modulus = 18.9
    isochrone_plotter(isochrone_points, mag_g, mag_i, distance_modulus, field, radius, PA)
    print('Finished with field', field)
