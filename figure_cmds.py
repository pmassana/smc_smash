import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from astropy.io import fits
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec


# Respectively u,g,r,i,z. Ex: u0 = u - 4.239*E(B-V)
extinction_coeff_sdss = np.array([4.239, 3.303, 2.285, 1.698, 1.263])
# Respectively u,g,r,i,z,Y. Ex: A_u = R_u*EBV_SFD98 ; u0 = u - A_u
extinction_coeff_decam = np.array([3.9631, 3.1863, 2.1401, 1.5690, 1.1957, 1.0476])

xbinsize = 0.02
ybinsize = 0.06
xbins = np.arange( start = -1, stop = 2+xbinsize, step = xbinsize)
ybins = np.arange( start = 14, stop = 24+ybinsize, step = ybinsize)


def read_fits (name):
    hdulist = fits.open(name)
    field_data = hdulist[1].data
    header = hdulist[1].header
    hdulist.close()
    return field_data


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
    plt.savefig('/home/pol/Documents/PhD/CMD_SMASH/Plots/CMD_%s_isochrone.pdf' % (field_number), bbox_inches='tight', dpi=100)
    plt.close()


#SMASH data
field_list = ['3', '13', '18', '139']
data_path = '/home/pol/PhD_DATA/SMASH_DATA/'
nrows = 2
ncols = 2
mag_g = [None]*4
mag_i = [None]*4
for i in range(4):
    data = read_fits(data_path+f'Field{field_list[i]}_allobj_stars.fits.gz')
    mag_g[i] = dereddener(data['G'], 1, data['EBV'], field_list[i])
    mag_i[i] = dereddener(data['I'], 3, data['EBV'], field_list[i])
    print('Data loaded for field ', field_list[i])


#Isochrone data
isochrones_data = np.loadtxt('/home/pol/PhD_DATA/ISOCHRONES/PARSEC_table.dat')

old_metallicity = -1.
young_metallicity = -0.4
log_age_old = 9.9
log_age_young = 8.
distance_modulus = 18.9

met_isochrones = isochrones_data[np.around(isochrones_data.T[1], 2)==old_metallicity]
old_isochrone = met_isochrones[np.around(met_isochrones.T[2], 2)==log_age_old]
met_isochrones = isochrones_data[np.around(isochrones_data.T[1], 2)==young_metallicity]
young_isochrone = met_isochrones[np.around(met_isochrones.T[2], 2)==log_age_young]

g_i_old = old_isochrone[:,12] - old_isochrone[:,14]
g_old = old_isochrone[:,12]
old_points = np.stack([g_i_old, g_old])
g_i_young = young_isochrone[:,12] - young_isochrone[:,14]
g_young = young_isochrone[:,12]
young_points = np.stack([g_i_young, g_young])


#Plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':18})

fig = plt.figure(1, figsize=(7,10.5))
gs = GridSpec(2, 1, height_ratios=[0.03,0.97])
gss = gs[1].subgridspec(2,2, hspace=0., wspace=0.)
colorbar_axis = fig.add_subplot(gs[0])
axis1 = fig.add_subplot(gss[0,0])
axis2 = fig.add_subplot(gss[0,1])
axis3 = fig.add_subplot(gss[1,0])
axis4 = fig.add_subplot(gss[1,1])
axes_list = [axis1, axis2, axis3, axis4]

print('Starting to plot fields...')

for (i, axs) in zip(range(4), axes_list):
    counts, xedges, yedges, imag = axs.hist2d(mag_g[i] - mag_i[i], mag_g[i], bins = [xbins, ybins], cmin = 1, zorder = 1, cmap='Greys',
        norm=colors.LogNorm(), label = 'SMASH', vmax=1e2)
    if i==0 or i==2:
        axs.scatter(young_points[0], young_points[1]+distance_modulus, c='blue', s=3.,
            label=f'{10**(log_age_young):.1f} yr isochrone')
    axs.scatter(old_points[0], old_points[1]+distance_modulus, c='red', s=3.,
        label=f'{10**(log_age_old):.1f} yr isochrone')

    axs.annotate(f'Field {field_list[i]}', (1.,15))
    axs.set_xlabel('$g-i$')
    axs.set_ylabel('$g$')
    axs.set_ylim(24,14)
    axs.set_xlim(-1,2)
    axs.minorticks_on()
    if i==2:
        axs.set_yticklabels([None, 16, 18, 20, 22, 24])
    if i==3:
        axs.set_xticklabels([None, 0, 1, 2])
    axs.label_outer()
    #axs.invert_yaxis()
    axs.tick_params(direction='inout', which='both')
    #axs.legend(loc='upper left')
    #print('Field', field_list[i], 'finished.')

plt.colorbar(imag, orientation = 'horizontal', cax = colorbar_axis)
colorbar_axis.xaxis.set_ticks_position('top')
fig.tight_layout()

plt.savefig('/home/pol/Documents/PhD/CMD_SMASH/Plots/CMD_isochrones_figure.pdf', rasterized=True)


