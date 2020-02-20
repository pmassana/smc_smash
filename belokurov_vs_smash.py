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


#plt.ion()

def phi12_rotmat(alpha,delta):
    '''
    Converts coordinates (alpha,delta) to ones defined by a rotation matrix R_phi12_radec, applied on the original coordinates

    Critical: All angles must be in degrees
    '''
    R_phi12_radec = np.array([[ 0.04842439,  0.30207466, -0.95205356],[-0.99433193,  0.10490506, -0.01728974],[ 0.09465244,  0.9474945 ,  0.30544245]])

    vec_radec = np.array([np.cos(alpha*np.pi/180.)*np.cos(delta*np.pi/180.),np.sin(alpha*np.pi/180.)*np.cos(delta*np.pi/180.),np.sin(delta*np.pi/180.)])

    vec_phi12 = np.zeros(np.shape(vec_radec))

    vec_phi12[0] = np.sum(R_phi12_radec[0][i]*vec_radec[i] for i in range(3))
    vec_phi12[1] = np.sum(R_phi12_radec[1][i]*vec_radec[i] for i in range(3))
    vec_phi12[2] = np.sum(R_phi12_radec[2][i]*vec_radec[i] for i in range(3))

    vec_phi12 = vec_phi12.T

    vec_phi12 = np.dot(R_phi12_radec,vec_radec).T

    phi1 = np.arctan2(vec_phi12[:,1],vec_phi12[:,0])*180./np.pi
    phi2 = np.arcsin(vec_phi12[:,2])*180./np.pi


    return [phi1,phi2]

# DON'T REUSE THIS FUNCTION.
def read_fits(name):
    hdulist = fits.open(name)
    field_data = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()
    return field_data

field_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','62','130','131','132','133','134','135','136','137','138','139','140','141','142','143','144','145','147','149','150','176','177','178','179','180','181']
inner_field_list = ['3','4','5','6','7','9','10','11','12','14','15','16','178']
outskirts_field_list = ['1','2','8','13','18','19','20','21','22','23']

SC_data = np.loadtxt('/user/HS128/pm00518/Documents/PhD/CMD_SMASH/SC_16Feb_2020_coordinates.dat').T
SCs, SCs_error = SC_data[3], SC_data[4]
field_number, field_radius = SC_data[0], SC_data[1]
ms_l, ms_b = SC_data[5], SC_data[6]

belokurov_data = read_fits('/scratch/CloudsInArms_data/lmc_den.fits')

#Plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':14})

fig = plt.figure(1, figsize=(9,4.8))
gs = GridSpec(1, 2)
gss = gs[0].subgridspec(2,1, height_ratios = [0.05,0.95], hspace=0.)

belokurov_axis = fig.add_subplot(gs[1])
colorbar_axis = fig.add_subplot(gss[0])
profile_axis = fig.add_subplot(gss[1])

smc_centre = SkyCoord.from_name('SMC')
lmc_centre = SkyCoord.from_name('LMC')
centres_msl, centres_msb = phi12_rotmat(np.array([smc_centre.ra.deg, lmc_centre.ra.deg]),np.array([smc_centre.dec.deg, lmc_centre.dec.deg]))

belokurov_axis.imshow(belokurov_data, extent=(-35,35,35,-35), cmap='Greys', norm=colors.LogNorm(vmax=140))
belokurov_axis.set_xlim(20,-26)
belokurov_axis.set_ylim(-30,22)
fourdeg_ring = patches.Circle((centres_msl[0],centres_msb[0]), radius=4., fill = False, linestyle='--', ec= 'w', lw = 1.)
eightdeg_ring = patches.Circle((centres_msl[0],centres_msb[0]), radius=8., fill = False, linestyle='--', ec ='k', lw = 1.)
belokurov_axis.add_patch(fourdeg_ring)
belokurov_axis.add_patch(eightdeg_ring)
scatter_fields = belokurov_axis.scatter(ms_l[field_radius < 15], ms_b[field_radius < 15], marker = 'H', s=130, linewidths=1., fc='none', ec='red')
for a, b, c in zip(field_number[field_radius < 15], ms_l[field_radius < 15], ms_b[field_radius < 15]):
    if str(a) not in inner_field_list:
        if a in [1,2,13,18,19,176,177]:
            belokurov_axis.annotate(str(int(a)),(b,c), size=8.5, ha='center', va='center',c='w')
        else:
            belokurov_axis.annotate(str(int(a)),(b,c), size=8, ha='center', va='center')

belokurov_axis.set_xlabel('MS$_{L} \, (\deg)$')
belokurov_axis.set_ylabel('MS$_{B} \, (\deg)$')
belokurov_axis.minorticks_on()
belokurov_axis.yaxis.tick_right()
belokurov_axis.yaxis.set_label_position('right')

color_circle = np.ones((256,3))*60
color_circle[:,1] = np.ones((256))*45
color_circle[:,2] = np.arange(0,360,360.0/256.0)
color_circle_rgb = cspace_convert(color_circle, "JCh","sRGB1")
cm2 = colors.ListedColormap(color_circle_rgb)

background = -7.
background_lower = -63.
background_upper = 25.3
plot = profile_axis.scatter(field_radius[field_radius < 15], SCs[field_radius < 15], c = SC_data[2][field_radius < 15], cmap=cm2, vmin = 0, vmax = 360, s = 50, zorder=10)
profile_axis.errorbar(field_radius[field_radius < 15], SCs[field_radius < 15], yerr = SCs_error[field_radius < 15], fmt='none', zorder=1, c='k', elinewith=0.7, alpha=0.7, capsize=2)
background_area = patches.Rectangle((0.0, background-abs(background_lower)), width=16., height=abs(background_lower)+abs(background_upper), alpha = 0.3, color = 'k', edgecolor = None, lw=0.)
profile_axis.add_patch(background_area)
profile_axis.hlines(background, xmin=0, xmax=25, linestyles='--', lw=1.)
for field_number, radius, sb in zip(SC_data[0], SC_data[1], SCs):
    if field_number in [8,17,23]:
        profile_axis.annotate(str(int(field_number)), (radius+0.2, sb), fontsize=11, ha='left')
    if field_number==142:
        profile_axis.annotate('142+\n 144', (radius+0.2, sb), fontsize=11, ha='left')

profile_axis.set_xlim(xmin=3.5, xmax=15.5)
profile_axis.set_ylim(ymin=-90., ymax=3e4)
profile_axis.set_yscale('symlog', linthreshy=100., linscaley=0.5, subsy=[2, 3, 4, 5, 6, 7, 8, 9])
#profile_axis.set_yscale('log')
profile_axis.set_xlabel('$r$ (deg)')
profile_axis.set_ylabel('Star Density (stars/deg$^{2}$)')
profile_axis.minorticks_on()

cb = plt.colorbar(plot, ax=profile_axis, cax=colorbar_axis, orientation ='horizontal')
colorbar_axis.set_xlabel('Position Angle')
colorbar_axis.xaxis.set_label_position('top')
colorbar_axis.xaxis.tick_top()
colorbar_axis.minorticks_on()
cb.set_ticks([  0., 90., 180., 270., 360.])

plt.tight_layout()

plt.savefig('/user/HS128/pm00518/Documents/PhD/CMD_SMASH/Plots/SMASH_vs_belokurov2018.pdf', overwrite=True)
