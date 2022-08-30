import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update({'font.size': 22})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 200
figsize = (7, 5)
legendsize = 5

def colorFader(mix): #fade from color c1 (at mix=0) to c2 (mix=0.5) to c3 (mix=1)
    c1 = 'yellow'
    c2 = 'red'
    c3 = 'blue'
    c1 = np.array(matplotlib.colors.to_rgb(c1))
    c2 = np.array(matplotlib.colors.to_rgb(c2))
    c3 = np.array(matplotlib.colors.to_rgb(c3))
    if mix < 0.5:
        return matplotlib.colors.to_hex( (0.5-mix)*c1 + 2*mix*c2       )
    else:
        return matplotlib.colors.to_hex( (2-2*mix)*c2 + (2*mix-1)*c3   )

###########################
plate_vel       = 0.0
slope_angle_deg = 20
x_slice         = 0.2  # 0.1987 + 0.19
frames          = np.arange(1,601,100)

# folder = "/media/blatny/harddrive4/larsie/conveyor_MCCmui_quad_ppc25_new_xi1_a" + str(slope_angle_deg) + "/"
# folder = "/media/blatny/harddrive4/larsie/conveyor_MCCmui_quad_ppc25_new_xi1_a" + str(slope_angle_deg) + "_test_b3/"
# folder = "/media/blatny/harddrive4/larsie/pbc_flip_muimcc_ppc20_L03/"
folder = "/media/blatny/harddrive4/larsie/pbc_flip_muidp_ppc20_L05/"
###########################

info = np.loadtxt(folder + "info.txt")
# end_frame = np.loadtxt(folder + "last_written.txt").astype(int)
# end_frame = 10

frame_dt  = 1.0 / float(info[1])
dx              = float(info[2])
Np              = int(info[7])
print("dx = ", dx)
print("Np = ", Np)

######################
### Bagnold prefactor
in_numb_ref    = 1e-3;
grain_diameter = 7e-4
rho            = 1500.0
rho_s          = 2450.0
mu_1           = 0.35
mu_2           = 0.7
slope_angle    = np.radians(slope_angle_deg)

fac_Q = in_numb_ref / (grain_diameter * np.sqrt(rho_s))
v_max_pred_prefac = (2.0/3.0) * fac_Q * ( np.tan(slope_angle) - mu_1 ) / ( mu_2 - np.tan(slope_angle) ) * np.sqrt(rho * 9.81 * np.cos(slope_angle))

# C_pred = np.sqrt(9.81) * 2.0/3.0 * in_numb_ref / grain_diameter * np.sqrt(rho/rho_s * np.cos(slope_angle)) * ( np.tan(slope_angle) - mu_1 ) / ( mu_2 - np.tan(slope_angle) )
# print("C_pred = ", C_pred)
# C_pred = np.sqrt(9.81) * 2.0/3.0 * in_numb_ref /(grain_diameter * np.sqrt(rho_s)) * ( np.tan(slope_angle) - mu_1 ) / ( mu_2 - np.tan(slope_angle) ) * np.sqrt(rho * np.cos(slope_angle))
# print("C_pred = ", C_pred)
# C_pred = 2.0/3.0 * in_numb_ref / grain_diameter * ( np.tan(slope_angle) - mu_1 ) / ( mu_2 - np.tan(slope_angle) ) * np.sqrt(rho/rho_s * 9.81 * np.cos(slope_angle))
# print("C_pred = ", C_pred)
######################

print("v_max_pred_prefac = ", v_max_pred_prefac)

cmap = plt.cm.get_cmap("rainbow", 1000)
norm = matplotlib.colors.Normalize(vmin=min(frames)*frame_dt, vmax=max(frames)*frame_dt)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

fig = plt.figure(figsize = figsize, dpi = dpi)
for frame in frames:
    ax = plt.subplot(111)
    print("frame = ", frame, " / ", frames[-1])
    data = np.loadtxt(folder + "out_grid_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0,1,3,4) )

    xg  = data[:,0]
    inds = np.argwhere(np.abs(xg - x_slice) < 0.4*dx).flatten();    print("    Num of nodes in y-dir = ", len(inds))
    yg  = data[inds,1]

    vxg = data[inds,2] + plate_vel # NB adding the velocity of the ground
    vyg = data[inds,3]

    vg = np.sqrt(vxg**2 + vyg**2)

    inds_sorted = yg.argsort()
    yg = yg[inds_sorted]
    vg = vg[inds_sorted]

    yg  =  yg[2:-2]
    vg  =  vg[2:-2]
    vxg = vxg[2:-2]

    for i in range(0,90):
        if (vxg[-1] < 1e-5+plate_vel):
            yg  =  yg[0:-1]
            vg  =  vg[0:-1]
            vxg = vxg[0:-1]

    height = yg[-1]

    v_max_pred = v_max_pred_prefac * height**(3.0/2.0)
    print("    Predicted v_max       = ", v_max_pred)

    def bagnold(y, v_max_prefac): # velocity as a function of y
        return v_max_prefac * height**(3.0/2.0) * ( 1 - (1 - y/height)**(3.0/2.0) )

    fit, pcov = curve_fit(bagnold, yg, vxg);
    print("    Fitted v_max          = ", fit[0] * height**(3.0/2.0) )

    yg_equi = np.linspace(np.min(yg), np.max(yg), 1000)
    # ax.plot(bagnold(yg_equi, fit[0]), yg_equi, 'k-')
    # ax.plot(bagnold(yg_equi, v_max_pred), yg_equi, 'k--')

    #ax.plot(vxg, yg, '.-', color=colorFader( (frame-np.min(frames)) / (np.max(frames)-np.min(frames)) ), label=frame)
    ax.plot(vxg, yg, '.-', color=cmap(norm(frame*frame_dt)), label=frame)


# plt.xlabel(r'$|v|$')
plt.xlabel(r'$v_x (m/s)$')
plt.ylabel(r'$y (m)$', rotation=0, labelpad=10)
plt.grid()

cbaxes = inset_axes(ax, width="30%", height="3%", loc='upper left')
cbar = plt.colorbar(sm, cax=cbaxes, orientation='horizontal', ticks=[round(np.min(frames)*frame_dt)+1, round(np.max(frame*frame_dt))-1])
cbar.set_label(r'$t (s)$',rotation=0, labelpad=5, fontsize=22)

# plt.legend(prop={'size': legendsize})
# plt.show()

plt.savefig(folder + "flow_profile.png", bbox_inches = 'tight')
