import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update({'font.size': 18})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 200
figsize = (13, 7)
legendsize = 10

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

################################################

plate_vel       = 0.0
slope_angle_deg = 28

# name = "feeder_ramp_nolid_a"+str(slope_angle_deg)
name = "feeder_ramp_nolid_a"+str(slope_angle_deg)+"_dx3"

folder    = "/media/blatny/harddrive4/larsie/" + name + "/"
frames = np.arange(50, 250+1, 50)
# frames = np.array([250])
x_slices = np.linspace(0.2, 8, num=300)

height_fixed = -1 # put negative if not fixed

in_numb_ref    = 0.279
grain_diameter = 0.7e-3
rho            = 1550
rho_s          = 2500
theta_1        = np.radians(20.9)
theta_2        = np.radians(32.76)
slope_angle = np.radians(slope_angle_deg)

################################################

fac_Q = in_numb_ref / (grain_diameter * np.sqrt(rho_s))
v_max_pred_prefac = (2.0/3.0) * fac_Q * ( np.tan(slope_angle) - np.tan(theta_1) ) / ( np.tan(theta_2) - np.tan(slope_angle) ) * np.sqrt(rho * 9.81 * np.cos(slope_angle))
v_max_pred_prefac = max(0, v_max_pred_prefac)

info = np.loadtxt(folder + "info.txt")

frame_dt  = 1.0 / float(info[1])
dx              = float(info[2])
Np              = int(info[7])
print("dx = ", dx)
print("Np = ", Np)
print("v_max_pred_prefac = ", v_max_pred_prefac)

if (height_fixed >= 0):
    v_max_pred = v_max_pred_prefac * height_fixed**(3.0/2.0)
    print("-----------------------------------------")
    print("Predicted v_max       = ", v_max_pred)
    print("-----------------------------------------")

################################################

cmap = plt.cm.get_cmap("rainbow", 1000)
norm = matplotlib.colors.Normalize(vmin=min(frames)*frame_dt, vmax=max(frames)*frame_dt)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

vmax_list      = np.zeros((len(x_slices), len(frames)))
vmax_pred_list = np.zeros((len(x_slices), len(frames)))
height_list    = np.zeros((len(x_slices), len(frames)))

index_f = 0
for frame in frames:

    print("Frame = ", frame, " / ", frames[-1])
    data = np.loadtxt(folder + "out_grid_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0,1,3,4) )

    index_x = 0
    for x_slice in x_slices:
        xg  = data[:,0]

        inds = np.argwhere(np.abs(xg - x_slice) < 0.5*dx).flatten();
        print("    Num of nodes in y-dir = ", len(inds))
        if len(inds) == 0:
            continue

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

        for i in range(0,1000):
            try:
                if (vxg[-1] < 1e-5+plate_vel):
                    yg  =  yg[0:-1]
                    vg  =  vg[0:-1]
                    vxg = vxg[0:-1]
            except IndexError:
                yg  = np.zeros(1)
                vg  = np.zeros(1)
                vxg = np.zeros(1)

        if (height_fixed < 0):
            try:
                height = yg[-1] - 1.5*dx
            except IndexError:
                height = 0
        else:
            height = height_fixed
        print("    Height                = ", height)

        v_max_pred = v_max_pred_prefac * height**(3.0/2.0)
        print("    Predicted v_max       = ", v_max_pred)

        try:
            vmax_list[index_x, index_f]      = vxg[-1]
            vmax_pred_list[index_x, index_f] = v_max_pred
            height_list[index_x, index_f]    = height
        except IndexError:
            print("Lars IndexError")

        index_x += 1
    index_f += 1



fig, (ax1, ax2) = plt.subplots(2, figsize=figsize, dpi=dpi)
fig.suptitle(name.replace("_", " "))

for f in range(0, len(frames)):
    ax1.plot(x_slices, vmax_list[:,f], '-',   color=cmap(norm(frames[f]*frame_dt)), label=r'$t = $ ' + str(round(frames[f]*frame_dt,3)) + r' s')
    # plt.plot(x_slices, vmax_pred_list[:,f], '-', color='k')
    ax2.plot(x_slices, height_list[:,f], '-', color=cmap(norm(frames[f]*frame_dt)) )

plt.xlabel(r'$x$ (m)')
ax1.set_ylabel(r'$v_x^{max}$ (m/s)', rotation=0, labelpad=45)
ax2.set_ylabel(r'$h$ (m)', rotation=0, labelpad=35)
ax1.grid()
ax2.grid()
ax1.set_ylim(bottom=0)
ax2.set_ylim(bottom=0, top=0.05)
ax1.set_xlim(left=0)
ax2.set_xlim(left=0)

ax1.legend(prop={'size': legendsize})
# plt.show()

plt.savefig(folder + "flow_vmax.png", bbox_inches = 'tight')
plt.close()
