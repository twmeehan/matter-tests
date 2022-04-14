import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 12})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 200
figsize = (5, 5.4)
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

folder = "/media/blatny/harddrive4/larsie/perzynaMCC/"
x_slice = 1.5
frames = np.arange(150,601,100)

###########################

info = np.loadtxt(folder + "info.txt")
# end_frame = np.loadtxt(folder + "last_written.txt").astype(int)
# end_frame = 10

frame_dt  = 1.0 / float(info[1])
dx              = float(info[2])
Np              = int(info[7])


fig = plt.figure(figsize = figsize, dpi = dpi)

for frame in frames:
    print("frame = ", frame, " / ", frames[-1])
    data = np.loadtxt(folder + "out_grid_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0,1,3,4) )
    xg  = data[:,0]

    inds = np.argwhere(np.abs(xg - x_slice) < 0.1*dx).flatten()
    # print(inds)

    yg  = data[inds,1]
    vxg = data[inds,2]
    vyg = data[inds,3]
    vg = np.sqrt(vxg**2 + vyg**2)

    inds_sorted = yg.argsort()
    yg = yg[inds_sorted]
    vg = vg[inds_sorted]

    yg = yg[2:-2]
    vg = vg[2:-2]

    if (vg[-1] < 1e-10):
        yg = yg[0:-1]
        vg = vg[0:-1]
    if (vg[-1] < 1e-10):
        yg = yg[0:-1]
        vg = vg[0:-1]

    height = yg[-1]
    def bagnold(y, C): # velocity as a function of y
        exponent = 3.0/2.0
        return C * ( height**exponent - (height-y)**exponent )
    fit, pcov = curve_fit(bagnold, yg, vg); print("C = ", fit[0])
    yg_equi = np.linspace(np.min(yg), np.max(yg), 1000)
    plt.plot(bagnold(yg_equi, fit[0]), yg_equi, 'k-')

    plt.plot(vg, yg, '.-', color=colorFader( (frame-np.min(frames)) / (np.max(frames)-np.min(frames)) ), label=frame)

plt.xlabel(r'$|v|$')
plt.ylabel(r'$y$', rotation=0, labelpad=10)
# plt.legend(prop={'size': legendsize})
plt.grid()
plt.show()
