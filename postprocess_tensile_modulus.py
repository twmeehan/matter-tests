import numpy as np
import sys
import math
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 22})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 100
figsize = (11, 6)
# from matplotlib.font_manager import FontProperties
# fontP = FontProperties()
# fontP.set_size('x-small')

folder_name = "/home/blatny/repos/larsiempm/build/dumps/threedim-basetest_E1e5_short/"
file_name = 'out_part_frame_'
Youngs = 1e5
num_frames = 21

epsilon = np.zeros(num_frames)
tauyy   = np.zeros(num_frames)

L0 = 1 # will be changed

for i in range(0, num_frames):
    full_path = folder_name + file_name + str(i) + ".csv"
    data = np.loadtxt(full_path, skiprows=1, delimiter=',')
    y = data[:,1]

    if i==0:
        L0 = np.max(y)-np.min(y)

    L = np.max(y)-np.min(y)

    print("L = ", L)

    epsilon[i] = np.log( L / L0 )

    sum = 0
    count = 0
    N = data.shape[0]
    for j in range(0, N):
        if y[j] < 0.95 and y[j] > 0.05:
            sum += data[j,15]
            count += 1
    tauyy[i] = sum / count

    print("tauyy = ", tauyy[i])



plt.figure(figsize=figsize, dpi=dpi)
plt.plot(epsilon, tauyy/1e3,'bo', label="MPM")
plt.plot((0,np.max(epsilon)),(0,np.max(epsilon)*Youngs/1e3),'k-', label="Theoretical")
plt.xlabel(r'$\epsilon_{yy}$ [-]')
plt.ylabel(r'$\tau_{yy}$ [kPa]')
#plt.xticks([0,0.05,0.1,0.15,0.2])
#plt.xlim([0,0.2])
#plt.ylim([0,2000])
plt.locator_params(nbins=6)
plt.legend()
plt.savefig("tensile_modulus.eps", bbox_inches = 'tight')
plt.close()
