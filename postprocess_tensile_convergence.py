import numpy as np
import sys
import math
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 22})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 100
figsize = (11, 6)

sim = np.array([1,  2,   3,  4, 5])
N   = np.array([40, 50, 60, 70, 80])

dx = 1.0 / N
E = np.zeros(len(sim))

tauxx   = np.zeros(len(sim))
tauxy   = np.zeros(len(sim))

for n in range(0, len(sim)):

    folder_name = "/home/blatny/repos/larsiempm/build/dumps/threedim-conv-" + str(sim[n]) + "/"
    file_name = 'out_part_frame_'
    Youngs = 1e5
    num_frames = 21
    L0 = 1 # will be overwritten

    epsilon = np.zeros(num_frames)
    tauyy   = np.zeros(num_frames)

    for i in range(0, num_frames):
        full_path = folder_name + file_name + str(i) + ".csv"
        data = np.loadtxt(full_path, skiprows=1, delimiter=',')
        y = data[:,1]

        if i==0:
            L0 = np.max(y)-np.min(y)

        L = np.max(y)-np.min(y)
        epsilon[i] = np.log( L / L0 )

        sum = 0
        count = 0
        for j in range(0, data.shape[0]):
            if y[j] < 0.95 and y[j] > 0.05:
                sum += data[j,15]
                count += 1
        tauyy[i] = sum / count

        if i == (num_frames-1):
            tauxx[n] = np.max(np.abs(data[:,11]))
            tauxy[n] = np.max(np.abs(data[:,12]))

    E[n] = np.polyfit(epsilon, tauyy, deg=1)[0]

plt.figure(figsize=figsize, dpi=dpi)
plt.plot(dx, np.abs(E-Youngs)/Youngs,'bo')
plt.xlabel(r'$\Delta x$ [m]')
plt.ylabel(r'$|E_{m}-E|/E$ [-]')
#plt.yticks([0, 0.0005, 0.001, 0.0015])
#plt.legend()
plt.savefig("tensile_convergence.eps", bbox_inches = 'tight')
plt.close()

plt.figure(figsize=figsize, dpi=dpi)
plt.plot(dx, tauxx,'bo', label=r'$\tau_{xx}$')
plt.plot(dx, tauxy,'bo', label=r'$\tau_{xy}$')
plt.xlabel(r'$\Delta x$ [m]')
plt.ylabel(r'$\tau_{xx/xy}$ [Pa]')
plt.legend()
plt.savefig("tensile_convergence_xx_xy.eps", bbox_inches = 'tight')
plt.close()
