import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 200
figsize = (6, 6)
#################################

basename = "/media/blatny/harddrive4/larsie/3d_ql_anal_rho300_"


names = ["xi0.3"]
vel = -0.2
Ly = 2.0

#################################

eps_ax_list = []
mean_p_list = []
mean_q_list = []
for name in names:

    info = np.loadtxt(basename + name + "/info.txt")
    Np = np.loadtxt(basename + name + "/out_part_frame_0.csv", delimiter=",", skiprows=1).shape[0]

    end_frame = 400 # int(info[0])
    frame_dt  = float(info[1])
    final_time = frame_dt * end_frame
    dx = float(info[2])

    frames = np.arange(0, end_frame+1) # last argument not included
    eps_ax = vel * frames*frame_dt / Ly

    mean_p = np.zeros(end_frame+1)
    mean_q = np.zeros(end_frame+1)

    for f in frames:
        pq_data = np.loadtxt(basename + name + "/out_pq_frame_" + str(f) + ".csv", delimiter=",")
        mean_p[f] = pq_data[0]
        mean_q[f] = pq_data[1]

    eps_ax_list.append(eps_ax)
    mean_p_list.append(mean_p)
    mean_q_list.append(mean_q)


plt.figure(figsize = figsize, dpi = dpi)

plt.subplot(211)
plt.title(basename.replace("_", " "))
for n in range(0, len(names)):
    plt.plot(-eps_ax_list[n], mean_p_list[n]/1e3, '.-',  label = names[n].replace("_", " "))
plt.ylabel(r'$\langle p \rangle$ [kPa]')
plt.legend()
plt.grid()

plt.subplot(212)
for n in range(0, len(names)):
    plt.plot(-eps_ax_list[n], mean_q_list[n]/1e3, '.-', label = names[n].replace("_", " "))
plt.xlabel(r'$\epsilon_y$ [-]')
plt.ylabel(r'$\langle q \rangle$ [kPa]')
plt.legend()
plt.grid()

plt.show()
