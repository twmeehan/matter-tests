import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 12})
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
dpi = 200
figsize = (5, 5.4)
legendsize = 5
#################################

#basename = "/media/blatny/harddrive4/larsie/3d_ql_anal_"
#names = ["rho300_xi0.3", "soft_rho300_xi0.3", "soft_monforte", "soft_monforte_v3"]

basename = "/media/blatny/harddrive4/larsie/2d_ql_anal_soft_alt"
names = ["0", "1", "3", "3_m100", "3_m100_xisoft0.1", "3_m100_xisoft0.3"]

vel = -0.2
Ly = 2.0
vmin_factor = 25;
load_factor = 75;

v_min = vel / vmin_factor

#################################

# while(True):

eps_ax_list = []
mean_p_list = []
mean_q_list = []
for name in names:

    info = np.loadtxt(basename + name + "/info.txt")
    last_written = np.loadtxt(basename + name + "/last_written.txt").astype(int)

    end_frame = last_written
    frame_dt  = 1.0 / float(info[1])
    final_time = frame_dt * end_frame
    dx = float(info[2])
    frames = np.arange(0, end_frame+1) # last argument not included

    eps_ax = np.zeros_like(frames).astype(float)
    for n in range(0,len(frames)):
        if n < load_factor:
            eps_ax[n] =   v_min * frame_dt*n                                              / Ly
        else:
            eps_ax[n] = ( v_min * frame_dt*load_factor + vel * frame_dt*(n-load_factor) ) / Ly
    # print(eps_ax)

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
plt.legend(prop={'size': legendsize})
plt.grid()

# plt.xlim([-0.0001,0.008])
# plt.ylim([-1,70])

plt.subplot(212)
for n in range(0, len(names)):
    plt.plot(-eps_ax_list[n], mean_q_list[n]/1e3, '.-', label = names[n].replace("_", " "))
plt.xlabel(r'$\epsilon_y$ [-]')
plt.ylabel(r'$\langle q \rangle$ [kPa]')
plt.legend(prop={'size': legendsize})
plt.grid()

# plt.xlim([-0.0001,0.008])
# plt.ylim([-1,70])

plt.show()

# plt.show(block=False)
# plt.pause(5)
# plt.close()
