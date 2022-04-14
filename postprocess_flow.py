import numpy as np
import matplotlib.pyplot as plt
import sys
import multiprocessing as mp

plt.rcParams.update({'font.size': 12})
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
dpi = 200
figsize = (5, 5.4)
legendsize = 5

###########################

folders = np.array([
# "perzynaMCC/",
# "perzynaMCCnew/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.4_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.5_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.7_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.8_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M1.2_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M2.4_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a10_f0.3_M0.6_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a20_f0.3_M0.6_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a40_f0.3_M0.6_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi0.001/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi0.01/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi0.05/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi0.1/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi1/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi2/",
# "/media/blatny/harddrive4/larsie/mcchard_a30_f0.3_M0.6_xi5/",
"/media/blatny/harddrive4/larsie/mcchard_a20_f0.3_M0.6_xi0.05/",
"/media/blatny/harddrive4/larsie/mcchard_a20_f0.3_M0.6_xi0.05_r0.001/",
])

###########################

def get_xmax(folder):
    info = np.loadtxt(folder + "info.txt")
    end_frame = np.loadtxt(folder + "last_written.txt").astype(int)
    # end_frame = 10

    frame_dt  = 1.0 / float(info[1])
    dx              = float(info[2])
    Np              = int(info[7])

    frames = np.arange(0, end_frame+1)
    time = frame_dt * frames

    xp_max = np.zeros(end_frame+1)
    for frame in frames:
        print("frame = ", frame, " / ", frames[-1])
        # data = np.loadtxt(folder + "out_part_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0,1,3,4) )
        # xp  = data[:,0]
        # yp  = data[:,1]
        # vxp = data[:,2]
        # vyp = data[:,3]
        # vp = np.sqrt(vxp**2 + vyp**2)
        # xp_max[frame] = np.max(xp)
        xp_max[frame] = np.max( np.loadtxt(folder + "out_part_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0)) )
    return (xp_max, time, frames)


procs = len(folders)
pool = mp.Pool(processes=procs)
outputs = [pool.apply_async( get_xmax, args=(folder,) ) for folder in folders]
xp_max_list = [out.get()[0] for out in outputs]
time_list   = [out.get()[1] for out in outputs]
frames_list = [out.get()[2] for out in outputs]


# xp_max_list = []
# time_list = []
# frames_list = []
#
# for n in range(0, len(folders)):
#     folder = folders[n]
#     print(folder)
#
#     info = np.loadtxt(folder + "info.txt")
#     end_frame = np.loadtxt(folder + "last_written.txt").astype(int)
#     # end_frame = 10
#
#     frame_dt  = 1.0 / float(info[1])
#     dx              = float(info[2])
#     Np              = int(info[7])
#
#     frames = np.arange(0, end_frame+1)
#     time = frame_dt * frames
#
#     xp_max = np.zeros(end_frame+1)
#     for frame in frames:
#         print("frame = ", frame, " / ", frames[-1])
#         # data = np.loadtxt(folder + "out_part_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0,1,3,4) )
#         # xp  = data[:,0]
#         # yp  = data[:,1]
#         # vxp = data[:,2]
#         # vyp = data[:,3]
#         # vp = np.sqrt(vxp**2 + vyp**2)
#         # xp_max[frame] = np.max(xp)
#         xp_max[frame] = np.max( np.loadtxt(folder + "out_part_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0)) )
#
#     xp_max_list.append(xp_max)
#     time_list.append(time)
#     frames_list.append(frames)

fig = plt.figure(figsize = figsize, dpi = dpi)
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
for n in range(0, len(folders)):
    ax1.plot(time_list[n], xp_max_list[n], '.-', label=folders[n])
    ax2.plot(frames_list[n], xp_max_list[n], '.-')
ax1.set_xlabel(r'Time')
ax2.set_xlabel(r'Frame')
ax1.set_ylabel(r'$x_{max}$', rotation=0, labelpad=10)
ax1.legend(prop={'size': legendsize})
plt.grid()
plt.show()
