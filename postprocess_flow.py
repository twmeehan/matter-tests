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

directory = "/media/blatny/harddrive4/larsie/"
names = [
# "slopestart_mccmui_a13/",
# "slopestart_mccmui_a15/",
# "slopestart_mccmui_a17/",
# "slopestart_mccmui_a19/",
# "slopestart_mccmui_a21/",
# "slopestart_mccmui_a23/",
# "slopestart_mccmui_a25/",
# "test_pmccmuihardv2_xi50_b0/",
# "test_pmccmuihardv2_xi50_b01/",
# "test_pmccmuihardv2_xi50_b02/",
# "test_pmccmuihardv2_xi50_b04/",
# "test_pmccmuihardv2_xi50_b04_a50/",
# "test_pmccmuihardv2_xi50_b06/",
# "test_pmccmuihardv2_xi50_b06_a60/",
# "feeder_h005_a10/",
# "feeder_h005_a12/",
# "feeder_h005_a14/",
# "feeder_h005_a16/",
# "feeder_h005_a18/",
# "feeder_h005_a20/",
# "feeder_h005_a22/",
# "feeder_h005_a24/",
# "feeder_h005_a25/",
# "feeder_h005_a27/",
# "feeder_h005_a29/",
# "feeder_h005_a30/",
# "feeder_h005_a35/",
"feeder_h005_a0/",
]

folders = []
for i in range(0, len(names)):
    folders.append(directory + names[i])

print("USING THE FOLLOWING SIMS:")
print(folders)

###########################

def get_xmax(folder):
    info = np.loadtxt(folder + "info.txt")
    end_frame = np.loadtxt(folder + "last_written.txt").astype(int)
    # end_frame = 10

    frame_dt  = 1.0 / float(info[1])
    dx              = float(info[2])
    Np              = int(info[7])

    frames = np.arange(1, end_frame+1)
    time = frame_dt * frames

    xp_max = np.zeros(len(frames))
    for frame in frames:
        print("frame = ", frame, " / ", frames[-1])
        ### ALT 1:
        # xp_max[frame] = np.max( np.loadtxt(folder + "out_part_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0)) )
        ### ALT 2:
        xp_max[frame-1] = np.max( np.loadtxt(folder + "out_grid_frame_"+str(frame)+".csv", delimiter=",", skiprows=1, usecols=(0)) ) - 2.5*dx
    return (xp_max, time, frames)


procs = len(folders)
pool = mp.Pool(processes=procs)
outputs = [pool.apply_async( get_xmax, args=(folder,) ) for folder in folders]
xp_max_list = [out.get()[0] for out in outputs]
time_list   = [out.get()[1] for out in outputs]
frames_list = [out.get()[2] for out in outputs]


fig = plt.figure(figsize = figsize, dpi = dpi)
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
for n in range(0, len(names)):
    ax1.plot(time_list[n], xp_max_list[n], '.-', label=names[n])
    ax2.plot(frames_list[n], xp_max_list[n], '.-')
ax1.set_xlabel(r'Time')
ax2.set_xlabel(r'Frame')
ax1.set_ylabel(r'$x_{max}$', rotation=0, labelpad=10)
ax1.legend(prop={'size': legendsize})
plt.grid()
plt.show()
