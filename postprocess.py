import numpy as np
import matplotlib.pyplot as plt
figsize = (15, 12)
#################################

name = "elastic"

#################################

info = np.loadtxt("dumps/" + name + "/info.txt")

end_frame = int(info[0])
frame_dt  = float(info[1])
final_time = frame_dt * end_frame

dx = float(info[2])

mu   = float(info[3])
lamb = float(info[4])
K = mu + lamb # in 2D
print("Expected bulk modulus K = ", K, " Pa")

frames = np.arange(0, end_frame+1) # last argument not included
max_y          = np.zeros(end_frame+1)
eps_pl_dev     = np.zeros(end_frame+1)
regularization = np.zeros(end_frame+1)
pressure       = np.zeros(end_frame+1)

for f in frames:
    particle_data = np.loadtxt("dumps/" + name + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)
    max_y[f]          = np.max(  particle_data[:,1] )
    eps_pl_dev[f]     = np.mean( particle_data[:,6] )
    regularization[f] = np.mean( particle_data[:,7] )
    pressure[f]       = np.mean( particle_data[:,8] )


#################################
L = 1.0 - dx*0.25
eps_axial = (L - max_y) / L
#################################
E_fit = np.polyfit(eps_axial, pressure, deg=1)[0]
print("Measured elastic modulus = ", E_fit, " Pa")
#################################

plt.figure(figsize = figsize)

plt.subplot(211)
plt.title(name)
plt.plot(eps_axial, pressure, 'k-*')
plt.plot(eps_axial, E_fit*eps_axial, 'r-', label="Fit: E = " + str(np.format_float_scientific(E_fit, precision=2)))
plt.legend()
plt.xlabel("eps_axial [-]")
plt.ylabel("avg(p) [Pa]")

plt.subplot(212)
plt.plot(eps_axial, eps_pl_dev, 'k-*')
plt.xlabel("eps_axial [-]")
plt.ylabel("avg(eps_pl_dev) [-]")
plt.show()
