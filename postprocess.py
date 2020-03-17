import numpy as np
import matplotlib.pyplot as plt
figsize = (18, 12)
#################################

name = "elastic_wip"

#################################

info = np.loadtxt("dumps/" + name + "/info.txt")
Np = np.loadtxt("dumps/" + name + "/out_part_frame_0.csv", delimiter=",", skiprows=1).shape[0]

end_frame = int(info[0])
frame_dt  = float(info[1])
final_time = frame_dt * end_frame
dx = float(info[2])
mu   = float(info[3])
lamb = float(info[4])

K = mu + lamb # in 2D
E = mu * (3*lamb + 2*mu) / (lamb + mu)
print("Theoretical bulk modulus    K = ", K, " Pa")
print("Theoretical Young's modulus E = ", E, " Pa")

frames = np.arange(0, end_frame+1) # last argument not included

mean_eps_pl_dev     = np.zeros(end_frame+1)
mean_regularization = np.zeros(end_frame+1)
mean_pressure       = np.zeros(end_frame+1)

max_tau_yy     = np.zeros(end_frame+1)
max_y          = np.zeros(end_frame+1)
max_pressure   = np.zeros(end_frame+1)

tau_yy  = np.zeros((end_frame+1, Np))
Fe_yx   = np.zeros((end_frame+1, Np))
Fe_yy   = np.zeros((end_frame+1, Np))

for f in frames:
    particle_data = np.loadtxt("dumps/" + name + "/out_part_frame_" + str(f) + ".csv", delimiter=",", skiprows=1)

    mean_eps_pl_dev[f]     = np.mean( particle_data[:,6] )
    mean_regularization[f] = np.mean( particle_data[:,8] )
    mean_pressure[f]       = np.mean( particle_data[:,9] )

    max_y[f]          = np.max(  particle_data[:,1] )
    max_pressure[f]   = np.max(  particle_data[:,9] )
    max_tau_yy[f]     = np.max(  particle_data[:,14] )

    tau_yy[f,:] = particle_data[:,14]
    Fe_yx[f,:]  = particle_data[:,17]
    Fe_yy[f,:]  = particle_data[:,18]

hencky_yy = 0.5 * np.log( Fe_yx**2 + Fe_yy**2 );

#################################
# These params will depend on case
L = 1.0 - dx*0.25
macro_logstrain_yy = 2*np.log(max_y / L)
eps_axial = (L - max_y) / L
#################################

plt.figure(figsize = figsize)

plt.subplot(221)
plt.title(name)

plt.plot(macro_logstrain_yy, max_tau_yy, 'k-*')
modulus_fit = np.polyfit(macro_logstrain_yy, max_tau_yy, deg=1)[0]
plt.plot(macro_logstrain_yy, modulus_fit_2*macro_logstrain_yy, 'r-', label="Fit: modulus = " + str(np.format_float_scientific(modulus_fit, precision=2)))
plt.legend()
plt.xlabel("Macroscopic logarithmic strain [-]")
plt.ylabel("Max(tau_yy) [Pa]")

plt.subplot(222)
p = 5
plt.plot(hencky_yy[:, p], tau_yy[:, p], 'k*')
modulus_fit = np.polyfit(hencky_yy[:, p], tau_yy[:, p], deg=1)[0]
plt.plot(hencky_yy[:, p], modulus_fit_3*hencky_yy[:, p], 'r-', label="Fit: modulus = " + str(np.format_float_scientific(modulus_fit, precision=2)))
plt.legend()
plt.xlabel("hencky_p_yy [-]")
plt.ylabel("tau_p_yy [Pa]")

plt.subplot(223)
plt.plot(eps_axial, mean_pressure, 'k-*')
modulus_fit = np.polyfit(eps_axial, mean_pressure, deg=1)[0]
plt.plot(eps_axial, modulus_fit_1*eps_axial, 'r-', label="Fit: modulus = " + str(np.format_float_scientific(modulus_fit, precision=2)))
plt.legend()
plt.xlabel("Engineering axial strain [-]")
plt.ylabel("Mean pressure [Pa]")

plt.subplot(224)
plt.plot(eps_axial, mean_eps_pl_dev, 'k-*')
plt.xlabel("Engineering axial strain [-]")
plt.ylabel("Mean plastic dev strain [-]")
plt.show()
