import numpy as np
import matplotlib.pyplot as plt
figsize = (15, 6)
plt.rcParams.update({'font.size': 14})
#################################
# Parameters
name = "micro-m65-mc9-phi026-IC-v2-pureelastic"
velocity = np.array([-0.0005, 0.0005, -0.0005, 0.0005, -0.0005, 0.0005])
L_rve    = 1.0
phi = 0.2653618570778334
#################################
# Read Info
info = np.loadtxt("dumps/" + name + "/info.txt")
Np = np.loadtxt("dumps/" + name + "/out_part_frame_0.csv", delimiter=",", skiprows=1).shape[0]
end_frame = int(info[0])
frame_dt  = float(info[1])
final_time = frame_dt * end_frame
dx = float(info[2])

print("Np        = ", Np)
print("end_frame = ", end_frame)
print("frame_dt  = ", frame_dt)
print("dx        = ", dx)

mu   = float(info[3])
lamb = float(info[4])
E = mu * (3*lamb + 2*mu) / (lamb + mu)
print("Ice Young's modulus E = ", E, " Pa")

#################################
# Get volumetric strain
frame = np.arange(0, end_frame+1) # end value not included in np.arange
time = (frame - 0.5) * frame_dt
time[0] = 0
#              T1 / A1 / B1  -  T2 / A2 / B2
eps_ax = time * (velocity[0] - velocity[1]) / L_rve
eps_sa = time * (velocity[2] - velocity[3]) / L_rve
eps_sb = time * (velocity[4] - velocity[5]) / L_rve
eps_vol = eps_ax + eps_sa + eps_sb

p = np.zeros(len(frame))
q = np.zeros(len(frame))
for f in frame:
    pq = np.loadtxt("dumps/" + name + "/out_pq_frame_" + str(f) + ".csv", delimiter=",")
    p[f] = pq[0]
    q[f] = pq[1]

p *= phi
q *= phi
# NB phi should adjusted according to J_ice and eps_vol!

plt.figure(figsize=figsize)
plt.plot(eps_vol, p, "k-*")
plt.xlabel("eps_vol [-]")
plt.ylabel("p [Pa]")
plt.savefig("dumps/" + name + "/fig_strain_p.png", bbox_inches = 'tight')

plt.figure(figsize=figsize)
plt.plot(eps_vol, q, "k-*")
plt.xlabel("eps_vol [-]")
plt.ylabel("q [Pa]")
plt.savefig("dumps/" + name + "/fig_strain_q.png", bbox_inches = 'tight')

plt.figure(figsize=figsize)
plt.plot(p, q, "k-*")
plt.xlabel("p [Pa]")
plt.ylabel("q [Pa]")
plt.savefig("dumps/" + name + "/fig_p_q.png", bbox_inches = 'tight')
