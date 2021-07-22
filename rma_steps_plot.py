import sys
import numpy as np
import matplotlib.pyplot as plt

if int(sys.argv[1]) == 0: # QUAD ANAL
    name = "test_rma_quad_anal"
elif int(sys.argv[1]) == 1: # QUAD ITER
    name = "test_rma_quad_iter"
elif int(sys.argv[1]) == 2: # MCC
    name = "test_rma_mcc"
else:
    print("Unknown model")

steps = np.loadtxt("dumps/" + name + "/rma_steps.txt")
info = np.loadtxt("dumps/" + name + "/plastic_info.txt")

M    = info[0]
p0   = info[1]
beta = info[2]

p_step = steps[:,1]
q_step = steps[:,2]

p_final = p_step[-1]
q_final = q_step[-1]

tol = 1e-3
p_true = np.linspace(-beta*p0+tol, p0-tol, 1000000)

p_equi = np.linspace(0.9*p_final, 1.1*p_final, 100)

if int(sys.argv[1]) == 0 or int(sys.argv[1]) == 1: # QUAD
    q_true = 2.0*M / (2*beta+1.0) * (p0-p_true)*(beta*p0+p_true) / p0
    dqdp = 2.0*M / (2*beta+1.0) * (1-beta-2*p_final/p0)
elif int(sys.argv[1]) == 2: # MCC
    q_true = M * np.sqrt( (p0-p_true)*(beta*p0+p_true) / (2*beta+1.0) )
    dqdp = M / np.sqrt((2*beta+1.0)) * (p0-beta*p0-2*p_final) / 2 / np.sqrt((p0-p_final)*(p_final+beta*p0))

intercept = q_final - dqdp * p_final
tangent = dqdp * p_equi + intercept

normal_slope = -1.0 / dqdp
intercept = q_final - normal_slope * p_final
normal = normal_slope * p_equi + intercept

plt.figure()
plt.axhline(y=0, color='k', linestyle='--')
plt.plot(p_true, q_true, 'k-', label="Yield surface")
plt.plot(p_step, q_step, 'r-*', label="Return mapping steps")
plt.plot([p_step[0], p_step[-1]], [q_step[0], q_step[-1]], 'r-.', label="Return mapping projection")
plt.plot(p_equi, normal, 'b-', label="Normal at projected point")
plt.plot(p_equi, tangent, 'g-', label="Tangent at projected point")
plt.axis('scaled')
plt.xlabel("p")
plt.ylabel("q")
plt.legend(loc="lower center")
if int(sys.argv[1]) == 0: # QUAD
    plt.title("Quad Anal")
elif int(sys.argv[1]) == 1: # MCC
    plt.title("Quad Iter")
elif int(sys.argv[1]) == 2: # MCC
    plt.title("MCC")
plt.show()
