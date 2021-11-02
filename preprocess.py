import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from bridson import poisson_disc_samples

Lx = float(sys.argv[1])
Ly = float(sys.argv[2])
r = float(sys.argv[3]) # smaller r means more samples!
k = 50

folder = "samples/sample_Lx" + str(Lx) + "_Ly" + str(Ly) + "_r" + str(r) + "/"
os.makedirs(os.path.dirname(folder), exist_ok=True)
np.savetxt(folder + "Lx.txt", np.array([Lx]).astype(float), fmt='%f')
np.savetxt(folder + "Ly.txt", np.array([Ly]).astype(float), fmt='%f')

seed = 42
np.random.seed(seed=seed)
samples = poisson_disc_samples(width=Lx, height=Ly, r=r, k=k, random=np.random.random)

Np = len(samples)
print("Number of particles created: ", Np)

x = np.zeros(Np)
y = np.zeros(Np)
for i in range(0, len(samples)):
    x[i] = samples[i][0]
    y[i] = samples[i][1]

np.savetxt(folder + "pos.txt", samples )
np.savetxt(folder + "num.txt", np.array([Np]).astype(int), fmt='%i')

plt.figure()
plt.title("Number of particles: " + str(Np))
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x,y,'k.')
plt.axis('scaled')

plt.show()
