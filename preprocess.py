import numpy as np
import sys
import matplotlib.pyplot as plt
from bridson import poisson_disc_samples

w = float(sys.argv[1])
h = float(sys.argv[2])
r = float(sys.argv[3]) # smaller r means more samples!
k = 50

seed = 42
np.random.seed(seed=seed)
samples = poisson_disc_samples(width=w, height=h, r=r, k=k, random=np.random.random)

Np = len(samples)
print("Number of particles created: ", Np)

x = np.zeros(Np)
y = np.zeros(Np)
for i in range(0, len(samples)):
    x[i] = samples[i][0]
    y[i] = samples[i][1]

filename = "samples/sample_w" + str(w) + "_h" + str(h) + "_r" + str(r)
np.savetxt(filename + "_pos.txt", samples )
np.savetxt(filename + "_num.txt", np.array([Np]).astype(int), fmt='%i')

plt.figure()
plt.title("Number of particles: " + str(Np))
plt.xlabel("x")
plt.ylabel("y")
plt.plot(x,y,'k.')
plt.axis('scaled')

plt.show()
