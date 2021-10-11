from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
import math

x, y, z = [], [], []
with open("C:/Users/jjg52/Desktop/crds/amber_enzymes/MTAN_EC/M06_corrs.txt", 'r') as f:
    lines = f.readlines()

for line in lines:
    split = line.split()
    x.append(float(split[0]))
    y.append(float(split[1]))
    z.append(float(split[2]))

xn = np.array(x)
yn = np.array(y)
points = np.array((xn, yn)).T
values = z
print(len(z))
c = interpolate.interp2d(x, y, z, kind='cubic')

x, y, z = [], [], []
with open("C:/Users/jjg52/Desktop/crds/amber_enzymes/MTAN_EC/full.txt", 'r') as f:
    lines = f.readlines()

for line in lines:
    split = line.split()
    x.append(float(split[0]))
    y.append(float(split[1]))
    z.append(float(split[2]))

xin = np.array(x)
yin = np.array(y)
xi = np.array((xin, yin)).T
print(len(points))
g = interpolate.griddata(points, values, xi, method='cubic')
# with open("C:/Users/jjg52/Desktop/crds/MTAN/corrections/corrected_RHF", 'w+') as f:
#     k = 0
#     for i in range(len(x)):
#         f.write("%20.10f%20.10f%20.10f\n" % (x[i], y[i], (z[i] + c(x[i], y[i]))))
#         k += 1
#     f.write("\n")

k = 0
for i in xi:
    if not math.isnan(g[k]):
        print("%20.10f%20.10f%20.10f" % (i[0], i[1], z[k] + g[k]))
    k += 1
#print(g)