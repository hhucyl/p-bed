import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

prefix = "/media/user/9EAEE48CAEE45DF1/cyl_temp/p-bed-data/1e4/"
prefix = "../../post/"
prefix = prefix + "test_pbed_i1_"
ppy =  131
R = 8
nu = 5e-4

name = prefix + str(999).zfill(4) + ".h5"
f = h5.File(name)
Nx = int(f['Nx'][0])
Ny = int(f['Ny'][0])
print("Nx ", Nx, "Ny ", Ny)
ga = (np.array(f['Gamma'])).reshape((Ny,Nx))
plt.pcolor(ga)
plt.savefig('Pos.png',dpi=500)
plt.clf()
f.close()
