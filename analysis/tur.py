import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

prefix = "/media/user/9EAEE48CAEE45DF1/cyl_temp/p-bed-data/1e4/"
prefix = "../../post/"
prefix = prefix + "turbulence_"
ppy =  131
R = 8
num = np.arange(0,999+1)
nu = 5e-4

name = prefix + str(num[0]).zfill(4) + ".h5"
f = h5.File(name)
Nx = int(f['Nx'][0])
Ny = int(f['Ny'][0])
print("Nx ", Nx, "Ny ", Ny)
ga = (np.array(f['Gamma'])).reshape((Ny,Nx))
ua = (np.array(f['Vave'])[0:-2:3]).reshape((Ny,Nx))
va = (np.array(f['Vave'])[1:-1:3]).reshape((Ny,Nx))
uaa = np.average(ua,axis=1)
vaa = np.average(va,axis=1)
plt.plot(uaa)
plt.savefig('uaa.png',dpi=500)
plt.clf()

us = np.zeros((Ny,Nx))+np.nan
vs = np.zeros((Ny,Nx))+np.nan
for i in range(len(uaa)):
	us[i,:] = ua[i,:]-uaa[i]
	vs[i,:] = va[i,:]-vaa[i]
uvs = us*vs
uvsy = np.average(uvs,axis=1)
plt.plot(uvsy)
plt.savefig('uvsy.png',dpi=500)
plt.clf()

uh = np.zeros((Ny,Nx,len(num)))+np.nan
vh = np.zeros((Ny,Nx,len(num)))+np.nan
uv = np.zeros((Ny,Nx,len(num)))+np.nan

for i in range(len(num)):
	name = prefix + str(num[i]).zfill(4) + ".h5"
	print("start process ", name)
	f = h5.File(name)
	
	vhx = (np.array(f['Vhas'])[0:-2:3]).reshape((Ny,Nx))
	vhy = (np.array(f['Vhas'])[1:-1:3]).reshape((Ny,Nx))
	ga[np.where(ga>1.0)] = 1.0
	uh[:,:,i] = vhx*(1-ga)
	vh[:,:,i] = vhy*(1-ga)
	uv[:,:,i] = vhx*vhy*(1-ga)
	# name = 't' + str(num[i])+'.png'
	# print(np.shape(ga))
	# plt.savefig(name,dpi=500)

uva = np.average(uv,axis=2)
plt.pcolor(uva)
plt.savefig('uva.png',dpi=500)
plt.clf()
uvhasy = np.average(uva,axis=1)
plt.plot(uvhasy)
plt.savefig('uvhasy.png',dpi=500)
plt.clf()

uu = np.sqrt(-min(uvhasy))
print("uu based on reynolds stress %.4e"%uu)
fi = (float(ppy)*float(Nx)-60.0*np.pi*R**2)/(float(ppy)*float(Nx))
k = 5.6e-3*fi**3/(1-fi)**2*(2*R)**2
Rek = uu*np.sqrt(k)/nu
print('porosity ', fi, 'k ',k)
print("Rek ",Rek)

du = np.zeros(Ny)
du[0:Ny-2] = uaa[1:Ny-1] - uaa[0:Ny-2]
du[Ny-1] = 0
tau = -uvhasy-uvsy+nu*du;
plt.plot(tau)
plt.savefig('tau.png',dpi=500)
plt.clf()
uu1 = np.sqrt(np.max(tau/1.0))
Rek1 = uu1*np.sqrt(k)/nu
print("uu1 based on tau %.4e"%uu1)
print("Rek1 ",Rek1)

f.close()

with open("Rek_1e4.txt","w") as f:
	f.write(str(uu)+"\t uu based on reynolds stress\n")
	f.write(str(uu1)+"\t uu1 based on tau\n")
	f.write(str(fi)+"\t porosity\n")
	f.write(str(k)+"\t k\n")
	f.write(str(Rek)+ "\t Rek\n")
	f.write(str(Rek1)+ "\t Rek1\n")


