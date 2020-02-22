import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

prefix = "/media/user/9EAEE48CAEE45DF1/cyl_temp/p-bed-data/1e4/"
prefix = "../../post/"
prefix = prefix + "test_pbed_r1_"
ppy =  131
R = 8
num = np.arange(0,200+1)
Dm = 1.15e-6; 

name = prefix + str(num[0]).zfill(4) + ".h5"
print(name)
f = h5.File(name)
PX = np.array(f['RWPposition'])
px = PX[0:-2:3]
py = PX[1:-1:3]
Ga = np.array(f['Gamma'])
Nx = int(f['Nx'][0])
Ny = int(f['Ny'][0])
ga = Ga.reshape((Ny,Nx))
ga[np.where(ga>1.0)] = 1.0
gga = np.sum(ga,axis=1)
# print(gga)
print(ppy-1,gga[ppy-1])
print(ppy,gga[ppy])


layer = np.array([ppy-6*(R+2),ppy-4*(R+2),ppy-2*(R+2),ppy])
layer = np.array([66,88,110,131]);
M0 = np.zeros((len(layer))) + np.nan
kkk = []
for i in range(len(layer)):
	kkk.append(np.where(py<=layer[i]))
	M0[i] = np.size(np.where(py<=layer[i]))
colors = cm.rainbow(np.linspace(0,1,len(kkk)))

print('initial num ',M0)

nn = np.zeros(len(num)) + np.nan
M = np.zeros((len(num),len(layer))) + np.nan
for i in range(len(num)):
	nn[i] = i
	name = prefix + str(num[i]).zfill(4) + ".h5"
	print("start process ", name)
	f = h5.File(name)
	PX = np.array(f['RWPposition'])
	px = PX[0:-2:3]
	py = PX[1:-1:3]
	
	for ii in range(len(layer)):
		M[i,ii] = float(np.size(np.where(py[kkk[ii]]<=layer[ii])))/M0[ii]

k1 = np.zeros(len(layer)) + np.nan
Vs = np.zeros(len(layer)) + np.nan
Vt = np.zeros(len(layer)) + np.nan
Vv = np.zeros(len(layer)) + np.nan
As = np.zeros(len(layer)) + np.nan
 
for ii,c in zip(range(len(layer)),colors):
	plt.plot(np.sqrt(nn),M[:,ii],'o',color=c,label=layer[ii])
	Z = np.polyfit(np.sqrt(nn),M[:,ii],1)
	k1[ii] = -Z[0]
	plt.plot(np.sqrt(nn),Z[0]*np.sqrt(nn)+Z[1],color=c)
	Vs[ii] = np.sum(np.sum(ga[1:layer[ii],:]))
	Vt[ii] = (layer[ii]*1.0-1.0)*Nx*1.0
	Vv[ii] = Vt[ii] - Vs[ii]
	As[ii] = Nx - np.sum(ga[layer[ii],:]) 


plt.legend()
plt.savefig('M.png',dpi=500)
plt.clf()

Deff = np.pi*(0.5*Vv/As*k1)**2/Dm
print('Deff/Dm ', Deff)
plt.plot(Deff,layer,'*')
plt.ylabel('Depth')
plt.xlabel('Deff/Dm')

plt.savefig('Deff',dpi=500)
plt.clf()

f.close()
