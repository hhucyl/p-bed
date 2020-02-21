import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# prefix = "/media/user/9EAEE48CAEE45DF1/cyl_temp/p-bed-data/1e4/"
prefix = "../"
prefix = prefix + "test_pbed_r1_"
ppy =  131
R = 10
num = np.arange(0,100+1)

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
gga = np.sum(ga,axis=1)
# print(gga)
print(ppy-1,gga[ppy-1])
print(ppy,gga[ppy])


layer = np.array([66,88,110,ppy])
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
	
	for ii,c in zip(range(len(layer)),colors):
		M[i,ii] = float(np.size(np.where(py[kkk[ii]]<=layer[ii])))/M0[ii]
		plt.plot(np.sqrt(nn),M[:,ii],'*',color=c,label=layer[ii])
	plt.legend()

	name = 'r' + str(num[i])+'.png'
	plt.savefig(name,dpi=500)
	plt.clf()
	print("generate ",name)
print M
f.close()
