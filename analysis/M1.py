import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.cm as cm

prefix = "../"
prefix = prefix + "test_pbed_r1_"
ppy =  131
R = 10
num = np.arange(0,1+1)

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


layer = np.array([ppy-8*R,ppy-6*R,ppy-4*R,ppy-2*R,ppy])

kkk = [np.size(np.where(py<l)) for l in layer]
colors = cm.rainbow(np.linspace(0,1,len(kkk)))

print('initial num ',kkk)

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
	kkk1 = [np.size(np.where(py<l)) for l in layer]
	ax = plt.subplot(111)
	for ii,c in zip(range(len(layer)),colors):
		M[i,ii] = float(kkk1[ii])/float(kkk[ii])
		ax.plot(np.sqrt(nn),M[:,ii],'*',color=c,label=layer[ii])
	plt.legend()

	name = 'r' + str(num[i])+'.png'
	plt.savefig(name,dpi=500)
	plt.clf()
	print("generate ",name)
print M
f.close()
