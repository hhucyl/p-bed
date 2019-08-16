import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys

prefix = "/home/pzhang/chen/p-bed/"
prefix = prefix + ""
prefix = prefix + "test_pbed_r_"
num = np.arange(0,1+1)
Py = 5
Px = 10
R = 10
nn = []
Vx = []
Vbx = []
Vpx = []

name = prefix + str(num[0]).zfill(4) + ".h5"
f = h5.File(name)
PX = np.array(f['RWPposition'])
px = PX[0:-2:3]
py = PX[1:-1:3]
pos = np.array(f['Pposition'])
ppy = int(np.max(pos[1:-1:3])+10)
print(ppy)
layer = np.array([22,43,65,86,ppy])

kkk = []
kkk1 = []
kkk.append(np.where(py<layer[0]))
kkk1.append(np.where(py<layer[0]))
for i in np.arange(1,np.size(layer)):
	# print(layer[i-1],' ',layer[i])
	k1 = np.where(py>layer[i-1])
	k2 = np.where(py<layer[i])
	kkk.append(np.intersect1d(k1,k2))
	kkk1.append(k2)
# print(kkk)
colors = cm.rainbow(np.linspace(0,1,len(kkk)))

V = np.zeros(len(num))+np.nan
Vbx = np.zeros(len(num))+np.nan
nn = np.zeros(len(num))+np.nan
M = np.zeros((len(num),len(layer)))+np.nan
for i in range(len(num)):
	nn[i] = i
	name = prefix + str(num[i]).zfill(4) + ".h5"
	f = h5.File(name)
	Nx = np.array(f['Nx'])
	Ny = np.array(f['Ny'])
	Ga = np.array(f['Gamma'])
	ga = Ga.reshape((Ny[0],Nx[0]))

	plt.figure(figsize=(10,10))

	ax = plt.subplot(3,1,1)
	ax.pcolor(ga)
	PX = np.array(f['RWPposition'])
	px = PX[0:-2:3]
	py = PX[1:-1:3]
	for k,c in zip(kkk,colors):
		ax.scatter(px[k],py[k],s=0.5e-1,color=c)
	# plt.hold(False)
	ax.set_xlim(0,Nx[0])
	ax.set_ylim(0,Ny[0])
	plt.xlabel('X')
	plt.ylabel('Y')			
	
	ax = plt.subplot(3,4,5)
	vel = np.array(f['Velocity_0'])
	vx = np.reshape(vel[0:-2:3],(Ny[0],Nx[0]))
	vy = np.reshape(vel[1:-1:3],(Ny[0],Nx[0]))
	vvx = np.average(vx,axis=1)
	ax.plot(vvx[:-1],np.arange(Ny[0]-1))
	plt.xlabel(r'$u_x$')
	plt.ylabel(r'depth')

	ax = plt.subplot(3,4,6)
	half = int(Nx[0]/2)
	V[i] = np.sum(vx[:,half])
	ax.plot(nn,V,'*')
	plt.ylabel(r'$\sum u_x$')
	plt.xlabel(r't')

	ax = plt.subplot(3,4,7)
	pos = np.array(f['Pposition'])
	ppy = np.max(pos[1:-1:3])+10
	plt.plot(vvx[:int(ppy)],np.arange(np.size(vvx[:int(ppy)])))
	plt.xlabel(r'$u_x^{bed}$')
	plt.ylabel(r'depth')
	# plt.hold(True)
	
	ax = plt.subplot(3,4,8)
	Vbx[i] = vvx[:int(ppy)].sum()
	plt.plot(nn,Vbx,'*')
	plt.ylabel(r'$\sum u_x^{bed}$')
	plt.xlabel(r't')
	
	ax = plt.subplot(3,1,3)
	for ii,c in zip(range(len(layer)),colors):
		M[i,ii] = np.sum(np.where(py[kkk1[ii]]<layer[ii]))
		ax.plot(np.sqrt(nn),M[:,ii],'*',color=c)
	plt.ylabel(r'$M$')
	plt.xlabel(r'sqrt(t)')

	plt.tight_layout()
	name = 'r' + str(num[i])+'.png'
	print(name)
	plt.savefig(name,dpi=500)
	plt.clf()
	plt.close('all')

f = h5.File('M.h5','w')
f.create_dataset('M',data=M)
f.close()


# plt.show()
