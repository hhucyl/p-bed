import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys

prefix = "/home/pzhang/chen/p-bed/"
prefix = prefix + ""
prefix = prefix + "test_pbed_r_"
num = np.arange(0,999+1)
Py = 5
Px = 10
R = 10
nn = []
Vx = []
Vbx = []
Vpx = []
for i in range(len(num)):
	name = prefix + str(num[i]).zfill(4) + ".h5"
	f = h5.File(name)
	Nx = np.array(f['Nx'])
	Ny = np.array(f['Ny'])
	vel = np.array(f['Velocity_0'])
	vx = vel[0:-2:3]
	vx = vx.reshape((Ny[0],Nx[0]))
	vvx = np.abs(np.average(vx,axis=1))
	
	plt.subplot(2,2,1)
	plt.plot(vvx,np.arange(Ny[0]))
	plt.xlabel(r'$u_x$')
	plt.ylabel(r'depth')
	# plt.hold(True)

	plt.subplot(2,2,2)
	nn.append(num[i])
	vx = np.array(vx)
	idx = int(Nx[0]/2)
	vvvx = vx[:,idx]
	Vx.append(np.abs(vvvx.sum()))
	plt.plot(nn,Vx,'*')
	plt.ylabel(r'$\sum u_x$')
	plt.xlabel(r'num')
	# plt.hold(False)

	
	ax = plt.subplot(2,2,3)
	pos = np.array(f['Pposition'])
	py = np.max(pos[1:-1:3])+10
	plt.plot(vvx[:int(py)],np.arange(np.size(vvx[:int(py)])))
	# ax.set_ylim(np.min(vvx[:int(py)]),np.max(vvx[:int(py)]))
	plt.xlabel(r'$u_x^{bed}$')
	plt.ylabel(r'depth')
	# plt.hold(True)
	
	ax = plt.subplot(2,2,4)
	Vbx.append(vvx[:int(py)].sum())
	plt.plot(nn,Vbx,'*')
	plt.ylabel(r'$\sum u_x^{bed}$')
	plt.xlabel(r'num')
	# plt.hold(False)
	
	plt.tight_layout()
	name = 'r' + str(num[i])+'.png'
	print(name)
	plt.savefig(name)
	plt.clf()



# plt.show()
