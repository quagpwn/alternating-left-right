#import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
#from matplotlib import rc


#sys.stdout = open('result.data', 'w')
#plt.style.use('mystyle')




def rotX(rot):
	r1 = np.array([1, 0, 0])
	r2 = np.array([0, np.cos(rot), -np.sin(rot)])
	r3 = np.array([0, np.sin(rot), np.cos(rot)])
	rotmat = np.array([r1, r2, r3])
	return rotmat
def rotY(rot):
	r1 = np.array([np.cos(rot), 0, np.sin(rot)])
	r2 = np.array([0, 1, 0])
	r3 = np.array([-np.sin(rot), 0, np.cos(rot)])
	rotmat = np.array([r1, r2, r3])
	return rotmat

def rotZ(rot):
	r1 = np.array([np.cos(rot), -np.sin(rot), 0])
	r2 = np.array([np.sin(rot), np.cos(rot), 0])
	r3 = np.array([0, 0, 1])
	rotmat = np.array([r1, r2, r3])
	return rotmat

def rotaxis(startmat,normalvec,angle):
	angle = angle*np.pi/180
	x = startmat[0]
	y = startmat[1]
	z = startmat[2]
	u = normalvec[0]
	v = normalvec[1]
	w = normalvec[2]

	k = -u * x - v * y - w * z
	X = -u*k * (1-np.cos(angle)) + x*np.cos(angle) + (v*z - w*y)*np.sin(angle)
	Y = -v*k * (1-np.cos(angle)) + y*np.cos(angle) + (w*x - u*z)*np.sin(angle)
	Z = -w*k * (1-np.cos(angle)) + z*np.cos(angle) + (u*y - v*x)*np.sin(angle)
	point = np.array([X,Y,Z])

	return point

def getDihedral(a,h,x, phi, psi):
#	rtd = 180/np.pi
	dtr = np.pi/180
	phi = phi*dtr
	psi = psi*dtr
	a = a*dtr
	h = h*dtr
	x = x*dtr
	A = np.cos(h)*np.cos(a) + np.sin(h)*np.sin(a)*np.cos(psi)
	B = np.cos(h)*np.sin(a) - np.sin(h)*np.cos(a)*np.cos(psi)

	theta = np.degrees(np.arccos(A*np.cos(x) + B*np.sin(x)*np.cos(phi) - np.sin(x)*np.sin(h)*np.sin(phi)*np.sin(psi)))

	L1_num = -B*np.sin(phi) - np.sin(h)*np.cos(phi)*np.sin(psi)
	L1_denom = A*np.sin(x) - B*np.cos(x)*np.cos(phi) + np.cos(x)*np.sin(h)*np.sin(phi)*np.sin(psi)

	L1 = np.degrees(np.arctan2(L1_num,L1_denom))
	if L1 >= (phi*180/np.pi):
		L1 = L1 - 180
	elif L1 <= (phi*180/np.pi):
	 	L1 = L1 + 180

	D = np.cos(x)*np.cos(a) + np.sin(x)*np.sin(a)*np.cos(phi)
	E = np.cos(x)*np.sin(a) - np.sin(x)*np.cos(a)*np.cos(phi)

	L2_num = -E*np.sin(psi) - np.sin(x)*np.cos(psi)*np.sin(phi)
	L2_denom = D*np.sin(h) - E*np.cos(h)*np.cos(psi) + np.cos(h)*np.sin(x)*np.sin(phi)*np.sin(psi)

	L2 = np.degrees(np.arctan2(L2_num,L2_denom))
	if L2 >= (psi*180/np.pi):
		L2 = L2 - 180
	elif L2 <= (psi*180/np.pi):
		L2 = L2 + 180

	dihedral = L1-L2 + 180
	if dihedral > 180:
		dihedral = dihedral - 360
	elif dihedral < -180:
		dihedral = dihedral + 360
	dihedral = abs(dihedral)

	d = 3.8*np.sin(theta*dtr)/np.cos(theta*dtr/2)

	#C3 = np.array([0,0,0])

	rot1 = (180-theta)/2 * dtr

	C32 = np.array([-3.8*np.cos(rot1),3.8*np.sin(rot1),0])
	C32 = C32/LA.norm(C32)

	C34 = np.array([3.8*np.cos(rot1),3.8*np.sin(rot1),0])
	C34 = C34/LA.norm(C34)

	C1 = np.array([-d,0,0])
	C5 = np.array([d,0,0])
	C31 = rotaxis(C1,C32,dihedral-180)
	C35 = rotaxis(C5,C34,180-dihedral)

	cost = np.dot(C31,C35)/(LA.norm(C31)*LA.norm(C35))
	cost = np.clip(cost, -1.0, 1.0)
	beta = np.degrees(np.arccos(cost))
	radius = 3.8*np.sin(theta*dtr)*np.sin(beta*dtr/2)/(np.sin(beta*dtr)*np.cos(theta*dtr/2))
	return dihedral,theta, L1, L2, beta, radius

# ------------------------------------------
X = []
Y = []
Z = []
Z2 = []
Z3 = []
Z4 = []
ticks = [-180, -120, -60, 0, 60, 120, 180]
increment = 5
for psi in range(-180, 181, increment):
	for phi in range(-180, 181, increment):
		dihed,theta,L1,L2,beta,radius = getDihedral(111,17.4,17.4,phi,psi)
		#alpha, eta, xi angles are N-CA-C, C2A-C1A-C, C0A-C1A-N

		X.append(phi)
		Y.append(psi)
		Z.append(dihed)
		Z2.append(theta)
		Z3.append(beta)
		Z4.append(radius)

df = pd.DataFrame({"phi":X, "psi":Y, "dihed":Z, "theta":Z2,"beta":Z3,"radius":Z4})

df_filtered_N12 = df[(abs(df.beta) >= 119) & (abs(df.beta) < 121)] #119 and 121
phi_N12 = df_filtered_N12.phi.values
psi_N12 = df_filtered_N12.psi.values

df_filtered_check = df[(df.phi == -75) & (df.psi == -25)]
print (df_filtered_check)

# ======================================

fig = plt.figure(figsize=(5,4.2))
ax1 = fig.add_subplot(2,2,2) #dihedral
ax2 = fig.add_subplot(2,2,1) #CA123 angle
ax3 = fig.add_subplot(2,2,3) #CA135 angle ~ betavature
ax4 = fig.add_subplot(2,2,4) #N linkage map

im1= ax1.scatter(X,Y,c=Z,s=5, cmap='jet',alpha=1)
ax1.set_xlabel(r'$\phi$', labelpad=0)
ax1.set_ylabel(r'$\psi$', labelpad=-8)
ax1.set_xticks(ticks)
ax1.set_yticks(ticks)
ax1.set_xlim(-180, 180)
ax1.set_ylim(-180, 180)
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(im1, cax=cax, orientation='vertical')
cbar.set_label(r'|Dihedral angle| ($|\gamma|$)', size = 'small', rotation = 270, labelpad = 10)

im2=ax2.scatter(X,Y,c=Z2,s=5, cmap='jet',alpha=1)
ax2.set_xlabel(r'$\phi$', labelpad=0)
ax2.set_ylabel(r'$\psi$', labelpad=-8)
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.set_xlim(-180, 180)
ax2.set_ylim(-180, 180)
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(im2, cax=cax, orientation='vertical')
cbar.set_label(r'CA1-CA2-CA3 angle ($\theta$)', size='small', rotation = 270, labelpad = 10)

im3 = ax3.scatter(X,Y,c=Z3,s=5, cmap='jet',alpha=1)
ax3.set_xlabel(r'$\phi$', labelpad=0)
ax3.set_ylabel(r'$\psi$', labelpad=-8)
ax3.set_xticks(ticks)
ax3.set_yticks(ticks)
ax3.set_xlim(-180, 180)
ax3.set_ylim(-180, 180)
divider = make_axes_locatable(ax3)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(im3, cax=cax, orientation='vertical')
cbar.set_label(r'CA1-CA3-CA5 Angle ($\beta$)', size = 'small', rotation = 270, labelpad = 10)

headers = ['Energy','Phi','Psi']
df_rama = pd.read_csv('LDalt10_edit.csv', names=headers)
df_rama = df_rama[(df_rama.Energy < 0)]
LD_phi = df_rama['Phi'].values
LD_psi = df_rama['Psi'].values
LD_energy = df_rama['Energy'].values
im4 = ax4.scatter(LD_phi,LD_psi,c=LD_energy,s=1, cmap='Greys_r',alpha=1)
ax4.scatter(phi_N12,psi_N12, c = 'blue', s=1, label = '12', alpha=1)
ax4.scatter(-150, 150, c='red', s=20)
ax4.scatter(-150, -100, c='red', s=20)
ax4.scatter(-75, -25, c='red', s=20)
ax4.scatter(-80, -140, c='red', s=20)
ax4.set_xlabel(r'$\phi$', labelpad=0)
ax4.set_ylabel(r'$\psi$', labelpad=-8)
ax4.set_xticks(ticks)
ax4.set_yticks(ticks)
ax4.set_xlim(-180, 180)
ax4.set_ylim(-180, 180)
divider = make_axes_locatable(ax4)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(im4, cax=cax, orientation='vertical')
cbar.set_label('Energy (kcal/mol)', size = 'small', rotation = 270, labelpad = 10)
ax4.grid()

plt.tight_layout()
plt.show()
