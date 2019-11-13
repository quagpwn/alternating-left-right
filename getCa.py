import sys
import numpy as np
from numpy import linalg as LA
import math
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from numpy.polynomial.polynomial import polyfit
from progress.bar import Bar

headers = ['Energy','Phi','Psi']
df_rama = pd.read_csv('LDalt10_edit.csv', names=headers)
df_rama = df_rama[(df_rama.Energy < 0)]

LD_phi = df_rama['Phi'].values
LD_psi = df_rama['Psi'].values
LD_energy = df_rama['Energy'].values

# print (LD_phi)
#sys.stdout = open('result.data', 'w')
plt.style.use('mystyle')

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
	rtd = 180/np.pi
	dtr = np.pi/180
	phi = phi*dtr
	psi = psi*dtr
	a = a*dtr
	h = h*dtr
	x = x*dtr
	A = np.cos(h)*np.cos(a) + np.sin(h)*np.sin(a)*np.cos(psi)
	B = np.cos(h)*np.sin(a) - np.sin(h)*np.cos(a)*np.cos(psi)

	ang = np.degrees(np.arccos(A*np.cos(x) + B*np.sin(x)*np.cos(phi) - np.sin(x)*np.sin(h)*np.sin(phi)*np.sin(psi)))


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


	dihedral = L1-L2 #L D hetero
	if dihedral <0:
		dihedral = dihedral + 360
	dihedral = abs(dihedral - 180)

	d = 3.8*np.sin(ang*dtr)/np.sin((180-ang)/2*dtr)

	C3 = np.array([0,0,0])

	rot1 = (180-ang)/2 * dtr

	C32 = np.array([-3.8*np.cos(rot1),3.8*np.sin(rot1),0])
	C32 = C32/LA.norm(C32)

	C34 = np.array([3.8*np.cos(rot1),3.8*np.sin(rot1),0])
	C34 = C34/LA.norm(C34)

	C1 = np.array([-d,0,0])
	C5 = np.array([d,0,0])
	C31 = rotaxis(C1,C32,180-dihedral)
	C35 = rotaxis(C5,C34,180+dihedral)

	cost = np.dot(C31,C35)/(LA.norm(C31)*LA.norm(C35))
	cost = np.clip(cost, -1.0, 1.0)
	cur = np.degrees(np.arccos(cost))

	getrad = (2*(1+np.cos(cur*dtr)))

	if getrad ==0:
		getrad = 0.000000000001

	radius = d/np.sqrt(getrad) #this is actually diameter
	radius2 = 3.8*np.sin(ang*dtr)*np.sin(cur*dtr/2)/(np.sin(cur*dtr)*np.cos(ang*dtr/2))
	print (radius, radius2)
	# if radius >=30:
	# 	radius = -1

	# print (phi*rtd, psi*rtd, dihedral, ang, L1, L2, cur, radius)


	return dihedral,ang, L1, L2, cur, radius
def getPolyang(theta, CAangle):

	t = theta*(np.pi/180)
	p = CAangle*(np.pi/180)
	#r = (np.sqrt(2)*3.8/2)/np.sin(t)
	r = 3.8*np.sqrt(1-np.cos(p))/np.sqrt(1-np.cos(2*t))

	z = np.sqrt(3.8**2 - ((r*np.cos(t) - r)**2 + (r*np.sin(t))**2))

	p1 = np.array([r, 0, 0])
	p2 = np.array([r*np.cos(t), r*np.sin(t), z])
	p3 = np.array([r*np.cos(2*t), r*np.sin(2*t), 0])
	p4 = np.array([r*np.cos(3*t), r*np.sin(3*t),z])

	v1 = p2 - p1
	v2 = p3 - p2
	v3 = p4 - p3

	n1 = np.cross(v1,v2)
	n2 = np.cross(v2,v3)

	o = (n1+n2) / LA.norm(n1+n2)
	curve = np.degrees(np.arccos(np.dot(n1,o) / (LA.norm(n1)*LA.norm(o))))

	angle = np.degrees(np.arccos(np.dot(n1,n2) / (LA.norm(n1)*LA.norm(n2))))

	return angle


phi, psi = -180, -180
X, Y, Z, Z2, Z3, Z4, Z5, Z6, xd, yd = [], [], [], [], [], [], [], [], [], []

increment = 5

ticks=[-180, -120, -60, 0, 60, 120, 180]
ticks2 = [0, 30, 60, 90, 120, 150, 180]

for resnum in range (6, 30, 2):
	polyangle = getPolyang(360/resnum, 90)
	xd.append(resnum)
	yd.append(polyangle)

# bar = Bar('Processing', max=100)
# for i in range(100):
# 	#do something
# 	bar.next()
# bar.finish()

for psi in range(-180, 181, increment):
	for phi in range(-180, 181, increment):
		dihed,ang,L1,L2,cur,radius = getDihedral(111,17.4,17.4,phi,psi)
		#alpha, eta, xi angles are N-CA-C, C2A-C1A-C, C0A-C1A-N

		X.append(phi)
		Y.append(psi)
		Z.append(dihed)
		Z2.append(ang)
		Z3.append(cur)
		Z4.append(L1)
		Z5.append(L2)
		Z6.append(radius)




df = pd.DataFrame({"phi":X, "psi":Y, "dihed":Z, "ang":Z2,
"cur":Z3, "L1":Z4, "L2":Z5, "Radius":Z6})

# print (df)

df_show = df[(df.phi == 150) & (df.psi == 100)]

# print (df_show.phi.values, df_show.psi.values,
# df_show.dihed.values,":dihedral", df_show.Radius.values,
# ":radius", df_show.cur.values,":curvature")
df_filtered_N6 = df[(abs(df.cur) >= 59) & (abs(df.cur) < 61)]
df_filtered_N8 = df[(abs(df.cur) >= 89) & (abs(df.cur) < 91)]
df_filtered_N10 = df[(abs(df.cur) >= 110) & (abs(df.cur) < 109)]
df_filtered_N12 = df[(abs(df.cur) >= 119) & (abs(df.cur) < 121)]
df_filtered_N14 = df[(abs(df.cur) >= 153) & (abs(df.cur) < 155)]
df_filtered_N16 = df[(abs(df.cur) >= 157) & (abs(df.cur) < 160)]
df_filtered_N18 = df[(abs(df.cur) >= 160) ]
# df_filtered_N20 = df[(abs(df.dihed) >= 154) & (abs(df.dihed) < 156)]



phi_N6 = df_filtered_N6.phi.values
psi_N6 = df_filtered_N6.psi.values

phi_N8 = df_filtered_N8.phi.values
psi_N8 = df_filtered_N8.psi.values

phi_N10 = df_filtered_N10.phi.values
psi_N10 = df_filtered_N10.psi.values

phi_N12 = df_filtered_N12.phi.values
psi_N12 = df_filtered_N12.psi.values

print (phi_N12, psi_N12)


phi_N14 = df_filtered_N14.phi.values
psi_N14 = df_filtered_N14.psi.values

phi_N16 = df_filtered_N16.phi.values
psi_N16 = df_filtered_N16.psi.values

phi_N18 = df_filtered_N18.phi.values
psi_N18 = df_filtered_N18.psi.values
#
phifit = np.array(phi_N12)
psifit = np.array(psi_N12)

fitx = []
fity = []
for i in phifit:
	for j in psifit:
		if i+5 < j and i+100 > j:
			fitx.append(i)
			fity.append(j)

# slope = psi_N12/phi_N12
# slope[phi_N12==0] = 0
# print (slope)
# phi_N20 = df_filtered_N20.phi.values
# psi_N20 = df_filtered_N20.psi.values

df_filtered_rad = df[(df.Radius <= 15)&(df.Radius >= 8)&(df.ang > 0) &(df.ang < 180)]
rad_filtered = df_filtered_rad.Radius.values
rad_phi = df_filtered_rad.phi.values
rad_psi = df_filtered_rad.psi.values





fig = plt.figure(figsize=(11,7.5))
ax1 = fig.add_subplot(3,4,1) #dihedral
ax2 = fig.add_subplot(3,4,9) #lambda 1
ax3 = fig.add_subplot(3,4,3) #Diameter
ax4 = fig.add_subplot(3,4,4) #angle
ax5 = fig.add_subplot(3,4,8) #angle hist
ax6 = fig.add_subplot(3,4,2) #Curvature
ax7 = fig.add_subplot(3,4,6) #Curv Hist
ax8 = fig.add_subplot(3,4,5) #res num in polygon
ax9 = fig.add_subplot(3,4,7) #Radius hist
ax10 = fig.add_subplot(3,4,11) #N linkage map
ax11 = fig.add_subplot(3,4,10) # lambda2
ax12 = fig.add_subplot(3,4,12) #LD ramachandran
# ax12 = fig.add_subplot(3,4,12) #LD ramachandran

im1= ax1.scatter(X,Y,c=Z,s=5, cmap='jet',alpha=1)
ax1.set_xlabel('phi')
ax1.set_ylabel('psi')
ax1.set_xticks(ticks)
ax1.set_yticks(ticks)
ax1.set_xlim(-180, 180)
ax1.set_ylim(-180, 180)
ax1.set_title('CA1~CA4 Dihedral')

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation='vertical')



ax2.scatter(X,Z4, s=1, alpha=1)
ax2.set_xlabel('phi')
ax2.set_ylabel('Lambda 1')
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.set_xlim(-180, 180)
ax2.set_ylim(-180, 180)
ax2.set_title('Lambda 1')

ax11.scatter(Y,Z5, s=1, cmap='Reds',alpha=1)
ax11.set_xlabel('psi')
ax11.set_ylabel('Lambda 2')
ax11.set_xticks(ticks)
ax11.set_yticks(ticks)
ax11.set_xlim(-180, 180)
ax11.set_ylim(-180, 180)
ax11.set_title('Lambda 2')

im4=ax4.scatter(X,Y,c=Z2,s=5, cmap='jet',alpha=1)
ax4.set_xlabel('phi')
ax4.set_ylabel('psi')
ax4.set_xticks(ticks)
ax4.set_yticks(ticks)
ax4.set_xlim(-180, 180)
ax4.set_ylim(-180, 180)
ax4.set_title('CA1-CA2-CA3 Angle')
divider = make_axes_locatable(ax4)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im4, cax=cax, orientation='vertical')

ax5.hist(Z2)
ax5.set_title('Angle Distribution')
ax5.set_xlabel('Angle')
ax5.set_ylabel('Count')
ax5.set_xticks(ticks2)


im6 = ax6.scatter(X,Y,c=Z3,s=5, cmap='jet',alpha=1)
ax6.set_xlabel('phi')
ax6.set_ylabel('psi')
ax6.set_xticks(ticks)
ax6.set_yticks(ticks)
ax6.set_xlim(-180, 180)
ax6.set_ylim(-180, 180)
ax6.set_title('CA1-CA3-CA5 Angle')
divider = make_axes_locatable(ax6)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im6, cax=cax, orientation='vertical')

ax7.hist(Z3)
ax7.set_xticks(ticks2)
ax7.set_xlabel('Angle')
ax7.set_title('CA1-CA3-CA5 Angle Distribution')
ax7.set_ylabel('Count')

ax8.scatter(xd, yd, s=5)
ax8.set_ylabel('Angle(Deg)')
ax8.set_xlabel('Number of residues in cyclic peptide')
ax8.set_xlim(5, 25)
ax8.set_ylim(80, 180)
ax8.xaxis.set_ticks([6, 8, 10, 12, 14, 16, 18, 20, 22, 24])
ax8.yaxis.set_ticks([90, 120, 150, 180])
ax8.grid()
ax8.set_title('CA1-CA4 dihedral')

im3=ax3.scatter(rad_phi,rad_psi,c=rad_filtered,s=5, cmap='jet',alpha=1)
ax3.set_xlabel('phi')
ax3.set_ylabel('psi')
# ax3.set_xticks([-180,-135,-90,-45,0])
# ax3.set_yticks([-180,-135,-90,-45,0])
ax3.set_xticks(ticks)
ax3.set_yticks(ticks)

ax3.set_xlim(-180, 180)
ax3.set_ylim(-180, 180)
ax3.set_title('Diameter')
divider = make_axes_locatable(ax3)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im3, cax=cax, orientation='vertical')

ax9.hist(rad_filtered)
ax9.set_xlabel('Diameter')
ax9.set_title('Radius distribution')
ax9.set_ylabel('Count')

ax10.scatter(LD_phi,LD_psi,c=LD_energy,s=1, cmap='Greys_r',alpha=1)
ax10.scatter(phi_N6,psi_N6, c = 'red', s=1, label = '6', alpha=1)
ax10.scatter(phi_N8,psi_N8, c = 'orange', s=1, label = '8', alpha=1)
ax10.scatter(phi_N10,psi_N10, c = 'green', s=1, label = '10', alpha=1)
ax10.scatter(phi_N12,psi_N12, c = 'blue', s=1, label = '12', alpha=1)
ax10.scatter(phi_N14,psi_N14, c = 'purple', s=1, label = '14', alpha=0.5)

#b, m = polyfit(fitx, fity, 1)
#ax10.scatter(fitx, fity)

# ax10.scatter(phi_N16,psi_N16, c = 'black', s=5, label = '16', alpha=0.5)
# ax10.scatter(phi_N18,psi_N18, c = 'grey', s=5, label = '18', alpha=0.5)
# ax10.scatter(phi_N20,psi_N20, c = 'blue', s=5, label = '20')
box = ax10.get_position()
ax10.set_position([box.x0, box.y0, box.width * 1, box.height])
ax10.legend(loc='center', bbox_to_anchor=(0.5, -0.5), ncol=7)

ax10.set_xlabel('phi')
ax10.set_ylabel('psi')
ax10.set_xticks(ticks)
ax10.set_yticks(ticks)
ax10.set_xlim(-180, 180)
ax10.set_ylim(-180, 180)
ax10.set_title('N residue cycle')

im12= ax12.scatter(LD_phi,LD_psi,c=LD_energy,s=3, cmap='Greys_r',alpha=1)
ax12.set_xlabel('phi')
ax12.set_ylabel('psi')
ax12.set_xticks(ticks)
ax12.set_yticks(ticks)
ax12.set_xlim(-180, 180)
ax12.set_ylim(-180, 180)
ax12.set_title('LD ramachandran')

divider = make_axes_locatable(ax12)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im12, cax=cax, orientation='vertical')



plt.tight_layout()
plt.show()
