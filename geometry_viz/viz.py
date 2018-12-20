from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
from matplotlib.mlab import PCA

np.random.seed(19680801)

def randrange(n,vmin, vmax):
	return (vmax- vmin)*np.random.rand(n) + vmin


fig = plt.figure()
ax = fig.add_subplot(111,projection = '3d')

n = 100

for c, m, zlow, zhigh in [('r','o', -50, -25),('b', '^', -30,-5)]:
	xs = randrange(n, 23, 32)
	ys = randrange(n, 0, 100)
	zs = randrange(n, zlow, zhigh)
	# ax.scatter(xs, ys, zs, c=c, marker=m)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

# plt.show()

def circle(r):
	x = np.empty((0,))
	y = np.empty((0,))

	for th in np.linspace(0,3.14*2, 100):
		x = np.append(x, r*np.cos(th) + r/20.0*np.random.rand()) 
		y = np.append(y, r*np.sin(th) + r/20.0*np.random.rand())  

	return x,y


def cylinder(r, l):
	x = np.empty((0,))
	y = np.empty((0,))
	z = np.empty((0,))

	for l in np.linspace(0, l, 30):
		for th in np.linspace(0,3.14*2, 50):
			x = np.append(x, r*np.cos(th) + r/10.0*np.random.rand()) 
			y = np.append(y, r*np.sin(th) + r/10.0*np.random.rand()) 
			z = np.append(z, l+np.random.rand())

	return np.array([x,y,z])


def rotx(th):
	return np.array([[1, 0,0],[0, np.cos(th), np.sin(th)],[0, -np.sin(th), np.cos(th)]])


def dist3d(cyl, dmin, dmax):

	new_cyl = np.empty((3,1))
	for i in range(cyl.shape[1]):
		norm = la.norm(cyl[:,i])
		if(norm<dmax and norm>dmin):
			new_cyl = np.hstack([new_cyl, cyl[None,:,i].T])
	return new_cyl



cyl1 = cylinder(10, 40)
cyl2 = rotx(np.deg2rad(60)) @ cylinder(10,40)+ np.array([[0,0, 40]]).T

cyl = np.hstack([cyl1,cyl2])


step = 15


colors = ['g', 'r', 'b', 'm','c', 'y' ]

for i in range(5):
	cyl_part = dist3d(cyl,i*step, (i+1)*step)
	ax.scatter(cyl_part[0,:],cyl_part[1,:], cyl_part[2,:], marker = '^', c=colors[i])
	mean = np.mean(cyl_part, axis = 1, keepdims = 1)
	ax.scatter(mean[0], mean[1], mean[2], color= 'k', linewidth = 10)

	cyl_n = cyl_part - mean
	w, v = la.eig(cyl_n @ cyl_n.T)
	direction = v[:, np.argmin(w)]
	ax.quiver(mean[0], mean[1], mean[2], 10*direction[0], 10*direction[1], 10*direction[2], linewidth= 7)

plt.show()