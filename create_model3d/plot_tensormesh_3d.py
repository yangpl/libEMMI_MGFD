from discretize import TensorMesh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits import mplot3d

from matplotlib import cm

nx = 100
ny = 100
nz = 100
dx = 200
dy = 200
dz = 40
x0 = -10000
y0 = -10000
z0 = -4000

hx = dx*np.ones(nx)
hy = dy*np.ones(ny)
hz = dz*np.ones(nz)
x1 = x0*np.ones(nx+1)
x2 = y0*np.ones(ny+1)
x3 = z0*np.ones(nz+1)
x1[1::] += np.cumsum(hx)
x2[1::] += np.cumsum(hy)
x3[1::] += np.cumsum(hz)


mesh = TensorMesh([hx, hy, hz], origin=(x1[0], x2[0], x3[0]))
print(mesh)

rho = np.fromfile('frho', dtype=np.float32)
rho = rho.reshape([nx, ny, nz], order='F')
rho = rho[:,:,::-1] #reverse the order of z-axis
mesh.plot_3d_slicer(rho,  xslice=2500, yslice=0, zslice=-1500,
                xlim=(x1[0], x1[nx-1]), ylim=(x2[0], x2[ny-1]), zlim=(x3[0], x3[nz-1]),
                pcolor_opts={'norm': LogNorm(vmin=0.1,vmax=100), 'cmap':'jet'})
mesh.plot_3d_slicer(rho,  xslice=-3500, yslice=1500, zslice=-2400,
                xlim=(x1[0], x1[nx-1]), ylim=(x2[0], x2[ny-1]), zlim=(x3[0], x3[nz-1]),
                pcolor_opts={'norm': LogNorm(vmin=0.1,vmax=100), 'cmap':'jet'})

                
rho = np.fromfile('frho_init', dtype=np.float32)
rho = rho.reshape([nx, ny, nz], order='F')
rho = rho[:,:,::-1] #reverse the order of z-axis
mesh.plot_3d_slicer(rho,  xslice=2500, yslice=0, zslice=-1500,
                xlim=(x1[0], x1[nx-1]), ylim=(x2[0], x2[ny-1]), zlim=(x3[0], x3[nz-1]),
                pcolor_opts={'norm': LogNorm(vmin=0.1,vmax=100), 'cmap':'jet'})


lenx = nx*dx
leny = ny*dy
x, y = np.meshgrid(x1, x2)
z = 900. + 100.*np.sin(2*np.pi*x/(3.*lenx)+np.pi/3)*np.sin(3.*np.pi*y/(2.*leny)-np.pi/3)

ax = plt.axes(projection ='3d')
ax.plot_surface(x,y,-z, cmap=cm.coolwarm)
ax.set_zlim(-1100,-700)
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
plt.show()
