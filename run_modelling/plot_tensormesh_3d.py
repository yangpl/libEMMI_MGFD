from discretize import TensorMesh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


nx = 200
ny = 200
nz = 100

hx = nx*[100,]
hy = ny*[100,]
hz = nz*[40,]

x0 = -10000
y0 = -10000
z0 = -4000
x1 = x0*np.ones(nx+1)
x2 = y0*np.ones(ny+1)
x3 = z0*np.ones(nz+1)
x1[1::] += np.cumsum(hx)
x2[1::] += np.cumsum(hy)
x3[1::] += np.cumsum(hz)


mesh = TensorMesh([hx, hy, hz], origin=(x1[0], x2[0], x3[0]))
print(mesh)

# rho = np.fromfile('frho11', dtype=np.float32)
# rho = rho.reshape([nx, ny, nz], order='F')
# rho = rho[:,:,::-1] #reverse the order of z-axis
# mesh.plot_3d_slicer(rho,  xslice=0, yslice=0, zslice=-1500,
#                 xlim=(x1[0], x1[nx-1]), ylim=(x2[0], x2[ny-1]), zlim=(x3[0], x3[nz-1]),
#                 grid=[4,4,3], pcolor_opts={'norm': LogNorm(vmin=0.3,vmax=100),'cmap':'rainbow'})


Rv = np.fromfile('Hz.bin', dtype=np.complex64)
Rv = Rv.reshape([nx, ny, nz], order='F')
Rv = Rv[:,:,::-1] #reverse the order of z-axis
mesh.plot_3d_slicer(np.abs(Rv),  xslice=0, yslice=0, zslice=-1550,
                xlim=(x1[0], x1[nx-1]), ylim=(x2[0], x2[ny-1]), zlim=(x3[0], x3[nz-1]),
                grid=[4,4,3], pcolor_opts={'norm': LogNorm(vmin=1e-6,vmax=1e-16),'cmap':'rainbow'})


