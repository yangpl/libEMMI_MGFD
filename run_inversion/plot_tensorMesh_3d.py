from discretize import TensorMesh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

nx = 100
ny = 100
nz = 100

hx = 200*np.ones(nx)
hy = 200*np.ones(ny)
hz = 40*np.ones(nz)

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

rho = np.fromfile('param_final_Rv', dtype=np.float32)
rho = rho.reshape([nx, ny, nz], order='F')
rho = rho[:,:,::-1] #reverse the order of z-axis
mesh.plot_3d_slicer(rho,  xslice=2000, yslice=0, zslice=-1500,
                xlim=(x1[0], x1[nx-1]), ylim=(x2[0], x2[ny-1]), zlim=(x3[0], x3[nz-1]),
                grid=[4,4,3], pcolor_opts={'norm': LogNorm(vmin=1e-1,vmax=1e2),'cmap':'jet'})

# Rv = np.fromfile('Ex.bin', dtype=np.complex64)
# Rv = Rv.reshape([nx, ny, nz], order='F')
# Rv = Rv[:,:,::-1] #reverse the order of z-axis
# mesh.plot_3d_slicer(np.abs(Rv),  xslice=0, yslice=0, zslice=-1100,
#                 xlim=(x1[0], x1[nx-1]), ylim=(x2[0], x2[ny-1]), zlim=(x3[0], x3[nz-1]),
#                 pcolor_opts={'norm': LogNorm(vmin=1e-17,vmax=1e-8),'cmap':'rainbow'})


