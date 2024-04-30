import emg3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

freq = 0.25
x0 = -10000
y0 = -10000
z0 = -4000
hx= 100*np.ones([200,])
hy= 100*np.ones([200,])
hz= 40*np.ones([101,])

grid = emg3d.TensorMesh([hx, hy, hz], origin=(x0, y0, z0))
print(grid)

res_x = np.loadtxt('fort.7', skiprows=0, unpack=True)
model = emg3d.Model(grid, res_x)
grid.plot_3d_slicer(
    model.property_x, zslice=-1900,
    xlim=(-10000, 10000), ylim=(-10000, 10000), zlim=(-4000, 40),
    pcolor_opts={'norm': LogNorm(vmin=0.3, vmax=200)}
)

grid_ext = emg3d.construct_mesh(
    frequency=freq,  # Frequenzy
    # props: 1st value irrelevant, because of `vector`.
    #        2nd value to calculate the boundary in xy
    #        3rd value to calculate the boundary to -z
    #        4th value to calculate the boundary to +z
    properties=[0.3, 5, 5, 1e6],
    # Computational boundary: By default 1 wavelength from the center
    # (In air all that makes no sense of course. However, usually
    # 50-100 km is fine for land to shallow marine. The deeper marine,
    # the less is required.)
    center=(0, 0, -500),
    vector=(grid.nodes_x, grid.nodes_y, grid.nodes_z),
)
print(grid_ext)# QC

# Interpolate the model to the extended grid (using volume averaging on log-scale).
model_ext = model.interpolate_to_grid(grid_ext)
grid_ext.plot_3d_slicer(
    model_ext.property_x, zslice=-1900,
    xlim=(-25000, 25000), ylim=(-15000, 15000), zlim=(-4000, 500),
    pcolor_opts={'norm': LogNorm(vmin=0.3, vmax=200)}
)


x, y, z, azimuth, dip, iTx = np.loadtxt('sources.txt', skiprows=1, unpack=True)
src = (x, y, -z, azimuth, dip)
efield = emg3d.solve_source(model_ext, src, freq, sslsolver=False, verb=4)
hfield = emg3d.get_magnetic_field(model_ext, efield)

# #display Ex
# grid_ext.plot_3d_slicer(
#         efield.fx.ravel('F'), view='abs', v_type='Fx',
#         xlim=(-10000, 10000), ylim=(-10000, 10000), zlim=(-4000, 40),
#         pcolor_opts={'norm': LogNorm(vmin=1e-16, vmax=1e-8)}
# )
# #display Hy
# grid_ext.plot_3d_slicer(
#         hfield.fx.ravel('F'), view='abs', v_type='Fy',
#         xlim=(-10000, 10000), ylim=(-10000, 10000), zlim=(-4000, 40),
#         pcolor_opts={'norm': LogNorm(vmin=1e-16, vmax=1e-8)}
# )


x, y, z, azimuth, dip, iRx = np.loadtxt('receivers.txt', skiprows=1, unpack=True)
rec = (x, y, -z, azimuth, dip)
Ex = efield.get_receiver(rec)
Hx = hfield.get_receiver(rec)

rec = (x, y, -z, azimuth+90, dip)
Ey = efield.get_receiver(rec)
Hy = hfield.get_receiver(rec)

f = open("emg3d_0001.txt",'w')
f.write("iTx 	 iRx    chrec  frequency/Hz 	 Real{E/H} Imag{E/H}\n")
for iRx, (dre, dim) in enumerate(zip(Ex.real, Ex.imag)):
    f.write('1 \t %d \t Ex \t %g \t %e \t %e\n'%(iRx+1, freq, dre, dim))
for iRx, (dre, dim) in enumerate(zip(Ey.real, Ey.imag)):
    f.write('1 \t %d \t Ey \t %g \t %e \t %e\n'%(iRx+1, freq, dre, dim))
for iRx, (dre, dim) in enumerate(zip(Hx.real, Hx.imag)):
    f.write('1 \t %d \t Hx \t %g \t %e \t %e\n'%(iRx+1, freq, dre, dim))
for iRx, (dre, dim) in enumerate(zip(Hy.real, Hy.imag)):
    f.write('1 \t %d \t Hy \t %g \t %e \t %e\n'%(iRx+1, freq, dre, dim))
f.close()

# plt.figure()
# plt.subplot(211)
# plt.plot(x/1e3, abs(E_rec))

# plt.yscale('log')
# plt.xlabel('x-coordinate (km)')
# plt.ylabel('$|E_x|$ (V/m)')

# plt.subplot(212)
# plt.plot(x/1e3, np.angle(E_rec, deg=True))
# plt.xlabel('x-coordinate (km)')
# plt.ylabel('$angle(E_x)$ (degree)')

# plt.tight_layout()
# plt.savefig('result_emg3d.png')
# plt.show()


