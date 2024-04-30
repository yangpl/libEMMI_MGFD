import matplotlib.pyplot as plt
import numpy as np
import cmath
import os

iTx = []
iRx = []
chrec = []
freq = []
rmse = []
x = []
y = []
z = []
f =open('rmse_0013.txt', 'r')
header = f.readline()
for line in f:
    line = line.strip()
    columns = line.split()
    iTx.append(int(columns[0]))
    iRx.append(int(columns[1]))
    chrec.append(columns[2])
    freq.append(float(columns[3]))
    rmse.append(float(columns[4]))
    x.append(float(columns[5]))
    y.append(float(columns[6]))
    z.append(float(columns[7]))
    
iTx = np.array(iTx)
iRx = np.array(iRx)
chrec = np.array(chrec)
freq = np.array(freq)
rmse = np.array(rmse)
x = np.array(x)
y = np.array(y)
z = np.array(z)


plt.figure(figsize=(11.5,9))
plt.subplot(221)
idx = (freq==0.25) & (chrec=='Ex')
plt.scatter(x[idx], y[idx], s=10, c=rmse[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('Ex, 0.25 Hz')
plt.colorbar()


plt.subplot(222)
idx = (freq==0.25) & (chrec=='Ey')
plt.scatter(x[idx], y[idx], s=10, c=rmse[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('Ey, 0.25 Hz')
plt.colorbar()


plt.subplot(223)
idx = (freq==0.25) & (chrec=='Hx')
plt.scatter(x[idx], y[idx], s=10, c=rmse[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('Hx, 0.25 Hz')
plt.colorbar()



plt.subplot(224)
idx = (freq==0.25) & (chrec=='Hy')
plt.scatter(x[idx], y[idx], s=10, c=rmse[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('Hy, 0.25 Hz')
plt.colorbar()

plt.tight_layout(pad=0.5)
plt.savefig('rmse_scatter.png', bbox_inches='tight')
plt.show()
