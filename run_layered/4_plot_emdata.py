import empymod
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------
x, y, z, azimuth, dip, irec=np.loadtxt('receivers.txt', skiprows=1, unpack=True) 

# Define column names and dtypes
col_names = ['iTx', 'iRx', 'chrec', 'frequency', 'Real', 'Imag']
col_dtypes = ['i4', 'i4', 'U10', 'f8', 'f8', 'f8']  # i4=int, U10=string up to 10 chars, f8=float

# Read the data
iTx, iRx, chrec, freq, dre, dim = np.genfromtxt('emf_0001.txt', 
                                                dtype=col_dtypes, 
                                                names=col_names,
                                                delimiter=None,  # automatic delimiter detection
                                                skip_header=1,   # skip the header line
                                                encoding='utf-8',
                                                unpack=True)
                     
dat = dre + dim*1j

#-----------------------------------------------ab=11,Ex                                                                                              
ref1d = empymod.dipole(ab=11,src=[0, 0, 950],
                       rec=[x, x*0, 1000],
                       depth=[0, 1000, 2300, 2400],
		       freqtime=[0.25,0.75,1.25],
                       res=[1e8, 0.3, 1.0, 50, 2],
                       aniso=[1,1,1,1.414,1],
                       verb=1)
ref1d = np.conj(ref1d)


plt.figure(figsize=(10,7))
plt.subplot(221)
idx = (freq==0.25) & (chrec=='Ex')
plt.plot(x, np.abs(dat[idx]), 'r--', markevery=50, marker='o', label='$E_x^{FD}$-0.25 Hz')
# idx = (ifreq==2) & (chrec=='Ex')
# plt.plot(x, np.abs(dat[idx]), 'g--', markevery=50, marker='d', label='$E_x^{FD}$-0.75 Hz')
# idx = (ifreq==3) & (chrec=='Ex')
# plt.plot(x, np.abs(dat[idx]), 'b--', markevery=50, marker='v', label='$E_x^{FD}$-1.25 Hz')

plt.plot(x, np.abs(ref1d[0,:]), 'r', label='$E_x^{ref}$-0.25 Hz')
# plt.plot(x, np.abs(ref1d[1,:]), 'g', label='$E_x^{ref}$-0.75 Hz')
# plt.plot(x, np.abs(ref1d[2,:]), 'b', label='$E_x^{ref}$-1.25 Hz')

plt.grid()

plt.yscale('log')
plt.ylabel('Amplitude (V/Am$^2$)')
plt.xlabel('Offset (m)')
plt.title('(a) Amplitude', fontweight='bold')
plt.legend()


plt.subplot(222)
idx = (freq==0.25) & (chrec=='Ex')
plt.plot(x, np.angle(dat[idx], deg=True), 'r--', markevery=50, marker='o', label='$E_x^{FD}$-0.25 Hz')
# idx = (ifreq==2) & (chrec=='Ex')
# plt.plot(x, np.angle(dat[idx], deg=True), 'g--', markevery=50, marker='d',  label='$E_x^{FD}$-0.75 Hz')
# idx = (ifreq==3) & (chrec=='Ex')
# plt.plot(x, np.angle(dat[idx], deg=True), 'b--', markevery=50, marker='v', label='$E_x^{FD}$-1.25 Hz')

plt.plot(x, np.angle(ref1d[0,:], deg=True), 'r', label='$E_x^{ref}$-0.25 Hz')
# plt.plot(x, np.angle(ref1d[1,:], deg=True), 'g', label='$E_x^{ref}$-0.75 Hz')
# plt.plot(x, np.angle(ref1d[2,:], deg=True), 'b', label='$E_x^{ref}$-1.25 Hz')

plt.grid()
plt.ylabel('Degree ($^o$)')
plt.xlabel('Offset (m)')
plt.title('(b) Phase', fontweight='bold')
plt.legend()


plt.subplot(223)
idx = (freq==0.25) & (chrec=='Ex')
err = np.abs(dat[idx])/np.abs(ref1d[0,:])-1
plt.plot(x, 100*err, 'r', markevery=50, marker='o', label='$E_x$-0.25 Hz')
# idx = (ifreq==2) & (chrec=='Ex')
# err = np.abs(dat[idx])/np.abs(ref1d[1,:])-1
# plt.plot(x, 100*err, 'g', markevery=50, marker='d', label='$E_x$-0.75 Hz')
# idx = (ifreq==3) & (chrec=='Ex')
# err = np.abs(dat[idx])/np.abs(ref1d[2,:])-1
# plt.plot(x, 100*err, 'b', markevery=50, marker='v', label='$E_x$-1.25 Hz')

plt.ylim([-10,10])
plt.grid()
plt.ylabel('$|E_x^{FD}|/|E_x^{ref}|$-1 (%)')
plt.xlabel('Offset (m)')
plt.title('(c) Amplitude error', fontweight='bold')
plt.legend()


plt.subplot(224)
idx = (freq==0.25) & (chrec=='Ex')
plt.plot(x, np.angle(dat[idx]/ref1d[0,:], deg=True), 'r', markevery=50, marker='o', label='$E_x$-0.25 Hz')
# idx = (ifreq==2) & (chrec=='Ex')
# plt.plot(x, np.angle(dat[idx]/ref1d[1,:], deg=True), 'g', markevery=50, marker='d', label='$E_x$-0.75 Hz')
# idx = (ifreq==3) & (chrec=='Ex')
# plt.plot(x, np.angle(dat[idx]/ref1d[2,:], deg=True), 'b', markevery=50, marker='v', label='$E_x$-1.25 Hz')
plt.ylim([-5,5])
plt.grid()
plt.ylabel('Degree ($^o$)')
plt.xlabel('Offset (m)')
plt.title('(d) Phase error', fontweight='bold')
plt.legend()

plt.tight_layout()
plt.savefig('comparison_ex.png')
plt.show()





