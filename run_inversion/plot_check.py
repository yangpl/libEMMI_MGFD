import empymod
import numpy as np
import matplotlib.pyplot as plt

#=====================================================================
#define a function to read a file based on given filename 
def read_emdata(filename):
    iTx = []
    iRx = []
    chrec = []
    freq = []
    dre = []
    dim = []
    
    f = open(filename, newline='', mode='r')
    header = f.readline()
    for line in f:
        thisline = line.strip()
        columns = thisline.split()
        iTx.append(int(columns[0]))
        iRx.append(int(columns[1]))
        chrec.append(columns[2])
        freq.append(float(columns[3]))
        dre.append(float(columns[4]))
        dim.append(float(columns[5]))
    iTx = np.array(iTx)
    iRx = np.array(iRx)
    chrec = np.array(chrec)
    freq = np.array(freq)
    dre = np.array(dre)
    dim = np.array(dim)
    dat = dre - dim*1j#note our Fourier transform convention: \partial_t <--> -i*omega

    return iTx, iRx, chrec, freq, dat
    
#===============================================================
iRx, x2, x3, azimuth, dip, iRx = np.loadtxt('receivers.txt', skiprows=1, unpack=True)

#-----------------------------------------------
iTx, iRx, chrec, freq, dat = read_emdata('emf_0001.txt')
idx = (freq==0.25) & (chrec=='Ex')

iTx, iRx2, chrec, freq, dat2 = read_emdata('syn_0001.txt')
idx2 = (freq==0.25) & (chrec=='Ex')


plt.figure(1,figsize=(12,8))
plt.subplot(221)
#idx = (freq==0.25) & (chrec=='Ex')
plt.plot(iRx[idx], np.abs(dat[idx]), 'r', label='$E_x^{obs}$-0.25 Hz')
plt.plot(iRx2[idx2], np.abs(dat2[idx2]), 'b', label='$E_x^{syn}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.yscale('log')
plt.ylabel('$|E_x|$ (V/m)')
plt.title('(a) Amplitude')
plt.legend()

plt.subplot(222)
#idx = (freq==0.25) & (chrec=='Ex')
plt.plot(iRx[idx], np.angle(dat[idx], deg=True), 'r', label='$E_x^{obs}$-0.25 Hz')
plt.plot(iRx2[idx2], np.angle(dat2[idx2], deg=True), 'b', label='$E_x^{syn}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.ylabel(r'$\angle E_x (^o)$')
plt.title('(b) Phase')
plt.legend()

plt.subplot(223)
#idx = (freq==0.25) & (chrec=='Ey')
plt.plot(iRx[idx], np.abs(dat[idx]), 'r', label='$E_y^{obs}$-0.25 Hz')
plt.plot(iRx2[idx2], np.abs(dat2[idx2]), 'b', label='$E_y^{syn}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.yscale('log')
plt.ylabel('$|E_y|$ (A/m)')
plt.title('(c) Amplitude')
plt.legend()

plt.subplot(224)
#idx = (freq==0.25) & (chrec=='Ey')
plt.plot(iRx[idx], np.angle(dat[idx], deg=True), 'r', label='$E_y^{obs}$-0.25 Hz')
plt.plot(iRx2[idx2], np.angle(dat2[idx2], deg=True), 'b', label='$E_y^{syn}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.ylabel(r'$\angle E_y (^o)$')
plt.title('(d) Phase')
plt.legend()


plt.tight_layout()
plt.savefig('amplitude_phase.png')
plt.show()

# #=========================================================
# plt.figure(2,figsize=(12,8))
# plt.subplot(221)
# idx = (freq==0.25) & (chrec=='Ex')
# rms = np.abs((dat[idx] - dat2[idx])/dat[idx])
# plt.plot(iRx[idx], rms, 'r', label='$E_x^{obs}$-0.25 Hz')


# plt.grid()
# plt.xlabel('Offset (m)')
# plt.yscale('log')
# plt.ylabel('$|E_x|$ (V/m)')
# plt.title('(a) Amplitude')
# plt.legend()

# plt.subplot(222)
# idx = (freq==0.25) & (chrec=='Ex')
# pha = np.angle((dat[idx] - dat2[idx])/dat[idx], deg=True)
# plt.plot(iRx[idx], pha, 'r', label='$E_x^{obs}$-0.25 Hz')


# plt.grid()
# plt.xlabel('Offset (m)')
# plt.ylabel(r'$\angle E_x (^o)$')
# plt.title('(b) Phase')
# plt.legend()

# plt.subplot(223)
# idx = (freq==0.25) & (chrec=='Ey')
# rms = np.abs((dat[idx] - dat2[idx])/dat[idx])
# plt.plot(iRx[idx], rms, 'r', label='$E_y^{obs}$-0.25 Hz')


# plt.grid()
# plt.xlabel('Offset (m)')
# plt.yscale('log')
# plt.ylabel('$|E_y|$ (A/m)')
# plt.title('(c) Amplitude')
# plt.legend()

# plt.subplot(224)
# idx = (freq==0.25) & (chrec=='Ey')
# pha = np.angle((dat[idx] - dat2[idx])/dat[idx], deg=True)

# plt.plot(iRx[idx], pha, 'r', label='$E_y^{obs}$-0.25 Hz')


# plt.grid()
# plt.xlabel('Offset (m)')
# plt.ylabel(r'$\angle E_y (^o)$')
# plt.title('(d) Phase')
# plt.legend()


# plt.tight_layout()
# plt.savefig('amplitude_phase.png')
# plt.show()

