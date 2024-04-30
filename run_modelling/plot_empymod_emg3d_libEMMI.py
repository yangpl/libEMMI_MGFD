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
x1, x2, x3, azimuth, dip, iRx = np.loadtxt('receivers.txt', skiprows=1, unpack=True)

#-----------------------------------------------
iTx, iRx, chrec, freq, dat = read_emdata('emf_0001.txt')
iTx, iRx, chrec, freq, dat2 = read_emdata('emg3d_0001.txt')
dat2 = np.conj(dat2)
#iTx, iRx, chrec, freq, dat3 = read_emdata('emf_0001_logcond.txt')

#-----------------------------------------------
offset = x1 #np.arange(-8000, 8001, 25)
ref1d11 = empymod.dipole(ab=11,src=[0, 0, 950],
                      rec=[offset, offset*0, 1000],
                      depth=[0, 1000, 1500, 1600],
                      freqtime=[0.25],
                      res=[1e8, 0.3, 1.5, 100, 2.5],
                      aniso=[1,1,1,1,1],
                      verb=1)
ref1d15 = empymod.dipole(ab=15,src=[0, 0, 950],
                      rec=[offset, offset*0, 1000],
                      depth=[0, 1000, 1500, 1600],
                      freqtime=[0.25],
                      res=[1e8, 0.3, 1.5, 100, 2.5],
                      aniso=[1,1,1,1,1],
                      verb=1)
pha = {'deg': True, 'lag': True, 'unwrap': False}


plt.figure(figsize=(12,8))
plt.subplot(221)
idx = (freq==0.25) & (chrec=='Ex')
plt.plot(x1, ref1d11.amp(), 'b', label='$E_x^{ref}$-0.25 Hz')
plt.plot(x1, np.abs(dat[idx]), 'r', label='$E_x^{libEMMI}$-0.25 Hz')
plt.plot(x1, np.abs(dat2[idx]), 'k', label='$E_x^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.abs(dat3[idx]), 'g', label='$E_x^{logcond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.yscale('log')
plt.ylabel('$|E_x|$ (V/m)')
plt.title('(a) Amplitude $E_x$')
plt.legend()

plt.subplot(222)
idx = (freq==0.25) & (chrec=='Ex')
plt.plot(x1, ref1d11.pha(**pha), 'b', label='$E_x^{ref}$-0.25 Hz')
plt.plot(x1, np.angle(dat[idx], deg=True), 'r', label='$E_x^{libEMMI}$-0.25 Hz')
plt.plot(x1, np.angle(dat2[idx], deg=True), 'k', label='$E_x^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.angle(dat3[idx], deg=True), 'g', label='$E_x^{logcond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.ylabel(r'$\angle E_x (^o)$')
plt.title('(b) Phase $E_x$')
plt.legend()

plt.subplot(223)
idx = (freq==0.25) & (chrec=='Hy')
plt.plot(x1, ref1d15.amp(), 'b', label='$H_y^{ref}$-0.25 Hz')
plt.plot(x1, np.abs(dat[idx]), 'r', label='$H_y^{libEMMI}$-0.25 Hz')
plt.plot(x1, np.abs(dat2[idx]), 'k', label='$H_y^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.abs(dat3[idx]), 'g', label='$H_y^{logcond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.yscale('log')
plt.ylabel('$|H_y|$ (A/m)')
plt.title('(c) Amplitude $H_y$')
plt.legend()

plt.subplot(224)
idx = (freq==0.25) & (chrec=='Hy')
plt.plot(x1, ref1d15.pha(**pha)+180, 'b', label='$H_y^{ref}$-0.25 Hz')
plt.plot(x1, np.angle(dat[idx], deg=True), 'r', label='$H_y^{libEMMI}$-0.25 Hz')
plt.plot(x1, np.angle(dat2[idx], deg=True), 'k', label='$H_y^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.angle(dat3[idx], deg=True), 'g', label='$H_y^{loncond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.ylabel(r'$\angle H_y (^o)$')
plt.title('(d) Phase $H_y$')
plt.legend()


plt.tight_layout()
plt.savefig('amplitude_phase.png')
plt.show()

