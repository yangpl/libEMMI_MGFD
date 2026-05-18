import os
import tempfile
from pathlib import Path

cache_dir = Path(tempfile.gettempdir()) / "libemmi_mgfd_plot_cache"
numba_cache_dir = cache_dir / "numba"
mpl_cache_dir = cache_dir / "matplotlib"
numba_cache_dir.mkdir(parents=True, exist_ok=True)
mpl_cache_dir.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(numba_cache_dir))
os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache_dir))

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


def select_data(iTx, iRx, chrec, freq, dat, tx, target_freq, component):
    idx = (iTx == tx) & np.isclose(freq, target_freq) & (chrec == component)
    order = np.argsort(iRx[idx])
    return iRx[idx][order], dat[idx][order]
    
#===============================================================
rx_x, rx_y, rx_z, azimuth, dip, rx_id = np.loadtxt('receivers.txt', skiprows=1, unpack=True)
rx_x_by_id = dict(zip(rx_id.astype(int), rx_x))

#-----------------------------------------------
iTx_emf, iRx_emf, chrec_emf, freq_emf, dat_emf = read_emdata('emf_0001.txt')
iTx_emg3d, iRx_emg3d, chrec_emg3d, freq_emg3d, dat_emg3d = read_emdata('emg3d_0001.txt')
dat_emg3d = np.conj(dat_emg3d)
#iTx, iRx, chrec, freq, dat3 = read_emdata('emf_0001_logcond.txt')

tx = 1
plot_freq = 0.25

iRx_ex_emf, ex_emf = select_data(iTx_emf, iRx_emf, chrec_emf, freq_emf, dat_emf, tx, plot_freq, 'Ex')
iRx_hy_emf, hy_emf = select_data(iTx_emf, iRx_emf, chrec_emf, freq_emf, dat_emf, tx, plot_freq, 'Hy')
iRx_ex_emg3d, ex_emg3d = select_data(iTx_emg3d, iRx_emg3d, chrec_emg3d, freq_emg3d, dat_emg3d, tx, plot_freq, 'Ex')
iRx_hy_emg3d, hy_emg3d = select_data(iTx_emg3d, iRx_emg3d, chrec_emg3d, freq_emg3d, dat_emg3d, tx, plot_freq, 'Hy')

x_ex_emf = np.array([rx_x_by_id[i] for i in iRx_ex_emf])
x_hy_emf = np.array([rx_x_by_id[i] for i in iRx_hy_emf])
x_ex_emg3d = np.array([rx_x_by_id[i] for i in iRx_ex_emg3d])
x_hy_emg3d = np.array([rx_x_by_id[i] for i in iRx_hy_emg3d])

#-----------------------------------------------
offset = rx_x #np.arange(-8000, 8001, 25)
ref1d11 = empymod.dipole(ab=11,src=[0, 0, 950],
                      rec=[offset, offset*0, 1000],
                      depth=[0, 1000, 1480, 1600],
                      freqtime=[plot_freq],
                      res=[1e8, 0.3, 1.5, 100, 2.5],
                      aniso=[1,1,1,1,1],
                      verb=1)
ref1d51 = empymod.dipole(ab=51,src=[0, 0, 950],
                      rec=[offset, offset*0, 1000],
                      depth=[0, 1000, 1480, 1600],
                      freqtime=[plot_freq],
                      res=[1e8, 0.3, 1.5, 100, 2.5],
                      aniso=[1,1,1,1,1],
                      verb=1)
pha = {'deg': True, 'lag': True, 'unwrap': False}


plt.figure(figsize=(12,8))
plt.subplot(221)
plt.plot(offset, ref1d11.amp(), 'b', label='$E_x^{ref}$-0.25 Hz')
plt.plot(x_ex_emf, np.abs(ex_emf), 'r', label='$E_x^{libEMMI}$-0.25 Hz')
plt.plot(x_ex_emg3d, np.abs(ex_emg3d), 'k', label='$E_x^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.abs(dat3[idx]), 'g', label='$E_x^{logcond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.yscale('log')
plt.ylabel('$|E_x|$ (V/m)')
plt.title('(a) Amplitude $E_x$')
plt.legend()

plt.subplot(222)
plt.plot(offset, ref1d11.pha(**pha), 'b', label='$E_x^{ref}$-0.25 Hz')
plt.plot(x_ex_emf, np.angle(ex_emf, deg=True), 'r', label='$E_x^{libEMMI}$-0.25 Hz')
plt.plot(x_ex_emg3d, np.angle(ex_emg3d, deg=True), 'k', label='$E_x^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.angle(dat3[idx], deg=True), 'g', label='$E_x^{logcond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.ylabel(r'$\angle E_x (^o)$')
plt.title('(b) Phase $E_x$')
plt.legend()

plt.subplot(223)
plt.plot(offset, ref1d51.amp(), 'b', label='$H_y^{ref}$-0.25 Hz')
plt.plot(x_hy_emf, np.abs(hy_emf), 'r', label='$H_y^{libEMMI}$-0.25 Hz')
plt.plot(x_hy_emg3d, np.abs(hy_emg3d), 'k', label='$H_y^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.abs(dat3[idx]), 'g', label='$H_y^{logcond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.yscale('log')
plt.ylabel('$|H_y|$ (A/m)')
plt.title('(c) Amplitude $H_y$')
plt.legend()

plt.subplot(224)
plt.plot(offset, ref1d51.pha(**pha), 'b', label='$H_y^{ref}$-0.25 Hz')
plt.plot(x_hy_emf, np.angle(hy_emf, deg=True), 'r', label='$H_y^{libEMMI}$-0.25 Hz')
plt.plot(x_hy_emg3d, np.angle(hy_emg3d, deg=True), 'k', label='$H_y^{emg3d}$-0.25 Hz')
#plt.plot(x1, np.angle(dat3[idx], deg=True), 'g', label='$H_y^{loncond}$-0.25 Hz')

plt.grid()
plt.xlabel('Offset (m)')
plt.ylabel(r'$\angle H_y (^o)$')
plt.title('(d) Phase $H_y$')
plt.legend()


plt.tight_layout()
plt.savefig('amplitude_phase.png')
if plt.get_backend().lower() != 'agg':
    plt.show()
