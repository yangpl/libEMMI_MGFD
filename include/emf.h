#ifndef emf_h
#define emf_h

typedef struct {
  int mode; //mode=0, modeling; mode=1, FWI; mode=2, FWI gradient only
  int verb;/* verbose display */
  int reciprocity;

  int nfreq;//number of frequencies
  float *freqs, *omegas;//a list of frequencies
  float rho_air;//resistivity in air
  float rho_water;//resistivity in sea water
  float lextend;//distance to extend survey domain on each side
  
  int nchsrc, nchrec;//number of active src/rec channels for E and H
  char **chsrc, **chrec;

  float x1min, x1max, x2min, x2max, x3min, x3max;//include extended modelling domain
  float *xx, *yy, *zz;

  int n1, n2, n3;
  float d1min, d2min, d3min;
  float *x1, *x2, *x3;
  float *x1s, *x2s, *x3s;
  float *d1, *d2, *d3;
  float *d1s, *d2s, *d3s;
  
  float ***sigma11, ***sigma22, ***sigma33, ***invmur;
  _Complex float ***dcal_fd, ***dobs_fd, ***dres_fd;
  _Complex float ****E1, ****E2, ****E3;//EM fields
  _Complex float ****H1, ****H2, ****H3;//EM fields
  
  float ***mur, ***rho11, ***rho22, ***rho33;
  float **bathy;
  float rhomin, rhomax, depthmax;
  
  float ox, oy, oz;
  float dx, dy, dz;
  int nx, ny, nz, nb;

  int addnoise;
  float amp_perc, delta_phi, delta_theta;//magnitude, azimuth and dip uncertainties  
  float noisefloorE, noisefloorH;

  float ***delta_emf;//uncertainties for EMF
  float ***rmse;
  float dsmute, drmute;
} emf_t; /* type of electromagnetic field (emf)  */

#endif
