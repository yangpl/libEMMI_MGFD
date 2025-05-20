/* inject source, extract EM data at receiver locations 
 * taking into account coordinate transformation between input grid (x,y,z) 
 * and acquisition grid(x',y',z')
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "acq.h"
#include "emf.h"
 

/*< find the index k in x[] such that x[k]<= val <x[k+1] >*/
int find_index(int n, float *x, float val);

/*< inject source >*/
void inject_source(acq_t *acq, emf_t *emf, complex *b, int ifreq)
{
  int ic, isrc, i, j, k;
  float w1, w2, w3, vol, phi, theta;
  complex src_vol, src_fd, Fx, Fy, Fz;

#undef id1
#undef id2
#undef id3
#define id1(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k))
#define id2(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + emf->n123pad)
#define id3(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + 2*emf->n123pad)

  src_fd = 1;
  memset(b, 0, 3*emf->n123pad*sizeof(complex));
  for(isrc=0; isrc<acq->nsrc; isrc++){
    phi = acq->src_azimuth[isrc];
    theta = acq->src_dip[isrc];
    Fx = 0;
    Fy = 0;
    Fz = 0;
    for(ic=0; ic<emf->nchsrc; ic++){
      if(strcmp(emf->chsrc[ic],"Ex")==0){
	Fx += src_fd*cos(phi);
	Fy += src_fd*sin(phi);
      }else if(strcmp(emf->chsrc[ic],"Ey")==0){
	Fx -= src_fd*sin(phi)*cos(theta);
	Fy += src_fd*cos(phi)*cos(theta);
	Fz += src_fd*sin(theta);
      }else if(strcmp(emf->chsrc[ic],"Ez")==0){
	Fx += src_fd*sin(phi)*sin(theta);
	Fy -= src_fd*cos(phi)*sin(theta);
	Fz += src_fd*cos(theta);
      }
    }

    //Jx(i+0.5,j,k)
    i = find_index(emf->n1pad, emf->x1s, acq->src_x1[isrc]);
    j = find_index(emf->n2pad, emf->x2, acq->src_x2[isrc]);
    k = find_index(emf->n3pad, emf->x3, acq->src_x3[isrc]);

    w1 = (acq->src_x1[isrc] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
    w2 = (acq->src_x2[isrc] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->src_x3[isrc] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);
    vol = (emf->x1s[i+1]- emf->x1s[i])*(emf->x2[j+1]- emf->x2[j])*(emf->x3[k+1]- emf->x3[k]);
    src_vol = I*emf->omegas[ifreq]*mu0*Fx/vol;

    b[id1(i,j,k)] += src_vol*(1.-w1)*(1.-w2)*(1.-w3);
    b[id1(i+1,j,k)] += src_vol*w1*(1.-w2)*(1.-w3);
    b[id1(i,j+1,k)] += src_vol*(1.-w1)*w2*(1.-w3);
    b[id1(i+1,j+1,k)] += src_vol*w1*w2*(1.-w3);
    b[id1(i,j,k+1)] += src_vol*(1.-w1)*(1.-w2)*w3;
    b[id1(i+1,j,k+1)] += src_vol*w1*(1.-w2)*w3;
    b[id1(i,j+1,k+1)] += src_vol*(1.-w1)*w2*w3;
    b[id1(i+1,j+1,k+1)] += src_vol*w1*w2*w3;


    //Jy(i,j+0.5,k)
    i = find_index(emf->n1pad, emf->x1, acq->src_x1[isrc]);
    j = find_index(emf->n2pad, emf->x2s, acq->src_x2[isrc]);
    k = find_index(emf->n3pad, emf->x3, acq->src_x3[isrc]);

    w1 = (acq->src_x1[isrc] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->src_x2[isrc] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
    w3 = (acq->src_x3[isrc] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);
    vol = (emf->x1[i+1]- emf->x1[i])*(emf->x2s[j+1]- emf->x2s[j])*(emf->x3[k+1]- emf->x3[k]);
    src_vol = I*emf->omegas[ifreq]*mu0*Fy/vol;

    b[id2(i,j,k)] += src_vol*(1.-w1)*(1.-w2)*(1.-w3);
    b[id2(i+1,j,k)] += src_vol*w1*(1.-w2)*(1.-w3);
    b[id2(i,j+1,k)] += src_vol*(1.-w1)*w2*(1.-w3);
    b[id2(i+1,j+1,k)] += src_vol*w1*w2*(1.-w3);
    b[id2(i,j,k+1)] += src_vol*(1.-w1)*(1.-w2)*w3;
    b[id2(i+1,j,k+1)] += src_vol*w1*(1.-w2)*w3;
    b[id2(i,j+1,k+1)] += src_vol*(1.-w1)*w2*w3;
    b[id2(i+1,j+1,k+1)] += src_vol*w1*w2*w3;

    //Jz(i,j,k+0.5)
    i = find_index(emf->n1pad, emf->x1, acq->src_x1[isrc]);
    j = find_index(emf->n2pad, emf->x2, acq->src_x2[isrc]);
    k = find_index(emf->n3pad, emf->x3s, acq->src_x3[isrc]);

    w1 = (acq->src_x1[isrc] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->src_x2[isrc] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->src_x3[isrc] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);

    vol = (emf->x1[i+1]- emf->x1[i])*(emf->x2[j+1]- emf->x2[j])*(emf->x3s[k+1]- emf->x3s[k]);
    src_vol = I*emf->omegas[ifreq]*mu0*Fz/vol;
	
    b[id3(i,j,k)] += src_vol*(1.-w1)*(1.-w2)*(1.-w3);
    b[id3(i+1,j,k)] += src_vol*w1*(1.-w2)*(1.-w3);
    b[id3(i,j+1,k)] += src_vol*(1.-w1)*w2*(1.-w3);
    b[id3(i+1,j+1,k)] += src_vol*w1*w2*(1.-w3);
    b[id3(i,j,k+1)] += src_vol*(1.-w1)*(1.-w2)*w3;
    b[id3(i+1,j,k+1)] += src_vol*w1*(1.-w2)*w3;
    b[id3(i,j+1,k+1)] += src_vol*(1.-w1)*w2*w3;
    b[id3(i+1,j+1,k+1)] += src_vol*w1*w2*w3;
  }

#undef id1
#undef id2
#undef id3
}


/*< inject adjoint source >*/
void inject_adj_source(acq_t *acq, emf_t *emf, complex *x, complex *b, int ifreq)
{
  int ic, irec, i, j, k;
  float w1, w2, w3, vol, phi, theta;
  complex rec_vol, Fx, Fy, Fz;

#undef id1
#undef id2
#undef id3
#define id1(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k))
#define id2(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + emf->n123pad)
#define id3(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + 2*emf->n123pad)

  memset(b, 0, 3*emf->n123pad*sizeof(complex));
  memset(x, 0, 3*emf->n123pad*sizeof(complex));
  for(irec=0; irec<acq->nrec; irec++){
    phi = acq->rec_azimuth[irec];
    theta = acq->rec_dip[irec];
    
    Fx = 0;
    Fy = 0;
    Fz = 0;
    for(ic=0; ic<emf->nchrec; ic++){
      if(strcmp(emf->chrec[ic],"Ex")==0){
	Fx += emf->dres_fd[ic][ifreq][irec]*cos(phi);
	Fy += emf->dres_fd[ic][ifreq][irec]*sin(phi);
      }else if(strcmp(emf->chrec[ic],"Ey")==0){
	Fx -= emf->dres_fd[ic][ifreq][irec]*sin(phi)*cos(theta);
	Fy += emf->dres_fd[ic][ifreq][irec]*cos(phi)*cos(theta);
	Fz += emf->dres_fd[ic][ifreq][irec]*sin(theta);
      }else if(strcmp(emf->chrec[ic],"Ez")==0){
	Fx += emf->dres_fd[ic][ifreq][irec]*sin(phi)*sin(theta);
	Fy -= emf->dres_fd[ic][ifreq][irec]*cos(phi)*sin(theta);
	Fz += emf->dres_fd[ic][ifreq][irec]*cos(theta);
      }
    }
    
    //Jx(i+0.5,j,k)
    i = find_index(emf->n1pad, emf->x1s, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
    w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);

    vol = (emf->x1s[i+1]- emf->x1s[i])*(emf->x2[j+1]- emf->x2[j])*(emf->x3[k+1]- emf->x3[k]);
    rec_vol = I*emf->omegas[ifreq]*mu0*Fx/vol;

    b[id1(i,j,k)] += rec_vol*(1.-w1)*(1.-w2)*(1.-w3);
    b[id1(i+1,j,k)] += rec_vol*w1*(1.-w2)*(1.-w3);
    b[id1(i,j+1,k)] += rec_vol*(1.-w1)*w2*(1.-w3);
    b[id1(i+1,j+1,k)] += rec_vol*w1*w2*(1.-w3);
    b[id1(i,j,k+1)] += rec_vol*(1.-w1)*(1.-w2)*w3;
    b[id1(i+1,j,k+1)] += rec_vol*w1*(1.-w2)*w3;
    b[id1(i,j+1,k+1)] += rec_vol*(1.-w1)*w2*w3;
    b[id1(i+1,j+1,k+1)] += rec_vol*w1*w2*w3;


    //Jy(i,j+0.5,k)
    i = find_index(emf->n1pad, emf->x1, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2s, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->rec_x2[irec] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
    w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);

    vol = (emf->x1[i+1]- emf->x1[i])*(emf->x2s[j+1]- emf->x2s[j])*(emf->x3[k+1]- emf->x3[k]);
    rec_vol = I*emf->omegas[ifreq]*mu0*Fy/vol;

    b[id2(i,j,k)] += rec_vol*(1.-w1)*(1.-w2)*(1.-w3);
    b[id2(i+1,j,k)] += rec_vol*w1*(1.-w2)*(1.-w3);
    b[id2(i,j+1,k)] += rec_vol*(1.-w1)*w2*(1.-w3);
    b[id2(i+1,j+1,k)] += rec_vol*w1*w2*(1.-w3);
    b[id2(i,j,k+1)] += rec_vol*(1.-w1)*(1.-w2)*w3;
    b[id2(i+1,j,k+1)] += rec_vol*w1*(1.-w2)*w3;
    b[id2(i,j+1,k+1)] += rec_vol*(1.-w1)*w2*w3;
    b[id2(i+1,j+1,k+1)] += rec_vol*w1*w2*w3;


    //Jz(i,j,k+0.5)
    i = find_index(emf->n1pad, emf->x1, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3s, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->rec_x3[irec] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);

    vol = (emf->x1[i+1]- emf->x1[i])*(emf->x2[j+1]- emf->x2[j])*(emf->x3s[k+1]- emf->x3s[k]);
    rec_vol = I*emf->omegas[ifreq]*mu0*Fy/vol;
	
    b[id3(i,j,k)] += rec_vol*(1.-w1)*(1.-w2)*(1.-w3);
    b[id3(i+1,j,k)] += rec_vol*w1*(1.-w2)*(1.-w3);
    b[id3(i,j+1,k)] += rec_vol*(1.-w1)*w2*(1.-w3);
    b[id3(i+1,j+1,k)] += rec_vol*w1*w2*(1.-w3);
    b[id3(i,j,k+1)] += rec_vol*(1.-w1)*(1.-w2)*w3;
    b[id3(i+1,j,k+1)] += rec_vol*w1*(1.-w2)*w3;
    b[id3(i,j+1,k+1)] += rec_vol*(1.-w1)*w2*w3;
    b[id3(i+1,j+1,k+1)] += rec_vol*w1*w2*w3;  

    //=======================================================
    Fx = 0;
    Fy = 0;
    Fz = 0;
    for(ic=0; ic<emf->nchrec; ic++){
      if(strcmp(emf->chrec[ic],"Hx")==0){
	Fx += emf->dres_fd[ic][ifreq][irec]*cos(phi);
	Fy += emf->dres_fd[ic][ifreq][irec]*sin(phi);
      }else if(strcmp(emf->chrec[ic],"Hy")==0){
	Fx -= emf->dres_fd[ic][ifreq][irec]*sin(phi)*cos(theta);
	Fy += emf->dres_fd[ic][ifreq][irec]*cos(phi)*cos(theta);
	Fz += emf->dres_fd[ic][ifreq][irec]*sin(theta);
      }else if(strcmp(emf->chrec[ic],"Hz")==0){
	Fx += emf->dres_fd[ic][ifreq][irec]*sin(phi)*sin(theta);
	Fy -= emf->dres_fd[ic][ifreq][irec]*cos(phi)*sin(theta);
	Fz += emf->dres_fd[ic][ifreq][irec]*cos(theta);
      }
    }
    
    //Mx(i,j+0.5,k+0.5)
    i = find_index(emf->n1pad, emf->x1, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2s, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3s, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->rec_x2[irec] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
    w3 = (acq->rec_x3[irec] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);

    vol = (emf->x1[i+1]- emf->x1[i])*(emf->x2s[j+1]- emf->x2s[j])*(emf->x3s[k+1]- emf->x3s[k]);
    rec_vol = Fx/vol;

    x[id1(i,j,k)] += rec_vol*(1.-w1)*(1.-w2)*(1.-w3);
    x[id1(i+1,j,k)] += rec_vol*w1*(1.-w2)*(1.-w3);
    x[id1(i,j+1,k)] += rec_vol*(1.-w1)*w2*(1.-w3);
    x[id1(i+1,j+1,k)] += rec_vol*w1*w2*(1.-w3);
    x[id1(i,j,k+1)] += rec_vol*(1.-w1)*(1.-w2)*w3;
    x[id1(i+1,j,k+1)] += rec_vol*w1*(1.-w2)*w3;
    x[id1(i,j+1,k+1)] += rec_vol*(1.-w1)*w2*w3;
    x[id1(i+1,j+1,k+1)] += rec_vol*w1*w2*w3;


    //My(i+0.5,j,k+0.5)
    i = find_index(emf->n1pad, emf->x1s, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3s, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
    w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->rec_x3[irec] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);

    vol = (emf->x1s[i+1]- emf->x1s[i])*(emf->x2[j+1]- emf->x2[j])*(emf->x3s[k+1]- emf->x3s[k]);
    rec_vol = Fy/vol;

    x[id2(i,j,k)] += rec_vol*(1.-w1)*(1.-w2)*(1.-w3);
    x[id2(i+1,j,k)] += rec_vol*w1*(1.-w2)*(1.-w3);
    x[id2(i,j+1,k)] += rec_vol*(1.-w1)*w2*(1.-w3);
    x[id2(i+1,j+1,k)] += rec_vol*w1*w2*(1.-w3);
    x[id2(i,j,k+1)] += rec_vol*(1.-w1)*(1.-w2)*w3;
    x[id2(i+1,j,k+1)] += rec_vol*w1*(1.-w2)*w3;
    x[id2(i,j+1,k+1)] += rec_vol*(1.-w1)*w2*w3;
    x[id2(i+1,j+1,k+1)] += rec_vol*w1*w2*w3;


    //Mz(i,j,k+0.5)
    i = find_index(emf->n1pad, emf->x1s, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2s, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
    w2 = (acq->rec_x2[irec] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
    w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);

    vol = (emf->x1s[i+1]- emf->x1s[i])*(emf->x2s[j+1]- emf->x2s[j])*(emf->x3[k+1]- emf->x3[k]);
    rec_vol = Fy/vol;
	
    x[id3(i,j,k)] += rec_vol*(1.-w1)*(1.-w2)*(1.-w3);
    x[id3(i+1,j,k)] += rec_vol*w1*(1.-w2)*(1.-w3);
    x[id3(i,j+1,k)] += rec_vol*(1.-w1)*w2*(1.-w3);
    x[id3(i+1,j+1,k)] += rec_vol*w1*w2*(1.-w3);
    x[id3(i,j,k+1)] += rec_vol*(1.-w1)*(1.-w2)*w3;
    x[id3(i+1,j,k+1)] += rec_vol*w1*(1.-w2)*w3;
    x[id3(i,j+1,k+1)] += rec_vol*(1.-w1)*w2*w3;
    x[id3(i+1,j+1,k+1)] += rec_vol*w1*w2*w3;  
  }
  for(k=1; k<emf->n3pad; k++){
    for(j=1; j<emf->n2pad; j++){
      for(i=1; i<emf->n1pad; i++){
	//here we assume mur=1
	b[id1(i,j,k)] += (x[id3(i,j,k)]-x[id3(i,j-1,k)])/emf->d2[j] - (x[id2(i,j,k)]-x[id2(i,j,k-1)])/emf->d3[k];
	b[id2(i,j,k)] += (x[id1(i,j,k)]-x[id1(i,j,k-1)])/emf->d3[k] - (x[id3(i,j,k)]-x[id3(i-1,j,k)])/emf->d1[i];
	b[id3(i,j,k)] += (x[id2(i,j,k)]-x[id2(i-1,j,k)])/emf->d1[i] - (x[id1(i,j,k)]-x[id1(i,j-1,k)])/emf->d2[j];
      }
    }
  }
#undef id1
#undef id2
#undef id3

}

/*< extract EM data at receiver locations >*/
void extract_emf_data(acq_t *acq, emf_t *emf, complex *x, complex *b, int ifreq)
{
  int i, j, k, ic, irec;
  float w1, w2, w3, phi, theta;
  complex Fx, Fy, Fz;

#undef id1
#undef id2
#undef id3
#define id1(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k))
#define id2(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + emf->n123pad)
#define id3(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + 2*emf->n123pad)

  for(irec=0; irec<acq->nrec; irec++){
    phi = acq->rec_azimuth[irec];
    theta = acq->rec_dip[irec];

    //=======================================================
    //Ex(i+0.5,j,k)
    i = find_index(emf->n1pad, emf->x1s, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
    w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);

    Fx = 0;
    Fx += x[id1(i,j,k)]*(1.-w1)*(1.-w2)*(1.-w3);
    Fx += x[id1(i+1,j,k)]*w1*(1.-w2)*(1.-w3);
    Fx += x[id1(i,j+1,k)]*(1.-w1)*w2*(1.-w3);
    Fx += x[id1(i+1,j+1,k)]*w1*w2*(1.-w3);
    Fx += x[id1(i,j,k+1)]*(1.-w1)*(1.-w2)*w3;
    Fx += x[id1(i+1,j,k+1)]*w1*(1.-w2)*w3;
    Fx += x[id1(i,j+1,k+1)]*(1.-w1)*w2*w3;
    Fx += x[id1(i+1,j+1,k+1)]*w1*w2*w3;

    //Ey(i,j+0.5,k)
    i = find_index(emf->n1pad, emf->x1, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2s, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->rec_x2[irec] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
    w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);

    Fy = 0;
    Fy += x[id2(i,j,k)]*(1.-w1)*(1.-w2)*(1.-w3);
    Fy += x[id2(i+1,j,k)]*w1*(1.-w2)*(1.-w3);
    Fy += x[id2(i,j+1,k)]*(1.-w1)*w2*(1.-w3);
    Fy += x[id2(i+1,j+1,k)]*w1*w2*(1.-w3);
    Fy += x[id2(i,j,k+1)]*(1.-w1)*(1.-w2)*w3;
    Fy += x[id2(i+1,j,k+1)]*w1*(1.-w2)*w3;
    Fy += x[id2(i,j+1,k+1)]*(1.-w1)*w2*w3;
    Fy += x[id2(i+1,j+1,k+1)]*w1*w2*w3;

    //Ez(i,j,k+0.5)
    i = find_index(emf->n1pad, emf->x1, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3s, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->rec_x3[irec] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);

    //interpolation over Jz which is continuous in the presence of interface
    Fz = 0;
    Fz += emf->sigma33[k][j][i]*x[id3(i,j,k)]*(1.-w1)*(1.-w2)*(1.-w3);
    Fz += emf->sigma33[k][j][i+1]*x[id3(i+1,j,k)]*w1*(1.-w2)*(1.-w3);
    Fz += emf->sigma33[k][j+1][i]*x[id3(i,j+1,k)]*(1.-w1)*w2*(1.-w3);
    Fz += emf->sigma33[k][j+1][i+1]*x[id3(i+1,j+1,k)]*w1*w2*(1.-w3);
    Fz += emf->sigma33[k+1][j][i]*x[id3(i,j,k+1)]*(1.-w1)*(1.-w2)*w3;
    Fz += emf->sigma33[k+1][j][i+1]*x[id3(i+1,j,k+1)]*w1*(1.-w2)*w3;
    Fz += emf->sigma33[k+1][j+1][i]*x[id3(i,j+1,k+1)]*(1.-w1)*w2*w3;
    Fz += emf->sigma33[k+1][j+1][i+1]*x[id3(i+1,j+1,k+1)]*w1*w2*w3;
    Fz *= emf->rho_water;//convert Jz to Ez: Ez=Jz/Fzigma_water

    //------------------------------------------------------------------------------
    //extract Ex,Ey,Ez at acquisition coordinate (x',y',z') by rotation from (x,y,z)
    //------------------------------------------------------------------------------
    for(ic=0; ic<emf->nchrec; ic++){
      if(strcmp(emf->chrec[ic],"Ex")==0){
	emf->dcal_fd[ic][ifreq][irec] = cos(phi)*Fx + sin(phi)*Fy;
      }else if(strcmp(emf->chrec[ic],"Ey")==0){
	emf->dcal_fd[ic][ifreq][irec] = -sin(phi)*cos(theta)*Fx + cos(phi)*cos(theta)*Fy - sin(theta)*Fz;
      }else if(strcmp(emf->chrec[ic],"Ez")==0){
	emf->dcal_fd[ic][ifreq][irec] = sin(phi)*sin(theta)*Fx - cos(phi)*sin(theta)*Fy + cos(theta)*Fz;
      }
    }//end for ic

    //=======================================================
    //Hx(i,j+0.5,k+0.5)
    i = find_index(emf->n1pad, emf->x1, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2s, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3s, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1[i])/(emf->x1[i+1]- emf->x1[i]);
    w2 = (acq->rec_x2[irec] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
    w3 = (acq->rec_x3[irec] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);

    Fx = 0;
    Fx += b[id1(i,j,k)]*(1.-w1)*(1.-w2)*(1.-w3);
    Fx += b[id1(i+1,j,k)]*w1*(1.-w2)*(1.-w3);
    Fx += b[id1(i,j+1,k)]*(1.-w1)*w2*(1.-w3);
    Fx += b[id1(i+1,j+1,k)]*w1*w2*(1.-w3);
    Fx += b[id1(i,j,k+1)]*(1.-w1)*(1.-w2)*w3;
    Fx += b[id1(i+1,j,k+1)]*w1*(1.-w2)*w3;
    Fx += b[id1(i,j+1,k+1)]*(1.-w1)*w2*w3;
    Fx += b[id1(i+1,j+1,k+1)]*w1*w2*w3;

    //Hy(i+0.5,j,k+0.5)
    i = find_index(emf->n1pad, emf->x1s, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3s, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
    w2 = (acq->rec_x2[irec] - emf->x2[j])/(emf->x2[j+1]- emf->x2[j]);
    w3 = (acq->rec_x3[irec] - emf->x3s[k])/(emf->x3s[k+1]- emf->x3s[k]);

    Fy = 0;
    Fy += b[id2(i,j,k)]*(1.-w1)*(1.-w2)*(1.-w3);
    Fy += b[id2(i+1,j,k)]*w1*(1.-w2)*(1.-w3);
    Fy += b[id2(i,j+1,k)]*(1.-w1)*w2*(1.-w3);
    Fy += b[id2(i+1,j+1,k)]*w1*w2*(1.-w3);
    Fy += b[id2(i,j,k+1)]*(1.-w1)*(1.-w2)*w3;
    Fy += b[id2(i+1,j,k+1)]*w1*(1.-w2)*w3;
    Fy += b[id2(i,j+1,k+1)]*(1.-w1)*w2*w3;
    Fy += b[id2(i+1,j+1,k+1)]*w1*w2*w3;

    //Hz(i+0.5,j+0.5,k)
    i = find_index(emf->n1pad, emf->x1s, acq->rec_x1[irec]);
    j = find_index(emf->n2pad, emf->x2s, acq->rec_x2[irec]);
    k = find_index(emf->n3pad, emf->x3, acq->rec_x3[irec]);

    w1 = (acq->rec_x1[irec] - emf->x1s[i])/(emf->x1s[i+1]- emf->x1s[i]);
    w2 = (acq->rec_x2[irec] - emf->x2s[j])/(emf->x2s[j+1]- emf->x2s[j]);
    w3 = (acq->rec_x3[irec] - emf->x3[k])/(emf->x3[k+1]- emf->x3[k]);

    Fz = 0;
    Fz += b[id3(i,j,k)]*(1.-w1)*(1.-w2)*(1.-w3);
    Fz += b[id3(i+1,j,k)]*w1*(1.-w2)*(1.-w3);
    Fz += b[id3(i,j+1,k)]*(1.-w1)*w2*(1.-w3);
    Fz += b[id3(i+1,j+1,k)]*w1*w2*(1.-w3);
    Fz += b[id3(i,j,k+1)]*(1.-w1)*(1.-w2)*w3;
    Fz += b[id3(i+1,j,k+1)]*w1*(1.-w2)*w3;
    Fz += b[id3(i,j+1,k+1)]*(1.-w1)*w2*w3;
    Fz += b[id3(i+1,j+1,k+1)]*w1*w2*w3;

    //------------------------------------------------------------------------------
    //extract Hx,Hy,Hz at acquisition coordinate (x',y',z') by rotation from (x,y,z)
    //------------------------------------------------------------------------------
    for(ic=0; ic<emf->nchrec; ic++){
      if(strcmp(emf->chrec[ic],"Hx")==0){
	emf->dcal_fd[ic][ifreq][irec] = cos(phi)*Fx + sin(phi)*Fy;
      }else if(strcmp(emf->chrec[ic],"Hy")==0){
	emf->dcal_fd[ic][ifreq][irec] = -sin(phi)*cos(theta)*Fx + cos(phi)*cos(theta)*Fy - sin(theta)*Fz;
      }else if(strcmp(emf->chrec[ic],"Hz")==0){
	emf->dcal_fd[ic][ifreq][irec] = sin(phi)*sin(theta)*Fx - cos(phi)*sin(theta)*Fy + cos(theta)*Fz;
      }
    }//end for ic

  }//end for irec

#undef id1
#undef id2
#undef id3
}

/*< extract EMF in the computing grid without rotation to build the inversion gradient >*/
void extract_electric_field(acq_t *acq, emf_t *emf, complex *x, int ifreq)
{
  int i, j, k;
  float *w1, *w2, *w3;
  int *i1, *i2, *i3;
  float *w1s, *w2s, *w3s;
  int *i1s, *i2s, *i3s;
  float x1, x2, x3;
  float wx, wy, wz;
  complex s;
  float sigma;
  int ii, jj, kk;

#undef id1
#undef id2
#undef id3
#define id1(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k))
#define id2(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + emf->n123pad)
#define id3(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + 2*emf->n123pad)

  w1 = alloc1float(emf->nx);
  w2 = alloc1float(emf->ny);
  w3 = alloc1float(emf->nz);
  i1 = alloc1int(emf->nx);
  i2 = alloc1int(emf->ny);
  i3 = alloc1int(emf->nz);
  w1s = alloc1float(emf->nx);
  w2s = alloc1float(emf->ny);
  w3s = alloc1float(emf->nz);
  i1s = alloc1int(emf->nx);
  i2s = alloc1int(emf->ny);
  i3s = alloc1int(emf->nz);

  for(i=0; i<emf->nx; i++){
    x1 = emf->ox + (i+0.5)*emf->dx;
    
    ii = find_index(emf->n1pad, emf->x1, x1);
    w1[i] = (x1 - emf->x1[ii])/(emf->x1[ii+1]- emf->x1[ii]);
    i1[i] = ii;

    ii = find_index(emf->n1pad, emf->x1s, x1);
    w1s[i] = (x1 - emf->x1s[ii])/(emf->x1s[ii+1]- emf->x1s[ii]);
    i1s[i] = ii;
  }
  for(j=0; j<emf->ny; j++){
    x2 = emf->oy + (j+0.5)*emf->dy;

    jj = find_index(emf->n2pad, emf->x2, x2);
    w2[j] = (x2 - emf->x2[jj])/(emf->x2[jj+1]- emf->x2[jj]);
    i2[j] = jj;

    jj = find_index(emf->n2pad, emf->x2s, x2);
    w2s[j] = (x2 - emf->x2s[jj])/(emf->x2s[jj+1]- emf->x2s[jj]);
    i2s[j] = jj;
  }
  for(k=0; k<emf->nz; k++){
    x3 = emf->oz + (k+0.5)*emf->dz;
    
    kk = find_index(emf->n3pad, emf->x3, x3);
    w3[k] = (x3 - emf->x3[kk])/(emf->x3[kk+1]- emf->x3[kk]);
    i3[k] = kk;

    kk = find_index(emf->n3pad, emf->x3s, x3);
    w3s[k] = (x3 - emf->x3s[kk])/(emf->x3s[kk+1]- emf->x3s[kk]);
    i3s[k] = kk;
  }

  //--------------------------Ex----------------------------
  for(k=0; k<emf->nz; k++){
    kk = i3[k];
    wz = w3[k];
    for(j=0; j<emf->ny; j++){
      jj = i2[j];
      wy = w2[j];
      for(i=0; i<emf->nx; i++){
  	ii = i1s[i];
  	wx = w1s[i];
	
	s = 0;
	s += x[id1(ii,jj,kk)]*(1.-wx)*(1.-wy)*(1.-wz);
	s += x[id1(ii+1,jj,kk)]*wx*(1.-wy)*(1.-wz);
	s += x[id1(ii,jj+1,kk)]*(1.-wx)*wy*(1.-wz);
	s += x[id1(ii+1,jj+1,kk)]*wx*wy*(1.-wz);
	s += x[id1(ii,jj,kk+1)]*(1.-wx)*(1.-wy)*wz;
	s += x[id1(ii+1,jj,kk+1)]*wx*(1.-wy)*wz;
	s += x[id1(ii,jj+1,kk+1)]*(1.-wx)*wy*wz;
	s += x[id1(ii+1,jj+1,kk+1)]*wx*wy*wz;
  	emf->E1[ifreq][k][j][i] = s;
      }
    }
  }

  //-------------------------Ey-----------------------------
  for(k=0; k<emf->nz; k++){
    kk = i3[k];
    wz = w3[k];
    for(j=0; j<emf->ny; j++){
      jj = i2s[j];
      wy = w2s[j];
      for(i=0; i<emf->nx; i++){
  	ii = i1[i];
  	wx = w1[i];
	
	s = 0;
	s += x[id2(ii,jj,kk)]*(1.-wx)*(1.-wy)*(1.-wz);
	s += x[id2(ii+1,jj,kk)]*wx*(1.-wy)*(1.-wz);
	s += x[id2(ii,jj+1,kk)]*(1.-wx)*wy*(1.-wz);
	s += x[id2(ii+1,jj+1,kk)]*wx*wy*(1.-wz);
	s += x[id2(ii,jj,kk+1)]*(1.-wx)*(1.-wy)*wz;
	s += x[id2(ii+1,jj,kk+1)]*wx*(1.-wy)*wz;
	s += x[id2(ii,jj+1,kk+1)]*(1.-wx)*wy*wz;
	s += x[id2(ii+1,jj+1,kk+1)]*wx*wy*wz;
  	emf->E2[ifreq][k][j][i] = s;
      }
    }
  }

  //---------------------------Ez-----------------------------
  for(k=0; k<emf->nz; k++){
    kk = i3s[k];
    wz = w3s[k];
    for(j=0; j<emf->ny; j++){
      jj = i2[j];
      wy = w2[j];
      for(i=0; i<emf->nx; i++){
  	ii = i1[i];
  	wx = w1[i];

	//we interpolate Jz=sigma*Ez as it is continuous along z
	s = 0;
	sigma = 0.25*(emf->sigma33[kk][jj][ii] + emf->sigma33[kk][jj-1][ii] + emf->sigma33[kk][jj][ii-1] + emf->sigma33[kk][jj-1][ii-1]);
	s += sigma*x[id3(ii,jj,kk)]*(1.-wx)*(1.-wy)*(1.-wz);
	sigma = 0.25*(emf->sigma33[kk][jj][ii+1] + emf->sigma33[kk][jj-1][ii+1] + emf->sigma33[kk][jj][ii] + emf->sigma33[kk][jj-1][ii]);
	s += sigma*x[id3(ii+1,jj,kk)]*wx*(1.-wy)*(1.-wz);
	sigma = 0.25*(emf->sigma33[kk][jj+1][ii] + emf->sigma33[kk][jj][ii] + emf->sigma33[kk][jj+1][ii-1] + emf->sigma33[kk][jj][ii-1]);
	s += sigma*x[id3(ii,jj+1,kk)]*(1.-wx)*wy*(1.-wz);
	sigma = 0.25*(emf->sigma33[kk][jj+1][ii+1] + emf->sigma33[kk][jj][ii+1] + emf->sigma33[kk][jj+1][ii] + emf->sigma33[kk][jj][ii]);
	s += sigma*x[id3(ii+1,jj+1,kk)]*wx*wy*(1.-wz);
	sigma = 0.25*(emf->sigma33[kk+1][jj][ii] + emf->sigma33[kk+1][jj-1][ii] + emf->sigma33[kk+1][jj][ii-1] + emf->sigma33[kk+1][jj-1][ii-1]);
	s += sigma*x[id3(ii,jj,kk+1)]*(1.-wx)*(1.-wy)*wz;
	sigma = 0.25*(emf->sigma33[kk+1][jj][ii+1] + emf->sigma33[kk+1][jj-1][ii+1] + emf->sigma33[kk+1][jj][ii] + emf->sigma33[kk+1][jj-1][ii]);
	s += sigma*x[id3(ii+1,jj,kk+1)]*wx*(1.-wy)*wz;
	sigma = 0.25*(emf->sigma33[kk+1][jj+1][ii] + emf->sigma33[kk+1][jj][ii] + emf->sigma33[kk+1][jj+1][ii-1] + emf->sigma33[kk+1][jj][ii-1]);
	s += sigma*x[id3(ii,jj+1,kk+1)]*(1.-wx)*wy*wz;
	sigma = 0.25*(emf->sigma33[kk+1][jj+1][ii+1] + emf->sigma33[kk+1][jj][ii+1] + emf->sigma33[kk+1][jj+1][ii] + emf->sigma33[kk+1][jj][ii]);
	s += sigma*x[id3(ii+1,jj+1,kk+1)]*wx*wy*wz;
  	emf->E3[ifreq][k][j][i] = s*emf->rho33[k][j][i];//convert Jz to Ez
      }
    }
  }
  
  free1float(w1);
  free1float(w2);
  free1float(w3);
  free1int(i1);
  free1int(i2);
  free1int(i3);
  free1float(w1s);
  free1float(w2s);
  free1float(w3s);
  free1int(i1s);
  free1int(i2s);
  free1int(i3s);

#undef id1
#undef id2
#undef id3
}

/*< extract EMF in the computing grid without rotation to build the inversion gradient >*/
void extract_magnetic_field(acq_t *acq, emf_t *emf, complex *b, int ifreq)
{
  int i, j, k;
  float *w1, *w2, *w3;
  int *i1, *i2, *i3;
  float *w1s, *w2s, *w3s;
  int *i1s, *i2s, *i3s;
  float x1, x2, x3;
  float wx, wy, wz;
  complex s;
  int ii, jj, kk;

#undef id1
#undef id2
#undef id3
#define id1(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k))
#define id2(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + emf->n123pad)
#define id3(i,j,k)  (i+emf->n1pad*(j + emf->n2pad*k) + 2*emf->n123pad)

  w1 = alloc1float(emf->nx);
  w2 = alloc1float(emf->ny);
  w3 = alloc1float(emf->nz);
  i1 = alloc1int(emf->nx);
  i2 = alloc1int(emf->ny);
  i3 = alloc1int(emf->nz);
  w1s = alloc1float(emf->nx);
  w2s = alloc1float(emf->ny);
  w3s = alloc1float(emf->nz);
  i1s = alloc1int(emf->nx);
  i2s = alloc1int(emf->ny);
  i3s = alloc1int(emf->nz);

  for(i=0; i<emf->nx; i++){
    x1 = emf->ox + (i+0.5)*emf->dx;
    
    ii = find_index(emf->n1pad, emf->x1, x1);
    w1[i] = (x1 - emf->x1[ii])/(emf->x1[ii+1]- emf->x1[ii]);
    i1[i] = ii;

    ii = find_index(emf->n1pad, emf->x1s, x1);
    w1s[i] = (x1 - emf->x1s[ii])/(emf->x1s[ii+1]- emf->x1s[ii]);
    i1s[i] = ii;
  }
  for(j=0; j<emf->ny; j++){
    x2 = emf->oy + (j+0.5)*emf->dy;

    jj = find_index(emf->n2pad, emf->x2, x2);
    w2[j] = (x2 - emf->x2[jj])/(emf->x2[jj+1]- emf->x2[jj]);
    i2[j] = jj;

    jj = find_index(emf->n2pad, emf->x2s, x2);
    w2s[j] = (x2 - emf->x2s[jj])/(emf->x2s[jj+1]- emf->x2s[jj]);
    i2s[j] = jj;
  }
  for(k=0; k<emf->nz; k++){
    x3 = emf->oz + (k+0.5)*emf->dz;
    
    kk = find_index(emf->n3pad, emf->x3, x3);
    w3[k] = (x3 - emf->x3[kk])/(emf->x3[kk+1]- emf->x3[kk]);
    i3[k] = kk;

    kk = find_index(emf->n3pad, emf->x3s, x3);
    w3s[k] = (x3 - emf->x3s[kk])/(emf->x3s[kk+1]- emf->x3s[kk]);
    i3s[k] = kk;
  }
  
  //--------------------------Hx------------------------------    
  for(k=0; k<emf->nz; k++){
    kk = i3s[k];
    wz = w3s[k];
    for(j=0; j<emf->ny; j++){
      jj = i2s[j];
      wy = w2s[j];
      for(i=0; i<emf->nx; i++){
  	ii = i1[i];
  	wx = w1[i];
	
	s = 0;
	s += b[id1(ii,jj,kk)]*(1.-wx)*(1.-wy)*(1.-wz);
	s += b[id1(ii+1,jj,kk)]*wx*(1.-wy)*(1.-wz);
	s += b[id1(ii,jj+1,kk)]*(1.-wx)*wy*(1.-wz);
	s += b[id1(ii+1,jj+1,kk)]*wx*wy*(1.-wz);
	s += b[id1(ii,jj,kk+1)]*(1.-wx)*(1.-wy)*wz;
	s += b[id1(ii+1,jj,kk+1)]*wx*(1.-wy)*wz;
	s += b[id1(ii,jj+1,kk+1)]*(1.-wx)*wy*wz;
	s += b[id1(ii+1,jj+1,kk+1)]*wx*wy*wz;
  	emf->H1[ifreq][k][j][i] = s;
      }
    }
  }

  //-------------------------Hy-------------------------------
  for(k=0; k<emf->nz; k++){
    kk = i3s[k];
    wz = w3s[k];
    for(j=0; j<emf->ny; j++){
      jj = i2[j];
      wy = w2[j];
      for(i=0; i<emf->nx; i++){
  	ii = i1s[i];
  	wx = w1s[i];
	
	s = 0;
	s += b[id2(ii,jj,kk)]*(1.-wx)*(1.-wy)*(1.-wz);
	s += b[id2(ii+1,jj,kk)]*wx*(1.-wy)*(1.-wz);
	s += b[id2(ii,jj+1,kk)]*(1.-wx)*wy*(1.-wz);
	s += b[id2(ii+1,jj+1,kk)]*wx*wy*(1.-wz);
	s += b[id2(ii,jj,kk+1)]*(1.-wx)*(1.-wy)*wz;
	s += b[id2(ii+1,jj,kk+1)]*wx*(1.-wy)*wz;
	s += b[id2(ii,jj+1,kk+1)]*(1.-wx)*wy*wz;
	s += b[id2(ii+1,jj+1,kk+1)]*wx*wy*wz;
  	emf->H2[ifreq][k][j][i] = s;
      }
    }
  }

  //---------------------------Hz-----------------------------
  for(k=0; k<emf->nz; k++){
    kk = i3[k];
    wz = w3[k];
    for(j=0; j<emf->ny; j++){
      jj = i2s[j];
      wy = w2s[j];
      for(i=0; i<emf->nx; i++){
  	ii = i1s[i];
  	wx = w1s[i];

	s = 0;
	s += b[id3(ii,jj,kk)]*(1.-wx)*(1.-wy)*(1.-wz);
	s += b[id3(ii+1,jj,kk)]*wx*(1.-wy)*(1.-wz);
	s += b[id3(ii,jj+1,kk)]*(1.-wx)*wy*(1.-wz);
	s += b[id3(ii+1,jj+1,kk)]*wx*wy*(1.-wz);
	s += b[id3(ii,jj,kk+1)]*(1.-wx)*(1.-wy)*wz;
	s += b[id3(ii+1,jj,kk+1)]*wx*(1.-wy)*wz;
	s += b[id3(ii,jj+1,kk+1)]*(1.-wx)*wy*wz;
	s += b[id3(ii+1,jj+1,kk+1)]*wx*wy*wz;
  	emf->H3[ifreq][k][j][i] = s;
      }
    }
  }
  
  free1float(w1);
  free1float(w2);
  free1float(w3);
  free1int(i1);
  free1int(i2);
  free1int(i3);
  free1float(w1s);
  free1float(w2s);
  free1float(w3s);
  free1int(i1s);
  free1int(i2s);
  free1int(i3s);

#undef id1
#undef id2
#undef id3
}
