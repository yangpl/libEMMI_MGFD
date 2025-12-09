/* Multigrid method for 3D CSEM modelling in the frequency domain
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include "emf.h"
#include "acq.h"

void regrid_init(acq_t *acq, emf_t *emf);
void regrid_free(emf_t *emf);

void gmg_init(emf_t *emf, int ifreq);
void gmg_free();
void gmg_apply(int n, complex *b, complex *x);

void inject_source(acq_t *acq, emf_t *emf, complex *b, int ifreq);
void extract_emf_data(acq_t *acq, emf_t *emf, complex *x, complex *b, int ifreq);
void extract_electric_field(acq_t *acq, emf_t *emf, complex *x, int ifreq);
void extract_magnetic_field(acq_t *acq, emf_t *emf, complex *b, int ifreq);

void write_data(acq_t *acq, emf_t *emf, char *fname, _Complex float ***dcal_fd);

void do_modelling(acq_t *acq, emf_t *emf)
{
  int ifreq, n, icheck;
  complex *x, *b;
  FILE *fp;
  char fname[sizeof("emf_0000.txt")];

  if(!getparint("icheck", &icheck)) icheck = 0;//1=check out EM field
  if(icheck){
    emf->E1 = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
    emf->E2 = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
    emf->E3 = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
    emf->H1 = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
    emf->H2 = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
    emf->H3 = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
  }
  emf->dcal_fd = alloc3complexf(acq->nrec, emf->nfreq, emf->nchrec);
  n = 3*(emf->n1+1)*(emf->n2+1)*(emf->n3+1);
  x = alloc1complex(n);//vector E=(Ex,Ey,Ez)^T
  b = alloc1complex(n);//vector Js=(Jx,Jy,Jz)^T
  memset(b, 0, n*sizeof(complex));//initialize x=0
  regrid_init(acq, emf);  
  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
    gmg_init(emf, ifreq);

    inject_source(acq, emf, b, ifreq);//initialize b=i*omega*mu*Js
    memset(x, 0, n*sizeof(complex));//initialize x=0
    
    gmg_apply(n, b, x);//after multigrid convergence, x=(Ex,Ey,Ez)^T
    extract_emf_data(acq, emf, x, b, ifreq);
    if(icheck){
      extract_electric_field(acq, emf, x, ifreq);
      extract_magnetic_field(acq, emf, b, ifreq);
    }
    
    gmg_free();
  }
  free1complex(x);
  free1complex(b);
  regrid_free(emf);
    
  sprintf(fname, "emf_%04d.txt", acq->shot_idx[iproc]);
  write_data(acq, emf, fname, emf->dcal_fd);
  
  if(icheck){
    if(emf->verb){
      n = emf->nfreq*emf->nx*emf->ny*emf->nz;
      fp = fopen("Ex.bin", "wb");
      fwrite(&emf->E1[0][0][0][0], n*sizeof(_Complex float), 1, fp);
      fclose(fp);
      fp = fopen("Ey.bin", "wb");
      fwrite(&emf->E2[0][0][0][0], n*sizeof(_Complex float), 1, fp);
      fclose(fp);
      fp = fopen("Ez.bin", "wb");
      fwrite(&emf->E3[0][0][0][0], n*sizeof(_Complex float), 1, fp);
      fclose(fp);
      fp = fopen("Hx.bin", "wb");
      fwrite(&emf->H1[0][0][0][0], n*sizeof(_Complex float), 1, fp);
      fclose(fp);
      fp = fopen("Hy.bin", "wb");
      fwrite(&emf->H2[0][0][0][0], n*sizeof(_Complex float), 1, fp);
      fclose(fp);
      fp = fopen("Hz.bin", "wb");
      fwrite(&emf->H3[0][0][0][0], n*sizeof(_Complex float), 1, fp);
      fclose(fp);
    }

    free4complexf(emf->E1);
    free4complexf(emf->E2);
    free4complexf(emf->E3);
    free4complexf(emf->H1);
    free4complexf(emf->H2);
    free4complexf(emf->H3);
  }
  free3complexf(emf->dcal_fd);  
}

