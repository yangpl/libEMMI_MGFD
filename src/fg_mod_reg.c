/* add model regularization term for function and gradient evaluation
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2021 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "acq.h"
#include "emf.h"
#include "fwi.h"
 

float regularization_tikhonov(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);
float regularization_tv(float *x, float *g, int n1, int n2, int n3, float d1, float d2, float d3);

/*< model regularization (assume fwi gradient has been computed and stored in grad) >*/
void fg_mod_reg(acq_t *acq, emf_t *emf, fwi_t *fwi, float *x, float *g)
{
  FILE *fp;
  int i1, i2, i3, i, k, nxyz;
  float gamma1, gamma2, beta, fcost_mod;
  float *gm, *xm;

  nxyz = emf->nx*emf->ny*emf->nz;
  beta = pow(0.8, fwi->iter);//beta is a cooling factor
  fwi->fcost_mod = 0;
  xm = alloc1float(fwi->n);
  gm = alloc1float(fwi->n);
  for(k=0; k<fwi->npar; k++){
    for(i3=0; i3<emf->nz; i3++){
      for(i2=0; i2<emf->ny; i2++){
	for(i1=0; i1<emf->nx; i1++){
	  i = i1 + emf->nx*(i2 + emf->ny*(i3+emf->nz*k));
	  xm[i] = (emf->oz+(i3+0.5)*emf->dz<emf->bathy[i2][i1]+emf->dz)?0:(x[i]-fwi->xref[i]);
	}
      }
    }
  }

  if(fwi->gamma1>0) {
    memset(gm, 0, fwi->n*sizeof(float));
    gamma1 = fwi->gamma1*beta;
    for(k=0; k<fwi->npar; k++){
      fcost_mod = regularization_tikhonov(&xm[k*nxyz], &gm[k*nxyz], emf->nx, emf->ny, emf->nz, emf->dx, emf->dy, emf->dz);
      fwi->fcost_mod += 0.5*gamma1*fcost_mod;
      for(i3= 0; i3<emf->nz; i3++){
	for(i2 = 0; i2<emf->ny; i2++){
	  for(i1 = 0; i1<emf->nx; i1++){
	    i = i1 + emf->nx*(i2 + emf->ny*(i3+emf->nz*k));
	    g[i] += 0.5*gamma1*gm[i];
	  }
	}
      }
    }//end for k    
    if(emf->verb){
      fp = fopen("gradient_tikhonov_Rh", "wb");
      fwrite(gm, nxyz*sizeof(float), 1, fp);
      fclose(fp);
      fp = fopen("gradient_tikhonov_Rv", "wb");
      fwrite(gm+nxyz, nxyz*sizeof(float), 1, fp);
      fclose(fp);
    }
  }

  if(fwi->gamma2>0){
    memset(gm, 0, fwi->n*sizeof(float));
    gamma2 = fwi->gamma2*beta;
    for(k=0; k<fwi->npar; k++){
      fcost_mod = regularization_tv(&xm[k*nxyz], &gm[k*nxyz], emf->nx, emf->ny, emf->nz, emf->dx, emf->dy, emf->dz);
      fwi->fcost_mod += 0.5*gamma2*fcost_mod;
      for(i3= 0; i3<emf->nz; i3++){
	for(i2 = 0; i2<emf->ny; i2++){
	  for(i1 = 0; i1<emf->nx; i1++){
	    i = i1 + emf->nx*(i2 + emf->ny*(i3+emf->nz*k));
	    g[i] += 0.5*gamma2*gm[i];
	  }
	}
      }
    }//end for k    
    if(emf->verb){
      fp = fopen("gradient_tv_Rh", "wb");
      fwrite(gm, nxyz*sizeof(float), 1, fp);
      fclose(fp);
      fp = fopen("gradient_tv_Rv", "wb");
      fwrite(gm+nxyz, nxyz*sizeof(float), 1, fp);
      fclose(fp);
    }
  }

  free1float(xm);
  free1float(gm);
}


/*< model regularization (assume fwi gradient has been computed and stored in grad) >*/
void Hv_mod_reg(emf_t *emf, fwi_t *fwi, float *r, float *Hv)
{
  int i1, i2, i3, k0, kp1, km1;
  float gamma1, t1, t2, t3;
    
  float beta = pow(0.8, fwi->iter);//beta is a cooling factor
  float _d1 = 1.;///(d1*d1);
  float _d2 = 1.;///(d2*d2);
  float _d3 = 1.;///(d3*d3);
  float c_h = 1;//large coefficient to penalize horizontal changes
  float c_v = 0.03;//small coefficients to panalize vertical changes
  gamma1 = fwi->gamma1*beta;
  
  for(i3=1; i3<emf->nz-1; i3++){
    for(i2=1; i2<emf->ny-1; i2++){
      for(i1=1; i1<emf->nx-1; i1++){

	k0 = i1 + emf->nx*(i2 + emf->ny*i3);
	if(i3>emf->bathy[i2][i1]/emf->dz){
	  km1 = k0-1;
	  kp1 = k0+1;
	  t1 = -(r[km1] -2.0*r[k0]+r[kp1])*_d1;

	  km1 = k0-emf->nx;
	  kp1 = k0+emf->nx;
	  t2 = -(r[km1] -2.0*r[k0]+r[kp1])*_d2;

	  km1 = k0-emf->nx*emf->ny;
	  kp1 = k0+emf->nx*emf->ny;
	  t3 = -(r[km1] -2.0*r[k0]+r[kp1])*_d3;

	  Hv[k0] += (c_h*(t1+t2)+c_v*t3)*gamma1;
	}
      }
    }
  }

}
