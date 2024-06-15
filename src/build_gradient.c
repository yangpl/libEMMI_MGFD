/* build inversion gradient for frequency domain CSEM inversion
 *----------------------------------------------------------------------
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *----------------------------------------------------------------------*/
#include "cstd.h"
#include "acq.h"
#include "emf.h"
#include "fwi.h"

void build_gradient(acq_t *acq, emf_t *emf, fwi_t *fwi,
		    _Complex float ****E1f, _Complex float ****E2f, _Complex float ****E3f,
		    _Complex float ****E1a, _Complex float ****E2a, _Complex float ****E3a)
{
  int i1, i2, i3, k, isrc, irec, ifreq;
  int i1min, i1max, i2min, i2max, i3min, i3max;
  float s;

  s = emf->dx*emf->dy*emf->dz;//cell volume
  for(k=0; k<fwi->npar; k++){
    for(i3=0; i3<emf->nz; i3++){
      for(i2=0; i2<emf->ny; i2++){
	for(i1=0; i1<emf->nx; i1++){
	  fwi->grad[k][i3][i2][i1] = 0;
	  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
	    if(fwi->idxpar[k]==1) {
	      fwi->grad[k][i3][i2][i1] += creal(E1f[ifreq][i3][i2][i1]*E1a[ifreq][i3][i2][i1]);
	      fwi->grad[k][i3][i2][i1] += creal(E2f[ifreq][i3][i2][i1]*E2a[ifreq][i3][i2][i1]);
	    }
	    if(fwi->idxpar[k]==2){
	      fwi->grad[k][i3][i2][i1] += creal(E3f[ifreq][i3][i2][i1]*E3a[ifreq][i3][i2][i1]);
	    }
	  }//end for ifreq
	  //convert from dJ/d(sigma) to dJ/dln(rho)
	  if(fwi->idxpar[k]==1) fwi->grad[k][i3][i2][i1] *= s/emf->rho11[i3][i2][i1];
	  if(fwi->idxpar[k]==2) fwi->grad[k][i3][i2][i1] *= s/emf->rho33[i3][i2][i1];

	  //mute gradient above bathymetry 
	  if(emf->oz + (i3+0.5)*emf->dz<= emf->bathy[i2][i1]+emf->dz) fwi->grad[k][i3][i2][i1] = 0.;

	  if(fwi->preco==1){
	    fwi->phess[k][i3][i2][i1] = 0.;
	    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
	      if(fwi->idxpar[k]==1) {
		fwi->phess[k][i3][i2][i1] += cabs(E1f[ifreq][i3][i2][i1])*cabs(E1f[ifreq][i3][i2][i1]);
		fwi->phess[k][i3][i2][i1] += cabs(E2f[ifreq][i3][i2][i1])*cabs(E2f[ifreq][i3][i2][i1]);
	      }
	      if(fwi->idxpar[k]==2){
		fwi->phess[k][i3][i2][i1] += cabs(E3f[ifreq][i3][i2][i1])*cabs(E3f[ifreq][i3][i2][i1]);
	      }
	    }//end for ifreq
	    //convert from dJ/d(sigma) to dJ/dln(rho)
	    if(fwi->idxpar[k]==1) fwi->phess[k][i3][i2][i1] *= s/emf->rho11[i3][i2][i1];
	    if(fwi->idxpar[k]==2) fwi->phess[k][i3][i2][i1] *= s/emf->rho33[i3][i2][i1];
	  }
	}//end for i1
      }//end for i2
    }//end for i3
  }//end for k

  //mute points in the vicinity of source, by default nsrc=1
  for(isrc=0; isrc<acq->nsrc; isrc++){
    i1min = (acq->src_x1[isrc]-emf->dsmute-emf->ox)/emf->dx;
    i1max = (acq->src_x1[isrc]+emf->dsmute-emf->ox)/emf->dx;
    i2min = (acq->src_x2[isrc]-emf->dsmute-emf->oy)/emf->dy;
    i2max = (acq->src_x2[isrc]+emf->dsmute-emf->oy)/emf->dy;
    i3min = (acq->src_x3[isrc]-emf->dsmute-emf->oz)/emf->dz;
    i3max = (acq->src_x3[isrc]+emf->dsmute-emf->oz)/emf->dz;
    i1min = MIN(MAX(i1min, 0), emf->nx-1);
    i1max = MIN(MAX(i1max, 0), emf->nx-1);
    i2min = MIN(MAX(i2min, 0), emf->ny-1);
    i2max = MIN(MAX(i2max, 0), emf->ny-1);
    i3min = MIN(MAX(i3min, 0), emf->nz-1);
    i3max = MIN(MAX(i3max, 0), emf->nz-1);
    for(k=0; k<fwi->npar; k++){
      for(i3=i3min; i3<=i3max; i3++){
	for(i2=i2min; i2<=i2max; i2++){
	  for(i1=i1min; i1<=i1max; i1++){
	    fwi->grad[k][i3][i2][i1] = 0.;
	  }
	}
      }
    }
  }//end for isrc
  
  //mute gradient in the vicinity of receivers
  for(irec=0; irec<acq->nrec; irec++){
    i1min = (acq->rec_x1[irec]-emf->drmute-emf->ox)/emf->dx;
    i1max = (acq->rec_x1[irec]+emf->drmute-emf->ox)/emf->dx;
    i2min = (acq->rec_x2[irec]-emf->drmute-emf->oy)/emf->dy;
    i2max = (acq->rec_x2[irec]+emf->drmute-emf->oy)/emf->dy;
    i3min = (acq->rec_x3[irec]-emf->drmute-emf->oz)/emf->dz;
    i3max = (acq->rec_x3[irec]+emf->drmute-emf->oz)/emf->dz;
    i1min = MIN(MAX(i1min, 0), emf->nx-1);
    i1max = MIN(MAX(i1max, 0), emf->nx-1);
    i2min = MIN(MAX(i2min, 0), emf->ny-1);
    i2max = MIN(MAX(i2max, 0), emf->ny-1);
    i3min = MIN(MAX(i3min, 0), emf->nz-1);
    i3max = MIN(MAX(i3max, 0), emf->nz-1);
    for(k=0; k<fwi->npar; k++){
      for(i3=i3min; i3<=i3max; i3++){
	for(i2=i2min; i2<=i2max; i2++){
	  for(i1=i1min; i1<=i1max; i1++){
	    fwi->grad[k][i3][i2][i1] = 0.;
	  }
	}
      }
    }
  }//end for irec

}
