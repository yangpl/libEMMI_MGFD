/* preconditioner by cascading depth preconditioning with triangle smoothing
 *----------------------------------------------------------------------
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include "fwi.h"
 

void triangle_smoothing(float ***mod, int n1, int n2, int n3, int r1, int r2, int r3, int repeat);

/*------------------------------------------------------*/
void precondition(emf_t *emf, fwi_t *fwi, float *g)
/*< precondition vector x of size n >*/
{
  int i1, i2, i3, k, i;
  float skin_depth, z, s1, s2;

  if(fwi->preco==1){
    s1 = 0.;
    s2 = 0.;
    for(k=0; k<fwi->npar; k++){
      z = 0;
      for(i3=0; i3<emf->nz; i3++){
	for(i2=0; i2<emf->ny; i2++){
	  for(i1=0; i1<emf->nx; i1++){
	    //find the maximum value of pseudo-Hessian
	    if(emf->oz + (i3+0.5)*emf->dz>emf->bathy[i2][i1]) z = MAX(z, fwi->phess[k][i3][i2][i1]);
	  }
	}
      }
      
      for(i3=0; i3<emf->nz; i3++){
	for(i2=0; i2<emf->ny; i2++){
	  for(i1=0; i1<emf->nx; i1++){
	    i = i1 + emf->nx*(i2 + emf->ny*(i3 + emf->nz*k));

	    s1 += g[i]*g[i];
	    g[i] /= (fwi->phess[k][i3][i2][i1] + 1e-5*z);//regularize to avoid division by 0
	    s2 += g[i]*g[i];
	  }
	}
      }
    }//end for k

    z = sqrt(s1/s2);
    for(k=0; k<fwi->npar; k++){
      for(i3=0; i3<emf->nz; i3++){
	for(i2=0; i2<emf->ny; i2++){
	  for(i1=0; i1<emf->nx; i1++){
	    i = i1 + emf->nx*(i2 + emf->ny*(i3 + emf->nz*k));
	    g[i] *= z;//rescale such that the scaling factor works
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for k
  }else if(fwi->preco==2){
    //depth preconditioning, see Plessix 2008 Inverse Problems
    skin_depth = sqrt(2./(emf->omegas[0]*mu0));//consider rho_ref=1 Ohm-m
    for(k=0; k<fwi->npar; k++){
      for(i3=0; i3<emf->nz; i3++){
	for(i2=0; i2<emf->ny; i2++){
	  for(i1=0; i1<emf->nx; i1++){
	    i = i1 + emf->nx*(i2 + emf->ny*(i3 + emf->nz*k) );
	    z = emf->oz + (i3+0.5)*emf->dz;
	    if(z>emf->bathy[i2][i1]){
	      z -= emf->bathy[i2][i1];
	      g[i] *= exp(z/skin_depth);
	    }
	  }//end for i1
	}//end for i2
      }//end for i3
    }//end for k

  }else if(fwi->preco==3){
    /* Reference: Harlan 1995, regularization by model reparameterization,
     * We should always use smoothing to achieve faster convergence **/
    memcpy(&fwi->grad[0][0][0][0], g, fwi->n*sizeof(float));
    for(k=0; k<fwi->npar; k++)
      triangle_smoothing(fwi->grad[k], emf->nx, emf->ny, emf->nz, fwi->r1, fwi->r2, fwi->r3, fwi->repeat);
    memcpy(g, &fwi->grad[0][0][0][0], fwi->n*sizeof(float));
  }
  
}
