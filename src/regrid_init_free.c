/* regridding input parameters from uniform grid to non-uniform grid to
 * obtain effective medium by homogenization 
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
 
int find_index(int n, float *x, float val);
float create_nugrid(int n, float len, float dx, float *x);
void homogenization(float ***in, float ***out,
		    int nx, int ny, int nz, float *xx, float *yy, float *zz,
		    int n1, int n2, int n3, float *x1, float *x2, float *x3);

/*< generate staggered grid coordinate >*/
void generate_staggered_xs_dx(int n1, float *x1, float *x1s, float *d1, float *d1s)
{
  int i;
  
  for(i=0; i<n1; i++) {
    d1s[i] = x1[i+1] - x1[i];
    x1s[i] = 0.5*(x1[i] + x1[i+1]);
  }
  x1s[n1] = x1[n1] + 0.5*d1s[n1-1];
  d1s[n1] = x1s[n1]-x1[n1];

  for(i=1; i<=n1; i++) d1[i] = x1s[i] - x1s[i-1];
  d1[0] = x1s[0] - x1[0];
}

/*< compute effective medium by regridding over non-uniform grid >*/
void regrid_init(acq_t *acq, emf_t *emf)
{
  int i1, i2, i3;
  int j1, j2, j3;
  int i, j, k;
  float r, dist, *x;
  int nxpad, nypad, nzpad, np1, np2;
  float ***sig11, ***sig22, ***rho33;

  if(iproc==0) printf("------- regrid_init --------\n");
  emf->xx = alloc1float(emf->nx + 2*emf->nb + 1);
  emf->yy = alloc1float(emf->ny + 2*emf->nb + 1);
  emf->zz = alloc1float(emf->nz + 2*emf->nb + 1);
  for(i1=0; i1<=emf->nx; i1++) emf->xx[i1+emf->nb] = emf->ox + i1*emf->dx;
  for(i2=0; i2<=emf->ny; i2++) emf->yy[i2+emf->nb] = emf->oy + i2*emf->dy;
  for(i3=0; i3<=emf->nz; i3++) emf->zz[i3+emf->nb] = emf->oz + i3*emf->dz;

  //outer nonuniform grid exterior to uniform grid
  x = alloc1float(emf->nb+1);
  r = create_nugrid(emf->nb, emf->lextend, emf->dx, x);
  for(i1=0; i1<=emf->nb; i1++) {
    emf->xx[i1] = emf->ox - x[emf->nb-i1];
    emf->xx[i1+emf->nx+emf->nb] = emf->ox + emf->nx*emf->dx + x[i1];
  }
  r = create_nugrid(emf->nb, emf->lextend, emf->dy, x);
  for(i2=0; i2<=emf->nb; i2++) {
    emf->yy[i2] = emf->oy - x[emf->nb-i2];
    emf->yy[i2+emf->ny+emf->nb] = emf->oy + emf->ny*emf->dy + x[i2];
  }
  r = create_nugrid(emf->nb, emf->lextend, emf->dz, x);
  for(i3=0; i3<=emf->nb; i3++) {
    emf->zz[i3] = emf->oz - x[emf->nb-i3];
    emf->zz[i3+emf->nz+emf->nb] = emf->oz + emf->nz*emf->dz + x[i3];
  }
  free(x);

  emf->x1min = emf->xx[0];
  emf->x1max = emf->xx[emf->nx + 2*emf->nb];
  emf->x2min = emf->yy[0];
  emf->x2max = emf->yy[emf->ny + 2*emf->nb];
  emf->x3min = emf->zz[0];
  emf->x3max = emf->zz[emf->nz + 2*emf->nb];
  if(emf->verb){
    printf("[x1min, x1max]=[%g,%g]\n", emf->x1min, emf->x1max);
    printf("[x2min, x2max]=[%g,%g]\n", emf->x2min, emf->x2max);
    printf("[x3min, x3max]=[%g,%g]\n", emf->x3min, emf->x3max);
  }

  //extend the input model
  nxpad = emf->nx + 2*emf->nb;
  nypad = emf->ny + 2*emf->nb;
  nzpad = emf->nz + 2*emf->nb;
  sig11 = alloc3float(nxpad, nypad, nzpad);
  sig22 = alloc3float(nxpad, nypad, nzpad);
  rho33 = alloc3float(nxpad, nypad, nzpad);
  //fill interior part
  for(i3=0; i3<emf->nz; i3++){
    j3 = i3 + emf->nb;
    for(i2=0; i2<emf->ny; i2++){
      j2 = i2 + emf->nb;
      for(i1=0; i1<emf->nx; i1++){
	j1 = i1 + emf->nb;
	sig11[j3][j2][j1] = 1./emf->rho11[i3][i2][i1];
	sig22[j3][j2][j1] = 1./emf->rho22[i3][i2][i1];
	rho33[j3][j2][j1] = emf->rho33[i3][i2][i1];
      }
    }
  }
  //we then extend model on each side
  for(k=0; k<nzpad; k++){
    for(j=0; j<nypad; j++){
      for(i=0; i<emf->nb; i++){
	sig11[k][j][i] = sig11[k][j][emf->nb];
	sig11[k][j][nxpad-1-i] = sig11[k][j][nxpad-1-emf->nb];
	sig22[k][j][i] = sig22[k][j][emf->nb];
	sig22[k][j][nxpad-1-i] = sig22[k][j][nxpad-1-emf->nb];
	rho33[k][j][i] = rho33[k][j][emf->nb];
	rho33[k][j][nxpad-1-i] = rho33[k][j][nxpad-1-emf->nb];
      }
    }
  }
  for(k=0; k<nzpad; k++){
    for(j=0; j<emf->nb; j++){
      for(i=0; i<nxpad; i++){
	sig11[k][j][i] = sig11[k][emf->nb][i];
	sig11[k][nypad-1-j][i] = sig11[k][nypad-1-emf->nb][i];
	sig22[k][j][i] = sig22[k][emf->nb][i];
	sig22[k][nypad-1-j][i] = sig22[k][nypad-1-emf->nb][i];
	rho33[k][j][i] = rho33[k][emf->nb][i];
	rho33[k][nypad-1-j][i] = rho33[k][nypad-1-emf->nb][i];
      }
    }
  }
  for(k=0; k<emf->nb; k++){
    for(j=0; j<nypad; j++){
      for(i=0; i<nxpad; i++){
	sig11[k][j][i] = 1./emf->rho_air;//fill in with air
	sig11[nzpad-1-k][j][i] = sig11[nzpad-1-emf->nb][k][i];
	sig22[k][j][i] = 1./emf->rho_air;//fill in with air
	sig22[nzpad-1-k][j][i] = sig22[nzpad-1-emf->nb][k][i];
	rho33[k][j][i] = emf->rho_air;//fill in with air
	rho33[nzpad-1-k][j][i] = rho33[nzpad-1-emf->nb][k][i];
      }
    }
  }
    
  emf->x1 = alloc1float(emf->n1+1);
  emf->x2 = alloc1float(emf->n2+1);
  emf->x3 = alloc1float(emf->n3+1);
  emf->x1s = alloc1float(emf->n1+1);
  emf->x2s = alloc1float(emf->n2+1);
  emf->x3s = alloc1float(emf->n3+1);
  emf->d1 = alloc1float(emf->n1+1);
  emf->d2 = alloc1float(emf->n2+1);
  emf->d3 = alloc1float(emf->n3+1);
  emf->d1s = alloc1float(emf->n1+1);
  emf->d2s = alloc1float(emf->n2+1);
  emf->d3s = alloc1float(emf->n3+1);

  //power-law grid stretching
  //-----------------------------
  np1 = emf->n1/2;
  np2 = emf->n1 - np1;

  x = alloc1float(np1+1);
  dist = acq->src_x1[0] - emf->x1min;
  r = create_nugrid(np1, dist, emf->d1min, x);
  for(i1=0; i1<=np1; i1++) emf->x1[np1-i1] = acq->src_x1[0] - x[i1];
  free1float(x);
  if(emf->verb) printf("left stretching r=%g\n", r);
  
  x = alloc1float(np2+1);
  dist = emf->x1max - acq->src_x1[0];
  r = create_nugrid(np2, dist, emf->d1min, x);
  for(i1=0; i1<=np2; i1++) emf->x1[np1+i1] = acq->src_x1[0] + x[i1];
  free1float(x);
  if(emf->verb) printf("right stretching r=%g\n", r);

  //-----------------------------
  np1 = emf->n2/2;
  np2 = emf->n2 - np1;

  x = alloc1float(np1+1);
  dist = acq->src_x2[0] - emf->x2min;
  r = create_nugrid(np1, dist, emf->d2min, x);
  for(i2=0; i2<=np1; i2++) emf->x2[np1-i2] = acq->src_x2[0] - x[i2];
  free1float(x);
  if(emf->verb) printf("front stretching r=%g\n", r);

  x = alloc1float(np2+1);
  dist = emf->x2max - acq->src_x2[0];
  r = create_nugrid(np2, dist, emf->d2min, x);
  for(i2=0; i2<=np2; i2++) emf->x2[np1+i2] = acq->src_x2[0] + x[i2];
  free1float(x);
  if(emf->verb) printf("rear stretching r=%g\n", r);
  
  //---------------------------------------------
  np1 = find_index(emf->nx + 2*emf->nb, emf->zz, acq->src_x3[0]) + 5;
  np2 = emf->n3 - np1;

  for(i3=0; i3<=np1; i3++) emf->x3[i3] = emf->zz[i3];
  if(emf->verb) printf("top above the source: use %d grid from input model\n", np1);
    
  x = alloc1float(np2+1);
  dist = emf->x3max - emf->x3[np1];
  r = create_nugrid(np2, dist, emf->d3min, x);
  for(i3=0; i3<=np2; i3++) emf->x3[np1+i3] = emf->x3[np1] + x[i3];
  free1float(x);
  if(emf->verb) printf("bottom stretching r=%g\n", r);
  

  //this will be used for extraction of EM data at receiver locations
  generate_staggered_xs_dx(emf->n1, emf->x1, emf->x1s, emf->d1, emf->d1s);
  generate_staggered_xs_dx(emf->n2, emf->x2, emf->x2s, emf->d2, emf->d2s);
  generate_staggered_xs_dx(emf->n3, emf->x3, emf->x3s, emf->d3, emf->d3s);

  emf->sigma11 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->sigma22 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->sigma33 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->invmur = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->vol = alloc3float(emf->n1, emf->n2, emf->n3);

  homogenization(sig11, emf->sigma11, nxpad, nypad, nzpad, emf->xx, emf->yy, emf->zz, emf->n1, emf->n2, emf->n3, emf->x1, emf->x2, emf->x3);
  homogenization(sig22, emf->sigma22, nxpad, nypad, nzpad, emf->xx, emf->yy, emf->zz, emf->n1, emf->n2, emf->n3, emf->x1, emf->x2, emf->x3);
  homogenization(rho33, emf->sigma33, nxpad, nypad, nzpad, emf->xx, emf->yy, emf->zz, emf->n1, emf->n2, emf->n3, emf->x1, emf->x2, emf->x3);
  for(i3=0; i3<emf->n3; i3++){
    for(i2=0; i2<emf->n2; i2++){
      for(i1=0; i1<emf->n1; i1++){
	emf->sigma33[i3][i2][i1] = 1./emf->sigma33[i3][i2][i1];//convert from rho to sigma
	emf->invmur[i3][i2][i1] = 1.;
	emf->vol[i3][i2][i1] = emf->d1s[i1]*emf->d2s[i2]*emf->d3s[i3];//volume assigned at cell center
      }
    }
  }
  
}

/*< free variables in emf >*/
void regrid_free(emf_t *emf)
{
  free1float(emf->xx);
  free1float(emf->yy);
  free1float(emf->zz);
  
  free1float(emf->x1);
  free1float(emf->x2);
  free1float(emf->x3);

  free1float(emf->x1s);
  free1float(emf->x2s);
  free1float(emf->x3s);

  free1float(emf->d1);
  free1float(emf->d2);
  free1float(emf->d3);

  free1float(emf->d1s);
  free1float(emf->d2s);
  free1float(emf->d3s);

  free3float(emf->sigma11);
  free3float(emf->sigma22);
  free3float(emf->sigma33);
  free3float(emf->invmur);
  free3float(emf->vol);
}


