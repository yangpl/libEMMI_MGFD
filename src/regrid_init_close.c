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
#include "constants.h"

float create_nugrid(int n, float len, float dx, float *x, int istretch);
void homogenization(emf_t *emf, float ***sigma_in, float ***sigma_out, int flag);

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
void regrid_init(acq_t *acq, emf_t *emf, int ifreq)
{
  int i1, i2, i3, nz;
  int j1, j2, j3;
  float vmin, vmax, dist, dleft, dright, r;
  float *x;
  int nb1, nb2, nb3, nb, logcond;
  float rhox, rhoy, rhoz, sigma;
  float lx, ly, lz, lair;
  
  if(!getparfloat("rhox", &rhox)) rhox = 5;//parameter to extend boundary in x direction
  if(!getparfloat("rhoy", &rhoy)) rhoy = 5;//parameter to extend boundary in y direction
  if(!getparfloat("rhoz", &rhoz)) rhoz = 10;//parameter to extend boundary in z direction
  if(!getparfloat("lair", &lair)) lair = 100000;//extend distance up to 100 km in air
  if(!getparint("logcond", &logcond)) logcond = 0;//1=average over log(conductivity)

  //extending distance=6*skin_depth, field becomes exp(-6.)=0.2%
  sigma = 1./rhox;
  lx = 6.*sqrt(2./(sigma*emf->omegas[ifreq]*mu0));
  sigma = 1./rhoy;
  ly = 6.*sqrt(2./(sigma*emf->omegas[ifreq]*mu0));
  sigma = 1./rhoz;
  lz = 6.*sqrt(2./(sigma*emf->omegas[ifreq]*mu0));
  if(emf->verb) {
    printf("---------------- emf init, freq=%g Hz ----------------\n", emf->freqs[ifreq]);
    printf("logcond=%d (1=average over log(cond), 0=not)\n", logcond);
    printf("parameters to extend domain: [rhox, rhoy, rhoz]=[%g, %g, %g]\n", rhox, rhoy, rhoz);
    printf("extended boundary: [lx, ly, lz, lair]=[%g, %g, %g, %g]\n", lx, ly, lz, lair);
  }
  
  emf->x1 = alloc1float(emf->n1pad);
  emf->x2 = alloc1float(emf->n2pad);
  emf->x3 = alloc1float(emf->n3pad);
  emf->x1s = alloc1float(emf->n1pad);
  emf->x2s = alloc1float(emf->n2pad);
  emf->x3s = alloc1float(emf->n3pad);
  emf->d1 = alloc1float(emf->n1pad);
  emf->d2 = alloc1float(emf->n2pad);
  emf->d3 = alloc1float(emf->n3pad);
  emf->d1s = alloc1float(emf->n1pad);
  emf->d2s = alloc1float(emf->n2pad);
  emf->d3s = alloc1float(emf->n3pad);

  nb1 = emf->nb_air;
  nb2 = emf->nb_air;
  nb3 = emf->nb_air;
  if(!emf->istretch){//uniform grid without stretching
    nb = (emf->n1 - emf->nx)/2;
    for(i1=0; i1<=emf->nx; i1++) emf->x1[nb+i1] = emf->ox + i1*emf->dx;

    x = alloc1float(nb+1);//outer nonuniform grid exterior to uniform grid
    r = create_nugrid(nb, lx, emf->dx, x, 1);
    for(i1=0; i1<=nb; i1++){
      emf->x1[nb+emf->nx+i1] = emf->x1[nb+emf->nx] + x[i1];
      emf->x1[nb-i1] = emf->x1[nb] - x[i1];
    }
    free1float(x);
    emf->x1min = emf->x1[0];
    emf->x1max = emf->x1[emf->n1];
    nb1 = nb;
        
    //------------------------------------------
    nb = (emf->n2 - emf->ny)/2;
    for(i2=0; i2<=emf->ny; i2++) emf->x2[nb+i2] = emf->oy + i2*emf->dy;

    x = alloc1float(nb+1);//outer nonuniform grid exterior to uniform grid
    r = create_nugrid(nb, ly, emf->dy, x, 1);
    for(i2=0; i2<=nb; i2++){
      emf->x2[nb-i2] = emf->x2[nb] - x[i2];
      emf->x2[nb+emf->ny+i2] = emf->x2[nb+emf->ny] + x[i2];
    }
    free1float(x);
    emf->x2min = emf->x2[0];
    emf->x2max = emf->x2[emf->n2];
    nb2 = nb;
    
    //------------------------------------------
    nb = (emf->n3 - emf->nz)/2;
    for(i3=0; i3<=emf->nz; i3++) emf->x3[nb+i3] = emf->oz + i3*emf->dz;

    x = alloc1float(nb+1);//outer nonuniform grid exterior to uniform grid
    r = create_nugrid(nb, lz, emf->dz, x, 1);
    for(i3=0; i3<=nb; i3++){
      emf->x3[nb+emf->nz+i3] = emf->x3[nb+emf->nz] + x[i3];
      emf->x3[nb-i3] = emf->x3[nb] - x[i3];
    }
    if(emf->addair){
      r = create_nugrid(nb, lair, emf->dz, x, 1);
      for(i3=0; i3<=nb; i3++) emf->x3[nb-i3] = emf->x3[nb] - x[i3];
    }
    free1float(x);
    emf->x3min = emf->x3[0];
    emf->x3max = emf->x3[emf->n3];
    nb3 = nb;

  }else{//power-law grid stretching
    //--------------------------------------------
    x = alloc1float(emf->n1/2+1);//outer nonuniform grid exterior to uniform grid
    dleft = acq->src_x1[0] - emf->ox;
    dright = emf->ox + emf->nx*emf->dx - acq->src_x1[0];
    dist = MAX(dleft, dright) + lx;
    r = create_nugrid(emf->n1/2, dist, emf->d1min, x, 1);
    emf->x1[emf->n1pad/2] = acq->src_x1[0];
    for(i1=1; i1<=emf->n1/2; i1++){
      emf->x1[emf->n1pad/2-i1] = emf->x1[emf->n1pad/2] - x[i1];
      emf->x1[emf->n1pad/2+i1] = emf->x1[emf->n1pad/2] + x[i1];
    }
    free1float(x);
    emf->x1min = emf->x1[0];
    emf->x1max = emf->x1[emf->n1pad-1];
    if(emf->verb) printf("extended domain [x1min, x1max]=[%g, %g], stretching r=%g\n", emf->x1min, emf->x1max, r);

    //--------------------------------------------
    x = alloc1float(emf->n2/2+1);//outer nonuniform grid exterior to uniform grid
    dleft = acq->src_x2[0] - emf->oy;
    dright = emf->oy + emf->ny*emf->dy - acq->src_x2[0];
    dist = MAX(dleft, dright) + ly;
    r = create_nugrid(emf->n2/2, dist, emf->d2min, x, 1);
    emf->x2[emf->n2pad/2] = acq->src_x2[0];
    for(i2=1; i2<=emf->n2/2; i2++){
      emf->x2[emf->n2pad/2-i2] = emf->x2[emf->n2pad/2] - x[i2];
      emf->x2[emf->n2pad/2+i2] = emf->x2[emf->n2pad/2] + x[i2];
    }
    free1float(x);
    emf->x2min = emf->x2[0];
    emf->x2max = emf->x2[emf->n2pad-1];
    if(emf->verb) printf("extended domain [x2min, x2max]=[%g, %g], stretching r=%g\n", emf->x2min, emf->x2max, r);
  
    //---------------------------------------------
    dist = emf->addair?lair:lz;//100 km
    emf->x3min = emf->oz - dist;
    emf->x3max = emf->oz + emf->dz*emf->nz + lz;
    if(emf->verb) printf("extended domain [x3min, x3max]=[%g, %g]\n", emf->x3min, emf->x3max);

    x = alloc1float(emf->nb_air+1);
    r = create_nugrid(emf->nb_air, dist, emf->dz, x, emf->istretch);
    if(emf->verb) printf("                x3 top stretching r=%g (above sea surface)\n", r);
    for(i3=0; i3<=emf->nb_air; i3++) emf->x3[emf->nb_air-i3] = emf->oz - x[i3];
    free1float(x);

    dist = emf->x3max-emf->oz;
    nz = emf->n3pad -1 - emf->nb_air;//remaining nz points with coordinates unassigned
    x = alloc1float(nz + 1);
    r = create_nugrid(nz, dist, emf->d3min, x, emf->istretch);
    if(emf->verb) printf("                x3 bottom stretching r=%g (below sea surface)\n", r);
    for(i3=0; i3<=nz; i3++) emf->x3[i3+emf->nb_air] = emf->x3[emf->nb_air] + x[i3];
    free1float(x);
  }

  //this will be used for extraction of EM data at receiver locations
  generate_staggered_xs_dx(emf->n1, emf->x1, emf->x1s, emf->d1, emf->d1s);
  generate_staggered_xs_dx(emf->n2, emf->x2, emf->x2s, emf->d2, emf->d2s);
  generate_staggered_xs_dx(emf->n3, emf->x3, emf->x3s, emf->d3, emf->d3s);

  //--------------------------------------------------------------
  vmax = emf->rho11[0][0][0];
  vmin = emf->rho11[0][0][0];
  for(i3=0; i3<emf->nz; i3++){
    for(i2=0; i2<emf->ny; i2++){
      for(i1=0; i1<emf->nx; i1++){
	vmax = MAX(vmax, emf->rho11[i3][i2][i1]);
	vmin = MIN(vmin, emf->rho11[i3][i2][i1]);
	vmax = MAX(vmax, emf->rho22[i3][i2][i1]);
	vmin = MIN(vmin, emf->rho22[i3][i2][i1]);
	vmax = MAX(vmax, emf->rho33[i3][i2][i1]);
	vmin = MIN(vmin, emf->rho33[i3][i2][i1]);
      }
    }
  }
  if(emf->verb) printf("unextended model [rhomin, rhomax]=[%g, %g]\n", vmin, vmax);

  emf->sigma11 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->sigma22 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->sigma33 = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->invmur = alloc3float(emf->n1, emf->n2, emf->n3);
  emf->vol = alloc3float(emf->n1, emf->n2, emf->n3);

  if(emf->istretch){//stretched grid requires medium homogenization
    //--------------------------------------------------------------
    if(!logcond){ //average horizontal conductivity and vertical resistivity
      homogenization(emf, emf->rho11, emf->sigma11, 1);
      homogenization(emf, emf->rho22, emf->sigma22, 1);
      homogenization(emf, emf->rho33, emf->sigma33, 2);
    }else{//average over log(conductivity)
      homogenization(emf, emf->rho11, emf->sigma11, 3);
      homogenization(emf, emf->rho22, emf->sigma22, 3);
      homogenization(emf, emf->rho33, emf->sigma33, 3);
    }
  }else{//uniform grid requires grid extension
    for(i3=0; i3<emf->nz; i3++){
      j3 = i3+nb3;
      for(i2=0; i2<emf->ny; i2++){
	j2 = i2+nb2;
	for(i1=0; i1<emf->nx; i1++){
	  j1 = i1+nb1;
	  emf->sigma11[j3][j2][j1] = 1./emf->rho11[i3][i2][i1];
	  emf->sigma22[j3][j2][j1] = 1./emf->rho22[i3][i2][i1];
	  emf->sigma33[j3][j2][j1] = 1./emf->rho33[i3][i2][i1];
	}
      }
    }
    for(i3=0; i3<emf->n3; i3++){
      for(i2=0; i2<emf->n2; i2++){
	for(i1=0; i1<nb1; i1++){
	  emf->sigma11[i3][i2][i1] = emf->sigma11[i3][i2][nb1];
	  emf->sigma22[i3][i2][i1] = emf->sigma22[i3][i2][nb1];
	  emf->sigma33[i3][i2][i1] = emf->sigma33[i3][i2][nb1];
	  j1 = emf->n1 - 1 - i1;
	  emf->sigma11[i3][i2][j1] = emf->sigma11[i3][i2][emf->n1-1-nb1];
	  emf->sigma22[i3][i2][j1] = emf->sigma22[i3][i2][emf->n1-1-nb1];
	  emf->sigma33[i3][i2][j1] = emf->sigma33[i3][i2][emf->n1-1-nb1];	  
	}
      }
    }
    for(i3=0; i3<emf->n3; i3++){
      for(i2=0; i2<nb2; i2++){
	j2 = emf->n2 - 1 - i2;
	for(i1=0; i1<emf->n1; i1++){
	  emf->sigma11[i3][i2][i1] = emf->sigma11[i3][nb2][i1];
	  emf->sigma22[i3][i2][i1] = emf->sigma22[i3][nb2][i1];
	  emf->sigma33[i3][i2][i1] = emf->sigma33[i3][nb2][i1];
	  
	  emf->sigma11[i3][j2][i1] = emf->sigma11[i3][emf->n2-1-nb2][i1];
	  emf->sigma22[i3][j2][i1] = emf->sigma22[i3][emf->n2-1-nb2][i1];
	  emf->sigma33[i3][j2][i1] = emf->sigma33[i3][emf->n2-1-nb2][i1];
	}
      }
    }

    for(i3=0; i3<nb3; i3++){
      j3 = emf->n3 - 1 - i3;
      for(i2=0; i2<emf->n2; i2++){
	for(i1=0; i1<emf->n1; i1++){
	  emf->sigma11[i3][i2][i1] = (emf->addair)?1./emf->rho_air:emf->sigma11[nb3][i2][i1];
	  emf->sigma22[i3][i2][i1] = (emf->addair)?1./emf->rho_air:emf->sigma22[nb3][i2][i1];
	  emf->sigma33[i3][i2][i1] = (emf->addair)?1./emf->rho_air:emf->sigma33[nb3][i2][i1];
	  
	  emf->sigma11[j3][i2][i1] = emf->sigma11[emf->n3-1-nb3][i2][i1];
	  emf->sigma22[j3][i2][i1] = emf->sigma22[emf->n3-1-nb3][i2][i1];
	  emf->sigma33[j3][i2][i1] = emf->sigma33[emf->n3-1-nb3][i2][i1];
	}
      }
    }
    
  }//end if

  //------------------------------------------------------------------
  vmax = emf->sigma11[0][0][0];
  vmin = emf->sigma11[0][0][0];
  for(i3=0; i3<emf->n3; i3++){
    for(i2=0; i2<emf->n2; i2++){
      for(i1=0; i1<emf->n1; i1++){
	vmax = MAX(vmax, emf->sigma11[i3][i2][i1]);
	vmin = MIN(vmin, emf->sigma11[i3][i2][i1]);
	vmax = MAX(vmax, emf->sigma22[i3][i2][i1]);
	vmin = MIN(vmin, emf->sigma22[i3][i2][i1]);
 	vmax = MAX(vmax, emf->sigma33[i3][i2][i1]);
	vmin = MIN(vmin, emf->sigma33[i3][i2][i1]);
	emf->invmur[i3][i2][i1] = 1.;
	emf->vol[i3][i2][i1] = emf->d1s[i1]*emf->d2s[i2]*emf->d3s[i3];//volume assigned at cell center
      }
    }
  }
  if(emf->verb) printf("extended model   [rhomin, rhomax]=[%g, %g]\n", 1./vmax, 1./vmin);

}

/*< free variables in emf >*/
void regrid_free(emf_t *emf)
{
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


