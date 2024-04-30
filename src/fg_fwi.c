/* function and gradient evaluation for EM full waveform inversion (FWI)
 *--------------------------------------------------------------------
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *--------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include "acq.h"
#include "emf.h"
#include "fwi.h"
#include "constants.h"
#include "mpi_info.h"

acq_t *acq;
emf_t *emf;
fwi_t *fwi;

void read_data(acq_t *acq, emf_t *emf);
void write_rmse(acq_t *acq, emf_t *emf);
void write_data(acq_t *acq, emf_t *emf, char *fname, _Complex float ***dcal_fd);

float cal_uncertainty_noise(acq_t *acq, emf_t *emf);

void regrid_init(acq_t *acq, emf_t *emf, int ifreq);
void regrid_close(emf_t *emf);

void gmg_init(emf_t *emf, int ifreq);
void gmg_close();
void gmg_apply(int n, complex *b, complex *x);

void inject_source(acq_t *acq, emf_t *emf, complex *b, int ifreq);
void inject_adj_source(acq_t *acq, emf_t *emf, complex *x, complex *b, int ifreq);
void extract_emf_data(acq_t *acq, emf_t *emf, complex *x, complex *b, int ifreq);
void extract_electric_field(acq_t *acq, emf_t *emf, complex *x, int ifreq);
void extract_magnetic_field(acq_t *acq, emf_t *emf, complex *b, int ifreq);

void build_gradient(acq_t *acq, emf_t *emf, fwi_t *fwi,
		    _Complex float ****E1f, _Complex float ****E2f, _Complex float ****E3f,
		    _Complex float ****E1a, _Complex float ****E2a, _Complex float ****E3a);
void fg_mod_reg(acq_t *acq, emf_t *emf, fwi_t *fwi, float *x, float *g);

/*--------------------------------------------------------*/
void fg_fwi_init(acq_t *acq_, emf_t *emf_, fwi_t *fwi_)
{
  int istat, i1, i2, n;
  float tmp;
  char *fbathy;
  FILE *fp;
  char fname[sizeof("delta_0000.txt")];
  
  acq = acq_;
  emf = emf_;
  fwi = fwi_;

  if(!getparstring("fbathy", &fbathy)) err("Need fbathy= for CSEM inversion");
  if(!getparint("addnoise", &emf->addnoise)) emf->addnoise = 1;//1=addnoise;0=not
  if(!getparfloat("amp_perc", &emf->amp_perc)) emf->amp_perc = 0.03;//amplitude uncertainty
  if(!getparfloat("delta_phi", &emf->delta_phi)) emf->delta_phi = 0;//2.*PI/180.;//azimuth uncertainty
  if(!getparfloat("delta_theta", &emf->delta_theta)) emf->delta_theta = 0;//dip uncertainty
  if(!getparfloat("noisefloorE", &emf->noisefloorE)) emf->noisefloorE = 1e-15;
  if(!getparfloat("noisefloorH", &emf->noisefloorH)) emf->noisefloorH = 1e-13;
  if(!getparfloat("dsmute", &emf->dsmute)) emf->dsmute = 800;//radius to mute gradient
  if(!getparfloat("drmute", &emf->drmute)) emf->drmute = 200;//radius to mute gradient
  
  if(!getparfloat("gamma1", &fwi->gamma1)) fwi->gamma1 = 10.;//Tikhonov
  if(!getparfloat("gamma2", &fwi->gamma2)) fwi->gamma2 = 0.;//TV
  if(!getparint("r1", &fwi->r1)) fwi->r1 = 3;
  if(!getparint("r2", &fwi->r2)) fwi->r2 = 3;
  if(!getparint("r3", &fwi->r3)) fwi->r3 = 1;
  if(!getparint("repeat", &fwi->repeat)) fwi->repeat = 1;
  if(emf->verb) {
    printf("preco=%d (0=off; 1=Shin's preconditioner; 2=depth weighting; 3=triangle smoothing)\n", fwi->preco);
    printf("Amplitude uncertainty: %g\%%\n", emf->amp_perc*100.);
    printf("noisefloorE=%g, noisefloorH=%g\n", emf->noisefloorE, emf->noisefloorH);
    printf("Tikhonov regularization: gamma1=%g\n", fwi->gamma1);
    printf("TV regularization:       gamma2=%g\n", fwi->gamma2);
    printf("----------- fwi init ------------\n");
  }
  
  emf->bathy = alloc2float(emf->nx, emf->ny);
  emf->dobs_fd = alloc3complexf(acq->nrec, emf->nfreq, emf->nchrec);
  emf->dcal_fd = alloc3complexf(acq->nrec, emf->nfreq, emf->nchrec);
  emf->dres_fd = alloc3complexf(acq->nrec, emf->nfreq, emf->nchrec);
  emf->delta_emf = alloc3float(acq->nrec, emf->nfreq, emf->nchrec);
  emf->rmse = alloc3float(acq->nrec, emf->nfreq, emf->nchrec);
    
  //step 1: read bathymetry file as a function depth=z(x.y)
  fp = fopen(fbathy, "rb");
  if(fp==NULL) err("cannot open file fbathy=%s\n", fbathy);
  istat = fread(&emf->bathy[0][0], sizeof(float), emf->nx*emf->ny, fp);
  if(istat != emf->nx*emf->ny) err("size parameter does not match the file-bathy!");
  fclose(fp);

  emf->depthmax = emf->bathy[0][0];
  for(i2=0; i2<emf->ny; i2++){
    for(i1=0; i1<emf->nx; i1++){
      emf->depthmax = MAX(emf->bathy[i2][i1], emf->depthmax);
    }
  }
  if(emf->verb) printf("maximum water depth: depthmax=%g\n", emf->depthmax);
  if(emf->depthmax < acq->src_x3[0]) err("coordinate error: source placed below seabed!");

  //step 2: read observed EM data, add certain amount of Gaussian white noise
  read_data(acq, emf);/* initialize dobs_fd */
  tmp = cal_uncertainty_noise(acq, emf);
  MPI_Allreduce(&tmp, &fwi->mse, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);/*sum over all process*/
  if(emf->verb) {
    n = nproc*emf->nchrec*emf->nfreq*acq->nrec;
    printf("addnoise=%d, target rmse=%e\n", emf->addnoise, sqrt(fwi->mse/n));
  }

  sprintf(fname, "obs_%04d.txt", acq->shot_idx[iproc]);
  write_data(acq, emf, fname, emf->dobs_fd);

  fwi->xref = alloc1float(fwi->n);
  fwi->grad = alloc4float(emf->nx, emf->ny, emf->nz, fwi->npar);
  if(fwi->preco==1) fwi->phess = alloc4float(emf->nx, emf->ny, emf->nz, fwi->npar);
  fwi->iter = 0;//initialize iteration counter
  fwi->firstgrad = 1;//preset flag to distinguish first gradient evaluation
  fwi->alpha = 1;//initialize the scaling factor
}

void fg_fwi_close()
{
  free2float(emf->bathy);
  free3complexf(emf->dobs_fd);
  free3complexf(emf->dcal_fd);
  free3complexf(emf->dres_fd);
  free3float(emf->delta_emf);
  free3float(emf->rmse);
    
  free1float(fwi->xref);
  free4float(fwi->grad);
  if(fwi->preco==1) free4float(fwi->phess);
}

/*< misfit function and gradient evaluation of EM FWI >*/
float fg_fwi(float *xx, float *g)
{
  float fcost, s, s1, s2;
  int i1, i2, i3, j, k, n;
  int ifreq, ichrec, irec;
  _Complex float ****E1f, ****E2f, ****E3f;
  _Complex float ****E1a, ****E2a, ****E3a;
  complex *x, *b;
  FILE *fp;
  char fname[sizeof("rmse_0000.txt")];
  
  //convert from logorithmic parameters to physical parameters
  emf->rhomin = exp(xx[0]);
  emf->rhomax = exp(xx[0]);
  for(j=0; j<fwi->npar; j++){
    for(i3=0; i3<emf->nz; i3++){
      for(i2=0; i2<emf->ny; i2++){
	for(i1=0; i1<emf->nx; i1++){
	  k = i1 + emf->nx*(i2 + emf->ny*(i3 + emf->nz*j));
	  s = exp(xx[k]);//from log(R) to R
	  if(fwi->idxpar[j]==1){
	    emf->rho11[i3][i2][i1] = s;
	    emf->rho22[i3][i2][i1] = s;
	  }
	  if(fwi->idxpar[j]==2) emf->rho33[i3][i2][i1] = s;	  
	  emf->rhomin = MIN(emf->rhomin, s);
	  emf->rhomax = MAX(emf->rhomax, s);
	}
      }
    }
  }
  if(emf->verb) printf("fg_fwi iter=%d: [rhomin,rhomax]=[%g,%g]\n", fwi->iter, emf->rhomin, emf->rhomax);
  
  E1f = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
  E2f = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
  E3f = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
  E1a = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
  E2a = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
  E3a = alloc4complexf(emf->nx, emf->ny, emf->nz, emf->nfreq);
  
  if(emf->verb) printf("--------step 1: forward modelling------------\n");
  emf->E1 = E1f;
  emf->E2 = E2f;
  emf->E3 = E3f;
  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
    regrid_init(acq, emf, ifreq);
    gmg_init(emf, ifreq);

    n = 3*emf->n123pad;
    x = alloc1complex(n);//vector E=(Ex,Ey,Ez)^T
    b = alloc1complex(n);//vector Js=(Jx,Jy,Jz)^T
    inject_source(acq, emf, b, ifreq);//initialize b=i*omega*mu*Js
    memset(x, 0, n*sizeof(complex));//initialize x=0
    
    gmg_apply(n, b, x);//after multigrid convergence, x=(Ex,Ey,Ez)^T, b=(Hx,Hy,Hz)^T
    
    extract_emf_data(acq, emf, x, b, ifreq);
    extract_electric_field(acq, emf, x, ifreq);

    free1complex(x);
    free1complex(b);
    
    gmg_close();
    regrid_close(emf);
  }
  sprintf(fname, "syn_%04d.txt", acq->shot_idx[iproc]);
  write_data(acq, emf, fname, emf->dcal_fd);
  
  if(emf->verb) printf("------step 2: compute fcost + adjoint source---------\n");
  memset(&emf->rmse[0][0][0], 0, acq->nrec*emf->nfreq*emf->nchrec*sizeof(float));
  fcost = 0.;
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      s = 5.5*sqrt(2.*emf->rho_water/(emf->omegas[ifreq]*mu0));//5 skin depth
      for(irec=0; irec<acq->nrec; irec++){
	s1 = acq->rec_x1[irec] - acq->src_x1[0];
	s2 = acq->rec_x2[irec] - acq->src_x2[0];
	if(sqrt(s1*s1 + s2*s2)>s){//only consider data > 5 skin depth
	  emf->dres_fd[ichrec][ifreq][irec] = emf->dobs_fd[ichrec][ifreq][irec] - emf->dcal_fd[ichrec][ifreq][irec];
	  emf->dres_fd[ichrec][ifreq][irec] /= emf->delta_emf[ichrec][ifreq][irec];// W(dobs-Ru)
	  emf->rmse[ichrec][ifreq][irec] = cabs(emf->dres_fd[ichrec][ifreq][irec]);
	  fcost += emf->dres_fd[ichrec][ifreq][irec]*conj(emf->dres_fd[ichrec][ifreq][irec]);
	  emf->dres_fd[ichrec][ifreq][irec] /= emf->delta_emf[ichrec][ifreq][irec];// W^T W(dobs-Ru)
	  emf->dres_fd[ichrec][ifreq][irec] = conj(emf->dres_fd[ichrec][ifreq][irec]);//take conjugate	  
	}//end if
      }/* end for irec */
    }//end for ifreq
  }// end for ichrec 
  fcost *= 0.5;
  write_rmse(acq, emf);//write down RMSE for each source
  n = emf->nchrec*emf->nfreq*acq->nrec;//number of data points for currrent source/processor
  printf("isrc=%d, fcost=%g, rmse=%g\n", acq->shot_idx[iproc], fcost, sqrt(2.*fcost/n));
  ierr = MPI_Allreduce(&fcost, &fwi->fcost_dat, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);/*sum over all process*/
  n *= nproc;//collect number of data points from all sources
  fwi->mse = sqrt(2.*fwi->fcost_dat/n);
  if(emf->verb) printf("fcost_dat=%g, rmse=%g\n", fwi->fcost_dat, fwi->mse);
  
  if(emf->verb) printf("-----------step3: adjoint modelling------------\n");
  emf->E1 = E1a;
  emf->E2 = E2a;
  emf->E3 = E3a;
  for(ifreq=0; ifreq<emf->nfreq; ifreq++){
    regrid_init(acq, emf, ifreq);
    gmg_init(emf, ifreq);

    n = 3*emf->n123pad;
    x = alloc1complex(n);//vector E=(Ex,Ey,Ez)^T
    b = alloc1complex(n);//vector Js=(Jx,Jy,Jz)^T, Js=data residual
    inject_adj_source(acq, emf, x, b, ifreq);//initialize b=i*omega*mu*Js,
    memset(x, 0, n*sizeof(complex));//initialize x=0

    gmg_apply(n, b, x);//after multigrid convergence, x=(Ex,Ey,Ez)^T, b=(Hx,Hy,Hz)^T
    
    extract_emf_data(acq, emf, x, b, ifreq);
    extract_electric_field(acq, emf, x, ifreq);

    free1complex(x);
    free1complex(b);
    
    gmg_close();
    regrid_close(emf);
  }

  if(emf->verb) printf("----------step 4: compute gradient --------------\n");
  build_gradient(acq, emf, fwi, E1f, E2f, E3f, E1a, E2a, E3a);
  if(fwi->preco==1){
    ierr = MPI_Allreduce(&fwi->phess[0][0][0][0], g, fwi->n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);//sum over different processors
    memcpy(&fwi->phess[0][0][0][0], g, fwi->n*sizeof(float));
  }
  ierr = MPI_Allreduce(&fwi->grad[0][0][0][0], g, fwi->n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);//sum over different processors

  free4complexf(E1f);
  free4complexf(E2f);
  free4complexf(E3f);
  free4complexf(E1a);
  free4complexf(E2a);
  free4complexf(E3a);

  if(emf->verb){
    printf("-------- stage 5: model regularization ----------\n");
    n = emf->nx*emf->ny*emf->nz;
    fp = fopen("gradient_dat_Rv", "wb");
    fwrite(g, n*sizeof(float), 1, fp);
    fclose(fp);
    fp = fopen("gradient_dat_Rh", "wb");
    fwrite(g+n, n*sizeof(float), 1, fp);
    fclose(fp);
  }
  fg_mod_reg(acq, emf, fwi, xx, g);
  fwi->fcost = fwi->fcost_dat + fwi->fcost_mod;
  if(emf->verb){
    printf("fcost_dat=%g\n", fwi->fcost_dat);
    printf("fcost_mod=%g\n", fwi->fcost_mod);
    printf("fcost=%g\n", fwi->fcost);
  }

  /* l-BFGS algorithm is not scale invariant, we need to estimate a scaling factor
     alpha at 1st iteration and use it for FWI objective function in later iterations */
  if(emf->mode == 1 && fwi->firstgrad){
    s1 = 0;
    s2 = 0;
    for(k = 0; k<fwi->n; k++){
      s1 += fabs(xx[k]);
      s2 += fabs(g[k]);
    }
    fwi->alpha = 2e-3*s1/s2;
    if(emf->verb) printf("|x|_1=%g, |g|_1=%g, scaling=%g\n", s1, s2, fwi->alpha);
  }
  for(k=0; k<fwi->n; k++) g[k] *= fwi->alpha;
  if(fwi->firstgrad) fwi->firstgrad = 0;

  return fwi->fcost*fwi->alpha; 
}


/*----------------------------------------------------------*/
void Hv_fwi(float *x, float *v, float *Hv)
/*< compute Gauss-Newton Hessian vector product for EM FWI >*/
{

}

