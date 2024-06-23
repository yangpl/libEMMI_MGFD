/* full waveform inversion of the CSEM data
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
#include "opt.h"
#include "fwi.h"

float l2norm(int n, float *a);
void flipsign(int n, float *a, float *b);
void lbfgs_save(int n, float *x, float *g, float **sk, float **yk, opt_t *opt);
void lbfgs_update(int n, float *x, float *g, float **sk, float **yk, opt_t *opt);
void lbfgs_descent(int n, float *g, float *d, float **sk, float **yk, float *q, float *rho, float *alp, opt_t *opt);
bool lbfgs_descent1(int n, float *g, float *q, float *rho, float *alp, float **sk, float **yk, opt_t *opt);
void lbfgs_descent2(int n, float *g, float *q, float *rho, float *alp, float **sk, float **yk, opt_t *opt);
void boundx(float *x, int n, float *xmin, float *xmax);
void line_search(int n, //dimension of x
		 float *x, //input vector x
		 float *g, //gradient of misfit function
		 float *d, //descent direction
		 opt_fg fg, //subroutine to evaluation function and gradient
		 opt_t *opt); //pointer of l-BFGS optimization parameters
void cg_solve(int n, //dimension of x
	      float *x, //input vector x
	      float *g, //gradient of misfit function
	      float *d, //descent direction
	      opt_Hv Hv, //subroutine to evaluation function and gradient
	      opt_t *opt); //pointer of l-BFGS optimization parameters

void fg_fwi_init(acq_t *acq_, emf_t *emf_, fwi_t *fwi_);
void fg_fwi_close();
float fg_fwi(float *x, float *g);
void Hv_fwi(float *x, float *v, float *Hv);
void precondition(emf_t *emf, fwi_t *fwi, float *g);

void do_fwi(acq_t *acqui, emf_t *emf)
{
  int i, j, i1, i2, i3;
  float fcost;
  opt_t *opt; /* pointer for opt_t parameters */
  fwi_t *fwi;
  FILE *fp;
  char fname[sizeof "Rv_iteration_00"];

  fwi = (fwi_t*)malloc(sizeof(fwi_t));
  opt = (opt_t*)malloc(sizeof(opt_t));
  
  /*-------------------------------------------------------------------------*/
  if(!getparint("niter", &opt->niter)) opt->niter = 30;/* maximum number of iterations */
  if(!getparint("nls", &opt->nls)) opt->nls = 5;/* maximum number of line searches */
  if(!getparfloat("tol", &opt->tol)) opt->tol = 1e-6;/* convergence tolerance */
  if(!getparint("npair", &opt->npair)) opt->npair = 5; /* l-BFGS memory length */
  if(!getparfloat("c1", &opt->c1)) opt->c1 = 1e-4; /* Nocedal value for Wolfe condition */
  if(!getparfloat("c2", &opt->c2)) opt->c2 = 0.9;  /* Nocedal value for Wolfe condition */
  if(!getparfloat("alpha", &opt->alpha)) opt->alpha = 1.;  /* initial step length */
  if(!getparint("bound", &opt->bound)) opt->bound = 0;/* 1 = bound on, 0 = off */
  if(!getparint("method", &opt->method)) opt->method = 0;//0=lBFGS; 1=Gauss-Newton
  if(!getparint("ncg", &opt->ncg)) opt->ncg = 5;//Gauss-Newton inversion
  if(!getparint("preco", &opt->preco)) opt->preco = 1;//1=precondition; 0=no precondition
  opt->verb = emf->verb; /* output message only on process 0, others remain silent */
  
  /*-------------------------------------------------------------------------*/
  fwi->preco = opt->preco;
  if(!getparint("npar", &fwi->npar)) fwi->npar = 2;//number of inversion parameters
  fwi->idxpar = alloc1int(fwi->npar);
  getparint("idxpar", fwi->idxpar);/* indices of the inversion parameters */
  if(opt->bound) {
    if(!(j=countparval("minpar"))) 
      err("Need lower bound of parameters minpar= vector");
    if(j!=fwi->npar) err("npar=%d, must have length[idxpar]=npar", fwi->npar);
    if(!(j=countparval("maxpar"))) 
      err("Need upper bound of parameters maxpar= vector");
    if(j!=fwi->npar) err("npar=%d, must have length[idxpar]=npar", fwi->npar);

    fwi->minpar=alloc1float(fwi->npar);
    fwi->maxpar=alloc1float(fwi->npar);
    getparfloat("minpar", fwi->minpar);
    getparfloat("maxpar", fwi->maxpar);
  }
  if(opt->verb){
    printf("---------------------------------------\n");
    printf("npar=%d\n", fwi->npar);
    printf("**** method=%d (0=l-BFGS; 1=Gauss-Newton)\n", opt->method);
    printf("**** niter=%d\n", opt->niter);
    printf("**** nls=%d\n", opt->nls);
    printf("**** bound=%d\n", opt->bound);
    printf("**** ncg=%d\n", opt->ncg);
    printf("**** preco=%d\n", opt->preco);
    if(opt->bound) printf("**** [lb,ub]=[%g, %g]\n", fwi->minpar[0], fwi->maxpar[0]);
  }  

  fwi->niter = opt->niter;
  fwi->iter = 0; /* before start set iter=0 */
  fwi->n = emf->nx*emf->ny*emf->nz*fwi->npar;/* lenth of unknown vector */
  fg_fwi_init(acqui, emf, fwi);

  /*---------------------------------------------------------------------------*/
  opt->x = alloc1float(fwi->n);/* model vector */
  opt->g = alloc1float(fwi->n);/* gradient vector */
  opt->d = alloc1float(fwi->n);/* descent direction */
  if(opt->method==0){
    opt->sk = alloc2float(fwi->n, opt->npair);
    opt->yk = alloc2float(fwi->n, opt->npair);
    opt->q = alloc1float(fwi->n);
  }
  for(j=0; j<fwi->npar; j++){
    for(i3=0; i3<emf->nz; i3++){
      for(i2=0; i2<emf->ny; i2++){
	for(i1=0; i1<emf->nx; i1++){
	  i = i1 + emf->nx*(i2 + emf->ny*(i3 + emf->nz*j));
	  if(fwi->idxpar[j]==1) opt->x[i] = log(emf->rho11[i3][i2][i1]);
	  if(fwi->idxpar[j]==2) opt->x[i] = log(emf->rho33[i3][i2][i1]);
	}
      }
    }
  }
  memcpy(fwi->xref, opt->x, fwi->n*sizeof(float));

  if(opt->bound){
    opt->xmin = alloc1float(fwi->n);
    opt->xmax = alloc1float(fwi->n);
    for(j = 0; j<fwi->npar; j++){
      for(i3=0; i3<emf->nz; i3++){
	for(i2=0; i2<emf->ny; i2++){
	  for(i1=0; i1<emf->nx; i1++){
	    i = i1 + emf->nx*(i2 + emf->ny*(i3 + emf->nz*j));
	    if(emf->oz+(i3+0.5)*emf->dz>emf->bathy[i2][i1]+emf->dz){
	      //below bathymetry, specify xmin,xmax
	      /* note that the unknown x[:] =log(m[:])*/
	      opt->xmin[i] = log(fwi->minpar[j]);
	      opt->xmax[i] = log(fwi->maxpar[j]);
	    }else{//above the bathymetry, use input parameter, same xmin,xmax
	      //this allows min value set to be formation 1 Ohm-m directly
	      opt->xmin[i] = opt->x[i];
	      opt->xmax[i] = opt->x[i];
	    }
	  }
	}
      }
    }
  }
  
  /*---------------------------------------------------------------------------*/  
  fcost = fg_fwi(opt->x, opt->g);  
  if(emf->mode==1){/* FWI */
    opt->f0 = fcost;
    opt->fk = fcost;
    opt->igrad = 0;
    opt->kpair = 0;
    opt->ils = 0;
    if(opt->verb){
      opt->gk_norm = l2norm(fwi->n, opt->g);
      fp=fopen("iterate.txt", "w");
      fprintf(fp, "==========================================================\n");
      fprintf(fp, "l-BFGS memory length: %d\n", opt->npair);
      fprintf(fp, "Maximum number of iterations: %d\n", opt->niter);
      fprintf(fp, "Convergence tolerance: %3.2e\n", opt->tol);
      fprintf(fp, "maximum number of line search: %d\n", opt->nls);
      fprintf(fp, "initial step length: alpha=%g\n", opt->alpha);
      fprintf(fp, "==========================================================\n");
      fprintf(fp, "iter    fk       fk/f0      ||gk||    alpha     nls   ngrad\n");
      fclose(fp);
    }
    /* l-BFGS/Newton-CG optimization */
    for(opt->iter=0; opt->iter<opt->niter; opt->iter++){
      fwi->iter = opt->iter;
      if(opt->verb){
	printf("==========================================================\n");
	printf("# iter=%d  fk/f0=%g\n", opt->iter, opt->fk/opt->f0);
	opt->gk_norm=l2norm(fwi->n, opt->g);

	fp=fopen("iterate.txt", "a");
	fprintf(fp, "%3d   %3.2e  %3.2e   %3.2e  %3.2e  %3d  %4d\n", 
		opt->iter, opt->fk, opt->fk/opt->f0, opt->gk_norm, opt->alpha, opt->ils, opt->igrad);
	fclose(fp);

	if(opt->iter==0) fp=fopen("rmse_misfit.txt", "w");
	else             fp=fopen("rmse_misfit.txt", "a");
	fprintf(fp, "%d \t %g\n", fwi->iter, fwi->mse);
	fclose(fp);
      }

      if(opt->method==1){
  	cg_solve(fwi->n, opt->x, opt->g, opt->d, Hv_fwi, opt);/* solve Hv = -g */
      }else{
	memcpy(opt->q, opt->g, fwi->n*sizeof(float));
  	if(opt->iter==0){/* first iteration, no stored gradient */
	  if(opt->preco) precondition(emf, fwi, opt->q);
	  flipsign(fwi->n, opt->q, opt->d);//descent direction=-gradient
	}else{
	  lbfgs_update(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
	  opt->rho = alloc1float(opt->kpair);
	  opt->alp = alloc1float(opt->kpair);
	  if(opt->preco) {
	    lbfgs_descent1(fwi->n, opt->g, opt->q, opt->rho, opt->alp, opt->sk, opt->yk, opt);
	    precondition(emf, fwi, opt->q);
	    if(opt->loop1) lbfgs_descent2(fwi->n, opt->g, opt->q, opt->rho, opt->alp, opt->sk, opt->yk, opt);
	  }else
	    lbfgs_descent(fwi->n, opt->g, opt->d, opt->sk, opt->yk, opt->q, opt->rho, opt->alp, opt);
	  free1float(opt->rho);
	  free1float(opt->alp);
	} 
	lbfgs_save(fwi->n, opt->x, opt->g, opt->sk, opt->yk, opt);
      }
      line_search(fwi->n, opt->x, opt->g, opt->d, fg_fwi, opt);

      /* not break, then line search succeeds or descent direction accepted */
      if(opt->verb) {/* print out inverted physical parameter at each iteration */
	sprintf(fname, "Rv_iteration_%02d", opt->iter);
	fp=fopen(fname, "wb");
	fwrite(&emf->rho33[0][0][0], emf->nx*emf->ny*emf->nz*sizeof(float), 1, fp);
	fclose(fp);

	sprintf(fname, "Rh_iteration_%02d", opt->iter);
	fp=fopen(fname, "wb");
	fwrite(&emf->rho11[0][0][0], emf->nx*emf->ny*emf->nz*sizeof(float), 1, fp);
	fclose(fp);
      }
    } 
    if(opt->verb && opt->iter==opt->niter) {
      fp=fopen("iterate.txt", "a");
      fprintf(fp, "==>Maximum iteration number reached!\n");
      fclose(fp);
    }
  }

  fg_fwi_close();
  free1float(opt->x);
  free1float(opt->g);
  free1float(opt->d);  
  if(opt->method==0){
    free2float(opt->sk);
    free2float(opt->yk);
    free1float(opt->q);
  }
  if(opt->bound){
    free1float(opt->xmin);
    free1float(opt->xmax);
    free1float(fwi->minpar);
    free1float(fwi->maxpar);
  }
  free1int(fwi->idxpar);
  free(opt);
  free(fwi);
}

