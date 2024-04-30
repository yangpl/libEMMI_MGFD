#ifndef fwi_h
#define fwi_h

typedef struct {
  int npar;
  int n; /* total number of unknowns in fwi */

  int preco;//0=no precodition; 1=depth precondition; 2=triangle smoothing
  int r1, r2, r3, repeat;
  
  int *idxpar;
  float *minpar;
  float *maxpar;
  int iter;
  int niter;

  float *xref; //a reference model to compute regularization term
  float gamma1;//Tikhonov regularization parameter
  float gamma2;//TV regularization parameter

  int firstgrad;//1=first grad; 0= grad at later iterations
  float alpha; //step length
  float fcost_dat, fcost_mod, fcost;
  float *fcost_list;
  float ****grad, ****phess, ***Hv;//gradient, pseudo-Hessian and Hessian vector product

  float mse; // target mse misfit
} fwi_t;

#endif
