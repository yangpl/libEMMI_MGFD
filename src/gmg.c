/* Solving 3D diffusive Maxwell equation using geometric multigrid method
 * It uses finite integration method!
 *------------------------------------------------------------------------
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"
#include "constants.h"

int verb;
int cycleopt;//1=v cycle; 2=f cycle; 3=w cycle
int icycle, ncycle;//number of multigrid cycles
int isemicoarsen;//1=semi-coarsening; 0=no semi-coarsening
int v1;//number of pre-smoothing
int v2;//number of post-smoothing
int lmax;//maximum number of levels for switching between fine and coarse grid
double tol;//tolerance for convergence
double rnorm0, rnorm;//norm of residual vector r
double omega;//frequency omega=2*PI*freqx
emf_t *emf_;

typedef struct{
  int n1, n2, n3;//number of cells in x, y and z directions
  double *x1, *x2, *x3;//node coordinates along x, y and z
  double *x1s, *x2s, *x3s;//staggered node coordinates along x, y and z
  double *d1, *d2, *d3;//cell sizes centered at integer nodes
  double *d1s, *d2s, *d3s;//cell sizes centered at staggered nodes
  complex ****u, ****f, ****r;//r=f-Au, u=(Ex,Ey,Ez)^T, b=i*omega*mu*(Jx,Jy,Jz)
  double ***sigma11, ***sigma22, ***sigma33;//electrical conductivities
  double ***invmur;//inverse of relative magnetic permeability
  int sc[3];//semi-coarsening factor (2 or 1) in x, y, z directions
} gmg_t;
gmg_t *gmg;

//=============================================================
/*< LU decomposition of n*n matrix A=LU, stored in place >*/
void lu_solve(complex *A, complex *b, int n)
{
  int i, j, k;
  complex s;
  
  //LU decomposition: A=LU, L and U are stored in place
  for(k=0; k<n-1; k++){
    for(i=k+1; i<n; i++) A[i*n+k] /= A[k*n+k];

    for(i=k+1; i<n; i++)
      for(j=k+1; j<n; j++)
	A[i*n+j] -= A[i*n+k]*A[k*n+j];
  }

  //z = L^{-1}y
  for(i=0; i<n; i++){
    s = 0;
    for(j=0; j<i; j++) s += A[i*n+j]*b[j];
    b[i] -= s;
  }
  //x = U^{-1}z
  for(i=n-1; i>=0; i--){
    s = 0;
    for(j=i+1; j<n; j++) s += A[i*n+j]*b[j];
    b[i] = (b[i]-s)/A[i*n+i];
  }

}


//=============================================================
/*< ILU decomposition of n*n matrix A=LU >*/
void ilu0_solve(complex *A, complex *b, int n)
{
  int i, j, k;
  int jj, kk, jj_;
  complex s;

  //A stored in place at non-zero locations in Compressed Diagonal Storage format
  for(k=0; k<n-1; k++){
    for(i=k+1; i<MIN(n, k+6); i++){
      kk = k-i+5;
      A[i*11+kk] /= A[k*11+5];//a_ik \neq 0
    }
    
    for(i=k+1; i<MIN(n, k+6); i++){//a_ik \neq 0
      kk = k-i+5;
      for(j=k+1; j<MIN(n, k+6); j++){//a_kj \neq 0
	jj = j-i+5;
	jj_ = j-k+5;
	A[i*11+jj] -= A[i*11+kk]*A[k*11+jj_];
      }
    }
  }

  //z = L^{-1}y
  for(i=0; i<n; i++){
    s = 0;
    for(j=MAX(0,i-5); j<i; j++){
      jj = j-i+5;
      s += A[i*11+jj]*b[j];
    }
    b[i] -= s;
  }
  //x = U^{-1}z
  for(i=n-1; i>=0; i--){
    s = 0;
    for(j=i+1; j<MIN(n,i+6); j++){
      jj = j-i+5;
      s += A[i*11+jj]*b[j];
    }
    b[i] = (b[i]-s)/A[i*11+5];//a_ii, j=i-i+5
  }
}


/*< smoothing by pointwise lexicographic Gauss-Seidel relaxation sweep >*/
void gauss_seidel(gmg_t *gmg, int lev, int iter)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  int n1, n2, n3;
  complex t1, t2, t3, t4, Ax;
  double *d1s, *d2s, *d3s;
  complex ***Ex, ***Ey, ***Ez;
  complex ***bx, ***by, ***bz;
  double ***sigma11, ***sigma22, ***sigma33, ***invmur;
  double sigma11_, sigma22_, sigma33_;
  double m1, m2, m3, m4;
  complex **A, *b;
  
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  Ex = gmg[lev].u[0];
  Ey = gmg[lev].u[1];
  Ez = gmg[lev].u[2];
  bx = gmg[lev].f[0];
  by = gmg[lev].f[1];
  bz = gmg[lev].f[2];
  sigma11 = gmg[lev].sigma11;
  sigma22 = gmg[lev].sigma22;
  sigma33 = gmg[lev].sigma33;
  invmur = gmg[lev].invmur;
  
  A = alloc2complex(6, 6);
  b = alloc1complex(6);
  for(k=1; k<n3; k++){
    if(iter%2==1) k = n3 - k;//symmetric Gauss-Seidel
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=1; j<n2; j++){
      if(iter%2==1) j = n2 - j;//symmetric Gauss-Seidel
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=1; i<n1; i++){
	if(iter%2==1) i = n1 - i;//symmetric Gauss-Seidel
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	//initialize all unknowns to 0, the residuals are then right hand sides
	Ex[k][j][im1] = 0;
	Ex[k][j][i] = 0;
	Ey[k][jm1][i] = 0;
	Ey[k][j][i] = 0;
	Ez[km1][j][i] = 0;
	Ez[k][j][i] = 0;
	/*
	  sigma11_ = 0.25*(sigma11[k][j][i] + sigma11[k][jm1][i] + sigma11[km1][j][i] + sigma11[km1][jm1][i]);
	  sigma22_ = 0.25*(sigma22[k][j][i] + sigma22[k][j][im1] + sigma22[km1][j][i] + sigma22[km1][j][im1]);
	  sigma33_ = 0.25*(sigma33[k][j][i] + sigma33[k][j][im1] + sigma33[k][jm1][i] + sigma33[k][jm1][im1]);
	*/  
	//===================================================================
	//Ex(I-1,j,k)
	sigma11_ = 0.25*(sigma11[k][j][im1] + sigma11[k][jm1][im1] + sigma11[km1][j][im1] + sigma11[km1][jm1][im1]);
	t1 = (Ey[k][j][i] - Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1] - Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[k][jm1][i] - Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1] - Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[kp1][j][im1] - Ex[k][j][im1])/d3s[k] - (Ez[k][j][i] - Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[k][j][im1] - Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i] - Ez[km1][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	m1 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	m2 = 0.5*(invmur[k][jm1][im1] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	m3 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m4 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax =  ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][im1];
	b[0] = bx[k][j][im1] - Ax;
	A[0][0] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coeffient for Ex(I-1,j,k)
	A[0][1] = 0;//coefficient for Ex(I,j,k)
	A[0][2] = -m2/(d1s[im1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	A[0][3] = m1/(d1s[im1]*d2s[j]);//coefficient for Ey(i,J,k)
	A[0][4] = -m4/(d1s[im1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	A[0][5] = m3/(d1s[im1]*d3s[k]);//coefficient for Ez(i,j,K)
	//Ex(I,j,k)
	sigma11_ = 0.25*(sigma11[k][j][i] + sigma11[k][jm1][i] + sigma11[km1][j][i] + sigma11[km1][jm1][i]);
	t1 = (Ey[k][j][ip1] - Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i] - Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[k][jm1][ip1] - Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i] - Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[kp1][j][i] - Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1] - Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[k][j][i] - Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1] - Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	m1 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m2 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m4 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][i];
	b[1] = bx[k][j][i] - Ax;
	A[1][0] = 0;//coeffient for Ex(I-1,j,k)
	A[1][1] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coefficient for Ex(I,j,k)
	A[1][2] = m2/(d1s[i]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	A[1][3] = -m1/(d1s[i]*d2s[j]);//coefficient for Ey(i,J,k)
	A[1][4] = m4/(d1s[i]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	A[1][5] = -m3/(d1s[i]*d3s[k]);//coefficient for Ez(i,j,K)

	//Ey(i,J-1,k)
	sigma22_ = 0.25*(sigma22[k][jm1][i] + sigma22[k][jm1][im1] + sigma22[km1][jm1][i] + sigma22[km1][jm1][im1]);
	t1 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[k][jm1][ip1]-Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i]-Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[k][jm1][i]-Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1]-Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	m1 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	m2 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	m3 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m4 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][jm1][i];
	b[2] = by[k][jm1][i] - Ax;
	A[2][0] = -m4/(d2s[jm1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	A[2][1] = m3/(d2s[jm1]*d1s[i]);//coefficient for Ex(I,j,k)
	A[2][2] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J-1,k)
	A[2][3] = 0;//coefficient for Ey(i,J,k)
	A[2][4] = -m2/(d2s[jm1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	A[2][5] = m1/(d2s[jm1]*d3s[k]);//coefficient for Ez(i,j,K)
	//Ey(i,J,k)
	sigma22_ = 0.25*(sigma22[k][j][i] + sigma22[k][j][im1] + sigma22[km1][j][i] + sigma22[km1][j][im1]);
	t1 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[k][j][ip1]-Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i]-Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[k][j][i]-Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1]-Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	m1 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m2 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m3 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m4 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][j][i];
	b[3] = by[k][j][i] - Ax;
	A[3][0] = m4/(d2s[j]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	A[3][1] = -m3/(d2s[j]*d1s[i]);//coefficient for Ex(I,j,k)
	A[3][2] = 0;//coefficient for Ey(i,J-1,k)
	A[3][3] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J,k)
	A[3][4] = m2/(d2s[j]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	A[3][5] = -m1/(d2s[j]*d3s[k]);//coefficient for Ez(i,j,K)
	
	//Ez(i,j,K-1)
	sigma33_ = 0.25*(sigma33[km1][j][i] + sigma33[km1][j][im1] + sigma33[km1][jm1][i] + sigma33[km1][jm1][im1]);
	t1 = (Ex[k][j][i]-Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1]-Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[k][j][im1]-Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i]-Ez[km1][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	m1 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	m2 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	m3 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m4 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[km1][j][i];
	b[4] = bz[km1][j][i] - Ax;
	A[4][0] = -m2/(d3s[km1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	A[4][1] = m1/(d3s[km1]*d1s[i]);//coefficient for Ex(I,j,k)
	A[4][2] = -m4/(d3s[km1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	A[4][3] = m3/(d3s[km1]*d2s[j]);//coefficient for Ey(i,J,k)
	A[4][4] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K-1)
	A[4][5] = 0;//coefficient for Ez(i,j,K)
	//Ez(i,j,K)
	sigma33_ = 0.25*(sigma33[k][j][i] + sigma33[k][j][im1] + sigma33[k][jm1][i] + sigma33[k][jm1][im1]);
	t1 = (Ex[kp1][j][i]-Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1]-Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[kp1][j][im1]-Ex[k][j][im1])/d3s[k] - (Ez[k][j][i]-Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	m1 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m2 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m4 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[k][j][i];
	b[5] = bz[k][j][i] - Ax;
	A[5][0] = m2/(d3s[k]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	A[5][1] = -m1/(d3s[k]*d1s[i]);//coefficient for Ex(I,j,k)
	A[5][2] = m4/(d3s[k]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	A[5][3] = -m3/(d3s[k]*d2s[j]);//coefficient for Ey(i,J,k)
	A[5][4] = 0;//coefficient for Ez(i,j,K-1)
	A[5][5] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K)
		
	//solve Ax=b, A is 6*6, solution stored in b[]
	lu_solve(&A[0][0], b, 6);
	
	//asign b[] to Ex,Ey,Ez
	Ex[k][j][im1] = b[0];
	Ex[k][j][i] = b[1];
	Ey[k][jm1][i] = b[2];
	Ey[k][j][i] = b[3];
	Ez[km1][j][i] = b[4];
	Ez[k][j][i] = b[5];
      }
    }
  }

  free2complex(A);
  free1complex(b);
}


/*< relaxation/smoothing step for multigrid along line x >*/
void line_gauss_seidel_x(gmg_t *gmg, int lev, int iter)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  int n1, n2, n3;
  complex t1, t2, t3, t4;
  double *d1s, *d2s, *d3s;
  complex ***Ex, ***Ey, ***Ez;
  complex ***bx, ***by, ***bz;
  double ***sigma11, ***sigma22, ***sigma33, ***invmur;
  double sigma11_, sigma22_, sigma33_;
  double m1, m2, m3, m4;
  complex **A, *b, Ax;
  int row_ind, ii, nr;
    
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  Ex = gmg[lev].u[0];
  Ey = gmg[lev].u[1];
  Ez = gmg[lev].u[2];
  bx = gmg[lev].f[0];
  by = gmg[lev].f[1];
  bz = gmg[lev].f[2];
  sigma11 = gmg[lev].sigma11;
  sigma22 = gmg[lev].sigma22;
  sigma33 = gmg[lev].sigma33;
  invmur = gmg[lev].invmur;

  nr = 5*n1-4;//nr rows in total
  //matrix stored as CDS format, including 11 non-zeros in each row
  A = alloc2complex(11, nr);//CDS=compressed diagonal storage format
  b = alloc1complex(nr);
  memset(&A[0][0], 0, 11*nr*sizeof(complex));
  
  for(k=1; k<n3; k++){
    if(iter%2==1) k = n3 - k;
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=1; j<n2; j++){
      if(iter%2==1) j = n2 - j;
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      //initialize all unknowns to 0, the residuals are then right hand sides
      for(i=1; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	Ex[k][j][im1] = 0;
	Ey[k][jm1][i] = 0;
	Ey[k][j][i] = 0;
	Ez[km1][j][i] = 0;
	Ez[k][j][i] = 0;
	if(i==n1-1) Ex[k][j][i] = 0;
      }
            
      //initialize b and matrix A only at non-zero locations
      for(i=1; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);
	//starting row index and column index before i-th 6*6 block 
	//col_ind = 5*(i-1);
	row_ind = 5*(i-1);

	//Ex(I-1,j,k)
	ii = row_ind;
	sigma11_ = 0.25*(sigma11[k][j][im1] + sigma11[k][jm1][im1] + sigma11[km1][j][im1] + sigma11[km1][jm1][im1]);
	t1 = (Ey[k][j][i] - Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1] - Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[k][jm1][i] - Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1] - Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[kp1][j][im1] - Ex[k][j][im1])/d3s[k] - (Ez[k][j][i] - Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[k][j][im1] - Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i] - Ez[km1][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	m1 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	m2 = 0.5*(invmur[k][jm1][im1] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	m3 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m4 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax =  ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][im1];
	b[ii] = bx[k][j][im1] - Ax;
	//jj = col_ind;
	A[ii][5] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coeffient for Ex(I-1,j,k)
	//jj = col_ind+1;
	A[ii][6] = -m2/(d1s[im1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+2;
	A[ii][7] = m1/(d1s[im1]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+3;
	A[ii][8] = -m4/(d1s[im1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][9] = m3/(d1s[im1]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][10] = 0;//coefficient for Ex(I,j,k)

	//Ey(i,J-1,k)
	ii = row_ind+1;
	sigma22_ = 0.25*(sigma22[k][jm1][i] + sigma22[k][jm1][im1] + sigma22[km1][jm1][i] + sigma22[km1][jm1][im1]);
	t1 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[k][jm1][ip1]-Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i]-Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[k][jm1][i]-Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1]-Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	m1 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	m2 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	m3 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m4 = 0.5*(invmur[k][jm1][im1] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][jm1][i];
	b[ii] = by[k][jm1][i] - Ax;
	//jj = col_ind-4;
	A[ii][0] = (i>1)?(-m4/(d1s[im1]*d1s[im1])):0.;//coefficient for Ey(i-1,J-1,k)
	//jj = col_ind-3;
	A[ii][1] = 0;
	//jj = col_ind-2;
	A[ii][2] = 0;
	//jj = col_ind-1;
	A[ii][3] = 0;
	//jj = col_ind;
	A[ii][4] = -m4/(d2s[jm1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+1;
	A[ii][5] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J-1,k)
	//jj = col_ind+2;
	A[ii][6] = 0;//coefficient for Ey(i,J,k)
	//jj = col_ind+3;
	A[ii][7] = -m2/(d2s[jm1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][8] = m1/(d2s[jm1]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][9] = m3/(d2s[jm1]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+6;
	A[ii][10] = (i<n1-1)? (-m3/(d1s[i]*d1s[i])):0.;//coefficient for Ey(i+1,J-1,k)
	  
	//Ey(i,J,k)
	ii = row_ind+2;
	sigma22_ = 0.25*(sigma22[k][j][i] + sigma22[k][j][im1] + sigma22[km1][j][i] + sigma22[km1][j][im1]);
	t1 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[k][j][ip1]-Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i]-Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[k][j][i]-Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1]-Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	m1 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m2 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m3 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m4 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][j][i];
	b[ii] = by[k][j][i] - Ax;
	//jj = col_ind-3;
	A[ii][0] = (i>1)? (-m4/(d1s[im1]*d1s[im1])):0.;//coefficient for Ey(i-1,J,k)
	//jj = col_ind-2;
	A[ii][1] = 0;
	//jj = col_ind-1;
	A[ii][2] = 0;
	//jj = col_ind;
	A[ii][3] = m4/(d2s[j]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+1;
	A[ii][4] = 0;//coefficient for Ey(i,J-1,k)
	//jj = col_ind+2;
	A[ii][5] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J,k)
	//jj = col_ind+3;
	A[ii][6] = m2/(d2s[j]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][7] = -m1/(d2s[j]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][8] = -m3/(d2s[j]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+6;
	A[ii][9] = 0;
	//jj = col_ind+7;
	A[ii][10] = (i<n1-1)?(-m3/(d1s[i]*d1s[i])):0.;//coefficient for Ey(i+1,J,k)

	//Ez(i,j,K-1)
	ii = row_ind+3;
	sigma33_ = 0.25*(sigma33[km1][j][i] + sigma33[km1][j][im1] + sigma33[km1][jm1][i] + sigma33[km1][jm1][im1]);
	t1 = (Ex[k][j][i]-Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1]-Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[k][j][im1]-Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i]-Ez[km1][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	m1 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	m2 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	m3 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m4 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[km1][j][i];
	b[ii] = bz[km1][j][i] - Ax;
	//jj = col_ind-2;
	A[ii][0] = (i>1)?(-m2/(d1s[im1]*d1s[im1])):0.;//coefficient for Ez(i-1,j,K-1)
	//jj = col_ind-1;
	A[ii][1] = 0;
	//jj = col_ind;
	A[ii][2] = -m2/(d3s[km1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+1;
	A[ii][3] = -m4/(d3s[km1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+2;
	A[ii][4] = m3/(d3s[km1]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+3;
	A[ii][5] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][6] = 0;//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][7] = m1/(d3s[km1]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+6;
	A[ii][8] = 0;
	//jj = col_ind+7;
	A[ii][9] = 0;
	//jj = col_ind+8;
	A[ii][10] = (i<n1-1)?(-m1/(d1s[i]*d1s[i])):0.;//coefficient for Ez(i+1,j,K-1)

	//Ez(i,j,K)
	ii = row_ind+4;
	sigma33_ = 0.25*(sigma33[k][j][i] + sigma33[k][j][im1] + sigma33[k][jm1][i] + sigma33[k][jm1][im1]);
	t1 = (Ex[kp1][j][i]-Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1]-Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[kp1][j][im1]-Ex[k][j][im1])/d3s[k] - (Ez[k][j][i]-Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	m1 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m2 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m4 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[k][j][i];
	b[ii] = bz[k][j][i] - Ax;
	//jj = col_ind-1;
	A[ii][0] = (i>1)?(-m2/(d1s[im1]*d1s[im1])):0.;//coefficient for Ez(i-1,j,K)
	//jj = col_ind;
	A[ii][1] = m2/(d3s[k]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+1;
	A[ii][2] = m4/(d3s[k]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+2;
	A[ii][3] = -m3/(d3s[k]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+3;
	A[ii][4] = 0;//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][5] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][6] = -m1/(d3s[k]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+6;
	A[ii][7] = 0;
	//jj = col_ind+7;
	A[ii][8] = 0;
	//jj = col_ind+8;
	A[ii][9] = 0;
	//jj = col_ind+9;
	A[ii][10] = (i<n1-1)?(-m1/(d1s[i]*d1s[i])):0.;//coefficient for Ez(i+1,j,K)
	  
	//Ex(I,j,k)
	ii = row_ind+5;
	m1 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m2 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m4 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	//jj = col_ind;
	A[ii][0] = 0;//coeffient for Ex(I-1,j,k)
	//jj = col_ind+1;
	A[ii][1] = m2/(d1s[i]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+2;
	A[ii][2] = -m1/(d1s[i]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+3;
	A[ii][3] = m4/(d1s[i]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][4] = -m3/(d1s[i]*d3s[k]);//coefficient for Ez(i,j,K)
	if(i==n1-1) {
	  sigma11_ = 0.25*(sigma11[k][j][i] + sigma11[k][jm1][i] + sigma11[km1][j][i] + sigma11[km1][jm1][i]);
	  t1 = (Ey[k][j][ip1] - Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i] - Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t2 = (Ey[k][jm1][ip1] - Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i] - Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	  t3 = (Ex[kp1][j][i] - Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1] - Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t4 = (Ex[k][j][i] - Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1] - Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t1 *= m1;
	  t2 *= m2;
	  t3 *= m3;
	  t4 *= m4;
	  Ax = ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][i];
	  b[ii] = bx[k][j][i] - Ax;

	  //jj = col_ind+5;
	  A[ii][5] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coefficient for Ex(I,j,k)
	}
      }//end for i
      ilu0_solve(&A[0][0], b, nr);
      
      //copy the solution back to Ex, Ey, Ez
      for(i=1; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	row_ind = 5*(i-1);
	Ex[k][j][im1] = b[row_ind];
	Ey[k][jm1][i] = b[row_ind+1];
	Ey[k][j][i] = b[row_ind+2];
	Ez[km1][j][i] = b[row_ind+3];
	Ez[k][j][i] = b[row_ind+4];
	if(i==n1-1) Ex[k][j][i] = b[row_ind+5];
      }
      
    }
  }
  
  free2complex(A);
  free1complex(b);
}

/*< relaxation/smoothing step for multigrid along line y >*/
void line_gauss_seidel_y(gmg_t *gmg, int lev, int iter)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  int n1, n2, n3;
  complex t1, t2, t3, t4;
  double *d1s, *d2s, *d3s;
  complex ***Ex, ***Ey, ***Ez;
  complex ***bx, ***by, ***bz;
  double ***sigma11, ***sigma22, ***sigma33, ***invmur;
  double sigma11_, sigma22_, sigma33_;
  double m1, m2, m3, m4;
  complex **A, *b, Ax;
  int row_ind, ii, nr;

  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  Ex = gmg[lev].u[0];
  Ey = gmg[lev].u[1];
  Ez = gmg[lev].u[2];
  bx = gmg[lev].f[0];
  by = gmg[lev].f[1];
  bz = gmg[lev].f[2];
  sigma11 = gmg[lev].sigma11;
  sigma22 = gmg[lev].sigma22;
  sigma33 = gmg[lev].sigma33;
  invmur = gmg[lev].invmur;

  nr = 5*n2-4;//nr rows in total
  A = alloc2complex(11, nr);//every row has 11 non-zeros,
  b = alloc1complex(nr);
  
  for(k=1; k<n3; k++){
    if(iter%2==1) k = n3 - k;
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(i=1; i<n1; i++){
      if(iter%2==1) i = n1 - i;
      ip1 = MIN(i+1, n1);
      im1 = MAX(i-1, 0);
      //initialize all unknowns to 0, the residuals are then right hand sides
      for(j=1; j<n2; j++){
	jp1 = MIN(j+1, n1);
	jm1 = MAX(j-1, 0);

	Ey[k][jm1][i] = 0;
	Ex[k][j][im1] = 0;
	Ex[k][j][i] = 0;
	Ez[km1][j][i] = 0;
	Ez[k][j][i] = 0;
	if(j==n2-1) Ey[k][j][i] = 0;
      }

      //initialize b and matrix A only at non-zero locations
      for(j=1; j<n2; j++){
	jp1 = MIN(j+1, n1);
	jm1 = MAX(j-1, 0);

	row_ind = 5*(j-1);
	//col_ind = 5*(j-1);
	//Ey(i,J-1,k)
	ii = row_ind;
	sigma22_ = 0.25*(sigma22[k][jm1][i] + sigma22[k][jm1][im1] + sigma22[km1][jm1][i] + sigma22[km1][jm1][im1]);
	t1 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[k][jm1][ip1]-Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i]-Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey -\partial_y Ex
	t4 = (Ey[k][jm1][i]-Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1]-Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	m1 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	m2 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	m3 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m4 = 0.5*(invmur[k][jm1][im1] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][jm1][i];
	b[ii] = by[k][jm1][i] - Ax;
	//jj = col_ind;
	A[ii][5] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J-1,k)
	//jj = col_ind+1;
	A[ii][6] = -m4/(d2s[jm1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][7] = m3/(d2s[jm1]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][8] = -m2/(d2s[jm1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][9] = m1/(d2s[jm1]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][10] = 0;//coefficient for Ey(i,J,k)
	
	//Ex(I-1,j,k)
	ii = row_ind+1;
	sigma11_ = 0.25*(sigma11[k][j][im1] + sigma11[k][jm1][im1] + sigma11[km1][j][im1] + sigma11[km1][jm1][im1]);
	t1 = (Ey[k][j][i] - Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1] - Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[k][jm1][i] - Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1] - Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[kp1][j][im1] - Ex[k][j][im1])/d3s[k] - (Ez[k][j][i] - Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[k][j][im1] - Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i] - Ez[km1][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	m1 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	m2 = 0.5*(invmur[k][jm1][im1] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	m3 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m4 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax =  ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][im1];
	b[ii] = bx[k][j][im1] - Ax;
	//jj = col_ind-4;
	A[ii][0] = (j>1)?-m2/(d2s[jm1]*d2s[jm1]):0.;//coefficient for Ex(I-1,j-1,k)
	//jj = col_ind-3;
	A[ii][1] = 0;
	//jj = col_ind-2;
	A[ii][2] = 0;
	//jj = col_ind-1;
	A[ii][3] = 0;
	//jj = col_ind;
	A[ii][4] = -m2/(d1s[im1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+1;
	A[ii][5] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coeffient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][6] = 0;//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][7] = -m4/(d1s[im1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][8] = m3/(d1s[im1]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][9] = m1/(d1s[im1]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+6;
	A[ii][10] = (j<n2-1)?-m1/(d2s[j]*d2s[j]):0.;//coefficient for Ex(I-1,j+1,k)
	
	//Ex(I,j,k)
	ii = row_ind+2;
	sigma11_ = 0.25*(sigma11[k][j][i] + sigma11[k][jm1][i] + sigma11[km1][j][i] + sigma11[km1][jm1][i]);
	t1 = (Ey[k][j][ip1] - Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i] - Ex[k][j][i])/d2s[j];//\partial_x Ey -\partial_y Ex
	t2 = (Ey[k][jm1][ip1] - Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i] - Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[kp1][j][i] - Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1] - Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[k][j][i] - Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1] - Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	m1 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m2 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m4 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][i];
	b[ii] = bx[k][j][i] - Ax;
	//jj = col_ind-3;
	A[ii][0] = (j>1)?-m2/(d2s[jm1]*d2s[jm1]):0;//coefficient for Ex(I,j-1,k)
	//jj = col_ind-2;
	A[ii][1] = 0;
	//jj = col_ind-1;
	A[ii][2] = 0;
	//jj = col_ind;
	A[ii][3] = m2/(d1s[i]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+1;
	A[ii][4] = 0;//coeffient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][5] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][6] = m4/(d1s[i]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][7] = -m3/(d1s[i]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][8] = -m1/(d1s[i]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+6;
	A[ii][9] = 0;
	//jj = col_ind+7;
	A[ii][10] = (j<n2-1)?-m1/(d2s[j]*d2s[j]):0.;//coefficient for Ex(I,j+1,k)
		
	//Ez(i,j,K-1)
	ii = row_ind+3;
	sigma33_ = 0.25*(sigma33[km1][j][i] + sigma33[km1][j][im1] + sigma33[km1][jm1][i] + sigma33[km1][jm1][im1]);
	t1 = (Ex[k][j][i]-Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1]-Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[k][j][im1]-Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i]-Ez[km1][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	m1 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	m2 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	m3 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m4 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[km1][j][i];
	b[ii] = bz[km1][j][i] - Ax;
	//jj = col_ind-2;
	A[ii][0] = (j>1)?-m4/(d2s[jm1]*d2s[jm1]):0;//coefficient for Ez(i,j-1,K-1)
	//jj = col_ind-1;
	A[ii][1] = 0;
	//jj = col_ind;
	A[ii][2] = -m4/(d3s[km1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+1;
	A[ii][3] = -m2/(d3s[km1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][4] = m1/(d3s[km1]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][5] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][6] = 0;//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][7] = m3/(d3s[km1]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+6;
	A[ii][8] = 0;
	//jj = col_ind+7;
	A[ii][9] = 0;
	//jj = col_ind+8;
	A[ii][10] = (j<n2-1)?-m3/(d2s[j]*d2s[j]):0.;//coefficient for Ez(i,j+1,K-1)
	
	//Ez(i,j,K)
	ii = row_ind+4;
	sigma33_ = 0.25*(sigma33[k][j][i] + sigma33[k][j][im1] + sigma33[k][jm1][i] + sigma33[k][jm1][im1]);
	t1 = (Ex[kp1][j][i]-Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1]-Ez[k][j][i])/d1s[i];//\partial_z Ex -\partial_x Ez
	t2 = (Ex[kp1][j][im1]-Ex[k][j][im1])/d3s[k] - (Ez[k][j][i]-Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t3 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	m1 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m2 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m4 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[k][j][i];
	b[ii] = bz[k][j][i] - Ax;
	//jj = col_ind-1;
	A[ii][0] = (j>1)?-m4/(d2s[jm1]*d2s[jm1]):0.;//coefficient for Ez(i,j-1,K)
	//jj = col_ind;
	A[ii][1] = m4/(d3s[k]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+1;
	A[ii][2] = m2/(d3s[k]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][3] = -m1/(d3s[k]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][4] = 0;//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][5] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K)
	//jj = col_ind+5;
	A[ii][6] = -m3/(d3s[k]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+6;
	A[ii][7] = 0;
	//jj = col_ind+7;
	A[ii][8] = 0;
	//jj = col_ind+8;
	A[ii][9] = 0;
	//jj = col_ind+9;
	A[ii][10] = (j<n2-1)?-m3/(d2s[j]*d2s[j]):0.;//coefficient for Ez(i,j+1,K)
		
	//Ey(i,J,k)
	ii = row_ind+5;
	m1 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m2 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m3 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m4 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	//jj = col_ind;
	A[ii][0] = 0;//coefficient for Ey(i,J-1,k)
	//jj = col_ind+1;
	A[ii][1] = m4/(d2s[j]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][2] = -m3/(d2s[j]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][3] = m2/(d2s[j]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+4;
	A[ii][4] = -m1/(d2s[j]*d3s[k]);//coefficient for Ez(i,j,K)
	if(j==n2-1){
	  sigma22_ = 0.25*(sigma22[k][j][i] + sigma22[k][j][im1] + sigma22[km1][j][i] + sigma22[km1][j][im1]);
	  t1 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t2 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	  t3 = (Ey[k][j][ip1]-Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i]-Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t4 = (Ey[k][j][i]-Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1]-Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t1 *= m1;
	  t2 *= m2;
	  t3 *= m3;
	  t4 *= m4;
	  Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][j][i];
	  b[ii] = by[k][j][i] - Ax;

	  //jj = col_ind+5;
	  A[ii][5] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J,k)
	}
      }
      ilu0_solve(&A[0][0], b, nr);

      //copy the solution back to Ex, Ey, Ez
      for(j=1; j<n2; j++){
	jp1 = MIN(j+1, n3);
	jm1 = MAX(j-1, 0);

	row_ind = 5*(j-1);	
	Ey[k][jm1][i] = b[row_ind];
	Ex[k][j][im1] = b[row_ind+1];
	Ex[k][j][i] = b[row_ind+2];
	Ez[km1][j][i] = b[row_ind+3];
	Ez[k][j][i] = b[row_ind+4];
	if(j<n2) Ey[k][j][i] = b[row_ind+5];
      }

    }
  }

  free2complex(A);
  free1complex(b);
}


/*< relaxation/smoothing step for multigrid along line z >*/
void line_gauss_seidel_z(gmg_t *gmg, int lev, int iter)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  int n1, n2, n3;
  complex t1, t2, t3, t4;
  double *d1s, *d2s, *d3s;
  complex ***Ex, ***Ey, ***Ez;
  complex ***bx, ***by, ***bz;
  double ***sigma11, ***sigma22, ***sigma33, ***invmur;
  double sigma11_, sigma22_, sigma33_;
  double m1, m2, m3, m4;
  complex **A, *b, Ax;
  int row_ind, ii, nr;

  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  Ex = gmg[lev].u[0];
  Ey = gmg[lev].u[1];
  Ez = gmg[lev].u[2];
  bx = gmg[lev].f[0];
  by = gmg[lev].f[1];
  bz = gmg[lev].f[2];
  sigma11 = gmg[lev].sigma11;
  sigma22 = gmg[lev].sigma22;
  sigma33 = gmg[lev].sigma33;
  invmur = gmg[lev].invmur;
  
  nr = 5*n3-4;//nr rows in total
  A = alloc2complex(11, nr);//every row has 11 non-zeros,
  b = alloc1complex(nr);
  
  for(j=1; j<n2; j++){
    if(iter%2==1) j = n2 - j;
    jp1 = MIN(j+1, n2);
    jm1 = MAX(j-1, 0);
    for(i=1; i<n1; i++){
      if(iter%2==1) i = n1 - i;
      ip1 = MIN(i+1, n1);
      im1 = MAX(i-1, 0);
      //initialize all unknowns to 0, the residuals are then right hand sides
      for(k=1; k<n3; k++){
	kp1 = MIN(k+1, n3);
	km1 = MAX(k-1, 0);

	Ez[km1][j][i] = 0;
	Ex[k][j][im1] = 0;
	Ex[k][j][i] = 0;
	Ey[k][jm1][i] = 0;
	Ey[k][j][i] = 0;
	if(k==n3-1) Ez[k][j][i] = 0.;
      }

      //initialize b and matrix A only at non-zero locations
      for(k=1; k<n3; k++){
	kp1 = MIN(k+1, n3);
	km1 = MAX(k-1, 0);

	row_ind = 5*(k-1);
	//col_ind = 5*(k-1);

	//Ez(i,j,K-1)
	ii = row_ind;
	sigma33_ = 0.25*(sigma33[km1][j][i] + sigma33[km1][j][im1] + sigma33[km1][jm1][i] + sigma33[km1][jm1][im1]);
	t1 = (Ex[k][j][i]-Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1]-Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t2 = (Ex[k][j][im1]-Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i]-Ez[km1][j][im1])/d1s[im1];//\partial_z Ex -\partial_x Ez
	t3 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t4 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	m1 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	m2 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	m3 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m4 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[km1][j][i];
	b[ii] = bz[km1][j][i] - Ax;
	//jj = col_ind;
	A[ii][5] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K-1)
	//jj = col_ind+1;
	A[ii][6] = -m2/(d3s[km1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][7] = m1/(d3s[km1]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][8] = -m4/(d3s[km1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+4;
	A[ii][9] = m3/(d3s[km1]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+5;
	A[ii][10] = 0;//coefficient for Ez(i,j,K)

	//Ex(I-1,j,k)
	ii = row_ind+1;
	sigma11_ = 0.25*(sigma11[k][j][im1] + sigma11[k][jm1][im1] + sigma11[km1][j][im1] + sigma11[km1][jm1][im1]);
	t1 = (Ey[k][j][i] - Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1] - Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[k][jm1][i] - Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1] - Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[kp1][j][im1] - Ex[k][j][im1])/d3s[k] - (Ez[k][j][i] - Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[k][j][im1] - Ex[km1][j][im1])/d3s[km1] - (Ez[km1][j][i] - Ez[km1][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	m1 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	m2 = 0.5*(invmur[k][jm1][im1] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	m3 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m4 = 0.5*(invmur[km1][j][im1] + invmur[km1][jm1][im1]);//Hy(I-1,j,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax =  ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][im1];
	b[ii] = bx[k][j][im1] - Ax;
	//jj = col_ind-4;
	A[ii][0] = (k>1)?-m4/(d3s[km1]*d3s[km1]):0.;//coeffient for Ex(I-1,j,k-1)
	//jj = col_ind-3;
	A[ii][1] = 0;
	//jj = col_ind-2;
	A[ii][2] = 0;
	//jj = col_ind-1;
	A[ii][3] = 0;
	//jj = col_ind;
	A[ii][4] = -m4/(d1s[im1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+1;
	A[ii][5] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coeffient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][6] = 0;//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][7] = -m2/(d1s[im1]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+4;
	A[ii][8] = m1/(d1s[im1]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+5;
	A[ii][9] = m3/(d1s[im1]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+6;
	A[ii][10] = (k<n3-1)?-m3/(d3s[k]*d3s[k]):0;//coeffient for Ex(I-1,j,k+1)

	//Ex(I,j,k)
	ii = row_ind+2;
	sigma11_ = 0.25*(sigma11[k][j][i] + sigma11[k][jm1][i] + sigma11[km1][j][i] + sigma11[km1][jm1][i]);
	t1 = (Ey[k][j][ip1] - Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i] - Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	t2 = (Ey[k][jm1][ip1] - Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i] - Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t3 = (Ex[kp1][j][i] - Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1] - Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	t4 = (Ex[k][j][i] - Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1] - Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	m1 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m2 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m4 = 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][i];
	b[ii] = bx[k][j][i] - Ax;
	//jj = col_ind-3;
	A[ii][0] = (k>1)?-m4/(d3s[km1]*d3s[km1]):0;//coeffient for Ex(I,j,k-1)
	//jj = col_ind-2;
	A[ii][1] = 0;
	//jj = col_ind-1;
	A[ii][2] = 0;
	//jj = col_ind;
	A[ii][3] = m4/(d1s[i]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+1;
	A[ii][4] = 0;//coeffient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][5] = (m1/(d2s[j]*d2s[j]) + m2/(d2s[jm1]*d2s[jm1])) + (m3/(d3s[k]*d3s[k]) + m4/(d3s[km1]*d3s[km1])) - I*omega*mu0*sigma11_;//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][6] = m2/(d1s[i]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+4;
	A[ii][7] = -m1/(d1s[i]*d2s[j]);//coefficient for Ey(i,J,k)
	//jj = col_ind+5;
	A[ii][8] = -m3/(d1s[i]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+6;
	A[ii][9] = 0;
	//jj = col_ind+7;
	A[ii][10] = (k<n3-1)?-m3/(d3s[k]*d3s[k]):0.;//coeffient for Ex(I,j,k+1)

	//Ey(i,J-1,k)
	ii = row_ind+3;
	sigma22_ = 0.25*(sigma22[k][jm1][i] + sigma22[k][jm1][im1] + sigma22[km1][jm1][i] + sigma22[km1][jm1][im1]);
	t1 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[km1][j][i]-Ez[km1][jm1][i])/d2s[jm1] - (Ey[k][jm1][i]-Ey[km1][jm1][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[k][jm1][ip1]-Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i]-Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[k][jm1][i]-Ey[k][jm1][im1])/d1s[im1] - (Ex[k][j][im1]-Ex[k][jm1][im1])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	m1 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	m2 = 0.5*(invmur[km1][jm1][i] + invmur[km1][jm1][im1]);//Hx(i,J-1,K-1)
	m3 = 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	m4 = 0.5*(invmur[k][jm1][im1] + invmur[km1][jm1][im1]);//Hz(I-1,J-1,k)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][jm1][i];
	b[ii] = by[k][jm1][i] - Ax;
	//jj = col_ind-2;
	A[ii][0] = (k>1)?-m2/(d3s[km1]*d3s[km1]):0;//coeffient for Ey(i,J-1,k-1)
	//jj = col_ind-1;
	A[ii][1] = 0;
	//jj = col_ind;
	A[ii][2] = -m2/(d2s[jm1]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+1;
	A[ii][3] = -m4/(d2s[jm1]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][4] = m3/(d2s[jm1]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][5] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J-1,k)
	//jj = col_ind+4;
	A[ii][6] = 0;//coefficient for Ey(i,J,k)
	//jj = col_ind+5;
	A[ii][7] = m1/(d2s[jm1]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+6;
	A[ii][8] = 0;
	//jj = col_ind+7;
	A[ii][9] = 0;
	//jj = col_ind+8;
	A[ii][10] = (k<n3-1)?-m1/(d3s[k]*d3s[k]):0;//coeffient for Ey(i,J-1,k+1)
	
	//Ey(i,J,k)
	ii = row_ind+4;
	sigma22_ = 0.25*(sigma22[k][j][i] + sigma22[k][j][im1] + sigma22[km1][j][i] + sigma22[km1][j][im1]);
	t1 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	t2 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	t3 = (Ey[k][j][ip1]-Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i]-Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	t4 = (Ey[k][j][i]-Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1]-Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	m1 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m2 = 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	m3 = 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	m4 = 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	t1 *= m1;
	t2 *= m2;
	t3 *= m3;
	t4 *= m4;
	Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][j][i];
	b[ii] = by[k][j][i] - Ax;
	//jj = col_ind-1;
	A[ii][0] = (k>1)?-m2/(d3s[km1]*d3s[km1]):0;//coeffient for Ey(i,J,k-1)
	//jj = col_ind;
	A[ii][1] = m2/(d2s[j]*d3s[km1]);//coefficient for Ez(i,j,K-1)
	//jj = col_ind+1;
	A[ii][2] = m4/(d2s[j]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][3] = -m3/(d2s[j]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][4] = 0;//coefficient for Ey(i,J-1,k)
	//jj = col_ind+4;
	A[ii][5] = (m1/(d3s[k]*d3s[k]) + m2/(d3s[km1]*d3s[km1])) + (m3/(d1s[i]*d1s[i]) + m4/(d1s[im1]*d1s[im1])) - I*omega*mu0*sigma22_;//coefficient for Ey(i,J,k)
	//jj = col_ind+5;
	A[ii][6] = -m1/(d2s[j]*d3s[k]);//coefficient for Ez(i,j,K)
	//jj = col_ind+6;
	A[ii][7] = 0;
	//jj = col_ind+7;
	A[ii][8] = 0;
	//jj = col_ind+8;
	A[ii][9] = 0;
	//jj = col_ind+9;
	A[ii][10] = (k<n3-1)?-m1/(d3s[k]*d3s[k]):0;//coeffient for Ey(i,J,k+1)

	//Ez(i,j,K)
	ii = row_ind+5;
	m1 = 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	m2 = 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	m3 = 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	m4 = 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	//jj = col_ind;
	A[ii][0] = 0;//coefficient for Ez(i,j,K-1)
	//jj = col_ind+1;
	A[ii][1] = m2/(d3s[k]*d1s[im1]);//coefficient for Ex(I-1,j,k)
	//jj = col_ind+2;
	A[ii][2] = -m1/(d3s[k]*d1s[i]);//coefficient for Ex(I,j,k)
	//jj = col_ind+3;
	A[ii][3] = m4/(d3s[k]*d2s[jm1]);//coefficient for Ey(i,J-1,k)
	//jj = col_ind+4;
	A[ii][4] = -m3/(d3s[k]*d2s[j]);//coefficient for Ey(i,J,k)
	if(k==n3-1){
	  sigma33_ = 0.25*(sigma33[k][j][i] + sigma33[k][j][im1] + sigma33[k][jm1][i] + sigma33[k][jm1][im1]);
	  t1 = (Ex[kp1][j][i]-Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1]-Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t2 = (Ex[kp1][j][im1]-Ex[k][j][im1])/d3s[k] - (Ez[k][j][i]-Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	  t3 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t4 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t1 *= m1;
	  t2 *= m2;
	  t3 *= m3;
	  t4 *= m4;
	  Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[k][j][i];
	  b[ii] = bz[k][j][i] - Ax;

	  //jj = col_ind+5;
	  A[ii][5] = (m1/(d1s[i]*d1s[i]) + m2/(d1s[im1]*d1s[im1])) + (m3/(d2s[j]*d2s[j]) + m4/(d2s[jm1]*d2s[jm1])) - I*omega*mu0*sigma33_;//coefficient for Ez(i,j,K)
	}
      }
      ilu0_solve(&A[0][0], b, nr);

      //copy the solution back to Ex, Ey, Ez
      for(k=1; k<n3; k++){
	kp1 = MIN(k+1, n3);
	km1 = MAX(k-1, 0);

	row_ind = 5*(k-1);
	Ez[km1][j][i] = b[row_ind];
	Ex[k][j][im1] = b[row_ind+1];
	Ex[k][j][i] = b[row_ind+2];
	Ey[k][jm1][i] = b[row_ind+3];
	Ey[k][j][i] = b[row_ind+4];
	if(k==n3-1) Ez[k][j][i] = b[row_ind+5];
      }
    }
  }

  free2complex(A);
  free1complex(b);
}

/*< multigrid smoothing/relaxation >*/
void smoothing(gmg_t *gmg, int lev, int iter)
{
  if(isemicoarsen){
    //line Gauss-Seidel smoothing (LGS)=a special case of block Gauss-Seidel (BGS)
    //block Gauss-Seidel (BGS) also known as multiplicative Schwarz method
    if(gmg[lev+1].sc[0]==2) line_gauss_seidel_x(gmg, lev, iter);
    if(gmg[lev+1].sc[1]==2) line_gauss_seidel_y(gmg, lev, iter);
    if(gmg[lev+1].sc[2]==2) line_gauss_seidel_z(gmg, lev, iter);
  }else{
    //pointwise Gauss-Seidel sweeping 
    gauss_seidel(gmg, lev, iter);
  }
}

/*< compute residual r=f-Au >*/
void residual(gmg_t *gmg, int lev)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int im1, jm1, km1;
  int n1, n2, n3, n;
  complex t1, t2, t3, t4, Ax;
  double *d1s, *d2s, *d3s;
  complex ***Ex, ***Ey, ***Ez;
  complex ***bx, ***by, ***bz;
  complex ***rx, ***ry, ***rz;
  double ***sigma11, ***sigma22, ***sigma33, ***invmur;
  double sigma11_, sigma22_, sigma33_;
  
  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  Ex = gmg[lev].u[0];
  Ey = gmg[lev].u[1];
  Ez = gmg[lev].u[2];
  bx = gmg[lev].f[0];
  by = gmg[lev].f[1];
  bz = gmg[lev].f[2];
  sigma11 = gmg[lev].sigma11;
  sigma22 = gmg[lev].sigma22;
  sigma33 = gmg[lev].sigma33;
  invmur = gmg[lev].invmur;
  rx = gmg[lev].r[0];
  ry = gmg[lev].r[1];
  rz = gmg[lev].r[2];

  n = (n1+1)*(n2+1)*(n3+1);
  memcpy(&rx[0][0][0], &bx[0][0][0], n*sizeof(complex));
  memcpy(&ry[0][0][0], &by[0][0][0], n*sizeof(complex));
  memcpy(&rz[0][0][0], &bz[0][0][0], n*sizeof(complex));
  
  for(k=0; k<n3; k++){
    kp1 = MIN(k+1, n3);
    km1 = MAX(k-1, 0);
    for(j=0; j<n2; j++){
      jp1 = MIN(j+1, n2);
      jm1 = MAX(j-1, 0);
      for(i=0; i<n1; i++){
	ip1 = MIN(i+1, n1);
	im1 = MAX(i-1, 0);

	if(j>0 && k>0){
	  sigma11_ = 0.25*(sigma11[k][j][i] + sigma11[k][jm1][i] + sigma11[km1][j][i] + sigma11[km1][jm1][i]);
	  t1 = (Ey[k][j][ip1] - Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i] - Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t2 = (Ey[k][jm1][ip1] - Ey[k][jm1][i])/d1s[i] - (Ex[k][j][i] - Ex[k][jm1][i])/d2s[jm1];//\partial_x Ey - \partial_y Ex
	  t3 = (Ex[kp1][j][i] - Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1] - Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t4 = (Ex[k][j][i] - Ex[km1][j][i])/d3s[km1] - (Ez[km1][j][ip1] - Ez[km1][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t1 *= 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	  t2 *= 0.5*(invmur[k][jm1][i] + invmur[km1][jm1][i]);//Hz(I,J-1,k)
	  t3 *= 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	  t4 *= 0.5*(invmur[km1][j][i] + invmur[km1][jm1][i]);//Hy(I,j,K-1)
	  Ax = ((t1/d2s[j]-t2/d2s[jm1]) - (t3/d3s[k]-t4/d3s[km1])) - I*omega*mu0*sigma11_*Ex[k][j][i];
	  rx[k][j][i] -= Ax;
	}

	if(i>0 && k>0){
	  sigma22_ = 0.25*(sigma22[k][j][i] + sigma22[k][j][im1] + sigma22[km1][j][i] + sigma22[km1][j][im1]);
	  t1 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t2 = (Ez[km1][jp1][i]-Ez[km1][j][i])/d2s[j] - (Ey[k][j][i]-Ey[km1][j][i])/d3s[km1];//\partial_y Ez - \partial_z Ey
	  t3 = (Ey[k][j][ip1]-Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i]-Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t4 = (Ey[k][j][i]-Ey[k][j][im1])/d1s[im1] - (Ex[k][jp1][im1]-Ex[k][j][im1])/d2s[j];//\partial_x Ey - \partial_y Ex
	  t1 *= 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	  t2 *= 0.5*(invmur[km1][j][i] + invmur[km1][j][im1]);//Hx(i,J,K-1)
	  t3 *= 0.5*(invmur[k][j][i] + invmur[km1][j][i]);//Hz(I,J,k)
	  t4 *= 0.5*(invmur[k][j][im1] + invmur[km1][j][im1]);//Hz(I-1,J,k)
	  Ax = ((t1/d3s[k]-t2/d3s[km1]) - (t3/d1s[i]-t4/d1s[im1])) - I*omega*mu0*sigma22_*Ey[k][j][i];
	  ry[k][j][i] -= Ax;
	}

	if(i>0 && j>0){
	  sigma33_ = 0.25*(sigma33[k][j][i] + sigma33[k][j][im1] + sigma33[k][jm1][i] + sigma33[k][jm1][im1]);
	  t1 = (Ex[kp1][j][i]-Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1]-Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	  t2 = (Ex[kp1][j][im1]-Ex[k][j][im1])/d3s[k] - (Ez[k][j][i]-Ez[k][j][im1])/d1s[im1];//\partial_z Ex - \partial_x Ez
	  t3 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t4 = (Ez[k][j][i]-Ez[k][jm1][i])/d2s[jm1] - (Ey[kp1][jm1][i]-Ey[k][jm1][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	  t1 *= 0.5*(invmur[k][j][i] + invmur[k][jm1][i]);//Hy(I,j,K)
	  t2 *= 0.5*(invmur[k][j][im1] + invmur[k][jm1][im1]);//Hy(I-1,j,K)
	  t3 *= 0.5*(invmur[k][j][i] + invmur[k][j][im1]);//Hx(i,J,K)
	  t4 *= 0.5*(invmur[k][jm1][i] + invmur[k][jm1][im1]);//Hx(i,J-1,K)
	  Ax = ((t1/d1s[i]-t2/d1s[im1]) - (t3/d2s[j]-t4/d2s[jm1])) - I*omega*mu0*sigma33_*Ez[k][j][i];
	  rz[k][j][i] -= Ax;
	}
      }
    }
  }

}


/*< compute H and store it in gmg[0].f, must have lev=0>*/
void compute_H_from_E(gmg_t *gmg, int lev)
{
  int i, j, k;
  int ip1, jp1, kp1;
  int n1, n2, n3, n;
  complex t1, t2, t3;
  double *d1s, *d2s, *d3s;
  complex ***Ex, ***Ey, ***Ez;
  complex ***Hx, ***Hy, ***Hz;

  n1 = gmg[lev].n1;
  n2 = gmg[lev].n2;
  n3 = gmg[lev].n3;
  d1s = gmg[lev].d1s;
  d2s = gmg[lev].d2s;
  d3s = gmg[lev].d3s;
  Ex = gmg[lev].u[0];
  Ey = gmg[lev].u[1];
  Ez = gmg[lev].u[2];
  Hx = gmg[lev].f[0];
  Hy = gmg[lev].f[1];
  Hz = gmg[lev].f[2];

  n = (n1+1)*(n2+1)*(n3+1);
  memset(&Hx[0][0][0], 0, n*sizeof(complex));
  memset(&Hy[0][0][0], 0, n*sizeof(complex));
  memset(&Hz[0][0][0], 0, n*sizeof(complex));
  
  for(k=0; k<n3; k++){
    kp1 = MIN(k+1, n3);
    for(j=0; j<n2; j++){
      jp1 = MIN(j+1, n2);
      for(i=0; i<n1; i++){
	ip1 = MIN(i+1, n1);

	//here we assume mu=mu0, note invmur has been multiplied with cell volume
	t1 = (Ez[k][jp1][i]-Ez[k][j][i])/d2s[j] - (Ey[kp1][j][i]-Ey[k][j][i])/d3s[k];//\partial_y Ez - \partial_z Ey
	Hx[k][j][i] = t1/(I*omega*mu0);//Hx(i,J,K)
	t2 = (Ex[kp1][j][i]-Ex[k][j][i])/d3s[k] - (Ez[k][j][ip1]-Ez[k][j][i])/d1s[i];//\partial_z Ex - \partial_x Ez
	Hy[k][j][i] = t2/(I*omega*mu0);//Hy(I,j,K)
	t3 = (Ey[k][j][ip1]-Ey[k][j][i])/d1s[i] - (Ex[k][jp1][i]-Ex[k][j][i])/d2s[j];//\partial_x Ey - \partial_y Ex
	Hz[k][j][i] = t3/(I*omega*mu0);//Hz(I,J,k)
      }
    }
  }
  
}

/*< prolongation u from (lev+1) to lev-th grid: coarse to fine grid >*/
void prolongation(gmg_t *gmg, int lev)
{
  int i, j, k;
  int ii, jj, kk;
  int ip1, jp1, kp1;
  double w1l, w2l, w3l;
  double w1r, w2r, w3r;
  complex ***ux, ***uy, ***uz;
  complex ***ex, ***ey, ***ez;

  ex = gmg[lev+1].u[0];
  ey = gmg[lev+1].u[1];
  ez = gmg[lev+1].u[2];
  ux = gmg[lev].u[0];
  uy = gmg[lev].u[1];
  uz = gmg[lev].u[2];
  for(kk=0; kk<gmg[lev].n3; kk++){
    k = (gmg[lev+1].sc[2]==2)?(kk/2):kk;
    kp1 = k+1;
    w3r = (gmg[lev].x3[kk] - gmg[lev+1].x3[k])/(gmg[lev+1].x3[kp1] - gmg[lev+1].x3[k]);
    w3l = 1. - w3r;
    for(jj=0; jj<gmg[lev].n2; jj++){
      j = (gmg[lev+1].sc[1]==2)?(jj/2):jj;
      jp1 = j+1;
      w2r = (gmg[lev].x2[jj] - gmg[lev+1].x2[j])/(gmg[lev+1].x2[jp1] - gmg[lev+1].x2[j]);
      w2l = 1. - w2r;
      for(ii=0; ii<gmg[lev].n1; ii++){
	i = (gmg[lev+1].sc[0]==2)?(ii/2):ii;
	ip1 = i+1;
	w1r = (gmg[lev].x1[ii] - gmg[lev+1].x1[i])/(gmg[lev+1].x1[ip1] - gmg[lev+1].x1[i]);
	w1l = 1. - w1r;

	//ii can be 2*i or 2*i+1, we perform same for both cases, that is constant interpolation along x, bilinear interpolation along y and z
	if(jj>0 && kk>0) ux[kk][jj][ii] += w2l*(w3l*ex[k][j][i] + w3r*ex[kp1][j][i]) + w2r*(w3l*ex[k][jp1][i] + w3r*ex[kp1][jp1][i]);
	//jj can be 2*j or 2*j+1, we perform same for both cases, that is constant interpolation along y, bilinear interpolation along x and z
	if(ii>0 && kk>0) uy[kk][jj][ii] += w1l*(w3l*ey[k][j][i] + w3r*ey[kp1][j][i]) + w1r*(w3l*ey[k][j][ip1] + w3r*ey[kp1][j][ip1]);
	//kk can be 2*k or 2*k+1, we perform same for both cases, that is constant interpolation along z, bilinear interpolation along x and y
	if(ii>0 && jj>0) uz[kk][jj][ii] += w1l*(w2l*ez[k][j][i] + w2r*ez[k][jp1][i]) + w1r*(w2l*ez[k][j][ip1] + w2r*ez[k][jp1][ip1]);
      }
    }
  }
  
}


/*< restrict r from lev to (lev+1)-th lev: fine to coarse grid >*/
//restriction should be the exact adjoint of prolongation!
void restriction(gmg_t *gmg, int lev)
{
  int n1, n2, n3;
  int i, j, k;
  int ii, jj, kk;
  int ip1, jp1, kp1;
  double w1l, w2l, w3l;
  double w1r, w2r, w3r;
  complex ***fx, ***fy, ***fz;
  complex ***rx, ***ry, ***rz;
  
  n1 = gmg[lev+1].n1;
  n2 = gmg[lev+1].n2;
  n3 = gmg[lev+1].n3;
  fx = gmg[lev+1].f[0];
  fy = gmg[lev+1].f[1];
  fz = gmg[lev+1].f[2];
  rx = gmg[lev].r[0];
  ry = gmg[lev].r[1];
  rz = gmg[lev].r[2];

  memset(&fx[0][0][0], 0, (n1+1)*(n2+1)*(n3+1)*sizeof(complex));
  memset(&fy[0][0][0], 0, (n1+1)*(n2+1)*(n3+1)*sizeof(complex));
  memset(&fz[0][0][0], 0, (n1+1)*(n2+1)*(n3+1)*sizeof(complex));
  for(kk=0; kk<gmg[lev].n3; kk++){
    k = (gmg[lev+1].sc[2]==2)?(kk/2):kk;
    kp1 = k+1;
    w3r = (gmg[lev].x3[kk] - gmg[lev+1].x3[k])/(gmg[lev+1].x3[kp1] - gmg[lev+1].x3[k]);
    w3l = 1. - w3r;
    for(jj=0; jj<gmg[lev].n2; jj++){
      j = (gmg[lev+1].sc[1]==2)?(jj/2):jj;
      jp1 = j+1;
      w2r = (gmg[lev].x2[jj] - gmg[lev+1].x2[j])/(gmg[lev+1].x2[jp1] - gmg[lev+1].x2[j]);
      w2l = 1. - w2r;
      for(ii=0; ii<gmg[lev].n1; ii++){
	i = (gmg[lev+1].sc[0]==2)?(ii/2):ii;
	ip1 = i+1;
	w1r = (gmg[lev].x1[ii] - gmg[lev+1].x1[i])/(gmg[lev+1].x1[ip1] - gmg[lev+1].x1[i]);
	w1l = 1. - w1r;

	//adjoint of prolongation
	if(jj>0 && kk>0){
	  fx[k][j][i] += w2l*w3l*rx[kk][jj][ii];
	  fx[kp1][j][i] += w2l*w3r*rx[kk][jj][ii];
	  fx[k][jp1][i] += w2r*w3l*rx[kk][jj][ii];
	  fx[kp1][jp1][i] += w2r*w3r*rx[kk][jj][ii];
	}
	if(ii>0 && kk>0){
	  fy[k][j][i] += w1l*w3l*ry[kk][jj][ii];
	  fy[kp1][j][i] += w1l*w3r*ry[kk][jj][ii];
	  fy[k][j][ip1] += w1r*w3l*ry[kk][jj][ii];
	  fy[kp1][j][ip1] += w1r*w3r*ry[kk][jj][ii];
	}
	if(ii>0 && jj>0){
	  fz[k][j][i] += w1l*w2l*rz[kk][jj][ii];
	  fz[k][jp1][i] += w1l*w2r*rz[kk][jj][ii];
	  fz[k][j][ip1] += w1r*w2l*rz[kk][jj][ii];
	  fz[k][jp1][ip1] += w1r*w2r*rz[kk][jj][ii];
	}
      }
    }
  }
  
}

/*< compute dot product between two vectors s=<x,y> >*/
complex inner_product(int n, complex *x, complex *y)
{
  int i;
  complex s;

  for(i=0, s=0; i<n; i++) s += x[i]*conj(y[i]);
  return s;
}

/*< multigrid V-cycle >*/
void v_cycle(gmg_t *gmg, int lev)
{
  int n, i;

  for(i=0; i<v1; i++) smoothing(gmg, lev, i);//pre-smoothing of u based on u,f at lev-th level
  if(lev<lmax-1 && gmg[lev+1].sc[0]*gmg[lev+1].sc[1]*gmg[lev+1].sc[2]>1){
    residual(gmg, lev);//residual r=f-Au at lev-th lev
    if(cycleopt==1 && lev==0){//compute the norm of the residual
      n = 3*(gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);
      rnorm = sqrt(creal(inner_product(n, &gmg[lev].r[0][0][0][0], &gmg[lev].r[0][0][0][0])));
      if(verb) printf("icycle=%d rnorm=%e\n", icycle, rnorm);

      if(icycle==0) rnorm0 = rnorm;
      else if(rnorm<rnorm0*tol) { icycle=ncycle; return; }
    }
    restriction(gmg, lev);//restrict gmg[lev].r to gmg[lev+1].f 

    n = 3*(gmg[lev+1].n1+1)*(gmg[lev+1].n2+1)*(gmg[lev+1].n3+1);
    memset(&gmg[lev+1].u[0][0][0][0], 0, n*sizeof(complex));
    v_cycle(gmg, lev+1);// another v-cycle at (lev+1)-th level

    prolongation(gmg, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
  }
  //if lev==lmax-1, then nx=ny=2, grid size=3*3, only 1 point at the center is unknwn
  //direct solve (equivalent to smoothing at center point) by one post-smoothing will do the job  
  for(i=0; i<v2; i++) smoothing(gmg, lev, i);//post-smoothing
}

/*< multigrid F-cycle >*/
void f_cycle(gmg_t *gmg, int lev)
{
  int n;

  if(lev==lmax-1){//coarsest grid, direct solve or smoothing
    n = 3*(gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);	
    memset(&gmg[lev].u[0][0][0][0], 0, n*sizeof(complex));
  }else{
    residual(gmg, lev);//residual r=f-Au at lev-th level
    if(cycleopt==2 && lev==0){//compute the norm of the residual vector at 1st iteration
      n = 3*(gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);
      rnorm = sqrt(creal(inner_product(n, &gmg[lev].r[0][0][0][0], &gmg[lev].r[0][0][0][0])));
      if(verb) printf("icycle=%d rnorm=%e\n", icycle, rnorm);
      
      if(icycle==0) rnorm0 = rnorm;
      else if(rnorm<rnorm0*tol) { icycle=ncycle; return; }
    }
    restriction(gmg, lev);//restrict r at lev-th lev to f at (lev+1)-th level

    n = 3*(gmg[lev+1].n1+1)*(gmg[lev+1].n2+1)*(gmg[lev+1].n3+1);
    memset(&gmg[lev+1].u[0][0][0][0], 0, n*sizeof(complex));
    f_cycle(gmg, lev+1);
    prolongation(gmg, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
  }
  v_cycle(gmg, lev);// another v-cycle at lev-th level  
}

/*< multigrid W cycle>*/
void w_cycle(gmg_t *gmg, int lev)
{
  int n, i;

  for(i=0; i<v1; i++) smoothing(gmg, lev, i);//pre-smoothing of u based on u,f at lev-th level
  if(lev<lmax-1){
    residual(gmg, lev);//residual r=f-Au at lev-th lev
    if(cycleopt==3 && lev==0){//compute the norm of the residual
      n = 3*(gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);
      rnorm = sqrt(creal(inner_product(n, &gmg[lev].r[0][0][0][0], &gmg[lev].r[0][0][0][0])));
      if(verb) printf("icycle=%d rnorm=%e\n", icycle, rnorm);
      
      if(icycle==0) rnorm0 = rnorm;
      else if(rnorm<rnorm0*tol) { icycle=ncycle; return; }
    }
    restriction(gmg, lev);//restrict r at lev-th lev to gmg[lev+1].f 

    n = 3*(gmg[lev+1].n1+1)*(gmg[lev+1].n2+1)*(gmg[lev+1].n3+1);
    memset(&gmg[lev+1].u[0][0][0][0], 0, n*sizeof(complex));
    w_cycle(gmg, lev+1);// another w-cycle at (lev+1)-th level
    w_cycle(gmg, lev+1);// another w-cycle at (lev+1)-th level

    prolongation(gmg, lev);//interpolate r^h=gmg[lev+1].u to r^2h from (lev+1) to lev-th level
  }
  for(i=0; i<v2; i++) smoothing(gmg, lev, i);//post-smoothing
}



/*< use x1[] to derive x1s[], d1[], d1s[] >*/
void generate_xs_dx(int n1, double *x1, double *x1s, double *d1, double *d1s)
{
  int i;
  
  for(i=0; i<n1; i++) {
    d1s[i] = x1[i+1] - x1[i];
    x1s[i] = 0.5*(x1[i] + x1[i+1]);
  }
  x1s[n1] = x1[n1] + 0.5*d1s[n1-1];
  d1s[n1] = x1s[n1] - x1[n1];

  for(i=1; i<=n1; i++) d1[i] = x1s[i] - x1s[i-1];
  d1[0] = x1s[0] - x1[0];
}

//get the total depth of multigrid 
int get_depth(int n)
{
  int lmax, nn;
  lmax = 0;
  nn = n;
  while(nn%2==0&& nn>1){
    lmax++;
    nn = nn>>1;
  }
  if(nn>2) lmax++;
  return lmax;
}

void grid_init(gmg_t *gmg, int lev);

/*< initialize multigrid solver with different parameters >*/
void gmg_init(emf_t *emf, int ifreq)
{
  int lmax1, lmax2, lmax3;
  
  emf_ = emf;
  omega = emf->omegas[ifreq];
  verb = emf->verb; //verbosity display, 1=verbose, 0=not 
  if(!getparint("cycleopt", &cycleopt)) cycleopt = 2;/* 1=V cycle; 2=F cycle; 3=W cycle */
  if(!getparint("ncycle", &ncycle)) ncycle = 30;/* number of multigrid cycles */  
  if(!getparint("v1", &v1)) v1 = 1;/* number of pre-smoothing */
  if(!getparint("v2", &v2)) v2 = 1;/* number of post-smoothing */
  if(!getpardouble("tol", &tol)) tol = 1e-6;/* stopping criteria */
  if(!getparint("isemicoarsen", &isemicoarsen)) isemicoarsen = 1;/*1=semi-coarsening; 0=no semi-coarsening */
  if(!getparint("lmax", &lmax)) {
    lmax1 = get_depth(emf->n1);
    lmax2 = get_depth(emf->n2);
    lmax3 = get_depth(emf->n3);
    lmax = MAX(MAX(lmax1, lmax2), lmax3);
  }
  if(verb) {
    printf("------------------ GMG init -------------------\n");
    printf("cycleopt=%d (1=V cycle; 2=F cycle; 3=W cycle)\n", cycleopt);
    printf("number of V/F/W cycles, ncycle=%d\n", ncycle);
    printf("number of pre-smoothing, v1=%d\n", v1);
    printf("number of post-smoothing, v2=%d\n", v2);
    printf("number of grid levels, [lmax1,lmax2,lmax3]=[%d,%d,%d]\n", lmax1,lmax2,lmax3);
    printf("isemicoarsen=%d (1=semi-coarsening, 0=full coarsening)\n", isemicoarsen);
    printf("----------------------------------------------\n");
  }  
  gmg = (gmg_t*)malloc(lmax*sizeof(gmg_t));
  grid_init(gmg, 0);
}

/*< setup the grid and the conductivity at different multigrid level >*/
void grid_init(gmg_t *gmg, int lev)
{
  int i, j, k;
  int ii, jj, kk;
  int iip1, jjp1, kkp1;
  int n1, n2, n3, repeat;
  double val;

  if(lev==0){
    gmg[lev].sc[0] = 1;//coarsening factor x
    gmg[lev].sc[1] = 1;//coarsening factor y
    gmg[lev].sc[2] = 1;//coarsening factor z
    n1 = emf_->n1;
    n2 = emf_->n2;
    n3 = emf_->n3;
  }else{ 
    gmg[lev].sc[0] = (gmg[lev-1].n1%2==0&&gmg[lev-1].n1>2)?2:1;//coarsening factor x
    gmg[lev].sc[1] = (gmg[lev-1].n2%2==0&&gmg[lev-1].n2>2)?2:1;//coarsening factor y
    gmg[lev].sc[2] = (gmg[lev-1].n3%2==0&&gmg[lev-1].n3>2)?2:1;//coarsening factor z
    if(isemicoarsen) gmg[lev].sc[icycle%3] = 1;
  
    n1 = gmg[lev-1].n1/gmg[lev].sc[0];
    n2 = gmg[lev-1].n2/gmg[lev].sc[1];
    n3 = gmg[lev-1].n3/gmg[lev].sc[2];    
    //if(verb &&gmg[lev].sc[0]*gmg[lev].sc[1]*gmg[lev].sc[2]>1) printf("lev=%d, [n1,n2,n3]=[%d,%d,%d], [scx,scy,scz]=[%d,%d,%d]\n", lev, n1, n2, n3, gmg[lev].sc[0], gmg[lev].sc[1], gmg[lev].sc[2]); 
  }
  gmg[lev].n1 = n1;
  gmg[lev].n2 = n2;
  gmg[lev].n3 = n3;
  gmg[lev].x1 = alloc1double(n1+1);
  gmg[lev].x2 = alloc1double(n2+1);
  gmg[lev].x3 = alloc1double(n3+1);
  gmg[lev].x1s = alloc1double(n1+1);
  gmg[lev].x2s = alloc1double(n2+1);
  gmg[lev].x3s = alloc1double(n3+1);
  gmg[lev].d1 = alloc1double(n1+1);
  gmg[lev].d2 = alloc1double(n2+1);
  gmg[lev].d3 = alloc1double(n3+1);
  gmg[lev].d1s = alloc1double(n1+1);
  gmg[lev].d2s = alloc1double(n2+1);
  gmg[lev].d3s = alloc1double(n3+1);
  gmg[lev].sigma11 = alloc3double(n1, n2, n3);
  gmg[lev].sigma22 = alloc3double(n1, n2, n3);
  gmg[lev].sigma33 = alloc3double(n1, n2, n3);
  gmg[lev].invmur = alloc3double(n1, n2, n3);
  gmg[lev].u = alloc4complex(n1+1, n2+1, n3+1, 3);
  gmg[lev].f = alloc4complex(n1+1, n2+1, n3+1, 3);
  gmg[lev].r = alloc4complex(n1+1, n2+1, n3+1, 3);

  for(i=0; i<=n1; i++) gmg[lev].x1[i] = (lev==0)?emf_->x1[i]:gmg[lev-1].x1[gmg[lev].sc[0]*i];
  for(i=0; i<=n2; i++) gmg[lev].x2[i] = (lev==0)?emf_->x2[i]:gmg[lev-1].x2[gmg[lev].sc[1]*i];
  for(i=0; i<=n3; i++) gmg[lev].x3[i] = (lev==0)?emf_->x3[i]:gmg[lev-1].x3[gmg[lev].sc[2]*i];
  generate_xs_dx(n1, gmg[lev].x1, gmg[lev].x1s, gmg[lev].d1, gmg[lev].d1s);
  generate_xs_dx(n2, gmg[lev].x2, gmg[lev].x2s, gmg[lev].d2, gmg[lev].d2s);
  generate_xs_dx(n3, gmg[lev].x3, gmg[lev].x3s, gmg[lev].d3, gmg[lev].d3s);

  if(lev==0){
    for(k=0; k<n3; k++){
      for(j=0; j<n2; j++){
	for(i=0; i<n1; i++){
	  gmg[lev].sigma11[k][j][i] = emf_->sigma11[k][j][i];
	  gmg[lev].sigma22[k][j][i] = emf_->sigma22[k][j][i];
	  gmg[lev].sigma33[k][j][i] = emf_->sigma33[k][j][i];
	  gmg[lev].invmur[k][j][i] = emf_->invmur[k][j][i];
	}
      }
    }
  }else{//restrict model from fine grid to coarse grid
    repeat = 8/(gmg[lev].sc[0]*gmg[lev].sc[1]*gmg[lev].sc[2]);
    for(k=0; k<n3; k++){
      kk = (gmg[lev].sc[2]==2)?2*k:k;
      kkp1 = (gmg[lev].sc[2]==2)?(2*k+1):k;
      for(j=0; j<n2; j++){
	jj = (gmg[lev].sc[1]==2)?2*j:j;
	jjp1 = (gmg[lev].sc[1]==2)?(2*j+1):j;
	for(i=0; i<n1; i++){
	  ii = (gmg[lev].sc[0]==2)?2*i:i;
	  iip1 = (gmg[lev].sc[0]==2)?(2*i+1):i;
	  
	  val = gmg[lev-1].sigma11[kk][jj][ii] 
	    + gmg[lev-1].sigma11[kk][jj][iip1] 
	    + gmg[lev-1].sigma11[kk][jjp1][ii] 
	    + gmg[lev-1].sigma11[kk][jjp1][iip1] 
	    + gmg[lev-1].sigma11[kkp1][jj][ii] 
	    + gmg[lev-1].sigma11[kkp1][jj][iip1]
	    + gmg[lev-1].sigma11[kkp1][jjp1][ii]
	    + gmg[lev-1].sigma11[kkp1][jjp1][iip1];
	  gmg[lev].sigma11[k][j][i] = val/repeat;

	  val = gmg[lev-1].sigma22[kk][jj][ii]
	    + gmg[lev-1].sigma22[kk][jj][iip1]
	    + gmg[lev-1].sigma22[kk][jjp1][ii]
	    + gmg[lev-1].sigma22[kk][jjp1][iip1]
	    + gmg[lev-1].sigma22[kkp1][jj][ii]
	    + gmg[lev-1].sigma22[kkp1][jj][iip1]
	    + gmg[lev-1].sigma22[kkp1][jjp1][ii]
	    + gmg[lev-1].sigma22[kkp1][jjp1][iip1];
	  gmg[lev].sigma22[k][j][i] = val/repeat;

	  val = gmg[lev-1].sigma33[kk][jj][ii]
	    + gmg[lev-1].sigma33[kk][jj][iip1]
	    + gmg[lev-1].sigma33[kk][jjp1][ii]
	    + gmg[lev-1].sigma33[kk][jjp1][iip1]
	    + gmg[lev-1].sigma33[kkp1][jj][ii]
	    + gmg[lev-1].sigma33[kkp1][jj][iip1]
	    + gmg[lev-1].sigma33[kkp1][jjp1][ii]
	    + gmg[lev-1].sigma33[kkp1][jjp1][iip1];
	  gmg[lev].sigma33[k][j][i] = val/repeat;

	  val = gmg[lev-1].invmur[kk][jj][ii]
	    + gmg[lev-1].invmur[kk][jj][iip1]
	    + gmg[lev-1].invmur[kk][jjp1][ii]
	    + gmg[lev-1].invmur[kk][jjp1][iip1]
	    + gmg[lev-1].invmur[kkp1][jj][ii]
	    + gmg[lev-1].invmur[kkp1][jj][iip1]
	    + gmg[lev-1].invmur[kkp1][jjp1][ii]
	    + gmg[lev-1].invmur[kkp1][jjp1][iip1];
	  gmg[lev].invmur[k][j][i] = val/repeat;	  
	}
      }
    }
      
  }//end if

}

/*< free the variables and the memory at different multigrid levels >*/
void grid_close(gmg_t *gmg, int lev)
{
  free1double(gmg[lev].x1);
  free1double(gmg[lev].x2);
  free1double(gmg[lev].x3);
  free1double(gmg[lev].x1s);
  free1double(gmg[lev].x2s);
  free1double(gmg[lev].x3s);
  free1double(gmg[lev].d1);
  free1double(gmg[lev].d2);
  free1double(gmg[lev].d3);
  free1double(gmg[lev].d1s);
  free1double(gmg[lev].d2s);
  free1double(gmg[lev].d3s);
  free3double(gmg[lev].sigma11);
  free3double(gmg[lev].sigma22);
  free3double(gmg[lev].sigma33);
  free3double(gmg[lev].invmur);
  free4complex(gmg[lev].u);
  free4complex(gmg[lev].f);
  free4complex(gmg[lev].r);
}

/*< free the multigrid solver >*/
void gmg_close()
{
  grid_close(gmg, 0);
  free(gmg);
}


/*< apply the multigrid as a linear solver/preconditioner >*/
void gmg_apply(int n, complex *b, complex *x)
{
  int i, j, k, lev;
  double vol;

  memcpy(&gmg[0].f[0][0][0][0], b, n*sizeof(complex));
  memset(&gmg[0].u[0][0][0][0], 0, n*sizeof(complex));
  for(k=0; k<gmg[0].n3; k++){
    for(j=0; j<gmg[0].n2; j++){
      for(i=0; i<gmg[0].n1; i++){
	//multiply volume 
	vol = gmg[0].d1s[i]*gmg[0].d2s[j]*gmg[0].d3s[k];
	gmg[0].sigma11[k][j][i] *= vol;
	gmg[0].sigma22[k][j][i] *= vol;
	gmg[0].sigma33[k][j][i] *= vol;
	gmg[0].invmur[k][j][i] *= vol;
      }
    }
  }
  for(k=0; k<=gmg[0].n3; k++){
    for(j=0; j<=gmg[0].n2; j++){
      for(i=0; i<=gmg[0].n1; i++){
	//multiply volume on both left and right sides
	vol = gmg[0].d1s[i]*gmg[0].d2[j]*gmg[0].d3[k];
	gmg[0].f[0][k][j][i] *= vol;
	vol = gmg[0].d1[i]*gmg[0].d2s[j]*gmg[0].d3[k];
	gmg[0].f[1][k][j][i] *= vol;
	vol = gmg[0].d1[i]*gmg[0].d2[j]*gmg[0].d3s[k];
	gmg[0].f[2][k][j][i] *= vol;
      }
    }
  }

  for(icycle=0; icycle<ncycle; icycle++){
    for(lev=1; lev<lmax; lev++) grid_init(gmg, lev);
    if(cycleopt==1) v_cycle(gmg, 0);
    if(cycleopt==2) f_cycle(gmg, 0);
    if(cycleopt==3) w_cycle(gmg, 0);//not recommended
    for(lev=1; lev<lmax; lev++) grid_close(gmg, lev);

    /*
    // Gauss-Seidel iterations without v-cycle
    residual(gmg, lev);//residual r=f-Au at lev-th lev
    n = 3*(gmg[lev].n1+1)*(gmg[lev].n2+1)*(gmg[lev].n3+1);
    rnorm = sqrt(creal(inner_product(n, &gmg[lev].r[0][0][0][0], &gmg[lev].r[0][0][0][0])));
    printf("rnorm=%e\n", rnorm);
    smoothing(gmg, lev, 1);
    */
  }
  memcpy(x, &gmg[0].u[0][0][0][0], n*sizeof(complex));//copy E into x
  compute_H_from_E(gmg, 0);//compute H and store it in gmg[0].f
  memcpy(b, &gmg[0].f[0][0][0][0], n*sizeof(complex));//copy H into b
}
