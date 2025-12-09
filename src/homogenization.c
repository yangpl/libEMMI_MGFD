/* Medium homogenization by volume averaging 
 * (tensor product of conservative averaging in 1D, see doc/note_averaging.pdf)
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

//volumetric averaging from in[nz][ny][nx] to out[n3][n2][n1]
//xx[],yy[],zz[] are of length: nx+1,ny+1,nz+1
//x1[],x2[],x3[] are of length: n1+1,n2+1,n3+1
void homogenization(float ***in, float ***out,
		    int nx, int ny, int nz, float *xx, float *yy, float *zz,
		    int n1, int n2, int n3, float *x1, float *x2, float *x3)
{
  float a, b, s;
  float left, right, overlap;
  float *rho1, *rho2;
  int ix, iy, iz, i, j, k;
  
  rho1 = malloc(n1*ny*nz*sizeof(float));
  for(iz=0; iz<nz; iz++){
    for(iy=0; iy<ny; iy++){

      i = 0;
      for(j=0; j<n1; j++) {
	a = x1[j];
	b = x1[j+1];
	s = 0.0;

	/* advance i to skip source cells entirely left of this target cell */
	while(i<nx && xx[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<nx && xx[k]<b) {
	  overlap = 0.0;
	  left = (xx[k]>a)? xx[k]:a;
	  right = (xx[k+1]<b)? xx[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += in[iz][iy][k]*overlap;
	  k++;
	}
	rho1[j + n1*(iy + ny*iz)] = s/(b-a);
      }
      
    }//end for iy
  }//end for iz
  
  rho2 = malloc(n1*n2*nz*sizeof(float));
  for(iz=0; iz<nz; iz++){
    for(ix=0; ix<n1; ix++){

      i = 0;
      for(j=0; j<n2; j++) {
	a = x2[j];
	b = x2[j + 1];
	s = 0.0;

	/* advance i to skip source cells entirely left of this target cell */
	while(i<ny && yy[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<ny && yy[k]<b) {
	  overlap = 0.0;
	  left = (yy[k]>a)? yy[k]:a;
	  right = (yy[k+1]<b)? yy[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += rho1[ix + n1*(k + ny*iz)]*overlap;
	  k++;
	}
	rho2[ix + n1*(j + n2*iz)] = s/(b-a);
      }
      
    }//end for ix
  }//end for iz
  free(rho1);

  for(iy=0; iy<n2; iy++){
    for(ix=0; ix<n1; ix++){

      i = 0;
      for(j=0; j<n3; j++) {
	a = x3[j];
	b = x3[j + 1];
	s = 0.0;
	/* advance i to skip source cells entirely left of this target cell */
	while(i<nz && zz[i+1]<=a) i++;

	k = i;
	/* accumulate overlaps with source cells intersecting [a,b] */
	while (k<nz && zz[k]<b) {
	  overlap = 0.0;
	  left = (zz[k]>a)? zz[k]:a;
	  right = (zz[k+1]<b)? zz[k+1]:b;
	  overlap = right - left;
	  if(overlap>0.0) s += rho2[ix + n1*(iy + n2*k)]*overlap;
	  k++;
	}
	out[j][iy][ix] = s/(b-a);
      }
      
    }//end for ix
  }//end for iy
  free(rho2);

}

