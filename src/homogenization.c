/* Medium homogenization by volume averaging
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

//medium averaging from rho_in[nz][ny][nx] to sigma_out[n3][n2][n1]
//flag=1, average over sigma
//flag=2, average over rho
//flag=3, average over log(sigma)
void homogenization(emf_t *emf, float ***rho_in, float ***sigma_out, int flag)
{
  int i1, i2, i3;
  int j1, j2, j3;
  int k1, k2, k3;
  float s, inta, intb;

  float ***intx = alloc3float(emf->nx, emf->ny, emf->nz);
  //integral along x
  for(i3=0; i3<emf->nz; i3++){
    for(i2=0; i2<emf->ny; i2++){
      s = 0;
      for(i1=0; i1<emf->nx; i1++){
	if(flag==1) s += 1./rho_in[i3][i2][i1]*emf->dx;
	if(flag==2) s += rho_in[i3][i2][i1]*emf->dx;
	if(flag==3) s += log(1./rho_in[i3][i2][i1])*emf->dx;
	intx[i3][i2][i1] = s;
      }
    }
  }  
  float ***avex = alloc3float(emf->n1, emf->ny, emf->nz);  
  for(i3=0; i3<emf->nz; i3++){
    for(i2=0; i2<emf->ny; i2++){
      for(i1=0; i1<emf->n1; i1++){
	j1 = (emf->x1[i1+1] - emf->ox)/emf->dx;//integer part
	if(j1<=0){
	  inta = intx[i3][i2][0]*(emf->x1[i1+1] - emf->ox)/emf->dx;
	}else if(j1>emf->nx-1){
	  s = (emf->x1[i1+1]-emf->ox-emf->nx*emf->dx)/emf->dx;
	  inta = intx[i3][i2][emf->nx-1] + s*(intx[i3][i2][emf->nx-1]-intx[i3][i2][emf->nx-2]);
	}else{
	  s = (emf->x1[i1+1] - emf->ox - j1*emf->dx)/emf->dx;//fraction after integer part
	  inta = intx[i3][i2][j1-1] + s*(intx[i3][i2][j1]-intx[i3][i2][j1-1]);
	}

	k1 = (emf->x1[i1] - emf->ox)/emf->dx;//integer part
	if(k1<=0){
	  intb = intx[i3][i2][0]*(emf->x1[i1] - emf->ox)/emf->dx;
	}else if(k1>emf->nx-1){
	  s = (emf->x1[i1]-emf->ox-emf->nx*emf->dx)/emf->dx;
	  intb = intx[i3][i2][emf->nx-1] + s*(intx[i3][i2][emf->nx-1]-intx[i3][i2][emf->nx-2]);
	}else{
	  s = (emf->x1[i1] - emf->ox - k1*emf->dx)/emf->dx;//fraction after integer part
	  intb = intx[i3][i2][k1-1] + s*(intx[i3][i2][k1]-intx[i3][i2][k1-1]);
	}

	avex[i3][i2][i1] = (inta - intb)/(emf->x1s[i1+1]-emf->x1s[i1]);//volume average
      }//end for i1
    }//end for i2
  }//end for i3
  free3float(intx);

  float ***inty = alloc3float(emf->n1, emf->ny, emf->nz);
  //integral along y
  for(i3=0; i3<emf->nz; i3++){
    for(i1=0; i1<emf->n1; i1++){
      s = 0;
      for(i2=0; i2<emf->ny; i2++){
	s += avex[i3][i2][i1]*emf->dy;
	inty[i3][i2][i1] = s;
      }
    }
  }
  free3float(avex);

  float ***avey = alloc3float(emf->n1, emf->n2, emf->nz);
  for(i3=0; i3<emf->nz; i3++){
    for(i1=0; i1<emf->n1; i1++){
      for(i2=0; i2<emf->n2; i2++){
	j2 = (emf->x2[i2+1]-emf->oy)/emf->dy;//integer part
	if(j2<=0){
	  inta = inty[i3][0][i1]*(emf->x2[i2+1] - emf->oy)/emf->dy;
	}else if(j2>emf->ny-1){
	  s = (emf->x2[i2+1]-emf->oy-emf->ny*emf->dy)/emf->dy;
	  inta = inty[i3][emf->ny-1][i1] + s*(inty[i3][emf->ny-1][i1]-inty[i3][emf->ny-2][i1]);
	}else{
	  s = (emf->x2[i2+1] - emf->oy - j2*emf->dy)/emf->dy;//fraction after integer part
	  inta = inty[i3][j2-1][i1] + s*(inty[i3][j2][i1]-inty[i3][j2-1][i1]);
	}

	k2 = (emf->x2[i2]-emf->oy)/emf->dy;//integer part
	if(k2<=0){
	  intb = inty[i3][0][i1]*(emf->x2[i2] - emf->oy)/emf->dy;
	}else if(k2>emf->ny-1){
	  s = (emf->x2[i2]-emf->oy-emf->ny*emf->dy)/emf->dy;
	  intb = inty[i3][emf->ny-1][i1] + s*(inty[i3][emf->ny-1][i1]-inty[i3][emf->ny-2][i1]);
	}else{
	  s = (emf->x2[i2] - emf->oy - k2*emf->dy)/emf->dy;//fraction after integer part
	  intb = inty[i3][k2-1][i1] + s*(inty[i3][k2][i1]-inty[i3][k2-1][i1]);
	}

	avey[i3][i2][i1] = (inta - intb)/(emf->x2[i2+1] - emf->x2[i2]);//volume average
      }//end for i1
    }//end for i2
  }//end for i3
  free3float(inty);
  
  float ***intz = alloc3float(emf->n1, emf->n2, emf->nz);
  //integral along z
  for(i2=0; i2<emf->n2; i2++){
    for(i1=0; i1<emf->n1; i1++){
      s = 0;
      for(i3=0; i3<emf->nz; i3++){
	s += avey[i3][i2][i1]*emf->dz;
	intz[i3][i2][i1] = s;
      }
    }
  }
  free3float(avey);

  for(i2=0; i2<emf->n2; i2++){
    for(i1=0; i1<emf->n1; i1++){
      for(i3=0; i3<emf->n3; i3++){
	j3 = (emf->x3[i3+1]-emf->oz)/emf->dz;
	if(j3<=0){
	  inta = intz[0][i2][i1]*(emf->x3[i3+1] - emf->oz)/emf->dz;
	}else if(j3>emf->nz-1){
	  s = (emf->x3[i3+1]-emf->oz-emf->nz*emf->dz)/emf->dz;
	  inta = intz[emf->nz-1][i2][i1] + s*(intz[emf->nz-1][i2][i1]-intz[emf->nz-2][i2][i1]);
	}else{
	  s = (emf->x3[i3+1] - emf->oz - j3*emf->dz)/emf->dz;//fraction after integer part
	  inta = intz[j3-1][i2][i1] + s*(intz[j3][i2][i1]-intz[j3-1][i2][i1]);
	}

	k3 = (emf->x3[i3]-emf->oz)/emf->dz;
	if(k3<=0){
	  intb = intz[0][i2][i1]*(emf->x3[i3] - emf->oz)/emf->dz;
	}else if(k3>emf->nz-1){
	  s = (emf->x3[i3]-emf->oz-emf->nz*emf->dz)/emf->dz;
	  intb = intz[emf->nz-1][i2][i1] + s*(intz[emf->nz-1][i2][i1]-intz[emf->nz-2][i2][i1]);
	}else{
	  s = (emf->x3[i3] - emf->oz - k3*emf->dz)/emf->dz;//fraction after integer part
	  intb = intz[k3-1][i2][i1] + s*(intz[k3][i2][i1]-intz[k3-1][i2][i1]);
	}
	
	s = (inta - intb)/(emf->x3[i3+1] - emf->x3[i3]);//volume average
	if(emf->addair && emf->x3s[i3]<emf->oz){
	  if(flag==1) s = 1e-8;
	  if(flag==2) s = 1e8;
	  if(flag==3) s = log(1e-8);	  
	}
	if(flag==1) sigma_out[i3][i2][i1] = s;
	if(flag==2) sigma_out[i3][i2][i1] = 1./s;
	if(flag==3) sigma_out[i3][i2][i1] = exp(s);	
      }//end for i1
    }//end for i2
  }//end for i3

  free3float(intz);
}

