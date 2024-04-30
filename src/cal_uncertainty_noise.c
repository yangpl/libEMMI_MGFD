/* calculate uncertainty and add Gaussian white noise
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
#include "fwi.h"

float frannor(void);//random noise generator

/*<calculate uncertainty and add Gaussian white noise>*/
float cal_uncertainty_noise(acq_t *acq, emf_t *emf)
{
  int ifreq, irec, ichrec;
  float phi, mse, s1, s2, s3;
  _Complex float cgwn;
  _Complex float Ex, Ey, Hx, Hy;
  _Complex float Exc, Eyc, Hxc, Hyc;

  memset(&emf->delta_emf[0][0][0], 0, emf->nchrec*emf->nfreq*acq->nrec*sizeof(float));
  mse = 0;  
  for(irec=0; irec<acq->nrec; irec++){
    //compute azimuth angle between source towline and receiver direction
    phi = acq->rec_azimuth[irec]-acq->src_azimuth[0];
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      //=======================================================
      //step 1: compute uncertainties
      //=======================================================
      Ex = 0;//Ex before rotation
      Ey = 0;//Ey before rotation
      Hx = 0;//Hx before rotation
      Hy = 0;//Hy before rotation
      Exc = 0;//Ex after rotation
      Eyc = 0;//Ey after rotation
      Hxc = 0;//Hx after rotation
      Hyc = 0;//Hy after rotation
      for(ichrec=0; ichrec<emf->nchrec; ichrec++){
	if(strcmp(emf->chrec[ichrec],"Ex")==0){
	  Ex = emf->dobs_fd[ichrec][ifreq][irec];
	  Exc += Ex*cos(phi);
	  Eyc += Ex*sin(phi);
	}else if(strcmp(emf->chrec[ichrec],"Ey")==0){
	  Ey = emf->dobs_fd[ichrec][ifreq][irec];
	  Exc -= Ey*sin(phi);
	  Eyc += Ey*cos(phi);
	}else if(strcmp(emf->chrec[ichrec],"Hx")==0){
	  Hx = emf->dobs_fd[ichrec][ifreq][irec];
	  Hxc += Hx*cos(phi);
	  Hyc += Hx*sin(phi);
	}else if(strcmp(emf->chrec[ichrec],"Hy")==0){
	  Hy = emf->dobs_fd[ichrec][ifreq][irec];
	  Hxc -= Hy*sin(phi);
	  Hyc += Hy*cos(phi);
	}
      }
      for(ichrec=0; ichrec<emf->nchrec; ichrec++){
	if(strcmp(emf->chrec[ichrec],"Ex")==0){
	  s1 = cabs(emf->amp_perc*Ex)*cos(phi);
	  s2 = cabs(emf->amp_perc*Ey)*sin(phi);
	  s3 = emf->delta_phi*cabs(Eyc);
	  emf->delta_emf[ichrec][ifreq][irec] = sqrt(s1*s1 + s2*s2 + s3*s3 + emf->noisefloorE*emf->noisefloorE);
	}else if(strcmp(emf->chrec[ichrec],"Ey")==0){
	  s1 = cabs(emf->amp_perc*Ex)*sin(phi);
	  s2 = cabs(emf->amp_perc*Ey)*cos(phi);
	  s3 = emf->delta_phi*cabs(Exc);
	  emf->delta_emf[ichrec][ifreq][irec] = sqrt(s1*s1 + s2*s2 + s3*s3 + emf->noisefloorE*emf->noisefloorE);
	}else if(strcmp(emf->chrec[ichrec],"Hx")==0){
	  s1 = cabs(emf->amp_perc*Hx)*cos(phi);
	  s2 = cabs(emf->amp_perc*Hy)*sin(phi);
	  s3 = emf->delta_phi*cabs(Hyc);
	  emf->delta_emf[ichrec][ifreq][irec] = sqrt(s1*s1 + s2*s2 + s3*s3 + emf->noisefloorH*emf->noisefloorH);
	}else if(strcmp(emf->chrec[ichrec],"Hy")==0){
	  s1 = cabs(emf->amp_perc*Hx)*sin(phi);
	  s2 = cabs(emf->amp_perc*Hy)*cos(phi);
	  s3 = cabs(emf->delta_phi*Hxc);
	  emf->delta_emf[ichrec][ifreq][irec] = sqrt(s1*s1 + s2*s2 + s3*s3 + emf->noisefloorH*emf->noisefloorH);
	}
	
	//=======================================================
	//step 2: add Gaussian white noise
	//=======================================================
	if(emf->addnoise){
	  //frannor() compute Gaussian white noise (GWN) following N(0,1) distribution
	  cgwn = (frannor() + I*frannor())/sqrt(2.);//complex-valued GWN
	  emf->dobs_fd[ichrec][ifreq][irec] += emf->delta_emf[ichrec][ifreq][irec]*cgwn;
	  mse += conj(cgwn)*cgwn;//accumulate mse value
	}//end if
      }	//end for ichrec
    }//end for ifreq
  }//end for irec
  
  return mse;//return mse on each source
}
