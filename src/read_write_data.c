/* read and write CSEM data in ASCII format
 *----------------------------------------------------------------------
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *----------------------------------------------------------------------*/
#include <unistd.h>
#include "cstd.h"
#include "acq.h"
#include "emf.h"
 

void write_data(acq_t *acq, emf_t *emf, char *fname, _Complex float ***dcal_fd)
/*< write synthetic data according to shot/process index >*/
{
  FILE *fp;
  int isrc, irec, ichrec, ifreq;
  float dp_re, dp_im;

  fp=fopen(fname,"w");
  if(fp==NULL) err("error opening file for writing");
  fprintf(fp, "iTx 	iRx    chrec  frequency/Hz    Real{E/H}     Imag{E/H}\n");
  isrc = acq->shot_idx[iproc];//index starts from 1
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      for(irec=0; irec<acq->nrec; irec++){
	dp_re = creal(dcal_fd[ichrec][ifreq][irec]);
	dp_im = cimag(dcal_fd[ichrec][ifreq][irec]);
	fprintf(fp, "%d \t %d \t %s \t %g \t %e \t %e\n",
		isrc, irec+1, emf->chrec[ichrec], emf->freqs[ifreq], dp_re, dp_im);
      }
    }
  }
  fclose(fp);
}


void read_data(acq_t *acq, emf_t *emf)
/*< read observed data according to shot/process index >*/
{
  char fname[sizeof("emf_0000.txt")];
  int isrc, ichrec, irec, ifreq, iseof;
  float dp_re, dp_im, freq;
  char chrec[sizeof("Ex")];
  FILE *fp;

  sprintf(fname, "emf_%04d.txt", acq->shot_idx[iproc]);
  fp=fopen(fname,"r");
  if(fp==NULL) err("error opening file for reading");
  fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
  while(1){
    iseof=fscanf(fp, "%d \t %d \t %s \t %g \t %e \t %e\n",
		 &isrc, &irec, chrec, &freq, &dp_re, &dp_im);
    if(iseof==EOF){
      break;
    }else{
      if(acq->shot_idx[iproc]==isrc && irec<=acq->nrec){
	for(ichrec=0; ichrec<emf->nchrec; ichrec++){
	  if(strcmp(emf->chrec[ichrec], chrec)==0){
	    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
		if(fabs(emf->freqs[ifreq]-freq)<1e-4){
		  emf->dobs_fd[ichrec][ifreq][irec-1] = dp_re + I* dp_im;
		  break;
		}//end if
	    }//end for ifreq

	    break;
	  }//end if
	}//end for

      }//end if
    }//end if 
  }//end while
  fclose(fp);    
}


/*< write synthetic data according to shot/process index >*/
void write_rmse(acq_t *acq, emf_t *emf)
{
  int isrc, irec, ichrec, ifreq;
  FILE *fp;

  char fname[sizeof("rmse_0000.txt")];
  sprintf(fname, "rmse_%04d.txt", acq->shot_idx[iproc]);
  fp=fopen(fname,"w");
  if(fp==NULL) err("error opening file for writing");
  fprintf(fp, "iTx iRx chrec frequency/Hz RMSE x y z\n");
  isrc = acq->shot_idx[iproc];//index starts from 1
  for(ichrec=0; ichrec<emf->nchrec; ichrec++){
    for(ifreq=0; ifreq<emf->nfreq; ifreq++){
      for(irec=0; irec<acq->nrec; irec++){
	fprintf(fp, "%d %d %s %g %e %g %g %g\n", isrc, irec+1,
		emf->chrec[ichrec], emf->freqs[ifreq], emf->rmse[ichrec][ifreq][irec],
		acq->rec_x1[irec], acq->rec_x2[irec], acq->rec_x3[irec]);		
      }//end for irec
    }//end for ifreq
  }//end for ichrec
  fclose(fp);

}

