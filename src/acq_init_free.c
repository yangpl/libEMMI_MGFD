/* read acquisition file for source and receiver geometry
 *----------------------------------------------------------------
 *
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *---------------------------------------------------------------*/
#include "cstd.h"
#include "acq.h"
#include "emf.h"
 
 

/*< read acquisition file to initialize survey geometry >*/
void acq_init(acq_t *acq, emf_t * emf)
{
  static int nd = 5000;//maximum dimensions for the number of source and receiver
  float src_x1[nd], src_x2[nd], src_x3[nd], src_hd[nd], src_pit[nd];/* source receiver coordinates */
  float rec_x1[nd], rec_x2[nd], rec_x3[nd], rec_hd[nd], rec_pit[nd];/* source receiver coordinates */
  int rec_idx[nd];/* reciver index associated with current processor */
  float x, y, z, hd, pit;
  float x1min, x1max, x2min, x2max, x3min, x3max;
  int isrc, irec, iseof, idx, i, nsrc;
  char *fsrc, *frec, *fsrcrec;
  FILE *fp=NULL;

  if(!(getparstring("fsrc", &fsrc))) err("Need fsrc= ");
  /* file to specify all possible source locations */
  if(!(getparstring("frec", &frec))) err("Need frec= ");
  /* file to specify all possible receiver locations */
  if(!(getparstring("fsrcrec", &fsrcrec))) err("Need fsrcrec= ");
  /* file to specify how source and receiver are combined */

  if(!getparfloat("x1min", &acq->x1min)) acq->x1min = emf->ox;
  /* minimum limit of the survey in x direction */
  if(!getparfloat("x1max", &acq->x1max)) acq->x1max = emf->ox+emf->nx*emf->dx;
  /* maximum limit of the survey in x direction */
  if(!getparfloat("x2min", &acq->x2min)) acq->x2min = emf->oy;
  /* minimum limit of the survey in y direction */
  if(!getparfloat("x2max", &acq->x2max)) acq->x2max = emf->oy+emf->ny*emf->dy;
  /* maximum limit of the survey in y direction */
  if(!getparfloat("x3min", &acq->x3min)) acq->x3min = emf->oz;
  /* minimum limit of the survey in z direction */
  if(!getparfloat("x3max", &acq->x3max)) acq->x3max = emf->oz+emf->nz*emf->dz;
  /* maximum limit of the survey in z direction */
  
  if(emf->ox>acq->x1min || emf->ox + emf->dx*emf->nx < acq->x1max)
    err("x - sources/receivers from acquisition file are out of domain");
  if(emf->oy>acq->x2min || emf->oy + emf->dy*emf->ny < acq->x2max)
    err("y - sources/receivers from acquisition file are out of domain");
  if(emf->oz>acq->x3min || emf->oz + emf->dz*emf->nz < acq->x3max)
    err("z - sources/receivers from acquisition file are out of domain");
  
  if(!getparint("nsubsrc", &acq->nsubsrc)) acq->nsubsrc = 1;
  /* number of subpoints to represent one source location */
  if(!getparint("nsubrec", &acq->nsubrec)) acq->nsubrec = 1;
  /* number of subpoints to represent one receiver location */
  if(!getparfloat("lensrc", &acq->lensrc)) acq->lensrc = 275.;
  /* length of the source antenna, default = 1 m */
  if(!getparfloat("lenrec", &acq->lenrec)) acq->lenrec = 8.;
  /* length of the receiver antenna, default = 8 m */
  if(fabs(acq->x1max-acq->x1min-emf->dx*emf->nx)>1e-15)
    err("inconsistent input: x1max-x1min!=dx*nx");
  if(fabs(acq->x2max-acq->x2min-emf->dy*emf->ny)>1e-15)
    err("inconsistent input: x2max-x2min!=dy*ny");
  
  acq->shot_idx = alloc1int(nproc);
  nsrc = countparval("shots");
  if(nsrc>0){
    if( nsrc<nproc) err("nproc > number of shot indices! ");
    getparint("shots", acq->shot_idx);/* a list of source index separated by comma */
  }
  if(nsrc==0){
    for(i=0; i<nproc; i++) acq->shot_idx[i] = i+1;//index starts from 1
  }


  idx = acq->shot_idx[iproc];

  //===========================================================================
  //step 0. establish source-receiver connection table considering reciprocity
  //==========================================================================
  fp = fopen(fsrcrec,"r");
  if(fp==NULL) err("file fsrcrec= missing!");
  fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
  i = 0;
  while(1){
    /* (northing,easting,depth)=(y,x,z); azimuth = heading; dip=pitch */
    iseof=fscanf(fp,"%d %d", &isrc, &irec);
    if(iseof==EOF)
      break;
    else{
      if(emf->reciprocity){
  	if(irec==idx) {//irec-th common receiver gather
  	  rec_idx[i] = isrc;//the global source index associated with current receiver
  	  i++;
  	}
      }else{
  	if(isrc==idx) {//isrc-th source gather
  	  rec_idx[i] = irec;//the global receiver index associated with current source
  	  i++;
  	}
      }
    }
  }
  acq->nrec = i;
  fclose(fp);


  /*============================================*/
  /* step 1: read all possible source locations */
  /*============================================*/
  x1min = acq->x1min;
  x1max = acq->x1max;
  x2min = acq->x2min;
  x2max = acq->x2max;
  x3min = acq->x3min;
  x3max = acq->x3max;
  fp = fopen(fsrc,"r");
  if(fp==NULL) err("file fsrc= missing"); 
  fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
  isrc = 0;
  while(1){
    /* (northing,easting,depth)=(y,x,z);   azimuth = heading;  dip=pitch */
    iseof=fscanf(fp,"%f %f %f %f %f %d", &x, &y, &z, &hd, &pit, &idx);
    if(iseof==EOF)
      break;
    else{
      src_x1[isrc] = x;
      src_x2[isrc] = y;
      src_x3[isrc] = z;
      src_hd[isrc] = hd;
      src_pit[isrc] = pit;

      x1min = MIN(x1min, src_x1[isrc]);
      x1max = MAX(x1max, src_x1[isrc]);
      x2min = MIN(x2min, src_x2[isrc]);
      x2max = MAX(x2max, src_x2[isrc]);
      x3min = MIN(x3min, src_x3[isrc]);
      x3max = MAX(x3max, src_x3[isrc]);

      isrc++;
    }
  }
  acq->nsrc_total = isrc;
  fclose(fp);
  if(x1min<acq->x1min) err("source location: x<x1min");
  if(x2min<acq->x2min) err("source location: y<x2min");
  if(x3min<acq->x3min) err("source location: z<x3min");
  if(x1max>acq->x1max) err("source location: x>x1max");
  if(x2max>acq->x2max) err("source location: y>x2max");
  if(x3max>acq->x3max) err("source location: z>x3max");

    
  /*==============================================*/
  /* step 2: read all possible receiver locations */
  /*==============================================*/
  fp = fopen(frec,"r");
  if(fp==NULL) err("file frec= missing"); 
  fscanf(fp, "%*[^\n]\n");//skip a line at the beginning of the file
  irec = 0;
  while(1){
    /* (northing,easting,depth)=(y,x,z);   azimuth = heading;  dip=pitch */
    iseof=fscanf(fp,"%f %f %f %f %f %d", &x, &y, &z, &hd, &pit, &idx);
    if(iseof==EOF)
      break;
    else{
      rec_x1[irec] = x;
      rec_x2[irec] = y;
      rec_x3[irec] = z;
      rec_hd[irec] = hd;
      rec_pit[irec] = pit;
      
      x1min = MIN(x1min, rec_x1[irec]);
      x1max = MAX(x1max, rec_x1[irec]);
      x2min = MIN(x2min, rec_x2[irec]);
      x2max = MAX(x2max, rec_x2[irec]);
      x3min = MIN(x3min, rec_x3[irec]);
      x3max = MAX(x3max, rec_x3[irec]);
      
      irec++;
    }
  }
  acq->nrec_total = irec;
  fclose(fp);
  if(x1min<acq->x1min) err("receiver location: x<x1min");
  if(x2min<acq->x2min) err("receiver location: y<x2min");
  if(x3min<acq->x3min) err("receiver location: z<x3min");
  if(x1max>acq->x1max) err("receiver location: x>x1max");
  if(x2max>acq->x2max) err("receiver location: y>x2max");
  if(x3max>acq->x3max) err("receiver location: z>x3max");

  //-----------------------------------------------------------
  acq->nsrc = 1; /* assume 1 source per process by default */
  acq->src_x1 = alloc1float(acq->nsrc);
  acq->src_x2 = alloc1float(acq->nsrc);
  acq->src_x3 = alloc1float(acq->nsrc);
  acq->src_azimuth = alloc1float(acq->nsrc);
  acq->src_dip = alloc1float(acq->nsrc);
  for(isrc=0; isrc< acq->nsrc; ++isrc){
    i = acq->shot_idx[iproc]-1;//index starts from 1
    if(emf->reciprocity){
      acq->src_x1[isrc] = rec_x1[i];
      acq->src_x2[isrc] = rec_x2[i];
      acq->src_x3[isrc] = rec_x3[i];
      acq->src_azimuth[isrc] = PI*rec_hd[i]/180.;
      acq->src_dip[isrc] = PI*rec_pit[i]/180.;
    }else{
      acq->src_x1[isrc] = src_x1[i];
      acq->src_x2[isrc] = src_x2[i];
      acq->src_x3[isrc] = src_x3[i];
      acq->src_azimuth[isrc] = PI*src_hd[i]/180.;
      acq->src_dip[isrc] = PI*src_pit[i]/180.;
    }
  }/* end for isrc */

  //-----------------------------------------------------------
  acq->rec_x1 = alloc1float(acq->nrec);
  acq->rec_x2 = alloc1float(acq->nrec);
  acq->rec_x3 = alloc1float(acq->nrec);
  acq->rec_azimuth = alloc1float(acq->nrec);
  acq->rec_dip = alloc1float(acq->nrec);
  for(irec=0; irec<acq->nrec; ++irec){//we always have: acq->nrec <= acq->nrec_total
    //nrec < nrec_total if only inline data are used
    //idx=index of the receivers associated with current source or common receiver gather
    i = rec_idx[irec]-1; //index starts from 1
    if(emf->reciprocity){
      acq->rec_x1[irec] = src_x1[i];
      acq->rec_x2[irec] = src_x2[i];
      acq->rec_x3[irec] = src_x3[i];
      acq->rec_azimuth[irec] = PI*src_hd[i]/180.;
      acq->rec_dip[irec] = PI*src_pit[i]/180.;
    }else{
      acq->rec_x1[irec] = rec_x1[i];
      acq->rec_x2[irec] = rec_x2[i];
      acq->rec_x3[irec] = rec_x3[i];
      acq->rec_azimuth[irec] = PI*rec_hd[i]/180.;
      acq->rec_dip[irec] = PI*rec_pit[i]/180.;
    }
    
  }/* end for irec */
  printf("****** isrc=%d nrec=%d, [x,y,z,azimuth,dip]=[%g, %g, %g, %g, %g]\n",
	 acq->shot_idx[iproc], acq->nrec, acq->src_x1[0], acq->src_x2[0],
	 acq->src_x3[0], 180*acq->src_azimuth[0]/PI, 180.*acq->src_dip[0]/PI);
  if(emf->verb){
    printf("lensrc=%g, nsubsrc=%d\n", acq->lensrc, acq->nsubsrc);
    printf("lenrec=%g, nsubrec=%d\n", acq->lenrec, acq->nsubrec);
    printf("original domain [x1min, x1max]=[%g, %g]\n", acq->x1min, acq->x1max);
    printf("original domain [x2min, x2max]=[%g, %g]\n", acq->x2min, acq->x2max);
    printf("original domain [x3min, x3max]=[%g, %g]\n", acq->x3min, acq->x3max);
    printf("nsrc_total=%d\n", acq->nsrc_total);
    printf("nrec_total=%d\n", acq->nrec_total);
  }
}

/*< free the allocated variables for acquisition >*/
void acq_free(acq_t *acq)
{
  free1float(acq->src_x1);
  free1float(acq->src_x2);
  free1float(acq->src_x3);
  free1float(acq->src_azimuth);
  free1float(acq->src_dip);

  free1float(acq->rec_x1);
  free1float(acq->rec_x2);
  free1float(acq->rec_x3);
  free1float(acq->rec_azimuth);
  free1float(acq->rec_dip);

  free1int(acq->shot_idx);
}

