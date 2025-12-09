/* electromagnetic field (EMF) initialization
 *------------------------------------------------------------------------
 *
 * Copyright (c) 2020-2022 Harbin Institute of Technology. All rights reserved.
 * Author: Pengliang Yang 
 * Email: ypl.2100@gmail.com
 * Homepage: https://yangpl.wordpress.com
 *-----------------------------------------------------------------------*/
#include "cstd.h"
#include "emf.h"

int find_good_size(int n);
int cmpfunc(const void *a, const void *b) { return ( *(int*)a - *(int*)b ); }

void emf_init(emf_t *emf)
{
  int ifreq, ic, istat;
  FILE *fp;
  char *frho11, *frho22, *frho33;

  if(iproc==0) printf("---------- emf_init -----------\n");
  if(!(emf->nfreq=countparval("freqs"))) err("Need freqs= vector");
  emf->freqs = alloc1float(emf->nfreq);
  emf->omegas = alloc1float(emf->nfreq);
  getparfloat("freqs", emf->freqs);/* a list of frequencies separated by comma */
  qsort(emf->freqs, emf->nfreq, sizeof(float), cmpfunc);/*sort frequencies in ascending order*/
  for(ifreq = 0; ifreq<emf->nfreq; ++ifreq) {
    emf->omegas[ifreq] = 2.*PI*emf->freqs[ifreq];
    if(emf->verb) printf("freq[%d]=%g \n", ifreq, emf->freqs[ifreq]);
  }
  /*-------- read active source channels ---------*/
  if((emf->nchsrc = countparval("chsrc"))!= 0) {
    emf->chsrc = (char**)alloc1(emf->nchsrc, sizeof(void*));
    getparstringarray("chsrc", emf->chsrc);
    /* active source channels: Ex, Ey, Ez, Hx, Hy, Hz or their combinations */
  }else{
    emf->nchsrc = 1;
    emf->chsrc = (char**)alloc1(emf->nchsrc, sizeof(void*));
    emf->chsrc[0] = "Ex";
  }
  /*-------- read active receiver channels ---------*/
  if((emf->nchrec = countparval("chrec"))!= 0) {
    emf->chrec = (char**)alloc1(emf->nchrec, sizeof(void*));
    getparstringarray("chrec", emf->chrec);
    /* active receiver channels: Ex, Ey, Ez, Hx, Hy, Hz or their combinations */
  }else{
    emf->nchrec = 1;
    emf->chrec = (char**)alloc1(emf->nchrec, sizeof(void*));
    emf->chrec[0]  = "Ex";
  }
  if(emf->verb){
    printf("Active source channels:");
    for(ic=0; ic<emf->nchsrc; ++ic) printf(" %s", emf->chsrc[ic]);
    printf("\n");
    printf("Active recever channels:");
    for(ic=0; ic<emf->nchrec; ++ic) printf(" %s", emf->chrec[ic]);
    printf("\n");
  }
  if(!getparint("reciprocity", &emf->reciprocity)) emf->reciprocity = 0;
  //dimensions and grid spacing over densely sampled uniform grid
  if(!getparint("nx", &emf->nx)) err("Need nx= ");
  if(!getparint("ny", &emf->ny)) err("Need ny= "); 
  if(!getparint("nz", &emf->nz)) err("Need nz= "); 
  if(!getparfloat("dx", &emf->dx)) err("Need dx= ");
  if(!getparfloat("dy", &emf->dy)) err("Need dy= "); 
  if(!getparfloat("dz", &emf->dz)) err("Need dz= "); 
  if(!getparfloat("ox", &emf->ox)) err("Need ox= ");
  if(!getparfloat("oy", &emf->oy)) err("Need oy= ");
  if(!getparfloat("oz", &emf->oz)) emf->oz = 0;//default, seafloor depth=0 m
  if(!getparfloat("d1min", &emf->d1min)) emf->d1min = 100;
  if(!getparfloat("d2min", &emf->d2min)) emf->d2min = 100;
  if(!getparfloat("d3min", &emf->d3min)) emf->d3min = 20;
  if(!getparint("n1", &emf->n1)) emf->n1 = 128;//make sure it is good for multigrid
  if(!getparint("n2", &emf->n2)) emf->n2 = 128;//make sure it is good for multigrid
  if(!getparint("n3", &emf->n3)) emf->n3 = 128;//make sure it is good for multigrid
  if(emf->verb) {
    printf("[nx, ny, nz]=[%d, %d, %d]\n", emf->nx, emf->ny, emf->nz);
    printf("[dx, dy, dz]=[%g, %g, %g]\n", emf->dx, emf->dy, emf->dz);
    printf("[ox, oy, oz]=[%g, %g, %g]\n", emf->ox, emf->oy, emf->oz);
    printf("[d1min, d2min, d3min]=[%g, %g, %g]\n", emf->d1min, emf->d2min, emf->d3min);
    printf("GMG grid size: [n1, n2, n3]=[%d, %d, %d]\n", emf->n1, emf->n2, emf->n3);
  }
 
  //resistivities 
  if(!(getparstring("frho11", &frho11))) err("Need frho11= ");
  if(!(getparstring("frho22", &frho22))) err("Need frho22= ");
  if(!(getparstring("frho33", &frho33))) err("Need frho33= ");

  emf->rho11 = alloc3float(emf->nx, emf->ny, emf->nz);
  emf->rho22 = alloc3float(emf->nx, emf->ny, emf->nz);
  emf->rho33 = alloc3float(emf->nx, emf->ny, emf->nz);

  fp = fopen(frho11, "rb");
  if(fp==NULL) err("cannot open file frho11=%s\n", frho11);
  istat = fread(&emf->rho11[0][0][0], sizeof(float), emf->nx*emf->ny*emf->nz, fp);
  if(istat != emf->nx*emf->ny*emf->nz) 
    err("size not match frho11: file=%d nx*ny*nz=%d\n", istat, emf->nx*emf->ny*emf->nz);
  fclose(fp);

  fp = fopen(frho22, "rb");
  if(fp==NULL) err("cannot open file frho22=%s\n", frho22);
  istat = fread(&emf->rho22[0][0][0], sizeof(float), emf->nx*emf->ny*emf->nz, fp);
  if(istat != emf->nx*emf->ny*emf->nz) 
    err("size not match frho22: file=%d nx*ny*nz=%d\n", istat, emf->nx*emf->ny*emf->nz);
  fclose(fp);

  fp = fopen(frho33, "rb");
  if(fp==NULL) err("cannot open file frho33=%s\n", frho33);
  istat = fread(&emf->rho33[0][0][0], sizeof(float), emf->nx*emf->ny*emf->nz, fp);
  if(istat != emf->nx*emf->ny*emf->nz) 
    err("size not match frho33: file=%d nx*ny*nz=%d\n", istat, emf->nx*emf->ny*emf->nz);
  fclose(fp);
 
  if(!getparfloat("lextend", &emf->lextend)) emf->lextend = 50e3;//extend 50 km on each side
  if(!getparfloat("rho_air", &emf->rho_air)) emf->rho_air = 1e8;
  if(!getparint("nb", &emf->nb)) emf->nb = 10;//nb layers on top of sea for air
  int i1, i2, i3;
  float vmax = emf->rho11[0][0][0];
  float vmin = emf->rho11[0][0][0];
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
  emf->rho_water = vmin;
  if(emf->verb){
    printf("input model [rhomin, rhomax]=[%g, %g]\n", vmin, vmax);
    printf("rho_air=%e\n", emf->rho_air);
    printf("rho_water=%g\n", emf->rho_water);
    printf("extended distance on each side: lextend=%g\n", emf->lextend);
  }

}


void emf_free(emf_t *emf)
{
  free3float(emf->rho11);
  free3float(emf->rho22);
  free3float(emf->rho33);
}
