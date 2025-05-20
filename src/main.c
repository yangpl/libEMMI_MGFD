/* CSEM modelling and inversion in frequency domain
 *----------------------------------------------------------------------
 *   Copyright (c) 2020, Harbin Institute of Technology, China
 *   Author: Pengliang Yang
 *   E-mail: ypl.2100@gmail.com
 *   Homepage: https://yangpl.wordpress.com
 *----------------------------------------------------------------------*/
#include <mpi.h>
#include "cstd.h"
#include "emf.h"
#include "acq.h"
 

int iproc, nproc, ierr;

void emf_init(emf_t *emf);
void emf_free(emf_t *emf);

void acq_init(acq_t *acq, emf_t * emf);
void acq_free(acq_t *acq);

void do_modelling(acq_t *acq, emf_t *emf);
void do_fwi(acq_t *acq, emf_t *emf);

int main(int argc, char **argv)
{
  emf_t *emf;
  acq_t *acq;
  char current_time[128];
  time_t      t;
  struct tm*  ptm;

  MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
  initargs(argc,argv);

  acq = malloc(sizeof(acq_t));
  emf = malloc(sizeof(emf_t));
  if(!getparint("mode", &emf->mode)) emf->mode = 0; 
  if(!getparint("verb", &emf->verb)) emf->verb = (iproc==0)?1:0;
  if(emf->verb){
    t = time(NULL);
    ptm = localtime(&t);
    strftime(current_time, 128, "%d-%b-%Y %H:%M:%S", ptm);
    printf("  Current date and time: %s\n", current_time);
    printf("=====================================================\n");
    printf("    Welcome to EM modeling and inversion code        \n");
    printf("    -------------lib_EMMI_MGFD---------------        \n");
    printf("            Author: Pengliang Yang                   \n");
    printf("            E-mail: ypl.2100@gmail.com               \n");
    printf("    Copyright (c) 2023. All rights reserved.         \n");
    printf("=====================================================\n");
    if(emf->mode==0)
      printf("Task: Electromagnetic forward modeling \n");
    else if(emf->mode==1)
      printf("Task: Inversion for medium parameters \n");
    else 
      printf("Task: Compute CSEM inversion gradient\n");
    printf("=====================================================\n");
  }
  
  emf_init(emf);
  acq_init(acq, emf);

  if(emf->mode==0)
    do_modelling(acq, emf);
  else
    do_fwi(acq, emf);

  acq_free(acq);
  emf_free(emf);

  free(emf);
  free(acq);
  
  MPI_Finalize();
  return 0;
}

