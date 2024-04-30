#!/bin/bash

mpirun -n 25 ../bin/main mode=1 \
       istretch=0 \
       addair=1 \
       freqs=0.25,1,2.75 \
       chsrc=Ex \
       chrec=Ex ,Ey,Hx,Hy \
       nx=100 \
       ny=100 \
       nz=100 \
       dx=200 \
       dy=200 \
       dz=40 \
       ox=-10000 \
       oy=-10000 \
       oz=0 \
       fbathy=fbathy \
       frho11=frho_init \
       frho22=frho_init \
       frho33=frho_init \
       fsrc=sources.txt \
       frec=receivers.txt \
       fsrcrec=src_rec_table.txt \
       niter=30 \
       preco=2 \
       npar=2 \
       bound=1 \
       idxpar=1,2 \
       minpar=1.0,1.0 \
       maxpar=100.0,100.0 \
       dsmute=500 \
       drmute=250 \
       gamma1=100 \
       gamma2=0 >out.txt&
       

