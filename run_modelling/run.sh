#!/bin/bash

#NB: the parameter input follows the convention of comma separated values!

#To run this test, you first need:
#gfortran create_model.f90; 
#./a.out (it will create resistivity model)
#gfortran create_acquisition.f90
#./a.out (it will create source-receiver acqusition files)
#then launch the test by: bash run.sh

#To mute grid stretching for modelling using your own resistivity, set istretch=0;
#If you want to look at the memory and runtime on different size of the model, set istretch=1;
#In addition, you can set n1=100, n2=100, n3=100 (try to change them), the code will do 
#automatic regridding to find a proper size for GMG modelling.

#To output 3D EM fields after modelling, set icheck=1;
#You can visualize them also by python3 plot_tensormesh_3d.py!

#To model with mutiple frequencies, set freqs=0.25,0.75

#To record EM fields from multiple components, set chrec=Ex,Ey,Hx,Hy

#The azimuth and dip angles are set in source and receiver files, check out them!

mpirun -n 1 time -v ../bin/main mode=0 \
     icheck=1 \
     istretch=1 \
     freqs=0.25 ,0.75 \
     chsrc=Ex \
     chrec=Ex,Ey,Hx,Hy \
     nx=200 \
     ny=200 \
     nz=100 \
     dx=100 \
     dy=100 \
     dz=40 \
     ox=-10000 \
     oy=-10000 \
     oz=0 \
     n1=160 \
     n2=160 \
     n3=160 \
     frho11=frho \
     frho22=frho \
     frho33=frho \
     fsrc=sources.txt \
     frec=receivers.txt \
     fsrcrec=src_rec_table.txt \
     

