# 3D CSEM modelling and inversion using multigrid frequency-domain solver

Author: Pengliang Yang

Department of Mathematics, Harbin Institute of Technology, China

E-mail: ypl.2100@gmail.com

Note: This modelling code requires installation of MPI (either mpich or openmpi),
as it has a parallel computing idea in mind since its design, to make it work 
for inversion simultaneously with multiple sources.

----------------------------------------------------

Instructions to run:

1. go into /src to compile the code using makefile:
make

Then the compiled executable will be placed at /bin

2. test the code in the modelling template /run_modelling
bash run.sh

All relevant parameters should be specified in run.sh, while the binary input 
resistivity should be prepared in advance.

Here I created a small Fortran code to generate acquisition files
sources.txt, receivers.txt and a source-receiver connection table src_rec_table.txt.

The input resistivity was created by another Fortran code as binaries.
You should prepare your own!

3. Test and compare with emg3d. I assume you have looked at the emg3d website and 
have installed it properly. (Otherwise, please consult Dieter for the help)
 Then, you can try:
 
python3 emg3d_test.py
Here, the input resistivity is in ASCII format, which will be created by my Fortran code.
This aims to maintain consistency between the model you input for my code and emg3d.

4. reproduce the 3D CSEM inversion experiment in /run_inversion.
First you need to create observed data using the true resistivity model frho.
Modify run.sh to set "mode=0 frho11=frho frho22=frho frho33=frho" and run "bash run.sh".

Then you can run inversion from a crude initial model frho_init.
Modify run.sh to set "mode=1 frho11=frho_init frho22=frho_init frho33=frho_init" and run "bash run.sh".


## Acknowledgements

Pengliang Yang is indebted to Dieter Werthmüller for responsive feedback in last two years to develop the
geometrical multigrid code using C programming language. Without his help and the Python code emg3d,
it is not possible to complete libEMMI_MGFD to share with the community.

## Bibliography

1. Pengliang Yang and An Ping, libEMMI MGFD: A program of marine controlled-source electromagnetic
modelling and inversion using frequency-domain multigrid solver, Computer Physics Communications 2024 305:109327
[doi:10.1016/j.cpc.2024.109327](https://doi.org/10.1016/j.cpc.2024.109327)
