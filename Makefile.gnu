#--------------------------------------------------------------------------
#
# This Gmake file will compile the PARAMESH library and create a
# set of library files to which you can link. To use it, make sure
# it is in the PARAMESH root directory.  
# It works by running gmake on the Makefile.gnu files which appear
# in the headers, source and mpi_source sub-directories.
# To simply create these PARAMESH library files, type
#     gmake -f Makefile.gnu
# when in the PARAMESH root directory. The library files will
# then be found in a newly created sub-directory called libs.
#
# If you type 
#     gmake -f Makefile.gnu Tests
# it will also compile and link the test programs in the Tests
# sub-directory. There is a file called Makefile.gnu inside Tests
# which is used.
# 
# To compile and link application files in a sub-directory called
# User_applic you could type
#     gmake -f Makefile.gnu User_applic
# provided you copy Makefile.gnu from Tests to User_applic, and modify
# it appropriately.
#
#
# Written : Ernest Mamikonyan        April 2002.
#
#--------------------------------------------------------------------------
export ARFLAGS = rsv

export cur-dir := $(shell pwd)
# Set the location of the paramesh top directory
export paramesh_dir = $(cur-dir)

# Define the fortran compiler
#export FC = /fserver/kevino/mpich-3.2.1/mpich/bin/mpif90
#export CC = /fserver/kevino/mpich-3.2.1/mpich/bin/mpicc
#export FC = mpif90
#export CC = mpicc
#export FC=ftn
#export FC = f90
#export FC = mpif90 -fc=nagfor
#export FC = mpif90 -fc=ifort 
#export CC = mpicc -cc=icc
#export FC = /home/kevin/SW/MPICH/bin/mpif90 -fc=nvfortran
#export CC = /home/kevin/SW/MPICH/bin/mpicc  -cc=nvcc
export FC = mpif90
export CC = mpicc 

#-----------------------------------------------
 
# Set the desired compilation flags

#XLF compiler on mac

#export FFLAGS = -O2 -qintsize=4 -qrealsize=8 -qfixed=72 -I../headers 
#export FFLAGS = -O2 -qintsize=4 -qrealsize=8 -I../headers -qsuffix=f=F90:cpp=F90

# INTEL compiler
##export FFLAGS = -O2 -qopenmp -no-vec -i4 -r8 -ftrapuv -traceback -I../headers
##export CFLAGS = -O2 -I../headers -I/Users/olson/Documents/SW/mpich-gfortran-gcc/HDF5/hdf5_intel/include
#export FFLAGS = -O2 -i4 -r8 -I../headers 
#export CFLAGS = -O2 -I../headers -I/fserver/kevino/hdf5-1.6.5/hdf5/include 

# gfortran
#export FFLAGS = -O3 -fopenmp -fdefault-real-8 -fdefault-double-8 -I../headers -I../headers -I/fserver/kevino/hdf5-1.6.5/hdf5/include 
# WITH MPICH
#export FFLAGS = -O3 -fdefault-real-8 -I../headers -I../headers -I/home/kevin/SW/HDF5_1.6/include -fPIE
#export CFLAGS = -O2 -I../headers -I/usr/local/include -I/home/kevin/SW/HDF5_1.6/include -fPIE
# Nvidia Fortran
export FFLAGS = -O3 -r8 -I../headers -I../headers -I/home/kevin/SW/HDF5_1.6/include 
export CFLAGS = -O2 -I../headers -I/usr/local/include -I/home/kevin/SW/HDF5_1.6/include 
# WITH OPENMPI
#export FFLAGS = -O3 -fdefault-real-8 -I../headers -I../headers -I/home/kevino/Documents/SW/openmpi-gfortran-gcc/include -fPIE
#export CFLAGS = -O2 -I../headers -I/usr/local/include -I/home/kevino/Documents/SW/openmpi-gfortran-gcc/include -fPIE
#export FFLAGS = -O3 -fdefault-real-8 -I../headers -I../headers -I/home/kevino/Documents/SW/SPACK/spack/opt/spack/linux-pop20-skylake/gcc-9.3.0/openmpi-3.1.6-tod27p37bnahlwqy6nd3jgrblj5vm7k3/include -I/home/kevino/Documents/SW/mpich-gfortran-gcc/include -fPIE 
#export CFLAGS = -O2 -I./headers -I../headers -I/usr/local/include -I/home/kevino/Documents/SW/SPACK/spack/opt/spack/linux-pop20-skylake/gcc-9.3.0/openmpi-3.1.6-tod27p37bnahlwqy6nd3jgrblj5vm7k3/include -I/home/kevino/Documents/SW/mpich-gfortran-gcc/include -fPIE

# NAG compiler
#export FFLAGS = -O3 -r8 -I../headers -mismatch -maxcontin=300 -fpp -I../headers -I/home/kevino/Documents/SW/mpich-gfortran-gcc/include
#export CFLAGS = -O2 -I../headers
#-I/home/kevino/Documents/SW/mpich-gfortran-gcc/include

# CRAY compiler
#export FFLAGS = -O3 -h noomp -s real64 -e p -I../headers -I../headers -I/lus/dal/users/kevin/HDF5/include
#export CFLAGS = -O2 -I../headers -I/lus/dal/users/kevin/HDF5/include

#-----------------------------------------------

#%.o:%.F90
#	$(FC) -c $(FFLAGS) $<


# Additional libraries to link to. You do not need
# To add the shmem library. This is automatically added
# if you define SHMEM=1 below.

# FOR XLF on mac
# FOR MPICH
#export ADD_LIB = /home/kevin/SW/HDF5_1.6/lib/libhdf5.a -lz -lm
# OPENMPI
#export ADD_LIB = /home/kevino/Documents/SW/openmpi-gfortran-gcc/lib/libhdf5.a -lz -lm -ldl
# CRAY
#export ADD_LIB = /lus/dal/users/kevin/HDF5/lib/libhdf5.a -lz -lm 
# laptop Nvidia Fortran
export ADD_LIB = /home/kevin/SW/HDF5_1.6/nvidia/lib/libhdf5.a -lz -lm 
# laptop gfortran
#export ADD_LIB = /home/kevin/SW/HDF5_1.6/lib/libhdf5.a -lz -lm 

#-----------------------------------------------

# some compilers can generate make rules to stdout from the source files
# if you have such a compiler, provide the flags, otherwise comment it out
#export MY_CPP := gcc -E -MM -MG  # for the GNU C Preprocessor

#-----------------------------------------------

# SHMEM or MPI ?
# uncomment to use SHMEM
#export SHMEM = 1

#--------------------------------------------------------------------------


.PHONY: all
ifdef SHMEM
all: libs headers source
else
all: libs headers mpi_source source
endif

.PHONY: headers
headers:
	$(MAKE) -C $@ -f Makefile.gnu
	cp -f headers/libmodules.a libs

.PHONY: mpi_source
mpi_source: headers
	$(MAKE) -C $@ -f Makefile.gnu
	cp -f mpi_source/libmpi_paramesh.a libs

.PHONY: source
source: headers
	$(MAKE) -C $@ -f Makefile.gnu
	cp -f source/libparamesh.a libs

.PHONY: clean
clean:
	$(RM) -r *~ libs
	for dir in headers {mpi_,}source Tests User_applic; do \
	  $(MAKE) -C $$dir -f Makefile.gnu clean; \
	done

.PHONY: Tests
Tests: all
	$(MAKE) -C $@ -f Makefile.gnu

libs:
	mkdir $@

# An example target to match an application directory name other than Tests
# in which the users application files are located.
User_applic:	all
	$(MAKE) -C $@ -f Makefile.gnu
