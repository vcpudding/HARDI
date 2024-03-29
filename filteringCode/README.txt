README File for DT-MRI Filtering Code
-------------------------------------

This README is organized into 3 sections
1. Contents 
2. Compiling Instructions
3. Usage of programs.

Contents
---------

This folder contains the following:
+  Code 
+  Data
+  Readme

CODE:
This consists of

a. dwiFilter.cxx 
b. tensorFilter.cxx 
c. besselFunction.h
d. itkModules
e. namicModules
f. CMakeLists.txt

dwiFilter.cxx is the driver program for running the DWI space
filtering. 
tensorFilter.cxx is the driver program for running tensor Space
filtering.

besselFunction.h contains helper functions

itkModules (Directory)
Contains the implementation of various filters as
itkModules.

namicModules(Directory)
Contains code to compute Riemannian space tensor
statistics.

CMakeLists.txt is the CMake file used for building the
executables.


DATA:

Some sample data is provided for running the
filtering techniques.This can be found in the data 
directory. The noisy data has a rician
noise with std. deviation 10.

a. cleanDWI.nhdr
b. noisyDWI_10.nhdr
c. cleanTensor.nhdr
d. noisyTensor_10.nhdr
 

Compiling Instructions
----------------------

1. Untar the contents of this package to a directory
say filteringCode. 

2. Create a directory called build
 (either inside the filteringCode or at the same level
as filteringCode).

3. Before compiling ensure you have
+ cmake(version 2.0.5) installed.
+ A working itk build. (The code was tested with
  the ITK_CVS repository (www.itk.org) as of 
  15th May 2006. In case there are issues with the latest
  itk version, it would be better to get the most stable release 
  close to this period.


a) Go inside the build directory and run ccmake.
b) Set the correct path to itk Build directory.
c) Set build type to Release.
d) Configure and generate the Makefile.
e) Run cmake <path-to-CmakeLists.txt>
f) Run make in the build directory

Two executable programs should now be present in
your current directory:
1. dwiFilter
2. tensorFilter


Usage of Programs
-----------------
dwiFilter:
The executable <dwiFilter> performs filtering in the
Diffusion Weighted Image Space (DWI-Space) and when 
executed without arguments provides the usage.

dwiFilter (no Arguments)

USAGE:dwiFilter <arguments>
Arguments:
1. Input File Name
2. Output File Name
3. NumIterations
4. Conductance
5. TimeStep
6. Filter Type : (Simple Aniso-0,Chi Squared-1,Rician-2,Gaussian-3)
7. Sigma for bias correction
8. Lamda (Rician Correction Term)
9. Lamda (Gaussian Correction Term)

Argument Description:

<Input File Name> 
Name of the DWI file to be filtered. For example
<noisyDWI_10.nhdr> is a noisy DWI file provided
in the data directory. It was generated by adding 
synthetic Rician noise with a sigma=10 to a cleanDWI.nhdr

<Output File Name>
Name of the filtered DWI file. For example
<filteredDWI.nhdr>

<NumIterations>
Number of iterations you want to run the filter for.

<Conductance>
The value of the conductance term in anisotropic 
diffusion filtering (Ex: 1.0)
Note: Large Conductance will oversmooth the image
It is important to tune the conductance to obtain
best results.

<Time Step>
This determines the step size in the gradient
descent. It can be atmost 0.0625.

<Filter Type>
   Can Take 3 values:
   0 means perform simple anisotropic diffusion
   
   1 means perform Chi-Squared smoothing (square the image and
   perform anisotropic diffusion and then subtract the variance
   of the noise, and take square root. (The square of a Rice 
   distribution is a Chi Squared distribution with known bias 
   equal to the variance of the noise)
   (Refer:Max Likelihood Est. of Rician Ditribution Parameters. 
   Sijbers et. al)

   2 means Perform Rician bias correction filtering.
   (Refer: Rician Noise Removal in DT-MRI.)

   3 is same as 2 except use a Gaussian Attachment Term .

<Sigma>
Estimate of noise in the data.
This can be done by squaring the airvoxels
in the real data. The sum of square of all
the intensities in the air region should equal
2*variance of the noise in the data.
(Sijbers et. al)

<lamda1, lamda2>
The weights for the Rician and Gaussian 
attachment terms. 

Example Usage:
-------------
dwiFilter ../data/noisyDWI_10.nhdr filteredDWI.nhdr 1 1.0 0.0625 2 10 100 0

Filters the noisyDWI_10.nhdr for 1 iteration with a conductance of 1.0
timeStep 0.0625 using Rician filtering with a Rician attachement term
weight of 100. The estimate of noise in the input image is a sigma of 10
The filtered image is filteredDWI.nhdr.

-----------------------------------------------------------------------
tensorFilter:
The executable <tensorFilter> performs filtering in the
tensor Space  and when 
executed without arguments provides the usage.

tensorFilter (no Arguments)

Usage: tensorDiffuseTest <Arguments>
1. FilterType:(0-Euclidean, 1-Log Space,2-Riemannian)
2. numIterations:Iterations For Anisotropic Diffusion
3. timeStep:timeStep Used in Anisotropic Diffusion
4. conductance:Conductance used for Anisotropic Diffusion
5  Input (filename of input data)

Arguments 2,3,4 have the same meaning as described 
for dwiFilter.

Argument 1 describes the filter type
0: Euclidean Space filtering (tensors are treated as
   6-d vectors)
1: Log Space filtering (Fast and Simple Calculus on Tensors in the Log-Euclidean Framework. In J. Duncan and G. Gerig, editors, Proceedings of the 8th Int. Conf. on Medical Image Computing and Computer-Assisted Intervention - MICCAI 2005, Part I, volume 3749 of LNCS, Palm Springs, CA, USA, October 26-29, pages 115-122, 2005. Springer Verlag)

2. Riemannian Space Filtering(A Riemannian Framework for the Processing of Tensor-Valued Images. In Ole Fogh Olsen, Luc Florak, and Arjan Kuijper, editors, Deep Structure, Singularities, and Computer Vision (DSSCV), number 3753 of LNCS, pages 112-123, June 2005. Springer Verlag.)

Currently, the Riemannian filter adjustment for negative eigen-values
is hard-coded in the source file and might need to be tweaked to get
better filtering performance.

Argument 5 is the name of the noisyTensor input.

Example Usage:
--------------
tensorFilter 2 1 0.0625 1.0 ../data/noisyTensor_10.nhdr


Output:
Currently the filtered files are named as follows:

FilterType |  OutputFileName
---------------------------------------
 0         | EuclideanFiltered.nhdr
 1         |  LogSpaceFiltered.nhdr 
 2         |  SymmetricSpaceFiltered.nhdr
--------------------------------------------



