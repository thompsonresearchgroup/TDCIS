# TDCIS
A time-dependent CI solver that can use either an orthogonal or nonorthogonal determinant basis expansion.

## Features
+ The truncation level of the orthogonal CI expansion is arbitrary (singles only by default). A CAS expansion is also implemented. Nonorthogonal CI expansions currently require a single matrix file for each determinant in the expansion, which are typically obtained from a Gaussian calculation (a perl script to automatically generate nonorthogonal singles and doubles is included). 
+ Several different types of electric field pulses can be used:
   - Rectangular pulse
   - Delta pulse
   - Oscillating rectangular envelope
   - Oscillating Gaussian envelope
   - Chirped Gaussian envelope
+ Different pulse magnitudes, pulse lengths, pulse widths, chirp parameters, oscillation frequencies, and pulse onset/maximum time can be set using options.
+ Different inital states can be prepared.
+ Time correlation functions can be obtained from a specified initial time. 
+ Direct matrix elements can be computed without the requirement to store 2ERIs on disk (requires modified version of Gaussian).
+ Changes in the nuclear positions can be made using either the sudden or adiabatic approximation to propagate the wavefunction. A user script can be called to modify the nuclear step in each iteration. 

## Prerequisites
+ Required Software:
  - gfortran v8 or newer
  - MQCPack (https://github.com/MQCPack/mqcPack)
+ Optional Software:
  - Gaussian 16 or higher
  - gauopen (https://gaussian.com/interfacing/)
  
## Installing
To build:
1. Edit line 1 of makefile to point to the mqcPack installation directory.
2. Edit line 2 of makefile to determine whether blas/lapack or intel (mkl) libraries will be used.
3. Edit line 3 of makefile to select compiler (gfortran or pgfortran).
2. Type make.

## Running
See documentation folder or use the --help option for command line arguments used to specify options.

## Contact Information
For any questions, please contact Lee M. Thompson, University of Louisville, lee.thompson.1@louisville.edu.

## License
This project is licensed under the MIT license - see the LICENSE file for details.
