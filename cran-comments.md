## Test environments
* Windows 10 (local) R 4.0.2
* win-builder (devel)
* macOS 10.13.6 High Sierra, R-release, brew
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Fiammetta Menchetti <fiammetta.menchetti@gmail.com>'
  
* This is a new release.

## Resubmission #1

This is a resubmission. The following changes were made:

* Package title was shortened to 49 characters ("MBSTS Models for Causal Inference and Forecasting")
* T/F were replaces with TRUE/FALSE
* Object names such as "T" were replaced with "Tt"
* Places in examples and vignette where user's par() values were modified were then reset back
