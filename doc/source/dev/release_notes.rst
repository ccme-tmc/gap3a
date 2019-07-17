Release Notes
=============

Version 2c
----------

Oct 2013 - Jun 2016

New features: 

* interfacing with Wannier90 
* constrained RPA for Hubbard U
* support with arbitary number of local orbitals (LO) for each l-channel and for l upto to lomax
* add a python script ``gap2_w2kset``, which can be used to reset IPRINT to 1 in case.inc 
  so that it is not necessary to change the source of WIEN2k 

Fixes:

* "rotdef.f90" -- the bug related to the lower/higher case of F,B,C... is corrected (16 Oct 2014)

