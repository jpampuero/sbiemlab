SBIEMLAB
a Spectral Boundary Integral Equation Method 
for 2D mode III rupture dynamics in Matlab

Version 0.1


CONTENTS
---------
	I.   Directory contents
	II.  Overview
	III. Pre-requisites
	IV.  Getting started
	V.   Questions and bug reports
	VI.  Copyright and license


I. DIRECTORY CONTENTS
----------------------

. README		the current file
. GNU_GPL		GNU General Public License
. SBIEM*jpg	example output figures
. src/SBIEM.m	the SBIEM solver
. src/friction.m	slip weakening friction laws
. src/SBIEM_ex*.m	examples (1 to 3)
. src/speed.m	estimates front speed

The function SBIEM creates and maintains a directory called 'kernels', 
in which some numerical coefficient tables are stored to avoid 
recomputations during subsequent calls.


II. OVERVIEW
-------------

SBIEMLAB is a collection of Matlab functions and scripts that solve 
antiplane (mode III) rupture dynamic problems with slip-weakening friction, 
on a 1D fault embedded in a 2D homogeneous elastic unbounded medium.
This is one of the most elementary models of dynamic earthquake rupture. 

The problem is initially formulated as a Boundary Integral Equation (BIE),
with space-time convolutions describing the elastodynamic stress transfer
along the fault. In the Spectral version of the BIE (SBIEM),
the space convolutions in the stress transfer functionals are expressed 
in the spectral domain associated with the along-fault spatial wavenumber (k).
The elastodynamic kernels are obtained analytically.
The method is introduced in detail by 

  J. W. Morrisey and P. H. Geubelle (1997),
  "A numerical scheme for mode III dynamic fracture problems"
  Int. J. Num. Meth. Eng., 40 (7), 1181-1196.

and has been improved and extensively used by Alain Cochard, Nadia Lapusta, etc.  
SBIEMLAB implements the "velocity" formulation of the stress transfer functional.
Time truncation of time-convolutions is implemented in SBIEMLAB as in 

  N. Lapusta et al. (2000), "Elastodynamic analysis for slow 
  tectonic loading with spontaneous rupture episodes on faults 
  with rate- and state-dependent friction"
  J. Geoph. Res. 105 (B10), 23765-23789.

The friction law assumed in SBIEMLAB is slip weakening (linear or not).
The method is however largely independent of the nature of the friction law,
and can be easily adapted to rate-and-state friction (e.g. Lapusta et al, 2000).

SBIEMLAB is intended to introduce new users to the SBIEM or 
to introduce researchers and students to computational earthquake dynamics.
For serious simulations you should turn to an optimized Fortran 
implementation of the SBIEM or of the Spectral Element Method (SEM):
  + MDSBI, https://pangea.stanford.edu/~edunham/codes/codes.html
  + SEM2DPACK, https://github.com/jpampuero/sem2dpack


III. PRE-REQUISITES
--------------------

I assume you have some familiarity with computational earthquake dynamics.
If this is not the case yet, the following paper is a good starting point:

  S. M. Day, L. A. Dalguer, N. Lapusta and Y. Liu (2005)
  Comparison of finite difference and boundary integral solutions 
  to three-dimensional spontaneous rupture 
  J. Geophys. Res. 110 (B12), B12307.

In particular you should pay attention to the grid spacing required to 
achieve adequate numerical resolution of the rupture process zone, otherwise
you will obtain numerical noise, notably in the slip rate histories.
SBIEMLAB gives limited warning about numerical resolution. You are responsible
for setting the grid spacing according to the usual criterion: several grid points
within the process zone size, defined as the zone between the rupture front
and the end of weakening where slip = Dc (e.g. see reference above).

The SBIEM implemented here assumes periodicity along the fault, beyond the
limits of the modelled fault segment. Stresses propagating at the S wave speed 
can wrap around the limits and create unphysical or undesired perturbations 
if the fault segment is not set large enough. SBIEMLAB gives limited warning
about this situation, and you are responsible for setting the fault size accordingly.

SBIEMLAB assumes you are familiar with Matlab's "structures" (variables containing
named fields).


IV. GETTING STARTED
--------------------

1. Under Matlab, type "help SBIEM" to get a complete description of the solver's
input and output arguments.

2. Type "help friction" for a description of the friction law parameters.

3. Explore the commented example script SBIEM_ex1.m,
   it illustrates the typical usage of SBIEM and visualization of results

4. If needed, explore the remaining examples scripts SBIEM_ex*.m 


V. QUESTIONS AND BUG REPORTS
-----------------------------

For any problem running SBIEMLAB contact the author:

Jean-Paul Ampuero
ampuero@geoazur.unice.fr

https://jpampuero.github.io/


VI. COPYRIGHT AND LICENSE
--------------------------

Copyright (C) 2006-2007 Jean-Paul Ampuero

This software is freely available for scientific research purposes. 
If you use this software in writing scientific papers or reports
include proper attributions to its author, Jean-Paul Ampuero.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

