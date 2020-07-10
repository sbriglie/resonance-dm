# resonance-dm
Public release of the [Mathematica](https://www.wolfram.com/mathematica) module used to generate the spectra shown
in https://doi.org/10.1007/JHEP11(2015)171.
[![DOI](https://zenodo.org/badge/275853039.svg)](https://zenodo.org/badge/latestdoi/275853039)


## Usage
All the necessary functions and the parameters that are
independent of the considered scenario are defined in resonance.m, which can
be imported either in a Mathematica notebook or in a *.m script. Two such 
sample scripts are provided in equil.m and noneq.m, to be used to solve the
differential equations (3.19), (3.20) (or (3.21)) in the cases where the lepton
densities are or are not equilibrated.
The *.dat files contain tabulated grids of the active neutrino widths and of the
thermodynamical parameters.

### Grid choice (introduced July 2020)

A new momentum grid has been added in v 1.1. The new grid has 250
momentum bins, while the original one has 200. The difference is that the new grid has a finer spacing and lower momenta in the IR, stretching down to *k*/*T*&approx;0.002 at *T*=1 MeV. For *k/T*> 0.3 (at *T*=1 MeV) the two grids are identical. 
At the moment this new grid has been tested only on a limited number of cases and it is switched off by default. See the variable grid2020 at the beginning of equil.m and noneq.m if you
want to turn it on.

### Limitations

This code reflects the approximation &mdash; originally made in https://doi.org/10.1007/JHEP11(2015)171 &mdash; of considering equal the distributions of the two helicity states of the sterile neutrino. If all that matters is the sum of the two distributions and the total DM abundance, these are obtained from (twice) our distribution *f, as long as it remains sufficiently smaller than the equilibrium distribution n<sub>F</sub>*. https://doi.org/10.1007/JHEP07(2019)078 presents a set of equations and rates that properly account for the distribution of the two helicity state. From that one can indeed see that the condition *n<sub>F</sub>*>>*f* guarantees that the helicity asymmetry has a negligible impact on the helicity-symmetric distribution.
We further note that, for many parameter points, the condition 
*n<sub>F</sub>*>>*f* holds generally; the points for which *f*(*k*<sub>T</sub>) approaches *n<sub>F</sub>*(*k*<sub>T</sub>) the most are the infrared ones, which may in some cases thermalise. Hence, this limitation will at most affect the IR part of the distribution *f*, whose contribution to the DM abundance is however small. Thus, the impact on the total DM abundance should remain within the overall uncertainties of the approach of https://doi.org/10.1007/JHEP11(2015)171. 

## History and changelog
### History
Code originally written in the spring of 2015 and released in June 2016 on [Mikko Laine's website](http://www.laine.itp.unibe.ch/dmpheno/).
### Changelog
- v 1.1 July 2020: added a new momentum grid, with finer spacing and lower momenta in the IR. At the moment it has been tested only on a limited number of cases and it is switched off by default. See the variable grid2020 at the beginning of equil.m and noneq.m
- v 1.0 Jun-July 2020: new release on github/zenodo, code cleanup and simpler handling of masses different from 7.1 keV and final temperatures greater than 1 MeV.

