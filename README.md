# resonance-dm
Public release of the [Mathematica](https://www.wolfram.com/mathematica) module used to generate the spectra shown
in https://doi.org/10.1007/JHEP11(2015)171. 

## Usage
All the necessary functions and the parameters that are
independent of the considered scenario are defined in resonance.m, which can
be imported either in a Mathematica notebook or in a *.m script. Two such 
sample scripts are provided in equil.m and noneq.m, to be used to solve the
differential equations (3.19), (3.20) (or (3.21)) in the cases where the lepton
densities are or are not equilibrated.
The *.dat files contain tabulated grids of the active neutrino widths and of the
thermodynamical parameters.

## History and changelog
Code originally written in the spring of 2015 and released in June 2016 on [Mikko Laine's website](http://www.laine.itp.unibe.ch/dmpheno/).
Jun-July 2020: new release on github/zenodo, code cleanup and simpler handling of masses different from 7.1 keV and final temperatures greater than 1 MeV.
