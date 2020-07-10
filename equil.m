(*driver for the case of equilibrated lepton densities. the numbers here 
refer to case c of the paper, with the central value for \sin^2(\theta)*)

(*written by Jacopo Ghiglieri, May-June 2015

updated in July 2020 to include the possibility (largely untested
and disabled by default) of running with a newer momentum grid, stretching
farther into the infrared
*)

(*the new grid is disabled by default*)
grid2020=False


<<resonance.m
(*set the parameters*)

(*mass of the sterile neutrino, in MeV*)
mi =  7.1 10^(-3);

(*these are the Dirac masses, related to the mixing angle
by Eq. (4.1)*)
md3 = Sqrt[1.75 10^(-11) mi^2]; md2 = 0; md1 = 0;

(*initial and final temperatures for the evolution. They must obey
4000 MeV >= Tstart>Tend>= 1 MeV*)
Tstart=4000.; Tend=1.;

(*the initial conditions for the lepton flavor asymmetries Y_a=n_a/s*)

(*in cases a and c, the initial condition is set as *)
startconde = vecYzeroeq[1, Tstart, 11.60 10^(-6)];
startcondmu = vecYzeroeq[2, Tstart, 11.60 10^(-6)];
startcondtau = vecYzeroeq[3, Tstart, 11.60 10^(-6)];

(*
the vecYzero and vecYzeroasym functions translate a neutrino asymmetry Y_nu into
a lepton flavor asymmetry Y_a
*)

(*one can clearly also directly input their own values for the Y_a/s into
the startcond* variables*)

(*define the desired accuracy. Increase this if you see noise in the spectrum*)
myacc = 10;

(*define ff*)
ff[T_] := Evaluate@Array[Unique[][T] &, nbins];


(*solve the coupled 3.19 and 3.21. need very high accuracy to catch
the resonance correctly*)
solsmallY = NDSolve[Join[{Yl'[T] == -1/T(totfitYintktbackfasteq[1,T,Yl[T],ff[T]]+       totfitYintktbackfasteq[2, T, Yl[T], ff[T]] + 
          totfitYintktbackfasteq[3, T, Yl[T], ff[T]]), 
      Yl[Tstart] == 
       startconde  + 
        startcondmu + 
        startcondtau},
(*this are the nbins differential equations for each gridpoint of f_{k_T}*)  
     Table[ff'[T][[i]] == -1/T totfitYfeq[T, Yl[T], i, ff[T][[i]]], {i, 
       1, nbins}], Table[ff[Tstart][[i]] == 0, {i, 1, nbins}]], {Yl, 
     ff[T]}, {T, Tstart, Tend},AccuracyGoal->myacc];


(*get a table with k* (in MeV) and f_{k*} and its interpolator*)
final = Table[{Tend gamma[[(ntemp-1)*nbins + i, 2]] (heff[Tend]/heff[1.])^(1/3),
     solsmallY[[1, 2, 2, i]] /. {T -> Tend}}, {i, 1, nbins}];
finalint = Interpolation[final];
(*the abundance is simply given by 3.28*)

abund = 6950 heff[1.]/heff[Tend]/(2 Pi^2) mi/(7.1 10^-3) NIntegrate[k^2 finalint[k], 
{k, gamma[[(ntemp-1)*nbins+1, 2]] Tend (heff[Tend]/heff[1.])^(1/3), gamma[[ntemp*nbins, 2]] (heff[Tend]/heff[1.])^(1/3) Tend}]/Tend^3;

