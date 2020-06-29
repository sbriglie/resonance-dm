(*driver for the case of non-equilibrated lepton densities. the numbers here 
refer to case d of the paper, with the central value for \sin^2(\theta)*)

(*written by Jacopo Ghiglieri, May-June 2015*)


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

(*in cases b and d, the initial condition is set as *)
startconde = vecYzero[1, Tstart, 17.80 10^(-6)];
startcondmu = vecYzero[2, Tstart, 17.80 10^(-6)];
startcondtau = vecYzero[3, Tstart, 17.80 10^(-6)];
(*in cases e-j the function vecYzeroasym should
	be used, i.e.
startconde = vecYzeroasym[1, Tstart, Ynue,Ynumu,Ynutau];
startcondmu = vecYzeroasym[2, Tstart, Ynue,Ynumu,Ynutau];
startcondtau = vecYzeroasym[3, Tstart, Ynue,Ynumu,Ynutau];

the vecYzero and vecYzeroasym functions translate a neutrino asymmetry Y_nu into
a lepton flavor asymmetry Y_a
*)

(*one can clearly also directly input their own values for the Y_a/s into
the startcond* variables*)

(*define the desired accuracy. Increase this if you see noise in the spectrum*)
myacc = 10;

(*define ff*)
ff[T_] := Evaluate@Array[Unique[][T] &, 200];

(*solve the coupled 3.19 and 3.20. need very high accuracy to catch
the resonance correctly*)
solsmallY =NDSolve[Join[{Ye'[T] ==
	 -1/T totfitYintktbackfast[1, T, Ye[T], Ymu[T], Yt[T], ff[T]], 
Ymu'[T] == -1/T totfitYintktbackfast[2, T, Ye[T], Ymu[T], Yt[T], ff[T]], 
Yt'[T] == -1/T totfitYintktbackfast[3, T, Ye[T], Ymu[T], Yt[T],ff[T]], 
(*initial conditions. in cases e-j the function vecYzeroasym should
	be used*)
Ye[Tstart] == startconde, 
     Ymu[Tstart] == startcondmu, 
     Yt[Tstart] == startcondtau},
(*this are the 200 differential equations for each gridpoint of f_{k_T}*) 
    Table[ff'[T][[i]] == -1/T totfitYf[T, Ye[T], Ymu[T], Yt[T], i, 
        ff[T][[i]]], {i, 1, 200}], 
    Table[ff[Tstart][[i]] == 0, {i, 1, 200}]], {Ye, Ymu, Yt, ff[T]},
 {T, Tstart, Tend},AccuracyGoal->myacc]; 

(*get a table with k* (in MeV) and f_{k*} and its interpolator*)
final = Table[{Tend gamma[[19800 + i, 2]] (heff[Tend]/heff[1.])^(1/3),
     solsmallY[[1, 4, 2, i]] /. {T -> Tend}}, {i, 1, 200}];
finalint = Interpolation[final];
(*the abundance is simply given by 3.28*)

abund = 6950 heff[1.]/heff[Tend]/(2 Pi^2) mi/(7.1 10^-3) NIntegrate[k^2 finalint[k], 
{k, gamma[[19801, 2]] Tend (heff[Tend]/heff[1.])^(1/3), gamma[[20000, 2]] (heff[Tend]/heff[1.])^(1/3) Tend}]/Tend^3;
