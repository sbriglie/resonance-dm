(*driver for the case of equilibrated lepton densities. the numbers here 
refer to case c of the paper, with the central value for \sin^2(\theta)*)

(*written by Jacopo Ghiglieri, May-June 2015*)



<<resonance.m
(*set the parameters*)
md3 = Sqrt[1.75 10^(-11) mi^2]; md2 = 0; md1 = 0;
Tstart=4000; Tend=1;
startcondeq = 11.60;


(*define ff*)
ff[T_] := Evaluate@Array[Unique[][T] &, 200];


(*solve the coupled 3.19 and 3.21. need very high accuracy to catch
the resonance correctly*)
solsmallY = NDSolve[Join[{Yl'[T] == -1/T(totfitYintktbackfasteq[1,T,Yl[T],ff[T]]+       totfitYintktbackfasteq[2, T, Yl[T], ff[T]] + 
          totfitYintktbackfasteq[3, T, Yl[T], ff[T]]), 
      Yl[Tstart] == 
       vecYzeroeq[1, Tstart, startcondeq 10^(-6)] + 
        vecYzeroeq[2, Tstart, startcondeq 10^(-6)] + 
        vecYzeroeq[3, Tstart, startcondeq 10^(-6)]},
(*this are the 200 differential equations for each gridpoint of f_{k_T}*)  
     Table[ff'[T][[i]] == -1/T totfitYfeq[T, Yl[T], i, ff[T][[i]]], {i, 
       1, 200}], Table[ff[Tstart][[i]] == 0, {i, 1, 200}]], {Yl, 
     ff[T]}, {T, Tstart, Tend},AccuracyGoal->10];

(*get a table with k* and f_{k*} and its interpolator *)
final = Table[{gamma[[19800 + i, 2]], 
   solsmallY[[1, 2, 2, i]] /. {T -> 1}}, {i, 1, 200}];
finalint = Interpolation[final];
(*the abundance is simply given by 3.28 *)
abund=6950/(2 Pi^2) NIntegrate[k^2 finalint[k], {k, gamma[[19801, 2]], 
   gamma[[20000, 2]]}];

