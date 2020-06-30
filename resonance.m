(*Mathematica module to obtain the spectra in 1506.06752.
this solves the equations for the non-equilibrated lepton active flavours,
i.e. equations (3.19) and (3.20) in the paper.

written by Jacopo Ghiglieri, May-June 2015
*)

(*the basic idea is that, as explained in Sec. 3.3., we solve along a trajectory
of redshifting momentum. Hence the function f_{k_T} is defined on a grid of 
200 momenta k whose k/T ratio is different at each T*)


(*first set up the parameters. all mass scales are in MeV*)
(*quarks*)
mu = 2.3; md = 4.8; ms = 95; mc = 1275; mb = 4180; 
(*leptons*)
me = 0.510998928 ; mmu = 105.6583715; mt = 1776.82;
ML = {me, mmu, mt};
(*electroweak parameters*)
 sinw2 = 0.23126; cosw = Sqrt[1 - sinw2]; 
GF = 1.1663787 10^(-5) 10^(-6); mw = 80385;
(*dirac masses for the neutrinos. they have to be
specified later on, when solving numerically the 
equations*)
MD = {md1, md2,md3}; 
(*planck mass*)
mpl = 1.2209*10^22;
(*define the chi function, as in Eq. (2.29)*) 
chi[m_, T_] := 
 1/(Pi^2 T) NIntegrate[p^2 Exp[Sqrt[p^2 + m^2]/T]/
(Exp[Sqrt[p^2 + m^2]/T] + 1)^2, {p, 0, Infinity},MaxRecursion->100,Method->
"LocalAdaptive"];
(*define \hat b, as in Eq. (3.10)*) 
hatb[a_, T_] := 16/( Pi GF Sqrt[2] mw^2/Pi) (7 Pi^2/360 cosw^2 + 
    2/(4 Pi^2) NIntegrate[
      x^2/(Exp[Sqrt[x^2 + ML[[a]]^2/T^2]] + 1)/
        Sqrt[x^2 + ML[[a]]^2/T^2] (4/3 x^2 + ML[[a]]^2/T^2), {x, 0, 
       Infinity},MaxRecursion->100,Method->"LocalAdaptive"]); 
(*import several auxiliary datasets.
the hatIQ files are the active neutrino  widths for the 3 flavours,
as in Eq. (3.9), as a function of the temperature and of a redshifting
momentum grid k_T as in (3.17)*)
gamma = Import["hatIQ_M0007_alpha1.dat", "Table",
	HeaderLines->6];
gamma2 =Import["hatIQ_M0007_alpha2.dat", "Table",
	HeaderLines->6];
gamma3 = Import["hatIQ_M0007_alpha3.dat", "Table",
	HeaderLines->6];
(*this is the list of temperature grid points (in MeV)*)
Tlist = Import["Tlist", "List"]; 
(*this is a table of N_c,eff, used to model hadronic effects as in Ref. [27].
It contains several other useful thermodynamical quantities.*)
ncefftab = Import["Nceff_expanded.dat", "Table",
	HeaderLines->2];

(*construct an array of stored values of chi for all the masses at the
grid temperatures*)
chistore = {}; Do[chietemp = chi[me, Tlist[[i]]]; 
 chimutemp = chi[mmu, Tlist[[i]]]; chittemp = chi[mt, Tlist[[i]]]; 
 chibtemp = chi[mb, Tlist[[i]]]; 
 chidstemp = chi[md, Tlist[[i]]] + chi[ms, Tlist[[i]]]; 
 chiuctemp = chi[mu, Tlist[[i]]] + chi[mc, Tlist[[i]]]; 
 AppendTo[chistore, {Tlist[[i]], chietemp, chimutemp, chittemp, 
   chietemp + chimutemp + chittemp, chiuctemp, chidstemp, chibtemp, 
   chidstemp + chibtemp, chiuctemp + chidstemp + chibtemp, 
   chiuctemp (chidstemp + chibtemp)}], {i, Length[Tlist], 1, -1}];
(*store the \hat b values as well*)
hatbtable = {}; Do[
 AppendTo[hatbtable, {Tlist[[i]], hatb[1, Tlist[[i]]], 
   hatb[2, Tlist[[i]]], hatb[3, Tlist[[i]]]}], {i, Length[Tlist], 1, -1}];


(*construct a set of interpolators
nceff= N_c,eff, g_eff and h_eff (see http://www.laine.itp.unibe.ch/eos15/), 
the speed of sound squared cs2*)

nceff = Interpolation[ncefftab[[All, {1, 2}]]];
geff =  Interpolation[ncefftab[[All, {1, 3}]]]; 
heff =  Interpolation[ncefftab[[All, {1, 4}]]]; 
cs2 =  Interpolation[ncefftab[[All, {1, 5}]]]; 
(*construct three tables with T/MeV, k/MeV, \hat I_Q*)
gammafit = Table[{{gamma[[i, 1]], gamma[[i, 1]] gamma[[i, 2]]}, 
   gamma[[i, 3]]}, {i, 1, 20000}];
gammafit2 =  Table[{{gamma2[[i, 1]], gamma2[[i, 1]] gamma2[[i, 2]]}, 
   gamma2[[i, 3]]}, {i, 1, 20000}];
gammafit3 = Table[{{gamma3[[i, 1]], gamma3[[i, 1]] gamma3[[i, 2]]}, 
   gamma3[[i, 3]]}, {i, 1, 20000}];
(*construct a vector of two-variable interpolators \hat I_Q(T,k)*)
g3d ={ Interpolation[ gammafit],
	Interpolation[  gammafit2],
	Interpolation[  gammafit3]}; 
(*an interpolation for hat b*)
bint = {Interpolation[hatbtable[[All, {1, 2}]]], 
  Interpolation[hatbtable[[All, {1, 3}]]], 
  Interpolation[hatbtable[[All, {1, 4}]]]}; 
(*and one for chi*)
chifit = {Interpolation[
   chistore[[All, {1, 2}]]], Interpolation[chistore[[All, {1, 3}]]], 
  Interpolation[chistore[[All, {1, 4}]]], 
  Interpolation[chistore[[All, {1, 5}]]], 
  Interpolation[chistore[[All, {1, 6}]]], 
  Interpolation[chistore[[All, {1, 7}]]], 
  Interpolation[chistore[[All, {1, 8}]]], 
  Interpolation[chistore[[All, {1, 9}]]], 
  Interpolation[chistore[[All, {1, 10}]]], 
  Interpolation[chistore[[All, {1, 11}]]]};


(*start defining functions. recall that we are calling c what is 
actually -c in the paper*)

(*this is a helper functions that returns the lepton flavor chemical potentials,
the charged and neutral lepton densities for each flavor ( see (2.33-2.35))
and the hadronic contribution to c (second line of (3.15)*)
vecmufitYall[T_Real, Ye_Real, Ymu_Real, Yt_Real] := 
  Module[{s, tempstruct, muq, nue, mue, numu, mumu, nut, mut, 
    mvec, chie, chimu, chit, chiemt, chiuc, chidsb, chiucdsb, 
    chizero,hadrc},
(*entropy and all the relevant \chi*) 
s = 2 Pi^2 T^3/45 heff[T];
   chie = chifit[[1]][T]; chimu = chifit[[2]][T]; 
   chit = chifit[[3]][T]; chiemt = chifit[[4]][T]; 
   chiuc = chifit[[5]][T]; chidsb = chifit[[8]][T]; 
   chiucdsb = chifit[[9]][T]; chizero = T^2/6; 
   (*this is what multiplies the sum in (2.33)*)
   tempstruct = chiucdsb/(chiemt chiucdsb + nceff[T] chiuc chidsb);
   (*extract the chemical potentials by inverting these relations*) 
   mue = s (2*chimu*(2*chit + chizero)*
         tempstruct*(chimu*Ye - chie*Ymu) - (2*chimu + 
           chizero)*((2*chit + chizero - 2*chit^2*tempstruct)*Ye + 
           2*chie*chit*tempstruct*Yt))/(2*
         chimu^2*(2*chie + chizero)*(2*chit + chizero)*
         tempstruct + (2*chimu + 
           chizero)*(-((2*chie + chizero)*(2*chit + chizero)) + 
           2*(2*chie*chit*(chie + chit) + (chie^2 + chit^2)*chizero)*
            tempstruct)); 
   mumu = s (-((2*chie + chizero)*(2*chit + chizero)*Ymu) + 
        2*(2*chie*chit*(chie + chit) + (chie^2 + chit^2)*chizero)*
         tempstruct*Ymu - 
        2*chimu*tempstruct*(chit*chizero*Yt + 
           chie*(chizero*Ye + 2*chit*(Ye + Yt))))/(-((2*chie + 
             chizero)*(2*chimu + chizero)*(2*chit + chizero)) + 
        2*(4*chie*chimu*chit*(chie + chimu + chit) + 
           2*(chie^2*(chimu + chit) + chimu*chit*(chimu + chit) + 
              chie*(chimu^2 + chit^2))*
            chizero + (chie^2 + chimu^2 + chit^2)*chizero^2)*
         tempstruct); 
   mut = s (-2*chit*
         tempstruct*(chimu*chizero*Ymu + 
           chie*(chizero*Ye + 
              2*chimu*(Ye + Ymu))) + (-((2*chie + chizero)*(2*chimu + 
                chizero)) + 
           2*(2*chie*chimu*(chie + chimu) + (chie^2 + chimu^2)*
               chizero)*tempstruct)*
         Yt)/(-((2*chie + chizero)*(2*chimu + chizero)*(2*chit + 
             chizero)) + 
        2*(4*chie*chimu*chit*(chie + chimu + chit) + 
           2*(chie^2*(chimu + chit) + chimu*chit*(chimu + chit) + 
              chie*(chimu^2 + chit^2))*
            chizero + (chie^2 + chimu^2 + chit^2)*chizero^2)*
         tempstruct);
(*the charge chemical potential*)	 
   muq = tempstruct (chie mue + chimu mumu + chit mut);
   (*hadronic contribution to c*)
   hadrc= (chie mue + chimu mumu + chit mut)nceff[T] chiuc chidsb
         /(chiemt chiucdsb + nceff[T] chiuc chidsb);
   mvec = {mue, mumu, mut}; 
   Return[Join[
    Table[{mvec[[a]], 2 chifit[[a]][T] (mvec[[a]] - muq), 
      chizero mvec[[a]]}, {a, 1, 3}],{hadrc}]]];

(*this is a helper functions that returns the lepton flavor chemical potential,
the charged and neutral lepton densities for each flavor ( see (2.36))
and the hadronic contribution to c (second line of (3.15) in the case
of equilibrated lepton densities*)
vecmufitYalleq[ T_Real, Y_Real] := 
 Module[{s, tempstruct, mulmuq, nua, mul, nnu, nel, hadrc}, 
  s = 2 Pi^2 T^3/45 heff[T]; nua = Y s;
(*mulmuq is 2 (\mu_L-\mu_Q)*)
mulmuq= 2 nceff[T]  chifit[[5]][
     T] chifit[[8]][
      T]/(chifit[[4]][T] (chifit[[5]][T] + chifit[[8]][T]) + 
       nceff[T] chifit[[5]][T] chifit[[8]][T]);
(*tempstruct is 2 \chi_{e\mu \tau}(\mu_L-\mu_Q)*)
  tempstruct = mulmuq chifit[[4]][T];
(*\mu_L*)
  mul = nua/(3 T^2/6 + tempstruct); 
  (*neutrino densities*)
  nnu = T^2/6 mul;
  (*charged lepton densities*) 
  nel =  mulmuq mul{chifit[[1]][T],chifit[[2]][T],chifit[[3]][T]};
(*hadronic contribution to c*)
hadrc = nceff[T] chifit[[4]][T] chifit[[5]][T] 
chifit[[8]][T] mul/(chifit[[4]][T] (chifit[[5]][T] + chifit[[8]][T]) + 
       nceff[T] chifit[[5]][T] chifit[[8]][T]); 
  Return[{mul, {nel[[1]],nel[[2]],nel[[3]]}, nnu, hadrc}]];

(*implementation of 3.20 . function of flavor, temperature, the three Y_a and
    the values of f_{k_T} at each point of the redshifting grid of k_T
    (see (3.17))*)
totfitYintktbackfast[a_, T_Real, Ye_Real, Ymu_Real, Yt_Real, fvec_] :=\
      Module[{mua, e1, array, ffit, gammaG, b, c, rpp, rpm, num, store}, 
  (*start by storing the result of vecmufitYall*)
  store = vecmufitYall[T, Ye, Ymu, Yt]; mua = store[[a, 1]]; 
  (*sterile neutrino energy*)
  e1 = Sqrt[mi^2 + kt^2];
(*construct an array where the first column is k_T and the second f_{k_T}*)
  array = 
   Table[{T gamma[[i, 2]] (heff[T]/heff[Tstart])^(1/3), 
     fvec[[i]]}, {i, 1, Length[fvec]}];
(*create an interpolator of the above*)  
 ffit = Interpolation[array];
(*introduce the width and the matter potentials b and c (see (3.9-3.14))*)
  gammaG = GF^2 T^4 e1 g3d[[a]][T, kt];
  b = GF^2 T^4 e1 bint[[a]][T];
  c = -Sqrt[2] GF ( +store[[a, 3]] + store[[a, 2]] + store[[1, 3]] + 
      store[[2, 3]] + 
      store[[3, 
        3]] + (-1/2 + 2 sinw2) (store[[1, 2]] + store[[2, 2]] + 
         store[[3, 2]])+2(1-2 sinw2)store[[4]]);
(*now build the rates in (3.8)*)
 num = MD[[a]]^2 mi^2 gammaG;
  rpp = num/((mi^2 + 2 e1 (b + c) + (b + c)^2)^2 + e1^2 gammaG^2); 
  rpm = num/((mi^2 + 2 e1 (b - c) + (b - c)^2)^2 + e1^2 gammaG^2);
(*return the r.h.s. of 3.20. The SymbolicProcessing makes the NIntegrate
	significantly faster. Having interpolators
for the relatively smooth functions above allows us to integrate without 
being constrained to a pre-determined grid, so that the resonance can 
be handled well, provided the accuracy used later in NDSolve is high enough.
  The part between zero and the first grid point is handled as a simple
  triangle.*)
  Return[1/(3 2 Pi^2 ( cs2[T] Sqrt[8 Pi/3 Pi^2/30 T^4 geff[T]]/mpl 2 Pi^2/
          45 T^3 heff[T]))( 
NIntegrate[kt^2 ((-(ffit[kt] - 1/(Exp[(e1 + mua)/T] + 1) ) rpm + (ffit[kt] 
- 1/(Exp[(e1 - mua)/T] + 1) ) rpp)), {kt,gamma[[1, 2]] T
     (heff[T]/heff[Tstart])^(1/3), gamma[[200, 2]] T
     (heff[T]/heff[Tstart])^(1/3)},   MaxRecursion -> 200, 
     Method -> {"LocalAdaptive", "SymbolicProcessing" -> 0}]+
(1/2 kt^3 ((-(ffit[kt] - 1/(Exp[(e1 + mua)/T] + 1) ) rpm + (ffit[kt] 
- 1/(Exp[(e1 - mua)/T] + 1) ) rpp))/.{kt->gamma[[1, 2]] T
     (heff[T]/heff[Tstart])^(1/3)}))]];


(*implementation of the r.h.s of Eq. (3.21). Lepton densities are
equilibrated here. Very similar to the previous function*)
totfitYintktbackfasteq[a_, T_Real, Y_Real, fvec_] :=
  Module[{mua, e1, array, ffit, gammaG, b, c, rpp, rpm, num, store}, 
  store = vecmufitYalleq[T, Y]; mua = store[[ 1]]; 
  e1 = Sqrt[mi^2 + kt^2];
  array = 
   Table[{T gamma[[i, 2]] (heff[T]/heff[Tstart])^(1/3), 
     fvec[[i]]}, {i, 1, Length[fvec]}];
  ffit = Interpolation[array]; 
gammaG = GF^2 T^4 e1 g3d[[a]][T, kt];
  b = GF^2 T^4 e1 bint[[a]][T];
  c = -Sqrt[2] GF ( 4 store[[ 3]] + store[[ 2,a]] 
         + (-1/2 + 2 sinw2) (store[[ 2,1]] + store[[2, 2]] + 
         store[[2, 3]])+2(1-2 sinw2) store[[4]]);
 num = MD[[a]]^2 mi^2 gammaG;
  rpp = num/((mi^2 + 2 e1 (b + c) + (b + c)^2)^2 + e1^2 gammaG^2); 
  rpm = num/((mi^2 + 2 e1 (b - c) + (b - c)^2)^2 + e1^2 gammaG^2);
  Return[1/(3 2 Pi^2 ( 
        cs2[T] Sqrt[8 Pi/3 Pi^2/30 T^4 geff[T]]/mpl 2 Pi^2/
          45 T^3 heff[T])) (NIntegrate[ 
     kt^2 ((-(ffit[kt] - 1/(Exp[(e1 + mua)/T] + 1) ) rpm + (ffit[
             kt] - 1/(Exp[(e1 - mua)/T] + 1) ) rpp)), {kt, 
      gamma[[1, 2]] T (heff[T]/heff[Tstart])^(1/3), gamma[[200, 2]] T (heff[T]/heff[Tstart])^(1/3)}, 
     MaxRecursion -> 200, 
     Method -> {"LocalAdaptive", "SymbolicProcessing" -> 0}]+
(1/2  kt^3 ((-(ffit[kt] - 1/(Exp[(e1 + mua)/T] + 1) ) rpm + (ffit[
             kt] - 1/(Exp[(e1 - mua)/T] + 1) ) rpp))/.{kt-> 
gamma[[1, 2]] T (heff[T]/heff[Tstart])^(1/3)} ) )]];



(*implementation of the r.h.s. of 3.19. i labels the ith site of
     the redshifting k_T grid*)
totfitYf[T_, Ye_Real, Ymu_Real, Yt_Real, i_, f_] :=
Module[{e1, tzero, gammaG, b, c, kt, num, rpm, rpp, mua, store, 
  ctilde}, tzero = Tstart;
(*start again by storing all the helper functions*)
store = vecmufitYall[T, Ye, Ymu, Yt]; 
mua = {store[[1, 1]],store[[2,1]],store[[3,1]]};
(*get k_T*)
 kt = (heff[T]/heff[tzero])^(1/3) gamma[[i, 2]] T;
 e1 = Sqrt[mi^2 + kt^2];
(*width vector*)
 gammaG = GF^2 T^4 e1 {g3d[[1]][T, kt],g3d[[2]][T, kt],g3d[[3]][T, kt]};
(*b, as above*)
 b = {GF^2 T^4 e1 bint[[1]][T], GF^2 T^4 e1 bint[[2]][T], 
   GF^2 T^4 e1 bint[[3]][T]};
(*construct a vector for the three c*) 
 ctilde = (+store[[1, 3]] + store[[2, 3]] + 
    store[[3, 
      3]] + (-1/2 + 2 sinw2) (store[[1, 2]] + store[[2, 2]] + 
       store[[3, 2]])+2(1-2 sinw2)store[[4]]); 
 c = -Sqrt[2] GF {ctilde + store[[1, 3]] + store[[1, 2]], 
    ctilde + store[[2, 3]] + store[[2, 2]], 
    ctilde + store[[3, 3]] + store[[3, 2]]};
(*now get the three R *)
 num = mi^2  ; rpp = {num gammaG[[1]] MD[[1]]^2/((mi^2 + 
          2 e1 (b[[1]] + c[[1]]) + (b[[1]] + c[[1]])^2)^2 + 
       e1^2 gammaG[[1]]^2), 
   num gammaG[[2]] MD[[2]]^2/((mi^2 + 
          2 e1 (b[[2]] + c[[2]]) + (b[[2]] + c[[2]])^2)^2 + 
       e1^2 gammaG[[2]]^2), 
   num gammaG[[3]] MD[[3]]^2/((mi^2 + 
          2 e1 (b[[3]] + c[[3]]) + (b[[3]] + c[[3]])^2)^2 + 
       e1^2 gammaG[[3]]^2)};
 rpm = {num gammaG[[1]] MD[[1]]^2/((mi^2 + 
          2 e1 (b[[1]] - c[[1]]) + (b[[1]] - c[[1]])^2)^2 + 
       e1^2 gammaG[[1]]^2), 
   num gammaG[[2]] MD[[2]]^2/((mi^2 + 
          2 e1 (b[[2]] - c[[2]]) + (b[[2]] - c[[2]])^2)^2 + 
       e1^2 gammaG[[2]]^2), 
   num gammaG[[3]] MD[[3]]^2/((mi^2 + 
          2 e1 (b[[3]] - c[[3]]) + (b[[3]] - c[[3]])^2)^2 + 
       e1^2 gammaG[[3]]^2)};
(*now return the r.h.s. of 3.19*)
 Return[1/(6 cs2[T] Sqrt[8 Pi/3 Pi^2/30 T^4 geff[T]]/
       mpl) Sum[((1/(Exp[(e1 - mua[[a]])/T] + 1) - 
         f) rpp[[a]] + (1/(Exp[(e1 + mua[[a]])/T] + 1) - 
         f) rpm[[a]]), {a, 1, 3}]]] ;

(*implementation of the r.h.s. of Eq. (3.19) for the case of 
equilibrated lepton densities. Very similar to the previous function*)
totfitYfeq[T_, Y_Real, i_, f_] :=
Module[{e1, tzero, gammaG, b, c, kt, num, rpm, rpp, mua, store, 
  ctilde}, tzero = Tstart;
 store = vecmufitYalleq[T, Y]; mua = store[[1]];
 kt = (heff[T]/heff[tzero])^(1/3) gamma[[i, 2]] T;
 e1 = Sqrt[mi^2 + kt^2];
 gammaG = {GF^2 T^4 e1 g3d[[1]][T, kt],GF^2 T^4 e1 g3d[[2]][T, kt],
GF^2 T^4 e1 g3d[[3]][T, kt]};
 b = {GF^2 T^4 e1 bint[[1]][T], GF^2 T^4 e1 bint[[2]][T], 
   GF^2 T^4 e1 bint[[3]][T]};
ctilde = (4 store[[ 3]] + (-1/2 + 2 sinw2) (store[[2, 1]] + store[[2, 2]] + 
       store[[2, 3]])+2(1-2 sinw2) store[[4]]); 
 c = -Sqrt[2] GF {ctilde  + store[[2, 1]], 
    ctilde  + store[[2, 2]], 
    ctilde + store[[2, 3]]};
 num = mi^2  ;
 rpp = {num gammaG[[1]] MD[[1]]^2/((mi^2 + 
          2 e1 (b[[1]] + c[[1]]) + (b[[1]] + c[[1]])^2)^2 + 
       e1^2 gammaG[[1]]^2), 
   num gammaG[[2]] MD[[2]]^2/((mi^2 + 
          2 e1 (b[[2]] + c[[2]]) + (b[[2]] + c[[2]])^2)^2 + 
       e1^2 gammaG[[2]]^2), 
   num gammaG[[3]] MD[[3]]^2/((mi^2 + 
          2 e1 (b[[3]] + c[[3]]) + (b[[3]] + c[[3]])^2)^2 + 
       e1^2 gammaG[[3]]^2)};
 rpm = {num gammaG[[1]] MD[[1]]^2/((mi^2 + 
          2 e1 (b[[1]] - c[[1]]) + (b[[1]] - c[[1]])^2)^2 + 
       e1^2 gammaG[[1]]^2), 
   num gammaG[[2]] MD[[2]]^2/((mi^2 + 
          2 e1 (b[[2]] - c[[2]]) + (b[[2]] - c[[2]])^2)^2 + 
       e1^2 gammaG[[2]]^2), 
   num gammaG[[3]] MD[[3]]^2/((mi^2 + 
          2 e1 (b[[3]] - c[[3]]) + (b[[3]] - c[[3]])^2)^2 + 
       e1^2 gammaG[[3]]^2)};
 Return[1/(6 cs2[T] Sqrt[8 Pi/3 Pi^2/30 T^4 geff[T]]/
       mpl) Sum[((1/(Exp[(e1 - mua)/T] + 1) - 
         f) rpp[[a]] + (1/(Exp[(e1 + mua)/T] + 1) - 
         f) rpm[[a]]), {a, 1, 3}]]] ; 


(*a function to get the initial condition, i.e. to get Y_a as a function of
       the initial active neutrino number density n_{\nu}/s (equal for the 
	three flavours) at a given temperature. to be used in cases (b) and
	(d) of the paper*)
vecYzero[a_, T_Real, Ynu_Real] := 
  Module[{tempstruct, mul}, 
   tempstruct = (chifit[[5]][T] + 
       chifit[[8]][
        T])/(chifit[[4]][T] (chifit[[5]][T] + chifit[[8]][T]) +  
       nceff[T] chifit[[5]][T] chifit[[8]][T]); mul = Ynu/(T^2/6); 
   Return[(T^2/6 + 
       2 chifit[[a]][T] (1 - tempstruct chifit[[4]][T])) mul]];

(*an analogous function for different n_{\nu_a}/s. To be used in cases
   (e) through (j)*)
vecYzeroasym[a_, T_Real, Ynue_Real, Ynumu_Real, Ynut_Real] := 
  Module[{tempstruct,tempstruct2, muv}, 
   tempstruct = (chifit[[5]][T] + chifit[[8]][T])/
(chifit[[4]][T] (chifit[[5]][T] + chifit[[8]][T]) +  
       nceff[T] chifit[[5]][T] chifit[[8]][T]); 
  muv = {Ynue,Ynumu,Ynut}/(T^2/6);
tempstruct2=chifit[[1]][T]muv[[1]]+chifit[[2]][T]muv[[2]]+
chifit[[3]][T]muv[[3]];
Return[(T^2/6 muv[[a]]+ 2 chifit[[a]][T] (muv[[a]] - tempstruct tempstruct2))]];

(*function to get the initial conditions in the case where
the lepton densities are equilibrated. See Eq. (2.36). To be used in
cases (a) and (c).*)
vecYzeroeq[a_, T_Real, Ynu_Real] := 
  Module[{tempstruct, mul}, 
   tempstruct = 
    2 nceff[T] chifit[[a]][T] chifit[[5]][
      T] chifit[[8]][
       T]/(chifit[[4]][T] (chifit[[5]][T] + chifit[[8]][T]) + 
        nceff[T] chifit[[5]][T] chifit[[8]][T]); mul = Ynu/(T^2/6); 
   Return[(T^2/6 + tempstruct) mul]];
