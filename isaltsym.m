BeginPackage["isaltsym`", { "orbitsv`", "prrandom`"}]

IsAltSym::usage = "Test for Symmetric or Alternating Group.";
checkPermForLongPrimeCycle::usage = "helper function";

Begin["`Private`"]

verbose = 0;

IsAltSym[X_List, OptionsPattern[]] := 
	Module[{n, orbit, l, d, nmax, prstate,rc=False, eps},
		verbose = OptionValue[IsAltSym, Verbose];
		eps = OptionValue[IsAltSym, Eps];
		n = Max[Map[PermutationMax, X]];
		orbit = OrbitSV[1, X];
		l = Length[orbit[[1]]];
		If[n!=l, Return[False]];
		If[n<8, Return[]];
		d = If[n>16, 0.57, 0.34] * (Log[2]/Log[n]);
		nmax = Floor[N[-Log[eps]/d]] + 1;
		If[BitAnd[verbose, 16^^1]!=0, Print["nmax=", nmax]];
		prstate = PrInitialize [X, 11, 43];
		Do[If[checkPermForLongPrimeCycle[PrRandom[prstate], n], rc=True; Break[]], {nmax}];
		Return[rc]; 
	]
	
checkPermForLongPrimeCycle[p_Cycles, n_Integer] :=
	Module[{lencycles, t0, rc},
		lencycles = Map[Length, p[[1]]];
		t0 = Map[((#>n/2) && (#<n-2) && PrimeQ[#])&, lencycles];
		rc = Apply[Or, t0];
		If[BitAnd[verbose, 16^^2]!=0, Print[{p, rc}]];
		rc
	]

Options[IsAltSym] = {Verbose->0, Eps->10^-9};
	
End[]

EndPackage[]