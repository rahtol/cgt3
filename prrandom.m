BeginPackage["prrandom`"]

PrInitialize::usage = "initialize random element generator for group <X>";
PrRandom::usage = "generate random element";

Begin["`Private`"]

PrInitialize[X_List,r_Integer,n_Integer] :=
	Module[{l = Length[X], x, x0, prstate},
		x = Table[X[[Mod[i-1,l]+1]], {i, 1, r}];
		x0 = Cycles[{}];
		prstate = {x0, x}; 
		Do [PrRandom[prstate], {n}];
		prstate		
	]

PrRandom[prstate_Symbol] := 
	Module[{x0=prstate[[1]], x=prstate[[2]], r=Length[prstate[[2]]], s, t, k, e},
		s = RandomInteger[{1, r}];
		t = RandomInteger[{1, r-1}]; If[t==s, t = r];
		k = RandomInteger[{1, 2}];
		e = 2 RandomInteger[{0, 1}] -1;
		If[k == 1,
			x[[s]] = PermutationProduct[x[[s]], PermutationPower[x[[t]], e]];
			x0 = PermutationProduct[x0, x[[s]]], 
			x[[s]] = PermutationProduct[PermutationPower[x[[t]], e], x[[s]]];
			x0 = PermutationProduct[x[[s]], x0]
		];
		prstate = {x0, x}; 
		x0	
	]

SetAttributes[PrRandom, HoldFirst];

End[] (* End Private Context *)

EndPackage[]