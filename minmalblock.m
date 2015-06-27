BeginPackage["minmalblock`"]

Minimalblock::usage = "Find minimal block system with alpha being part of a block in group <X>";
(*
Rep::usage = "helper function";  
Merge::usage = "helper function";
*)  

Begin["`Private`"] 

Minimalblock[X_List, alpha_List] :=
	Module[{n, k, i, gamma, delta, p, c, q, nu, gammax, deltax},
		n = Max[Map[PermutationMax, X]];
		k = Length[alpha];
		p = Table [i, {i, 1, n}];
		c = Table [1, {n}];
		For[i=1, i<k, i++, p[[alpha[[i+1]]]]=alpha[[1]]];
		q = Table[alpha[[i+1]], {i,1,k-1}];
		c[[alpha[[1]]]] = k;
		While[Length[q] > 0,
			gamma = First[q];
			q = Rest[q];
			Map[Function[x,
					delta = Rep[gamma, p];
					gammax = PermutationReplace[gamma, x];
					deltax = PermutationReplace[delta, x];
					nu = Merge[gammax, deltax, c, p];
					If[nu>0, AppendTo[q, nu]]
				], X
			]
		];
		For[i=1, i<=n, i++, Rep[i, p]];
		p
	]
	
Rep[kappa_, p_] := 
	Module[{lambda, rho, mu},
		lambda = kappa;
		rho = p[[lambda]];
		While[lambda != rho, lambda = rho; rho = p[[lambda]];];
		(* now perform path compression *)
		mu = kappa; rho = p[[mu]];
		While[lambda != rho, p[[mu]] = lambda; mu = rho; rho = p[[mu]];];
		lambda
	]

Merge[kappa_, lambda_, c_, p_] := 
	Module[{fi, psi, mu, nu},
		fi = Rep[kappa, p];
		psi = Rep [lambda, p];
		If[fi == psi, 
			0,
			If[c[[fi]] >= c[[psi]], mu = fi; nu = psi, nu = fi; mu = psi];
			p[[nu]] = mu;
			c[[mu]] = c[[mu]] + c[[nu]];
			nu
		]
	]

SetAttributes[Rep, HoldRest];
SetAttributes[Merge, HoldAll];

Options[Minimalblock] = {Verbose->0};
Options[Rep] = {Verbose->0};
Options[Merge] = {Verbose->0};

End[]

EndPackage[]