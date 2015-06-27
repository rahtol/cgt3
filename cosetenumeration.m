BeginPackage["cosetenumeration`"]

CosetTable::usage = "Coset table data structure. Internal use.";
CreateInitialCosetTable::usage = "Create an initial coset table with only the 1-coset representing subgroup H of G. Allow growing up to M cosets."
DEFINE::usage = "Adjoins a new coset beta to coset table ct by setting beta=alpha^x. Returns beta.";
SCAN::usage = "Scan coset alpha under relator w from A*. Returns scan code from 1 to 5.";
CosetEnumerationR::usage = "Perform coset enumeration using HLT method.";
COMPRESS::usage = "Compress the coset table, i.e. make Omega equal to 1..n.";

addcoset::usage = "helper function";  
getp::usage = "helper function";  
setp::usage = "helper function";  
getchi::usage = "helper function";  
setchi::usage = "helper function";
inv::usage = "helper function";
coincidence::usage = "helper function";
rep::usage = "helper function";
merge::usage = "helper function";
SCANANDFILL::usage = "helper function";

Begin["`Private`"]

(*
	ct1 = CosetTable[{1}, {{0, 0, 0, 0}}, 100]  initial coset table with two generators, e.g. x,y
	[[1]]  p : [1..n] -> [1..n] where p[a]=a for all live cosets a and p[x]<x for merged cosets
	[[2]]  chi : [1..n] x A -> [0..n] defines action of coset alpha on x el A, 0=undefined
	[[3]]  M : integer defining the maximum number of cosets, i.e. n<=M
	storing n explicitly is not necessary, can be derived by e.g. n=Length[p]
	tau is omitted in implementation as described in the book
*)

(*
	It's a good idea to make the call-by-ref paramter the first one (HoldFirst)
	If possible use only one call-by-ref parameter
*)

addcoset[ct_Symbol] := Module[{n, k},
   Assert[Head[ct] === CosetTable];
   {n , k} = Dimensions[ct[[2]]];
   If[n < ct[[3]],
    AppendTo[ct[[1]], 0];
    AppendTo[ct[[2]], Table[0, {k}]];
    n+1 (* return the new coset *),
    Abort[]
    ]
   ];

getp[ct_CosetTable, alpha_] := ct[[1, alpha]];
setp[ct_Symbol, alpha_, beta_] := (ct[[1, alpha]] = beta;);
getchi[ct_CosetTable, alpha_, x_] := ct[[2, alpha, x]];
setchi[ct_Symbol, alpha_, x_, beta_] := (ct[[2, alpha, x]] = beta;);
inv[i_Integer] := If[EvenQ[i], i-1, i+1]; (* x[i] -> x[i]^-1 such that inv[inv[i]]==i holds for all i *) 

CreateInitialCosetTable [X_List, M_Integer:200] := CosetTable[{1}, {Table[0,{2*Length[X]}]}, M];

DEFINE[ct_Symbol, alpha_, x_] :=
	Module[{beta},
		If[(beta=getchi[ct,alpha,x])==0,
			beta = addcoset[ct];
			setp[ct, beta, beta];
			setchi[ct, alpha, x, beta];
			setchi[ct, beta, inv[x], alpha];
			beta,
			-beta
		]
	]

SCAN[ct_Symbol, alpha_, w_List] :=
	Module[{r=Length[w], i=1, j=Length[w], f=alpha, b=alpha, fxi, bxi,rc=0},
		(*Scan forwards *)
		While[i<=r && (fxi=getchi[ct, f, w[[i]]])!=0, f=fxi; i++];
		If[i>r,
			If[f!=alpha, coincidence[ct,f,alpha]; rc=1, rc=2],
			(* Forward scan was incomplete. Scan backwards *)
			While[j>=i && (bxi=getchi[ct, b, inv[w[[j]]]])!=0, b=bxi; j--];
			If[j<i, 
				coincidence[ct,f,b]; rc=3,
				If[j==i, 
					(* Deduction *)
					setchi[ct, f, w[[i]], b];
					setchi[ct, b, inv[w[[i]]], f];
					rc=4,
					(* j>i, scan is incomplete and yields no information *)
					rc=5
				]
			]
		];
		rc
	]

SCANANDFILL[ct_Symbol, alpha_, w_List] :=
	Module[{r=Length[w], i=1, j=Length[w], f=alpha, b=alpha, fxi, bxi,cont=1},
		While[cont!=0,
			cont=0; (* do it only once, except cont is modified, i.e. when scan is incomplete. *)
			(*Scan forwards. *)
			While[i<=r && (fxi=getchi[ct, f, w[[i]]])!=0, f=fxi; i++];
			If[i>r,
				If[f!=alpha, coincidence[ct,f,alpha]],
				(* Forward scan was incomplete. Scan backwards. *)
				While[j>=i && (bxi=getchi[ct, b, inv[w[[j]]]])!=0, b=bxi; j--];
				If[j<i,
					Assert[f!=b]; 
					coincidence[ct,f,b],
					If[j==i, 
						(* Deduction *)
						setchi[ct, f, w[[i]], b];
						setchi[ct, b, inv[w[[i]]], f],
						(* j>i, scan is incomplete => Make a new definition and continue *)
						DEFINE[ct, f, w[[i]]];
						cont=1; (* continue with this scan and current values of i and j *)
					]
				]
			]
		]
	]

coincidence[ct_Symbol, alpha_, beta_] :=
	Module[{n, k, gamma, q={}, x, delta, mu, nu, mux, nuinvx, gamma2},
	    {n , k} = Dimensions[ct[[2]]];
		gamma2 = merge[ct[[1]], alpha, beta];
		If[gamma2 != 0, AppendTo[q, gamma2]];
		While[Length[q] > 0,
			gamma = First[q];
			q = Rest[q];
			Assert[ct[[1,gamma]]<gamma]; (* not a live coset any more *)
			For[x=1, x<=k, x++,
				delta = getchi[ct, gamma, x];
				If[delta>0,
					setchi[ct, delta, inv[x], 0]; (* undefine delta^(x^-1) *)
					mu = rep[ct[[1]], gamma];
					nu = rep[ct[[1]], delta];
					If[(mux=getchi[ct,mu,x])!=0,
						gamma2 = merge[ct[[1]], nu, mux];
						If[gamma2 != 0, AppendTo[q, gamma2]],
						If[(nuinvx=getchi[ct,nu,inv[x]])!=0,
							gamma2 = merge[ct[[1]], mu, nuinvx];
							If[gamma2 != 0, AppendTo[q, gamma2]],
							setchi[ct, mu, x, nu];
							setchi[ct, nu, inv[x], mu]
						]
					]
				]
			]
		]				
	]

rep[p_, kappa_] :=
	Module[{lambda, rho, mu, tmp},
		lambda = kappa;
		rho = p[[lambda]];
		While[lambda != rho, lambda = rho; rho = p[[lambda]];];
		(* now perform path compression *)
		mu = kappa; rho = p[[mu]];
		While[lambda != rho, tmp = p[[mu]]; p[[mu]] = lambda; mu = rho; rho = tmp;];
		lambda
	]

merge[p_, kappa_, lambda_] :=
	Module[{fi, psi, mu, nu},
		fi = rep[p, kappa];
		psi = rep [p, lambda];
		If[fi != psi, 
			mu = Min[fi,psi];
			nu = Max[fi,psi];
			Assert[p[[mu]]==mu];
			Assert[p[[nu]]==nu];
			p[[nu]] = mu;
			nu,
			0
		]
	]

CosetEnumerationR[X_List, R_List, Y_List, M_Integer:200] :=
	Module[{n, k, ct, alpha, w, x},
		ct = CreateInitialCosetTable[X, M];
	    {n , k} = Dimensions[ct[[2]]];
		Map[SCANANDFILL[ct,1,#]&, Y];
		For[alpha=1, alpha<=(n=Length[ct[[1]]]), alpha++,
			For[w=1, (w<=Length[R]) && (ct[[1,alpha]]==alpha), w++,
				SCANANDFILL[ct, alpha, R[[w]]];
			];
			If[ct[[1,alpha]]==alpha, 
				For[x=1, x<=k, x++,
					If[getchi[ct, alpha, x] == 0, DEFINE[ct, alpha, x]]
				]
			]
		];
		ct
	]

COMPRESS[ctin_CosetTable (* call by value ! *)] :=
	Module[{ct=ctin, n, k, gamma=0, beta, alpha, x, i},  (* copy to local variable is really necessary *)
	    {n , k} = Dimensions[ct[[2]]];
		For[alpha=1, alpha<=n, alpha++, If [ct[[1,alpha]]==alpha,
			gamma++;
			If[gamma!=alpha,
				(* Replace alpha by gamma in coset table *)
				For[x=1, x<=k, x++,
					beta=getchi[ct, alpha, x];
					If[beta==alpha, beta=gamma];
					setchi[ct, gamma, x, beta];
					setchi[ct, beta, inv[x], gamma];
				]
			]
		]];
		n=gamma;
		CosetTable[Table[i, {i,1,n}], ct[[2, 1;;n, All]], n]
	]

SetAttributes[addcoset, HoldFirst];
SetAttributes[setp, HoldFirst];
SetAttributes[setchi, HoldFirst];
SetAttributes[coincidence, HoldFirst];
SetAttributes[rep, HoldFirst];
SetAttributes[merge, HoldFirst];
SetAttributes[DEFINE, HoldFirst];
SetAttributes[SCAN, HoldFirst];
SetAttributes[SCANANDFILL, HoldFirst];

End[]

EndPackage[]