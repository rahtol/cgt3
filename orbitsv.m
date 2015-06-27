BeginPackage["orbitsv`"]

OrbitSV::usage = "Calculates Orbit, Schreier-vector and Generators of Stabilizer of alpha with resprct to <X>";
UBeta::usage = "Calculates Coset representative of alpha in its stabilizer in <X> from alphas schreier vector";

Begin["`Private`"]

UBeta[betain_Integer, v_List, X_List] :=
    Module[ {beta, u, k},
        beta = betain;
        If[ v[[beta]] == 0,
            {},
            u = Cycles[{}];
            k = v[[beta]];
            While[k != -1,
                u = PermutationProduct[X[[k]], u];
                beta = PermutationReplace[beta, InversePermutation[X[[k]]]];
                k = v[[beta]]
            ];
            u
        ]
    ]
  
  
OrbitSV[alpha_, X_List] :=
    Module[ {r, n, i, sv, orbit, insymbols, beta, betaxi, ub, ubx, ubxi, y, Y},
        r = Length[X];
        n = Max[alpha, Max[Map[PermutationMax, X]]];
        sv = Table[0, {n}];  (* schreier vector of alpha^<X> *)
        sv[[alpha]] = -1;
        Y = {}; (* list of schreier generators of the stabilizer of alpha in <X> *)
        orbit = {alpha};
        insymbols = {alpha};
        While [Length [insymbols] > 0,
        	beta = First [insymbols];
        	insymbols = Rest [insymbols];
            ub = UBeta[beta, sv, X];
            For[i = 1, i <= r, i++,
         	    betaxi = PermutationReplace[beta, X[[i]]];
         	    ubx = PermutationProduct[ub, X[[i]]];
                ubxi = UBeta[betaxi, sv, X];
                If[ !MemberQ[orbit, betaxi],
                    AppendTo[orbit, betaxi];
                    AppendTo[insymbols, betaxi];
                    sv[[betaxi]] = i,
                    y = PermutationProduct[ubx, InversePermutation[ubxi]];
                    If[ !MemberQ[Y, y] && y!=Cycles[{}],
                        AppendTo[Y, y]
                    ];
                    If[OptionValue[Verbose] > 0,
                        Print[
                        	"beta=", beta, 
                        	", betaxi=",betaxi, 
                        	", i=", i, 
                        	", y=", y,
                        	", ub=", ub,
                        	", ubx=", ubx,
                        	", ubxi^-1=", InversePermutation[ubxi]
                        ]
                    ];
                ];
            ];
        ];
        {orbit, sv, DeleteDuplicates[Y]}
    ]

Options[UBeta] = {Verbose->0};
Options[OrbitSV] = {Verbose->0};

End[]

EndPackage[]