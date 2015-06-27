BeginPackage["bsgs`", { "orbitsv`"}]

SchreierSims::usage = "construct BSGS using schreiersims algorithmus";
Strip::usage = "reduce permutatation by BSGS";
exS::usage = "helper function";


Begin["`Private`"] 

exS[x_Cycles, B_, i_Integer] := Length[Intersection[PermutationSupport[x], B[[1 ;; i - 1]]]] == 0;
exS[S_List, B_List, i_Integer] := Select[S, exS[#, B, i] &];
exS[S_List, B_List] := Table[exS[S, B, i], {i, 1, Length[B] + 1}]

Strip[g_Cycles, B_List, S_List, sv_List] :=
    Module[ {h, beta, i, k, ui},
          (* sv[[i[[ is the Schreier Vektor for beta[i] in <S^(i)>=H(i) *)
        k = Length[B]; (* also Length[sv] *)
        h = g;
        For[i = 1, i <= k, i++,
            (* h fixes base points beta[1],...,beta[i-1] *)
            beta = PermutationReplace[B[[i]], h];
            If[ sv[[i, beta]] == 0,
                Return[{h, i}]
            ];
            ui = UBeta[beta, sv[[i]], S[[i]]];
            h = PermutationProduct[h, InversePermutation[ui]];
        ];
        {h, k + 1}
    ]

SchreierSims[Bin_List, Sin_List] :=
    Module[ {n, B, S, sv, i, beta, xi, x, ubeta, betax, ubetax, h, j, y, l},
        n = Max[Map[PermutationMax, Sin]];
    	
        B = Bin;
        Scan[Function[x,
            If[ Length[Intersection[PermutationSupport[x], B]] == 0,
                AppendTo[B, First[PermutationSupport[x]]]
            ]
        ], Sin];
        S = exS[Sin, B];
        sv = Table[PadRight[OrbitSV[B[[i]], S[[i]]][[2]], n], {i, 1, Length[B]}];
        i = Length[B];
        
        While[i >= 1,
            y = True;
            For[beta = 1, y && (beta <= Length[sv[[i]]]), beta++, 
	          	If[ sv[[i, beta]] != 0,
	              	ubeta = UBeta[beta, sv[[i]], S[[i]]];
	              	For[xi = 1, y && (xi <= Length[S[[i]]]), xi++,
		                x = S[[i, xi]];
		                betax = PermutationReplace[beta, x];
		                ubetax = UBeta[betax, sv[[i]], S[[i]]];
		                If[ PermutationProduct[ubeta, x] != ubetax,
		                    {h, j} = Strip[PermutationProduct[ubeta, x, InversePermutation[ubetax]], B, S, sv];
		                    If[ j <= Length[B],
		                        y = False,
		                        If[ h != Cycles[{}],
		                            y = False;
		                            AppendTo[B, First[PermutationSupport[h]]];
		                            AppendTo[S, {}];
		                            AppendTo[sv, {}];
		                        ]
		                    ];
		                    If[ ! y,
		                        For[l = i + 1, l <= j, l++,
		                            AppendTo[S[[l]], h];
		                            sv[[l]] = PadRight[OrbitSV[B[[l]], S[[l]]][[2]], n];
		                        ];
		                        i = j;
		                    ];
		                ];
		            ];
		        ]
		    ];
	        If[ y,
	            i = i - 1
	        ];
    	];
        {B, S, sv}
    ]

End[]

EndPackage[]