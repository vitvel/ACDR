(* This math-file contains code, which is used in the paper Velizhanin V.N. "Analytic continuation of harmonic sums: dispersion representation" *)
(* https://arxiv.org/abs/2205.15182 *)
(* It demands two files with databases: *)
(* 1. DiffRootUpToW6.m with the linear difference equations for nested harmonic sums up to weight 6 *)
(* 2. ACtablePolesUpToW7.m with the pole expression for the analytically continued nested harmonic sums up to weight 7 *)
(* Database for the pole expresson for the analytically continued nested harmonic sums up to weight 13 can be found here: *)
(* https://github.com/vitvel/ACHS *)
(* Total time and RAM for the evaluation of this file is about 600 sec. and slightly more then 2 Gb *)
(* File can be loaded from MATHEMATICA FrontEnd with <<"DispRepr.m" with full path, if necessary *)
(* The extended version of necessary databases up to weight 11 can be found on GitHub *)
(* https://github.com/vitvel/ACDR *)

(* Generate all nested harmonic sums up to weight "weight" *)

weight = 6;
tuphs = Table[Tuples[{1, 0, -1}, weightI] // (# //. {a1___, 0, a_ /; a != 0, a2___} :> {a1, Sign[a] (Abs[a] + 1), a2} ) & // DeleteCases[#, {___,0}] & // Sort[#] &,{weightI,1,weight}]//Flatten[#,1]&;
  
(* Load DifferenceRoot for weight "weight" *)
  
diffroots = Get["DiffRootUpToW"<>ToString[weight]<>".m"];

(* Generate numerical values for nested harmonic sums of weight "weight" for 2000 integers with WorkingPrecision -> "prec" *)

prec = 100;
NR = 2000;
Print["Generate HS up to weight 6 for 2000 integers. Time for evaluation is about 60 sec."];
Print["Time ",Timing[rectab = RecurrenceTable[{diffroots[[#]][[1, 2]][[1]], Sequence @@ Thread[diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 1]] == (diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 2]] //  N[#, prec + 100] &)]}, \[FormalY], {\[FormalN], 1, NR},WorkingPrecision -> prec] & /@ Range[Length[tuphs]];][[1]]];

(* Set values for nested harmonic sums hsHS[a1,a2,...,an,i] *)

Thread[Set[Evaluate[Flatten[Table[hsHS[Sequence @@ tuphs[[#]], i], {i, 1, NR}] & /@ Range[Length[tuphs]]]], Flatten[rectab]]];

Print["RAM ",MaxMemoryUsed[]];

Print["Check ",hsHS[-1,-1,-1,-1,-1,-1,NR-1]];

(* Load the pole expression for the analytically continued nested harmonic sums of weight "weight+1" *)

ACP = Get["ACtablePolesUpToW7.m"];
DACP = Dispatch[ACP];

(* Extract pole expression for HS[2, 1, 1, 1, 1, 1], as in Eq.(31) *)

PES=S[2, 1, 1, 1, 1, 1] /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &;

(* Produce dispersion representation *)

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M))-(PES/.w->(i)))/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

(* Evaluate dispersion representation at M=2 *)

Print["Time/Result for S[2,1,1,1,1,1][2] ",Timing[res=DispRepr[2,NR]//N[#,10]&]];

(* Nested harmonic sums for integers *)

HS[a_ /; (a != 0), jList__, n_?IntegerQ] := HS[a, jList, n] = Sum[(Sign[a]^i/i^Abs[a])*HS[jList, i], {i, 1, n}];
HS[a_ /; (a != 0), n_?IntegerQ] := HS[a, n] = Sum[Sign[a]^i/i^Abs[a], {i, 1, n}];

(* Compare with exact result *)

Print["Compare with exact result ",res-HS[2,1,1,1,1,1,2]];

(* The same for HS[-2, 1, 1, 1, 1, 1] *)

PES=S[-2, 1, 1, 1, 1, 1] /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &;

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M))-(PES/.w->(i)))/.signN->(-1)^i/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

Print["Time/Result for S[-2,1,1,1,1,1][2]",Timing[res=DispRepr[2,NR]//N[#,10]&]];

Print["Compare with exact result ",res-HS[-2,1,1,1,1,1,2]];

(* For analytic continuation from odd values we must change signN to (-signN) only in first brackets *)
(* For odd M we must sum first term or second term up to NN-1, therefore we mutipply first term by UnitStep[NN - i - 1] *)

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M)/.signN->-signN)*UnitStep[NN - i - 1]-(PES/.w->(i)))/.signN->(-1)^i/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

Print["Time/Result for S[-2,1,1,1,1,1][3]",Timing[res=DispRepr[3,NR]//N[#,10]&]];

Print["Compare with exact result ",res-HS[-2,1,1,1,1,1,3]];

(* For harmonic sums started with "1" we must extract S[1] with the help of SRemoveLeading1 from HarmonicSums-package by Jakob Ablinger https://risc.jku.at/sw/harmonicsums/ *)
(* For example we cosider S[1, 1, -2, 1, -2]. S1=PolyGamma[M+1]-PolyGamma[1] is S[1] *)

HSTest = (S1^2*S[-2, 1, -2])/2 + S1*(S[-3, 1, -2] + S[-2, 1, -3] + S[-2, 2, -2] - S[-2, 1, -2, 1] - 2*S[-2, 1, 1, -2]) - (3*S[-2, 1, 2, -2])/2 - S[-2, 2, -2, 1] -  (3*S[-2, 2, 1, -2])/2 + S[2, -2, 1, -2]/2 + S[-2, 1, -2, 1, 1] + 2*S[-2, 1, 1, -2, 1] + 3*S[-2, 1, 1, 1, -2]+S[-4, 1, -2]/2 + S[-3, 1, -3] + S[-3, 2, -2] + S[-2, 1, -4]/2 +S[-2, 2, -3] + S[-2, 3, -2]/2 - S[-3, 1, -2, 1] -  2*S[-3, 1, 1, -2] - S[-2, 1, -3, 1] - S[-2, 1, -2, 2]/2 - 2*S[-2, 1, 1, -3];

PES=HSTest /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &;

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M))-(PES/.w->(i)))/.signN->(-1)^i/.S1->PolyGamma[M+1]-PolyGamma[1]/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

Print["Time/Result for S[1,1,-2,1,-2][2]",Timing[res=DispRepr[2,NR]//N[#,10]&]];

Print["Compare with exact result ",res-HS[1,1,-2,1,-2,2]];

(* For analytic continuation from odd values we must change signN to (-signN) only in first brackets *)

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M)/.signN->-signN)*UnitStep[NN - i - 1]-(PES/.w->(i)))/.signN->(-1)^i/.S1->PolyGamma[M+1]-PolyGamma[1]/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

Print["Time/Result for S[1,1,-2,1,-2][3]",Timing[res=DispRepr[3,NR]//N[#,10]&]];

Print["Compare with exact result ",res-HS[1,1,-2,1,-2,3]];

(*----------------------------------------------------------------------------------------------------------------*)
(* In practice, we don't need to generate values for all nested harmonic sums, only for those, which enter into expression *)
(* This will save time and RAM *)
(* Recompute the same as above with 100000 terms *)

(* Generate all nested harmonic sums up to weight "weight" *)

weight = 6;
tuphs = Table[Tuples[{1, 0, -1}, weightI] // (# //. {a1___, 0, a_ /; a != 0, a2___} :> {a1, Sign[a] (Abs[a] + 1), a2} ) & // DeleteCases[#, {___,0}] & // Sort[#] &,{weightI,1,weight}]//Flatten[#,1]&;
  
(* Load DifferenceRoot up to weight "weight" *)
  
diffroots = Get["DiffRootUpToW"<>ToString[weight]<>".m"];

(* Load the pole expression for the analtically continued nested harmonic sums of weight "weight+1" *)

ACP = Get["ACtablePolesUpToW7.m"];
DACP = Dispatch[ACP];

(* Extract pole expression for HS[2, 1, 1, 1, 1, 1], as in Eq.(31) *)

PES=S[2, 1, 1, 1, 1, 1] /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &;

(* Obtain the list of variables *)

PESvar=PES//Variables//Cases[#,S[a__]:>{a}]&;

(* Find the position inside the list "tuphs" for nested harmonic sums, entering into expression for PES *)

tuphsC=(Position[tuphs,#]&/@PESvar)//Flatten;

(* Generate numerical values for necessary nested harmonic sums for 100000 integers with WorkingPrecision -> "prec" *)

prec = 100;
NR = 100000;
Print["Generate HS for S[2,1,1,1,1,1] for 100000 integers. Time for evaluation is about 30 sec."];
Print["Time ",Timing[rectab = RecurrenceTable[{diffroots[[#]][[1, 2]][[1]], Sequence @@ Thread[diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 1]] == (diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 2]] //  N[#, prec + 100] &)]}, \[FormalY], {\[FormalN], 1, NR},WorkingPrecision -> prec] & /@ tuphsC;][[1]]];

(* Clear hsHS *)

Clear[hsHS];

(* Set values for nested harmonic sums hsHS[a1,a2,...,an,i] *)

Thread[Set[Evaluate[Flatten[Table[hsHS[Sequence @@ tuphs[[#]], i], {i, 1, NR}] & /@ tuphsC]], Flatten[rectab]]];

Print["RAM ",MaxMemoryUsed[]];

(* Produce dispersion representation *)

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M))-(PES/.w->(i)))/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

(* Evaluate dispersion representation at M=2 *)

Print["Calculate S[2,1,1,1,1,1][2]. Time for evaluation is about 20 sec. "];
Print["Time/Result for S[2,1,1,1,1,1][2]",Timing[res=DispRepr[2,NR]//N[#,10]&]];

(* Compare with exact result *)

Print["Compare with exact result ",res-HS[2,1,1,1,1,1,2]];

(* The same for HS[-2, 1, 1, 1, 1, 1] *)

PES=S[-2, 1, 1, 1, 1, 1] /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &;

PESvar=PES//Variables//Cases[#,S[a__]:>{a}]&;

tuphsC=(Position[tuphs,#]&/@PESvar)//Flatten;

prec = 100;
NR = 100000;
Print["Generate HS for S[-2,1,1,1,1,1] for 100000 integers. Time for evaluation is about 30 sec."];
Print["Time ",Timing[rectab = RecurrenceTable[{diffroots[[#]][[1, 2]][[1]], Sequence @@ Thread[diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 1]] == (diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 2]] //  N[#, prec + 100] &)]}, \[FormalY], {\[FormalN], 1, NR},WorkingPrecision -> prec] & /@ tuphsC;][[1]]];

Clear[hsHS];

Thread[Set[Evaluate[Flatten[Table[hsHS[Sequence @@ tuphs[[#]], i], {i, 1, NR}] & /@ tuphsC]], Flatten[rectab]]];

Print["RAM ",MaxMemoryUsed[]];

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M)/.signN->-signN)*UnitStep[NN - i - 1]-(PES/.w->(i)))/.signN->(-1)^i/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

(* Evaluate dispersion representation at M=3 *)

Print["Calculate S[-2,1,1,1,1,1][2]. Time for evaluation is about 20 sec. "];
Print["Time/Result for S[-2,1,1,1,1,1][3]",Timing[res=DispRepr[3,NR]//N[#,10]&]];

(* Compare with exact result *)

Print["Compare with exact result ",res-HS[-2,1,1,1,1,1,3]];

(* The same for S[1, 1, -2, 1, -2] *)

HSTest = (S1^2*S[-2, 1, -2])/2 + S1*(S[-3, 1, -2] + S[-2, 1, -3] + S[-2, 2, -2] - S[-2, 1, -2, 1] - 2*S[-2, 1, 1, -2]) - (3*S[-2, 1, 2, -2])/2 - S[-2, 2, -2, 1] -  (3*S[-2, 2, 1, -2])/2 + S[2, -2, 1, -2]/2 + S[-2, 1, -2, 1, 1] + 2*S[-2, 1, 1, -2, 1] + 3*S[-2, 1, 1, 1, -2]+S[-4, 1, -2]/2 + S[-3, 1, -3] + S[-3, 2, -2] + S[-2, 1, -4]/2 +S[-2, 2, -3] + S[-2, 3, -2]/2 - S[-3, 1, -2, 1] -  2*S[-3, 1, 1, -2] - S[-2, 1, -3, 1] - S[-2, 1, -2, 2]/2 - 2*S[-2, 1, 1, -3];

PES=HSTest /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &;

PESvar=PES//Variables//Cases[#,S[a__]:>{a}]&;

tuphsC=(Position[tuphs,#]&/@PESvar)//Flatten;

prec = 100;
NR = 100000;
Print["Generate HS for S[1,1,-2,1,-2] for 100000 integers. Time for evaluation is about 90 sec."];
Print["Time ",Timing[rectab = RecurrenceTable[{diffroots[[#]][[1, 2]][[1]], Sequence @@ Thread[diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 1]] == (diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 2]] //  N[#, prec + 100] &)]}, \[FormalY], {\[FormalN], 1, NR},WorkingPrecision -> prec] & /@ tuphsC;][[1]]];

Clear[hsHS];

Thread[Set[Evaluate[Flatten[Table[hsHS[Sequence @@ tuphs[[#]], i], {i, 1, NR}] & /@ tuphsC]], Flatten[rectab]]];

Print["RAM ",MaxMemoryUsed[]];

DispRepr[M_,NN_] := Sum[((PES/.w->(i+M)/.signN->-signN)*UnitStep[NN - i - 1]-(PES/.w->(i)))/.signN->(-1)^i/.S1->PolyGamma[M+1]-PolyGamma[1]/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

(* Evaluate dispersion representation at M=3 *)

Print["Calculate S[1,1,-2,1,-2][3]. Time for evaluation is about 200 sec. "];
Print["Time/Result for S[1,1,-2,1,-2][3]",Timing[res=DispRepr[3,NR]//N[#,10]&]];

(* Compare with exact result *)

Print["Compare with exact result ",res-HS[1,1,-2,1,-2,3]];

(*----------------------------------------------------------------------------------------------------------------*)
(* Second order of BFKL Pomeron eigenvalues arXiv:1507.04010 for the conformal spin n=0 and BFKL intercept function for arbitrary conformal spin arXiv:1802.06908 *)
(* Eq.(4.15) from arXiv:1802.06908 and Eq.(4) from arXiv:1507.04010 multiplied by 2 due to Eq.(3) of arXiv:1507.04010 *)

BFKLC=16*S[-3] + 8*Pi^2*S[1]/3 + 16*S[3] - 32*S[-2, 1]/.S[1]->S1;
BFKLP=2*(4*ln2*Pi^2 - 6*z3 + 4*Pi^2*S[-1] + 4*Pi^2*S[1]/3 + 8*S[3] -  16*S[-2, 1])/.S[1]->S1//Expand;

(* Extract pole expressions for HS *)

PESBFKLC=BFKLC /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &
PESBFKLP=BFKLP /. S[i__] :> S[R[i], N] /. DACP /. signx -> 1 // Series[#, {w, 0, -1}] & // Normal // Expand // Collect[# /. S[R[a__], i] :> S[a], w, Factor] &

PESvar={PESBFKLC,PESBFKLP}//Variables//Cases[#,S[a__]:>{a}]&;

(* Find position nested harmonic sums entering into expressions for PES *)

tuphsC=(Position[tuphs,#]&/@PESvar)//Flatten;

(* Generate numerical values for nested harmonic sums of weight "weight" for 100000 integers with WorkingPrecision -> "prec" *)

prec = 100;
NR = 100000;
Print["Generate HS for NLO BFKL for 100000 integers ",Timing[rectab = RecurrenceTable[{diffroots[[#]][[1, 2]][[1]], Sequence @@ Thread[diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 1]] == (diffroots[[#]][[1, 2]][[2 ;; -1]][[All, 2]] //  N[#, prec + 100] &)]}, \[FormalY], {\[FormalN], 1, NR},WorkingPrecision -> prec] & /@ tuphsC;]];

(* Clear hsHS *)

Clear[hsHS];

(* Set values for nested harmonic sums hsHS[a1,a2,...,an,i] *)

Thread[Set[Evaluate[Flatten[Table[hsHS[Sequence @@ tuphs[[#]], i], {i, 1, NR}] & /@ tuphsC]], Flatten[rectab]]];

Print["RAM ",MaxMemoryUsed[]];

(* Produce dispersion representation *)

DispReprBFKLC[M_,NN_] := Sum[((PESBFKLC/.w->(i+M))-(PESBFKLC/.w->(i)))/.signN->(-1)^i/.S1->PolyGamma[M+1]-PolyGamma[1]/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;
DispReprBFKLP[M_,NN_] := Sum[((PESBFKLP/.w->(i+M))-(PESBFKLP/.w->(i)))/.signN->(-1)^i/.S1->PolyGamma[M+1]-PolyGamma[1]/.S[a__]:>hsHS[a,i-1]/.hsHS[__,0]:>0/.Table[ToExpression["z"<>ToString[k]]->Zeta[k],{k,1,13}],{i,1,NN}]//N[#,100]&;

(* Evaluate dispersion representations at M=-1/2 *)

Print["Calculate Difference at M=-1/2. Time for evaluation is about 30 sec. "];
Print["Difference at M=-1/2 ",Timing[DispReprBFKLC[-1/2,NR]-(DispReprBFKLP[-1/2,NR]+(BFKLP/.S1->0/.S[__]:>0/.ln2->Log[2]/.z3->Zeta[3]))//N[#,20]&]];


