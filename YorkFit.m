(* ::Package:: *)

(* 	York Fit by Karolis Misiunas karolis@misiunas.com. Adaptation from:

  http://mathematica.stackexchange.com/questions/13054/estimate-error-on-slope-of-linear-regression-given-data-with-associated-uncertai
  http://aapt.scitation.org/doi/pdf/10.1119/1.1632486


  Requirements: M10
  Tests: M11.2

	Versions:
	1.0 - initial release. 
  2.0 (2017-12-05) - updated version to York 2004 with options mimicking LinearModelFit

 *)

BeginPackage["YorkFit`"]

YorkFit::usage = "\!\(\*
StyleBox[\"YorkFit\",\n\
FontVariations->{\"Underline\"->True}]\)[{{x,y},...}, \"Errors\"-> \
err] York linear regression fit:
  http://aapt.scitation.org/doi/pdf/10.1119/1.1632486
  
  \"Errors\" can be supplied with:
    Automatic - assumes 5% error variation in X and Y
    None - no errors in x and y
    p_NumberQ - assumes the erros are p fraction of the data
    {pX_NumberQ, pY_NumberQ} - assumes the erros are {pX,pY} fraction of the data
    {errY1,errY2,...} - errors only in Y
    {{errX1,errY1},...} - errors in x and y
  
  Weights -> None
    Allows to manually specify weight matrix

  \"ErrorCorrelations\" -> None
    Correlations between X and Y errors

  IncludeConstantBasis -> True
    if False, this leads to fit y=b*x
    achieved by adding {0,0} point with very high weights 

  Output: <|\"fit\" -> best fit, \"a\" -> offset, \"b\" -> slope, \"\[Sigma]a\" -> error of offset, 
    \"\[Sigma]b\" -> error of slope |>";


LinearModelFitYork::usage = 
  "\!\(\*
StyleBox[\"LinearModelFitYork\",\n\
FontVariations->{\"Underline\"->True}]\)[{{x,y},...}, \"Errors\"-> \
err] will perform Krystek/Anton and York linear regression fit:
  http://iopscience.iop.org/0957-0233/18/11/025/
  http://www.nrcresearchpress.com/doi/abs/10.1139/p66-090#.VXYHm9SrT2Q
  http://aapt.scitation.org/doi/pdf/10.1119/1.1632486
  
  \"Errors\" can be supplied with:
    Automatic - assumes 5% error variation in X and Y
    None - no errors in x and y
    p_NumberQ - assumes the erros are p fraction of the data
    {pX_NumberQ, pY_NumberQ} - assumes the erros are {pX,pY} fraction of the data
    {errY1,errY2,...} - errors only in Y
    {{errX1,errY1},...} - errors in x and y
  
  Output: {fitted fn: y=A*x+B, {error in A, error in B} }";


LinearModelFitYork::errBadInput = 
  "Did not recognise \"Errors\" input consult docs";
LinearModelFitYork::errFormatData = 
  "Input data must be a Nx2 matrix";
LinearModelFitYork::errFormatErrors = 
  "\"Errors\"  must be a convertable to Nx2 matrix";
LinearModelFitYork::errErrorsNoZeros = 
  "\"Errors\" can't contain zero elements";

Options[LinearModelFitYork] = {
   "Errors" -> Automatic
   };

Options[YorkFit] = {
   "Errors" -> Automatic,
   Weights -> None,
   IncludeConstantBasis -> True,
   "ErrorCorrelations" -> None
   };

(* Implementations *)
Begin["`Private`"]

(*original function from stackexchange*)
ortlinfit[data_?MatrixQ, errs_?MatrixQ] := Module[
  {n = Length[data], c, ct, dk, dm, k, m, p, s, st, ul, vl, w, wt, xm, ym},
  (*yes, I know I could have used FindFit[] for this...*)
  {ct, st, k} = Flatten[MapAt[Normalize[{1, #}] &, 
      NArgMin[Norm[
        Function[{x, y}, y - \[FormalM] x - \[FormalK]] @@@ 
         data], {\[FormalM], \[FormalK]}], 1]];
   (*find orthogonal regression coefficients*){c, s, p} = 
    FindArgMin[{Total[(data.{-\[FormalS], \[FormalC]} - \
\[FormalP])^2/((errs^2).{\[FormalS]^2, \[FormalC]^2})], \[FormalC]^2 \
+ \[FormalS]^2 == 1}, {{\[FormalC], ct}, {\[FormalS], 
       st}, {\[FormalP], k/ct}}];
   (*slope and intercept*){m, k} = {s, p}/c;
   wt = 1/errs^2; w = (Times @@@ wt)/(wt.{1, m^2});
   {xm, ym} = w.data/Total[w];
   ul = data[[All, 1]] - xm; vl = data[[All, 2]] - ym;
   (*uncertainties in slope and intercept*)
   dm = w.(m ul - vl)^2/w.ul^2/(n - 2);
   dk = dm (w.data[[All, 1]]^2/Total[w]);
   {Function[\[FormalX], Evaluate[{m, k}.{\[FormalX], 1}]], 
    Sqrt[{dm, dk}]}] /; Dimensions[data] === Dimensions[errs]

(*testing cases*)
(*
data = {{0, 5.9}, {0.9, 5.4}, {1.8, 4.4}, {2.6, 4.6}, {3.3, 
    3.5}, {4.4, 3.7}, {5.2, 2.8}, {6.1, 2.8}, {6.5, 2.4}, {7.4, 1.5}};

errs = {{1000., 1.}, {1000., 1.8}, {500., 4.}, {800., 8.}, {200., 
     20.}, {80., 20.}, {60., 70.}, {20., 70.}, {1.8, 100.}, {1, 
     500.}} // Sqrt[1/#] &;

{lin, {sm, sk}} = ortlinfit[data, errs] 
*)

(* ====================================== Cosmetic method ================================== *)


LinearModelFitYork[data_?MatrixQ, opts : OptionsPattern[] ] := Module[
  {checkErrors, compute, er},
  
  compute[errs_] := If[   Dimensions[data][[2]] == 2,
    ortlinfit[ data, checkErrors[errs]],
    (*else*)  Message[LinearModelFitYork::errFormatData]
    ];
  
  checkErrors[errs_] := If[ Dimensions[errs] == Dimensions[data],
    If[ ! MemberQ[errs, 0.0, 2],
     errs,
      Throw@Message[LinearModelFitYork::errErrorsNoZeros]],
    Throw@ Message[LinearModelFitYork::errFormatErrors] ];
  
  Switch[ er = OptionValue["Errors"],
   Automatic, compute[Abs[data*0.05 + $MinMachineNumber]],
   None,  compute[Abs[data*0.0 + $MinMachineNumber]] ,
   _?NumberQ,  compute[Abs[data*er + $MinMachineNumber]],
   {_?NumberQ, _?NumericQ},   
    compute[
    Abs[Transpose@{data[[All, 1]]*er[[1]], 
        data[[All, 2]]*er[[2]]} + $MinMachineNumber]],
   _?VectorQ, 
              
   compute[Abs[
     Transpose@{data[[All, 1]]*0.0 + $MinMachineNumber, N@er}]],
   _?MatrixQ, compute[N@er],
   _, Message[LinearModelFitYork::errBadInput]
   ]
  ]





(* ====================================== New Implementation =============================== *)


(*reads better from Notebook interface *)
yorkFitRaw[Xi_, Yi_, \[Omega]Xi_, \[Omega]Yi_, ri_] := Module[
  {iterateParametersFor, fit, b0, 
   \[Alpha]i, Wi, Wsum, Xmean, Ymean, Ui, Vi, \[Beta]i, a, b, xi, 
   xmean, ui, \[Sigma]bSq, \[Sigma]aSq, 
   n = 1},
  
  (* Step 1 *)
  
  fit = LinearModelFit[Transpose[{Xi, Yi}], x, x, 
    Weights -> \[Omega]Yi];
  b0 = fit["BestFitParameters"][[2]];
  
  (* Step 2 *)
  (* done at input time*)
  
  iterateParametersFor[b_] := Module[{bNew},
    (* Step 3 *)
    \[Alpha]i = Sqrt[\[Omega]Xi*\[Omega]Yi];
    Wi = (\[Omega]Xi*\[Omega]Yi)/(\[Omega]Xi + b^2  \[Omega]Yi - 
        2 b  ri  \[Alpha]i);
    
    (* Step 4 *)
    Wsum = Total[Wi];
    Xmean = Total[Wi* Xi]/Wsum ;
    Ymean = Total[Wi* Yi]/Wsum ;
    
    Ui = Xi - Xmean;
    Vi = Yi - Ymean;
    
    \[Beta]i = 
     Wi * (Ui / \[Omega]Yi  + 
        b Vi / \[Omega]Xi - (b Ui - Vi) ri / \[Alpha]i );
    (*\[Beta]mean = Total[Wi* \[Beta]i]/Wsum;*)
    
    (* Step 5 *)
    
    bNew = Total[Wi * \[Beta]i * Vi] / Total[Wi*\[Beta]i*Ui];
    bNew
    ];
  
  (* Step 6 *)
  
  While[ n < 10000,
   b = iterateParametersFor[b0];
   If[Abs[b - b0] < 10^-12, Break[]];
   (*bList = Append[bList, Abs[b-b0]];*)
   b0 = b;
   n++
   ];
  If[n >= 10000 - 1 , 
   Print["Warning: York fit failed to converge in ", n, " iterations"]];
  
  (* Step 7 *)
  a = Ymean - b * Xmean;
  
  (* Step 8 *)
  xi = Xmean - \[Beta]i;
  (*yi = Ymean - b* \[Beta]i;*)
  
  (* Step 9 *) 
  xmean = Total[Wi* xi]/Wsum ;
  (*ymean = Total[Wi* yi]/Wsum ;*)
  ui = xi - xmean;
  (*vi = yi - ymean;*)
  (*S= Total[Wi (Yi -b Xi - a)^2];*)
  
  (* Step 10 *)
  \[Sigma]bSq = 1/Total[Wi ui^2];
  \[Sigma]aSq = 1/Wsum  + xmean^2 \[Sigma]bSq;
  
  (* return all*)
  
  Return[<|"a" -> a, "b" -> b, "\[Sigma]a" -> Sqrt[\[Sigma]aSq], 
    "\[Sigma]b" -> Sqrt[\[Sigma]bSq]|> ];
]




YorkFit[data_?MatrixQ, opts : OptionsPattern[] ] := Module[
  {length, r, errors, weights, result, er},
  
  (* Check data *)
  length = Length[data];
  If[Dimensions[data][[2]] != 2, Throw @ Message[LinearModelFitYork::errFormatData]];
  dataIn = data;

  (* get weights *)
  errors = Switch[ er = OptionValue["Errors"],
    Automatic, Abs[data*0.05 + $MinMachineNumber],
    None,  Abs[data*0.0 + $MinMachineNumber] ,
    _?NumberQ,  Abs[data*er + $MinMachineNumber],
    {_?NumberQ, _?NumericQ},  Abs[Transpose@{data[[All, 1]]*er[[1]],  data[[All, 2]]*er[[2]]} + $MinMachineNumber],
    _?VectorQ, Abs[ Transpose@{data[[All, 1]]*0.0 + $MinMachineNumber, N@er}], 
    _?MatrixQ, N@er,
    _, Message[LinearModelFitYork::errBadInput]
  ];

  weights = If[ OptionValue["Weights"] == None,
    (* use errors*)
    Transpose[{ errors[[All,1]]^(-2) , errors[[All,2]]^(-2)}]
    , (*else use provided weights*)
    OptionValue["Weights"]
  ];

  (* check errors *)
  If[ Dimensions[weights] != Dimensions[data], Throw @ Message[LinearModelFitYork::errFormatErrors] ];
  If[ MemberQ[weights, 0.0, 2], Throw@Message[LinearModelFitYork::errErrorsNoZeros]];


  (* get correlations *)
  r = If[ OptionValue["ErrorCorrelations"] === None,
        ConstantArray[0.0, length],
        OptionValue["ErrorCorrelations"]
      ];

  (* check correlations *)

  (*no offset term?*)
  If[ !OptionValue[IncludeConstantBasis],
    dataIn = Append[dataIn, {0,0}];
    weights = Append[weights, {10^100,10^100}];
    r = Append[r, 0];
  ];


  (* compute *)
  result = yorkFitRaw[
              dataIn[[All,1]],
              dataIn[[All,2]],
              weights[[All,1]],
              weights[[All,2]],
              r
           ];

  (* format output *)
  Join[result,<|
    "fit" -> Function[\[FormalX], Evaluate[{result["a"], result["b"]}.{1, \[FormalX]}]],
    "fitFunction" -> Function[\[FormalX], Evaluate[{"a", "b"}.{1, \[FormalX]}]],
    "BestFitParameters" -> {result["a"], result["b"]},
    "Errors" -> {result["\[Sigma]a"], result["\[Sigma]b"]}
    |>]
]





End[ ]

EndPackage[ ]



