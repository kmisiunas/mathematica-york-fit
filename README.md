# mathematica-york-fit
Linear model fit for data with X and Y errors using [York 2004 method](http://aapt.scitation.org/doi/pdf/10.1119/1.1632486)

## Usage

``` 
YorkFit[{{x,y},...}, \"Errors\"-> err] York linear regression fit:
  http://aapt.scitation.org/doi/pdf/10.1119/1.1632486
  
  "Errors" can be supplied with:
    Automatic - assumes 5% error variation in X and Y
    None - no errors in x and y
    p_NumberQ - assumes the erros are p fraction of the data
    {pX_NumberQ, pY_NumberQ} - assumes the erros are {pX,pY} fraction of the data
    {errY1,errY2,...} - errors only in Y
    {{errX1,errY1},...} - errors in x and y
  
  Weights -> None
    Allows to manually specify weight matrix

  "ErrorCorrelations" -> None
    Correlations between X and Y errors

  IncludeConstantBasis -> True
    if False, this leads to fit y=b*x
    achieved by adding {0,0} point with very high weights 

  Output: <|"fit" -> best fit, "a" -> offset, "b" -> slope, "\[Sigma]a" -> error of offset, 
    "\[Sigma]b" -> error of slope |>"
```
