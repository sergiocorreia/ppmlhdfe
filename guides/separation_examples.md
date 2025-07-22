# Examples of separation in Poisson, Logit, and Multinomial Logit models

- About `ppmlhdfe`: [Github Readme](https://github.com/sergiocorreia/ppmlhdfe/tree/master?tab=readme-ov-file#ppmlhdfe-poisson-pseudo-likelihood-regression-with-multiple-levels-of-fixed-effects) | [Working Paper](https://arxiv.org/abs/1903.01690) | [Stata Journal](https://doi.org/10.1177/1536867X20909691) | [Help File](http://scorreia.com/help/ppmlhdfe.html) | [Undocumented Options](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/undocumented.md)
- About Separation: [Working Paper](https://arxiv.org/abs/1903.01633) | [Primer](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_primer.md) | [Examples](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_primer.md) [Software Benchmarks](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_benchmarks.md) | [Further Reading](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/further_reading.md)
- Sections: [UPDATETHIS](#why-might-ml-estimates-not-exist) | [When?](#when-does-this-issue-occur) | [Solutions](#how-does-ppmlhdfe-actually-find-the-separated-observations) | [∞](#to-infinity-and-beyond)

*(These examples complement [Verifying the existence of maximum likelihood estimates for generalized linear models](https://arxiv.org/abs/1903.01633); please see the links above for related guides.)*

This section discusses several

## Logit / Logistic


### Agresti Eight-Point Example



### Agresti Ten-Point Example


### Endometrial Example

- Source: Agresti (2015) - Section 5.7.1 - Page 186
- For an earlier use case, see also Heinze and Schemper (2002)
- For an R implementation using the simplex method that also solves this example, see [Kosmidis and Schumacher (2021)](https://cran.r-project.org/web/packages/detectseparation/vignettes/separation.html).



## Multinomial Logit

To solve the Multinomial Logit model, we take advantage of the Multinomial-Poisson equivalence, also known as the Multinomial-Poisson "trick" (see [Palmgreen 1981](https://www.jstor.org/stable/2335606), [Baker 1994](https://www.jstor.org/stable/2348134), [Guimarães 2004](https://journals.sagepub.com/doi/abs/10.1177/1536867X0400400304), amongst others).



### Alligators ([Stata do-file](./code/example-aligators.do))

This example on the food of choice of alligators comes from Agresti (2002, Table 7.4). The corresponding regression discussed by Agresti has finite estimates, but [Kosmidis (2017)](https://cran.r-project.org/web/packages/brglm2/vignettes/multinomial.html) shows how a simple sparsification of the dataset, achieved by updating the frequency weight to 1/3 of their original values, leads to infinite estimates.

In particular:

```stata
. replace freq = round(freq/3)
(39 real changes made)

. mlogit foodchoice ib2.size ib4.lake [fw=freq], tol(1e-12) nolog

Multinomial logistic regression                         Number of obs =     69
                                                        LR chi2(16)   =  25.49
                                                        Prob > chi2   = 0.0616
Log likelihood = -80.956515                             Pseudo R2     = 0.1360

------------------------------------------------------------------------------
  foodchoice | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-------------+----------------------------------------------------------------
Fish         |  (base outcome)
-------------+----------------------------------------------------------------
// omitted
-------------+----------------------------------------------------------------
Reptile      |
        size |
      <=2.3  |  -1.056638    1.28096    -0.82   0.409    -3.567274    1.453998
             |
        lake |
    Hancock  |   15.60054   2003.427     0.01   0.994    -3911.043    3942.245
   Oklawaha  |   16.29274   2003.427     0.01   0.994    -3910.351    3942.937
   Trafford  |   16.42811   2003.427     0.01   0.993    -3910.216    3943.072
             |
// omitted
```

As we can see here, the Reptile category as infinite estimates for the coefficients of each lake (values of 15, etc. approach infinite once exponenciated).

To estimate it with our algorithm we will first appply the Logit-Poisson trick:

```stata
gen obs = _n // Record original obs.
expand 5 // Transform the dataset to long (make five copies as dependent variable has five outcomes)
bys obs: gen byte category = _n // Record each of the five categories
bys obs: gen byte y = foodchoice==category // New outcome variable
```

Then, we can apply our algorithm while clustering for hte original observations:

```stata
. ppmlhdfe y ib1.category##(ib2.size ib4.lake) [fw=freq], a(obs) nocons sep(relu)
note: 4 variables omitted because of collinearity: 1bn.size 1bn.lake 2bn.lake 3bn.lake
(ReLU method dropped 17 separated observations in 92 iterations)
note: 2 variables omitted because of collinearity: 3.category#3.lake 4.category#2.lake
Iteration 1:   deviance = 1.9461e+02  eps = .         iters = 1    tol = 1.0e-04  min(eta) =  -3.31  P   
// omitted
Converged in 6 iterations and 6 HDFE sub-iterations (tol = 1.0e-08)


HDFE PPML regression                              No. of obs      =        170
Absorbing 1 HDFE group                            Residual df     =        136
                                                  Wald chi2(12)   =      21.79
Deviance             =  89.23908325               Prob > chi2     =     0.0399
Log pseudolikelihood = -78.61954162               Pseudo R2       =     0.1139
----------------------------------------------------------------------------------------
                       |               Robust
                     y | Coefficient  std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
// omitted
         category#lake |
 Invertebrate#Hancock  |  -2.601252   1.150367    -2.26   0.024     -4.85593   -.3465742
      Reptile#Hancock  |  -.8821981   1.017438    -0.87   0.386     -2.87634    1.111944
         Bird#Hancock  |   .7881922   1.262154     0.62   0.532    -1.685585    3.261969
        Other#Hancock  |   .2364445   .9539814     0.25   0.804    -1.633325    2.106214
// omitted
```

From this output, we can see three things:

1. The ReLu method dropped 17 observations as they were perfectly separated.
2. After dropping these observations, three variables with infinite estimates are now omitted due to collinearity
3. Other estimates, which did not had infinite estimates and were previously reported incorrect values, now report the correct estimates (the same values one would obtain had the software package been able to correctly estimate the respective infinite coefficients).

<!-- 
We can also obtain Geyer's (2007) "direction of recessions":

```stata
. ppmlhdfe y ib1.category##(ib2.size ib4.lake) [fw=freq], a(obs) nocons tagsep(sep) zvar(z) r2

```
WE NEED TO ADJUST THE CODE SO WE DROP ALL ROWS OF THE OBS
-->

## Poisson

### Geyer's contingency table example

Below we reproduce Table 1, Example 2.3 of Geyer (2009):

```stata
import delimited using "http://www.stat.umn.edu/geyer/gdor/catrec.txt", delim(" ") clear
ppmlhdfe y i.(v*)#i.(v*)#i.(v*) , tagsep(sep) zvar(z) r2 // Get certificate of separation Z, and regress it against the Xs
* Code below is just to present a prettier output:
matrix b = e(b)
mata: vars = st_matrixcolstripe("b")
mata: directions = round(st_matrix("b"), 0.001)'
mata: idx = selectindex(directions)
mata: (vars, strofreal(directions))[idx, .]
```

<p align="center"><img src="./figures/primer_geyer.png" alt="screenshot-poisson" width="60%"/></p>

As we can see, we are able to recover Geyer's "direction of recession" by employing the IR algorithm, which has the added advantage of being easy to implement, and not requiring exact algebra routines.


### Seventeen Examples
