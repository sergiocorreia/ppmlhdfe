## (Appendix) PPML Specifics

Because `ppml` only performs some of the steps required to identify separated observations, it might fail to identify some separated observations. Further, because it uses some extra heuristics, it might incorrectly identify and drop observations that are not separated (and regressors as well).

### I. False positives: some non-separated observations are dropped

In this case, there are no separated observations, but `ppml` drops some observations and/or regressors.

```stata
clear
input double y byte(x1 x2 x3)
        0 0 0 1
        0 0 0 1
        0 0 0 2
2.2596662 1 1 2
2.4177196 2 2 2
0.9788354 3 2 1
2.6114680 4 2 1
end

* No observations are separated
ppmlhdfe y x*, tagsep(sep) zvar(z) r2
assert z == 0

poisson y x* // runs normally

ppml y x* // two regressors incorrectly omitted
ppml y x*, strict // three regressors and one observation incorrectly omitted
```

### II. False negatives: not all separated observations are omitted

In this case, there are four separated observations, but `ppml` does not drop them.

```stata
cls
clear
input double y byte(x1 x2 x3)
  0  2 3 4
  0 -5 0 1
  0  0 2 1
  0  0 0 3
  0 -2 0 1
  0  0 0 2
1.1  1 1 2
 .5  2 2 2
3.3  3 3 2
2.1  4 2 4
end

* Create certificate of separation
gen z = x2 + x3 - x1 - 2 // alternative: ppmlhdfe y x*, tagsep(sep) zvar(z) r2
assert z >=0 if y==0
assert z == 0 if y>0

poisson y x* if z==0 // works once we restrict the sample

ppml y x* // does not detect separated observations
ppml y x*, strict // drops x2, but does not detect the four separated observations
```
