# REGHDFE: Poisson pseudo-likelihood regression with multiple levels of fixed effects

- Current version: ` 2.0.0 28feb2018`
- Jump to: [`references`](#references) [`citation`](#citation) [`install`](#installation)

**ppmlhdfe** is a Stata package that implements Poisson pseudo-maximum likelihood regressions (PPML) with multi-way fixed effects, as described by [Correia, Guimarães, Zylkin (2019a)](http://scorreia.com/research/ppmlhdfe.pdf). The estimator employed is robust to statistical separation and convergence issues, due to the procedures developed in [Correia, Guimarães, Zylkin (2019b)](http://scorreia.com/research/separation.pdf).


## Citation

```bibtex
ADD BIBTEX FILES HERE
```


## References

Quick information on the command can be glanced from the [help file](help-ppmlhdfe.html)

For more detailed information:

- The [`ppmlhdfe` paper](http://scorreia.com/research/ppmlhdfe.pdf) contains more detailed information and examples
- The paper on [statistical separation](http://scorreia.com/research/separation.pdf) discusses the crucial step of solving the separation issue that can lead to incorrect convergence (or no convergence) in Poisson and other GLM models.

For introductory guides on separation, and to how `ppmlhdfe` internally address it, see any of the following documents

- [Separation primer](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_primer.md): a quick practical introduction to separation on Poisson models.
- [Separation benchmarks](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_benchmarks.md): using examples, shows how all common statistical packages are vulnerable to the problem.
- [Undocumented options](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/undocumented.md): advanced and undocumented options of `reghdfe`, that might be useful for advanced users.


## Installation

`ppmlhdfe` requires the latest versions of [`ftools`](https://github.com/sergiocorreia/ftools) and [`reghdfe`](https://github.com/sergiocorreia/reghdfe):

```stata
* Install ftools
cap ado uninstall ftools
net install ftools, from("https://raw.githubusercontent.com/sergiocorreia/ftools/master/src/")

* Install reghdfe
cap ado uninstall reghdfe
net install reghdfe, from("https://raw.githubusercontent.com/sergiocorreia/reghdfe/master/src/")

* Install ppmlhdfe
cap ado uninstall ppmlhdfe
net install ppmlhdfe, from("https://raw.githubusercontent.com/sergiocorreia/ppmlhdfe/master/src/")

* (optional) create compiled files
ftools, compile
reghdfe, compile

* Check versions
ppmlhdfe, version

* Test program
sysuse auto, clear
reghdfe price weight, a(turn)
ppmlhdfe price weight, a(turn)
```
