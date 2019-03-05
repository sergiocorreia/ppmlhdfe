# `ppmlhdfe`: Poisson pseudo-likelihood regression with multiple levels of fixed effects

- Current version: ` 2.0.1 03mar2019`
- Jump to: [`references`](#references) [`citation`](#citation) [`install`](#installation)

**ppmlhdfe** is a Stata package that implements Poisson pseudo-maximum likelihood regressions (PPML) with multi-way fixed effects, as described in [Correia, Guimarães, Zylkin (2019a)](http://scorreia.com/research/ppmlhdfe.pdf). The estimator employed is robust to statistical separation and convergence issues, due to the procedures developed in [Correia, Guimarães, Zylkin (2019b)](http://scorreia.com/research/separation.pdf).


## Citation

```bibtex
ADD BIBTEX FILES HERE
```


## References

Quick information on the command can be glanced from the [help file](http://scorreia.com/help/ppmlhdfe.html).

For detailed information:

- The [ppmlhdfe paper](http://scorreia.com/research/ppmlhdfe.pdf) explains the command in depth, provides examples, etc.
- The paper on [statistical separation](http://scorreia.com/research/separation.pdf) discusses the crucial step of solving the separation issue, that can otherwise lead to incorrect convergence (or no convergence) in Poisson and other GLM models.

For introductory guides on separation, and on how `ppmlhdfe` internally address it, see the following documents:

- [Separation primer](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_primer.md): a quick practical introduction to separation in Poisson models.
- [Separation benchmarks](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_benchmarks.md): shows how separation affects all common statistical packages.
- [Undocumented options](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/undocumented.md): this pages briefly lists otherwise undocumented options of `ppmlhdfe`, which might be useful for advanced users.


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
