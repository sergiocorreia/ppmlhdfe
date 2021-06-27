# `ppmlhdfe`: Poisson pseudo-likelihood regression with multiple levels of fixed effects

![GitHub release (latest by date)](https://img.shields.io/github/v/release/sergiocorreia/ppmlhdfe?label=last%20version)
![GitHub Release Date](https://img.shields.io/github/release-date/sergiocorreia/ppmlhdfe)
![GitHub commits since latest release (by date)](https://img.shields.io/github/commits-since/sergiocorreia/ppmlhdfe/latest)
![StataMin](https://img.shields.io/badge/stata-%3E%3D%2012.1-blue)
[![DOI](https://zenodo.org/badge/63449974.svg)](https://zenodo.org/badge/latestdoi/63449974)

- Jump to: [`citation`](#citation) [`references`](#references) [`install`](#installation)
- Also see: [`ppmlhdfe` Paper](http://scorreia.com/research/ppmlhdfe.pdf) | [Separation Paper](http://scorreia.com/research/separation.pdf) | [Help File](http://scorreia.com/help/ppmlhdfe.html) | [Separation Primer](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_primer.md) | [Separation Benchmarks](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/separation_benchmarks.md) | [Undocumented](https://github.com/sergiocorreia/ppmlhdfe/blob/master/guides/undocumented.md)


**ppmlhdfe** is a Stata package that implements Poisson pseudo-maximum likelihood regressions (PPML) with multi-way fixed effects, as described in [Correia, Guimarães, Zylkin (2019a)](http://scorreia.com/research/ppmlhdfe.pdf). The estimator employed is robust to statistical separation and convergence issues, due to the procedures developed in [Correia, Guimarães, Zylkin (2019b)](http://scorreia.com/research/separation.pdf).

## Recent updates

- **Version 2.3 27jun2021:** minor changes due to reghdfe's v6 update. Currently, ppmlhdfe is still using the code from reghdfe v5, which the new version ships with. A port is planned at some point in the future, but because some Mata functions changed their behavior, this needs to be done carefully.
- **Version 2.2 02aug2019:** major speedups due to improved IRLS acceleration (see [page 7 of the paper](https://arxiv.org/pdf/1903.01690.pdf)) and due to faster separation checks.
- **Version 2.1 04apr2019:** added experimental [step-halving](https://journal.r-project.org/archive/2011/RJ-2011-012/index.html). Not as useful for Poisson models as with other GLMs, so it's turned off by default. You can enable it by including the option `use_step_halving(1)`. Other options you can set are `step_halving_memory(0.9)` and `max_step_halving(2)` (default values in parenthesis).


## Citation

[(Download BibTex file here)](https://raw.githubusercontent.com/sergiocorreia/ppmlhdfe/master/ppmlhdfe.bib)

#### As text

<ul>
<li>
Sergio Correia, Paulo Guimarães, Thomas Zylkin: “Verifying the existence of maximum likelihood estimates for generalized linear models”, 2019; <a href='http://arxiv.org/abs/1903.01633'>arXiv:1903.01633</a>.
</li>
</ul>

<ul>
<li>
Sergio Correia, Paulo Guimarães, Thomas Zylkin: “ppmlhdfe: Fast Poisson Estimation with High-Dimensional Fixed Effects”, 2019; <a href='http://arxiv.org/abs/1903.01690'>arXiv:1903.01690</a>.
</li>
</ul>

#### As BibTex

```bibtex
@Misc{ExistenceGLM,
  Author = {Correia, Sergio and Guimar{\~a}es, Paulo and Zylkin, Thomas},
  Title = {Verifying the existence of maximum likelihood estimates for generalized linear models},
  Year = {2019},
  Eprint = {arXiv:1903.01633},
}

@Misc{ppmlhdfe,
  Author = {Correia, Sergio and Guimar{\~a}es, Paulo and Zylkin, Thomas},
  Title = {{ppmlhdfe: Fast Poisson Estimation with High-Dimensional Fixed Effects}},
  Year = {2019},
  Eprint = {arXiv:1903.01690},
}
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

`ppmlhdfe` requires the latest versions of [`ftools`](https://github.com/sergiocorreia/ftools) and [`reghdfe`](https://github.com/sergiocorreia/reghdfe).

To install stable versions from SSC:

```stata
cap ado uninstall ftools
cap ado uninstall reghdfe
cap ado uninstall ppmlhdfe

ssc install ftools
ssc install reghdfe
ssc install ppmlhdfe

clear all
ftools, compile
reghdfe, compile

* Test program
sysuse auto, clear
reghdfe price weight, a(turn)
ppmlhdfe price weight, a(turn)
```

To install the latest versions from Github:

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

* Create compiled files
ftools, compile
reghdfe, compile

* Check versions
ppmlhdfe, version

* Clear programs already in memory
program drop _all

* Test program
sysuse auto, clear
reghdfe price weight, a(turn)
ppmlhdfe price weight, a(turn)
```
