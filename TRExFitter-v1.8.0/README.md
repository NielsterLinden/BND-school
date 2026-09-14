# TRExFitter   [![build status](https://gitlab.cern.ch/TRExStats/TRExFitter/badges/master/pipeline.svg "build status")](https://gitlab.cern.ch/TRExStats/TRExFitter/commits/master) [![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14845712.svg "Zenodo DOI")](https://doi.org/10.5281/zenodo.14845712)

This package provides a framework to perform profile likelihood fits.
In addition to that, many convenient features are available.
`TRExFitter` was previously also known as `TtHFitter`.

## Getting TRExFitter

`TRExFitter` is part of the `StatAnalysis` environment.
This makes it quick and easy to start using `TRExFitter`.
All you need is the following:

```bash
setupATLAS
asetup StatAnalysis,0.7.2
```

and you are ready to go!

For more information on the various ways to set up and run `TRExFitter` and an overview of the tool, [check out our website](https://trexfitter-docs.web.cern.ch/) (restricted to ATLAS members)!

## Getting help with TRExFitter

- If you find a bug in `TRExFitter`, please open a [Gitlab issue](https://gitlab.cern.ch/TRExStats/TRExFitter/-/issues/new).
Remember to look through previous issues to see if the problem was already solved!

## Useful links

Here are a few important references to make use of:

* [TRExFitter website](https://trexfitter-docs.web.cern.ch/) collects information in a convenient format
* Make sure to read the [FAQ](https://trexfitter-docs.web.cern.ch/trexfitter-docs/latest/faq/) section and take a look at the [Fit Problems Tutorial twiki](https://twiki.cern.ch/twiki/bin/view/AtlasProtected/FitProblemsTutorial) which describes common issues (only for ATLAS members)
* For syntax highlighting of configuration files in VSCode, you can use the dedicated extension available [here](https://marketplace.visualstudio.com/items?itemName=LukaVomberg.trexfitter-config)

## Citation
If you used TRExFitter in your analysis, please consider citing its zenodo record [10.5281/zenodo.14845712](https://zenodo.org/records/14845712).

You can find an example bibtex entry below.
```
@software{TRExFitter,
  author       = {{ATLAS Collaboration}},
  title        = {TRExFitter},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.14845712},
  url          = {https://doi.org/10.5281/zenodo.14845712},
}
```

## TRExFitter package authors

* Mohamed Aly [mohamed.aly@cern.ch](mohamed.aly@cern.ch)
* Tomas Dado [tomas.dado@cern.ch](tomas.dado@cern.ch)
* Alexander Held [alexander.held@cern.ch](alexander.held@cern.ch)
* Michele Pinamonti [michele.pinamonti@gmail.com](michele.pinamonti@gmail.com)
* Louis-Guillaume Gagnon [louis.guillaume.gagnon@cern.ch](louis.guillaume.gagnon@cern.ch)
* Volker Andreas Austrup [volker.andreas.austrup@cern.ch](volker.andreas.austrup@cern.ch)

## Past contributors

* Loic Valery

## Contributing to TRExFitter

Contributions to TRExFitter are welcome.
Please have a look at [CONTRIBUTING.md](CONTRIBUTING.md) to get started.
