# BowshockPy

_A Python package for generating spectral channel maps of a jet-driven bowshock model_

[![PyPI - Version](https://badge.fury.io/py/bowshockpy.svg)](https://pypi.org/project/bowshockpy/)
[![PyPI - Python Version](https://img.shields.io/python/required-version-toml?tomlFilePath=https://raw.githubusercontent.com/gblazquez/bowshockpy/main/pyproject.toml)](https://pypi.org/project/bowshockpy/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://github.com/gblazquez/bowshockpy/blob/main/LICENSE)
[![Tests](https://github.com/gblazquez/bowshockpy/actions/workflows/tests.yml/badge.svg)](https://github.com/gblazquez/bowshockpy/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/gblazquez/bowshockpy/graph/badge.svg?token=EGA6WEJCYR)](https://codecov.io/gh/gblazquez/bowshockpy)
[![Documentation Status](https://app.readthedocs.org/projects/bowshockpy/badge/?version=latest)](https://bowshockpy.readthedocs.io/en/latest/)

`BowshockPy` is a Python package that generates synthetic spectral cubes, position-velocity diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in [Ostriker et al. (2001)](https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract) and [Tabone et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract). The software computes column density and is able to calculate the intensities of low-J rotational transitions of a linear molecule such as CO, providing mock observations of the emission that radio interferometers as ALMA are able to detect at millimeter wavelengths.

## Documentation

An extensive documentation on `BowshockPy` can be found [here](https://bowshockpy.readthedocs.io/en/latest/).

## Requirements

`BowshockPy` requires:

- Python3
- astropy
- matplotlib
- numpy
- photutils
- scipy

It has been tested in python versions 3.8-3.12, but it could work with other versions also.

## Installation

You can install `BowshockPy` from PyPI.

```bash
$ pip install bowshockpy
```

## How to use

There are two different ways to use `BowshockPy`:

1. Run it from the terminal specifying an input file: Use an [example of input file](https://github.com/gblazquez/bowshockpy/tree/main/examples) and modify the input parameters according your scientific goals. Then, run `BowshockPy` in your terminal:

```bash
$ bowshockpy -r inputfile.py
```

2. Import `BowshockPy` package in your Python code: We include an [notebook tutorial](https://github.com/gblazquez/bowshockpy/tree/main/examples/notebook_tutorial.ipynb) that shows how to use the main classes.

See the [documentation](https://bowshockpy.readthedocs.io/en/latest/) for more details on the usage of `BowshockPy`.

## Contributing

If you are interested in contributing, see [contributing](CONTRIBUTING.md)

## License

This project is licensed under the MIT License. For details see the [LICENSE](LICENSE).

## Citation

```tex
@software{gblazquez2025,
  author    = {Blazquez-Calero, Guillermo AND et al.},
  title     = {{BowshockPy}: A Python package for the generation of synthetic spectral channel maps of a jet-driven
bowshock model},
  year      = {2025},
  version   = {0.3.1},
  url       = {https://github.com/gblazquez/bowshockpy}
}
```
