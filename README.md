# BowshockPy

_A Python package for generating spectral channel maps of a jet-driven bowshock model_

`BowshockPy` is a Python package that generates synthetic spectral cubes, position-velocity diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in [Ostriker et al. (2001)](https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract) and [Tabone et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract). The software computes the intensities of low-J rotational transitions of of a linear molecule as CO, providing mock observations of the emission that radio interferometers as ALMA are able to detect at millimeter wavelengths.

<!--
 computes spectral channel maps of jet-driven bowshock model. The bowshock shell morphology and kinematics are determined from the momentum conservation in the interaction of jet material ejected sideways by an internal working surface and the ambient medium (or a surrounding disk wind moving in the jet axis direction). Well mixing between the jet and ambient material are assumed.
-->

## Documentation

An extensive documentation on `BowshockPy` can be found [here](https://bowshockpy.readthedocs.io/en/latest/)

## Requirements

`BowshockPy` requires:

- Python3
- astropy
- matplotlib
- numpy
- scipy

It has been tested in python versions 3.9-3.12, but it could work with other versions also.

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

2. Importing `BowshockPy` package in your Python code: We include an [notebook tutorial](https://github.com/gblazquez/bowshockpy/tree/main/examples/notebook_tutorial.ipynb) that shows how to use the main classes.

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
  version   = {0.2.6},
  url       = {https://github.com/gblazquez/bowshockpy}
}
```
