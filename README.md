# bowshockpy

*A Python package for generating spectral channel maps of a jet-driven bowshocks model*

``Bowshockpy`` is a Python package that generates synthetic spectral cubes, position-velocity diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in [Ostriker et al. (2001)](https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract) and [Tabone et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract).

<!--
 computes spectral channel maps of jet-driven bowshock model. The bowshock shell morphology and kinematics are determined from the momentum conservation in the interaction of jet material ejected sideways by an internal working surface and the ambient medium (or a surrounding disk wind moving in the jet axis direction). Well mixing between the jet and ambient material are assumed.
-->

## Documentation

An extensive documentation on ``bowshockpy`` can be found [here](https://bowshockpy.readthedocs.io/en/latest/)


## Requirements
bowshockpy requires:

* Python3 
* astropy
* matplotlib
* numpy
* scipy 

It has been tested with `python == 3.10`, but it could work with previous versions.


## Installation

You can install ``bowshockpy`` from PyPI. 

```bash
(.venv) $ pip install bowshockpy 
```

## Usage

The easiest way to use ``bowshockpy`` is using an input file that contains all the parameters of the model to be generated. You can tell ``bowshockpy`` to read a input file and generate your models either from terminal

```bash
(.venv) $ bowshockpy --read params.py 
```

If you want to use an example of an input file, you can print some examples. If you want to print example 1:

```bash
(.venv) $ bowshockpy --print-example 1
```

For a more flexible use, you use ``bowshockpy`` as a python module

```python
import bowshockpy as bp

bp.generate_bowshock("params.py")
```

See the [documentation](https://bowshockpy.readthedocs.io/en/latest/) for more details on the modular usage of bowshockpy.

## Contributing

If you are interested in contributing, see [contributing](CONTRIBUTING)

## License

This project is licensed under the MIT License. For details see the [LICENSE](LICENSE).


## Citation

```tex
@software{gblazquez2025,
  author    = {Blazquez-Calero, Guillermo AND et al.},
  title     = {{BowshockPy}: A Python package for the generation of synthetic spectral channel maps of a jet-driven
bowshock model},
  year      = {2025},
  version   = {0.1.0},
  url       = {https://github.com/gblazquez/bowshockpy}
}
```
