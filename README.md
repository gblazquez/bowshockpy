# bowshockpy

*A Python package for the generation of synthetic spectral channel maps of a jet-driven bowshock model*

This program computes spectral channel maps of jet-driven bowshock model. The bowshock shell morphology and kinematics are determined from the momentum conservation in the interaction of jet material ejected sideways by an internal working surface and the ambient medium (or a surrounding disk wind moving in the jet axis direction). Well mixing between the jet and ambient material are assumed.

## Documentation

An extensive documentation on BowshockPy can be found [here](https://bowshockpy.readthedocs.io/en/latest/)

## Requirements
bowshockpy requires:

* Python3 
* astropy
* matplotlib
* numpy
* scipy 

It has been tested with `python == 3.10`, but it could work with previous versions.

## Installation

You can install ``bowshockpy`` using ``pip``, either from PyPI or from the source repository. The will be automatically installed.

### Installation from PyPI

You can get ``bowshockpy`` from PyPI using pip:

```bash
(.venv) $ pip install bowshockpy 
```


### Installation from source repository

You can install ``bowshockpy`` from the source by cloning this repository:

```bash
$ git clone https://github.com/gblazquez/bowshockpy.git 
$ cd bowshockpy
$ pip install .
```

## Usage

There are two ways to use bowshockpy, either using a configuration file with the model parameters to be generated or importing bowshockpy as a package and run the model manually.


### From the command-line using an input file

The easiest way to use ``bowshockpy`` is using an input file that contains all the parameters of the model to be generated. You can tell ``bowshockpy`` to read a input file and generate your models either from terminal

```bash
(.venv) $ bowshockpy --read params.py 
```

If you want to use an example of an input file, you can print some examples. If you want to print example 1:

```bash
(.venv) $ bowshockpy --print_example 1
```

Then, you can modify the example file to your needs. 


### Importing bowshockpy module

This the most flexible way to use ``bowshockpy``. In order to import ``bowshockpy`` as a python package:

```python
import bowshockpy as bp
```

Using bowshockpy as a package allows you to either load the model parameters from an input file or to define the parameters in you script and create the bowshock model. The input file that contains all the model parameters, "params.py" can be read in the following manner. 

```python
bp.generate_bowshock("params.py")
```

If you would like to print an example of the input file

```python
bp.print_example("1")
```

See the [documentation](https://bowshockpy.readthedocs.io/en/latest/) for more details on the modular usage of bowshockpy

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
