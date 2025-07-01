# bowshockpy

*Spectral channel maps generator of a bowshock model*

This program computes spectral channel maps of jet-driven bowshock model. The bowshock shell morphology and kinematics are determined from the momentum conservation in the interaction of jet material ejected sideways by an internal working surface and the ambient medium (or a surrounding disk wind moving in the jet axis direction). Well mixing between the jet and ambient material are assumed.

## Requirements
bowshockpy requires:

* Python3 
* astropy
* matplotlib
* numpy
* scipy 

It has been tested with `python == 3.10`, but it could work with previous versions.

## Installation

Within a python evironment:

```bash
(pyenv)$ pip install -r requirements.txt
```

## How to use bowshockpy

- One should set the right parameters in bowshock_params.py
- Run generate_bowshock.py:
 ```bash
 (pyenv)$ python generate_bowshock.py
 ```
- Your model will be saved in model folder

