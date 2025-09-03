How to use
==========

There are two ways to use ``BowshockPy``:

1. *The quick way*: running it from the command-line, using an input file that contains all the input parameters. 
2. *The flexible way*: Importing ``BowshockPy`` as a package in your Python code.


From the command-line using an input file
------------------------------------------

The easiest way to use ``BowshockPy`` is using an input file that contains all the parameters of the model to be generated. This documentation provide four :doc:`examples of input files<../examples/examples_inputfile>`. Also, you can generate these examples of input file running ``BowshockPy`` with ``--print`` flag and specifying the example you are interested in. For example,

.. code-block:: console

  $ bowshockpy --print example1.py

will print example 1 to your working directory. Then, you can modify the parameters according to your scientific needs. See :doc:`parameters<inputparams>` for a description of the parameters.

After defining the parameters in your input file, you can generate your models using ``--read`` flag and providing the name of your customized input file

.. code-block:: console

  $ bowshockpy --read inputfile.py 



Importing ``BowshockPy`` package
--------------------------------

This the most flexible way to use ``BowshockPy``. The basic functioning can be performed by for classes, that can be directly imported from the package 

>>> from bowshockpy import BowshockModel, ObsModel, BowshockCube, CubeProcessing

In a nutshell, these are the main tasks of each class:

- **BowshockModel**: This class generates the analytic momentum-conserving bowshock model.
- **ObsModel**: Allows to project the kinematics and morphology of the bowshock model.
- **BowshockCube**: Computes the masses in each pixel and velocity channel of the spectral cube.
- **CubeProcessing**: Performs the radiative transfer, obtaining the intensities of the molecule rotational transition. It also allows to convolve the cube and obtain position-velocity diagrams and moment images.

This documentation includes a :doc:`notebook tutorial <../examples/notebook_tutorial>` that shows how to use these classes. For a more detailed description of how to use these classes and other functions available in ``BowshockPy`` package, see the :doc:`API Reference <../api/index>`.

..
    Using ``BowshockPy`` as a package allows you to either load the model parameters from an input file or to define the parameters in you script and create the bowshock model. The input file that contains all the model parameters, "params.py" can be read in the following manner. 
