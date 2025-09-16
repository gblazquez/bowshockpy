How to use
==========

There are two ways to use ``BowshockPy``:

1. Run it from the command-line, using an input file that contains all the input parameters (see :doc:`examples of input files<../examples/examples_inputfile>`). 
2. Import ``BowshockPy`` as a package in your Python code (see :doc:`notebook tutorial <../examples/notebook_tutorial>`).


From the command-line using an input file
------------------------------------------

You can run ``BowshockPy`` using an input file that contains all the parameters of the model to be generated. The documentation provides several :doc:`examples of input files<../examples/examples_inputfile>`, which can also be downloaded from the `GitHub repository <https://github.com/gblazquez/bowshockpy>`_. Alternatively, you can generate these examples of input files running ``BowshockPy`` with ``--print`` flag and specifying the example you are interested in. For example,

.. code-block:: console

  $ bowshockpy --print example1.py

will print example 1 to your working directory. Then, you can modify the parameters according to your scientific needs. See :doc:`Input file parameters section<inputparams>` for a description of the parameters.

After defining the parameters in your input file, you can generate your models using ``--read`` flag and providing the name of your customized input file

.. code-block:: console

  $ bowshockpy --read inputfile.py 



Importing ``BowshockPy`` package
--------------------------------

This the most flexible way to use ``BowshockPy``. The basic functioning can be performed by for classes, that can be directly imported from the package 

>>> from bowshockpy import BowshockModel, ObsModel, MassCube, CubeProcessing

In a nutshell, the main tasks of each class are the following:

- **BowshockModel**: This class generates the analytic momentum-conserving bowshock model.
- **ObsModel**: Allows to project the kinematics and morphology of the bowshock model.
- **MassCube**: Computes the masses in each pixel and velocity channel of the spectral cube.
- **CubeProcessing**: Computes the column densities and performs the radiative transfer, obtaining the intensities of a linear molecule rotational transition assuming the rigid rotor approximation (although custom models to derive the intensities from the column densities can be implemented by the user and passed to CubeProcessing class). It also allows to convolve the cube and obtain position-velocity diagrams and moment images.

This documentation includes a :doc:`notebook tutorial <../examples/notebook_tutorial>` that shows how to use these classes. For a more detailed description of how to use these classes and other functions available in ``BowshockPy`` package, see the :doc:`API Reference <../api/index>`.

..
    Using ``BowshockPy`` as a package allows you to either load the model parameters from an input file or to define the parameters in you script and create the bowshock model. The input file that contains all the model parameters, "params.py" can be read in the following manner. 
