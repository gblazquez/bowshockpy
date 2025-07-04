How to use
==========

There are two ways to use ``bowshockpy``, either using a configuration file with the model parameters to be generated or importing ``bowshockpy`` as a package and run the model manually.


From the command-line using an input file
-------------------------------------------

The easiest way to use ``bowshockpy`` is using an input file that contains all the parameters of the model to be generated. You can tell ``bowshockpy`` to read a input file and generate your models either from terminal

.. code-block:: console

  (.venv) $ bowshockpy --read params.py 

You can download some :doc:`examples<examples_inputfile>` of input files. See :doc:`parameters<inputparams>` for a description of the parameters.


Importing bowshockpy module
---------------------------------------------

This the most flexible way to use ``bowshockpy``. In order to import ``bowshockpy`` as a python package:

>>> import bowshockpy as bp

Using bowshockpy as a package allows you to either load the model parameters from an input file or to define the parameters in you script and create the bowshock model. The input file that contains all the model parameters, "params.py" can be read in the following manner. 

>>> bp.generate_bowshock("params.py")

If you would like to print an example of the input file

>>> bp.print_example("1")



