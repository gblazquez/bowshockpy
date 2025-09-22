

==========
BowshockPy
==========

|PyPI - Version| |Python - Version| |PyPI - Python license| |Tests| |codecov|
|Documentation Status| 

.. |PyPI - Version| image:: https://badge.fury.io/py/bowshockpy.svg
   :target: https://pypi.org/project/bowshockpy/
.. |Python - Version| image:: https://img.shields.io/python/required-version-toml?tomlFilePath=https://raw.githubusercontent.com/gblazquez/bowshockpy/main/pyproject.toml
   :target: https://pypi.org/project/bowshockpy/
.. |PyPI - Python license| image:: https://img.shields.io/pypi/l/ansicolortags.svg
   :target: https://github.com/gblazquez/bowshockpy/blob/main/LICENSE
.. |Tests| image:: https://github.com/gblazquez/bowshockpy/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/gblazquez/bowshockpy/actions/workflows/tests.yml
.. |codecov| image:: https://codecov.io/gh/gblazquez/bowshockpy/graph/badge.svg?token=EGA6WEJCYR
   :target: https://codecov.io/gh/gblazquez/bowshockpy
.. |Documentation Status| image:: https://app.readthedocs.org/projects/bowshockpy/badge/?version=latest
   :target: https://bowshockpy.readthedocs.io/en/latest/


``BowshockPy`` is an open-source Python package for generating synthetic spectral cubes, position-velocity diagrams, and moment images of an analytical momentum-conserving bowshock model driven by a protostellar jet. The software computes the column densities of the bowshock shell and is capable of estimating the intensities of low-J rotational transitions of linear molecules (such as CO), providing mock observations of the emission that radio interferometers as ALMA are able to detect at millimeter wavelengths.

..
   .. note::
   
      This project is under development.
   

Getting started
---------------

You can :doc:`install<user_guide/installation>` it from PyPI using pip: 

.. code-block:: console

  $ pip install bowshockpy 

There are two different ways to :doc:`use<user_guide/usage>` the software:

1. Run it from the terminal specifying an input file: Use an :doc:`example of input file<examples/examples_inputfile>` and modify the :doc:`input parameters<user_guide/inputparams>` according your scientific goals. Then, run ``BowshockPy`` in your terminal:

   .. code-block:: console
   
     $ bowshockpy -r inputfile.py 

2. Import ``BowshockPy`` package in your Python code: We include an :doc:`notebook tutorial<examples/notebook_tutorial>` that shows how to use the main classes. Although the notebook contains an explanation of their basic functioning, see the :doc:`API Reference<api/index>` for a detailed documentation.

A description of the outputs of ``BowshockPy`` can be found in the :doc:`output section<user_guide/outputs>`.

If you are interested in the physics behind ``BowshockPy`` an its workflow, see the :doc:`background section<user_guide/background>`.


Contributing
------------

We welcome contributions and issue reports to this project. If you are interested, please follow our `contributing guidelines <https://github.com/gblazquez/bowshockpy/blob/main/CONTRIBUTING.md>`_.


License
-------

This project is licensed under the `MIT License <https://github.com/gblazquez/bowshockpy/blob/main/LICENSE>`_. 


Citation
--------

If you use ``BowshockPy`` in your research work, please cite it as

.. code-block:: tex

  @software{gblazquez2025,
    author    = {Blazquez-Calero, Guillermo AND et al.},
    title     = {{BowshockPy}: A Python package for the generation of synthetic spectral channel maps of a jet-driven
  bowshock model},
    year      = {2025},
    version   = {1.0.0},
    url       = {https://github.com/gblazquez/bowshockpy}
  }

Also, please cite its :doc:`dependencies<user_guide/dependencies>`.


Table of contents
-----------------

.. toctree::
   :caption: User Guide
   :maxdepth: 2

   user_guide/installation
   user_guide/dependencies
   user_guide/usage
   user_guide/inputparams
   user_guide/outputs
   user_guide/background

.. toctree::
   :caption: Tutorials and Examples
   :maxdepth: 2
   
   examples/notebook_tutorial
   examples/examples_inputfile

.. toctree::
   :caption: Code documentation
   :maxdepth: 2

   api/index
 


