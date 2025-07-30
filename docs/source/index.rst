..
   .. toctree::
      :hidden:
      :caption: User Guide
      :maxdepth: 2
   
      user_guide/installation
      user_guide/dependencies
      user_guide/usage
      user_guide/inputparams
      user_guide/outputs
      user_guide/background
   
   .. toctree::
      :hidden:
      :caption: Examples
      :maxdepth: 2
      
      examples/example_notebook
      examples/examples_inputfile
   
   .. toctree::
      :hidden:
      :caption: Code documentation
      :maxdepth: 2

    api/index
 
..
   .. toctree::
      :hidden:
      :caption: Background
      :maxdepth: 1
      user_guide/background

..
   .. toctree::
      :hidden:
      :caption: Caption
      :maxdepth: 2
   
      user_guide/index
      examples/index
      api/index

..
  .. image:: images/logo.png
      :width: 800
      :alt: a description of the logo

..
   ``bowshockpy`` is a Python package that generates synthetic spectral cubes, PV diagrams, and moment images for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in `Ostriker et al. (2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract>`_ and `Tabone et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract>`_. 

==========
BowshockPy
==========

``bowshockpy`` is an open-source Python package for generating synthetic spectral cubes, PV diagrams, and moment images of an analytical momentum-conserving bowshock model driven by an protostellar jet. The software computes the intensities of low-J rotational transitions of the CO molecule, providing mock observations of the CO emission that radio telescopes as ALMA are able to detect at millimeter wavelengths.

..
   .. note::
   
      This project is under development.
   

Getting started
---------------

You can :doc:`install<user_guide/installation>` it from PyPI using pip: 

.. code-block:: console

  $ pip install bowshockpy 

There are two different ways to :doc:`use<user_guide/usage>` the software:

1. Run it from the terminal specifying an input file: Use an :doc:`example of input file<examples/examples_inputfile>` and modify the :doc:`input parameters<user_guide/inputparams>` according your scientific goals. Then, run ``bowshockpy`` in your terminal:

   .. code-block:: console
   
     $ bowshockpy -r inputfile.py 

2. Importing ``bowshockpy`` package in your Python code: We include an :doc:`example of a notebook<examples/example_notebook>` that shows how to use the main classes. Although the notebook contains an explanation of their basic functioning, see the :doc:`API Reference<api/index>` for a detailed documentation.

A description of the outputs of ``bowshockpy`` can be found in the :doc:`output section<user_guide/outputs>`.

If you are interested in the physics behind ``bowshockpy`` an its workflow, see the :doc:`background section<user_guide/background>`.


Contributing
------------

We welcome contributions and issue reports to this project. If you are intersted, please follow our `contributing guidelines <https://github.com/gblazquez/bowshockpy/blob/main/CONTRIBUTING.md>`_.


License
-------

This project is licensed under the `MIT License <https://github.com/gblazquez/bowshockpy/blob/main/LICENSE>`_. 


Citation
--------

If you use ``bowshockpy`` in your research work, please cite it as

.. code-block:: tex

  @software{gblazquez2025,
    author    = {Blazquez-Calero, Guillermo AND et al.},
    title     = {{BowshockPy}: A Python package for the generation of synthetic spectral channel maps of a jet-driven
  bowshock model},
    year      = {2025},
    version   = {0.2.2b},
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
   :caption: Examples
   :maxdepth: 2
   
   examples/example_notebook
   examples/examples_inputfile

.. toctree::
   :caption: Code documentation
   :maxdepth: 2

   api/index
 



..
   for a simple analytical jet-driven bowshock model, using the prescription for protostellar jets presented in `Ostriker et al. (2001) <https://ui.adsabs.harvard.edu/abs/2001ApJ...557..443O/abstract>`_ and `Tabone et al. (2018) <https://ui.adsabs.harvard.edu/abs/2018A%26A...614A.119T/abstract>`_. 

..
   .. grid:: 1 2 2 2
      :gutter: 2
   
      .. grid-item-card:: :octicon:`book;7em`
         :link: user_guide/index.html
         :class-card: frontpage-display-card
         :text-align: center
   
         User Guide
   
      .. grid-item-card:: :octicon:`telescope;7em`
         :link: examples/index.html
         :class-card: frontpage-display-card
         :text-align: center
   
         Examples
   
      .. grid-item-card:: :octicon:`list-unordered;7em`
         :link: api/index.html
         :class-card: frontpage-display-card
         :text-align: center
   
         API Reference
   
      .. grid-item-card:: :octicon:`mark-github;7em`
         :link: https://github.com/gblazquez/bowshockpy
         :class-card: frontpage-display-card
         :text-align: center
   
         GitHub
