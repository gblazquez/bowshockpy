How to install
==============

You can install ``bowshockpy`` using ``pip``, either from PyPI or from the source repository. The :doc:`dependencies<dependencies>` will be automatically installed.


Installation from PyPI
----------------------

You can get ``bowshockpy`` from PyPI using pip:

.. code-block:: console

   $ pip install bowshockpy 


Installation from source repository
-----------------------------------

You can install ``bowshockpy`` from the source by cloning the `GitHub repository <https://github.com/gblazquez/bowshockpy>`_:

.. code-block:: console

    $ git clone https://github.com/gblazquez/bowshockpy.git 
    $ cd bowshockpy
    $ pip install .


Tests
-----

``bowshockpy`` includes tests that check the proper functioning of the software. In order to run the tests, you should have `pytest <https://docs.pytest.org/en/stable/getting-started.html>`_ installed in your python environment. Then, you can run the tests by executing in the installation folder

.. code-block:: console

   $ pytest tests/

