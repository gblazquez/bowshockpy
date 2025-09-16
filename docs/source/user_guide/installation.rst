How to install
==============

You can install ``BowshockPy`` using ``pip``, either from PyPI or from the source repository. The :doc:`dependencies<dependencies>` will be automatically installed.

We recommend to use a clean Python environment in order to install ``BowshockPy`` and its dependencies. The Python environment can be created using ``venv``. For example, in order to create a Python environment called "env", run in your terminal

.. code-block:: console

   $ python3 -m venv env

You should activate the environment in order to install ``BowshockPy`` and its dependencies, and every time you would like to use the package. You can activate the Python environment by executing

.. code-block:: console

   $ source /path/to/new/virtual/environment/bin/activate

After activating your environment, you can proceed with the installation of ``BowshockPy``, either from PyPI or from the source repository.


Installation from PyPI
----------------------

You can get ``BowshockPy`` from PyPI using pip:

.. code-block:: console

   (env)$ pip install bowshockpy 


Installation from source repository
-----------------------------------

You can install ``BowshockPy`` from the source by cloning the `GitHub repository <https://github.com/gblazquez/BowshockPy>`_:

.. code-block:: console

    (env)$ git clone https://github.com/gblazquez/bowshockpy.git 
    (env)$ cd bowshockpy
    (env)$ pip install .


Tests
-----

``BowshockPy`` includes tests that check the proper functioning of the software. In order to run the tests, you should have `pytest <https://docs.pytest.org/en/stable/getting-started.html>`_ installed in your python environment. Then, you can run the tests by executing in the installation folder

.. code-block:: console

   (env)$ pytest tests/

