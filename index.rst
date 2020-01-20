.. scatt documentation master file, created by
   sphinx-quickstart on Fri Jul 12 11:40:30 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to scatt's documentation!
*********************************

.. toctree::
   :maxdepth: 4
   :caption: Contents:

Install scatt
=============

Use the yaml conda environment file to create a dedicated environment.

.. code-block:: bash

   conda env create -f conda-scatt.yml


Run scatt
=========

.. automodule:: scatt
    :members: main

Majatools PACKAGE
=================

Classes
-------

Aoi
+++
.. autoclass:: majatools.majatools.Aoi
    :members:

Context
+++++++
.. autoclass:: majatools.majatools.Context
    :members:

Image
+++++
.. autoclass:: majatools.majatools.Image
    :members:

Run
+++
.. autoclass:: majatools.majatools.Run
    :members:

Timeseries
++++++++++
.. autoclass:: majatools.majatools.Timeseries
    :members:

Modules
-------

diffmap
+++++++
.. automodule:: majatools.majatools.diffmap
    :members:

get_geodata
+++++++++++
.. automodule:: majatools.majatools.get_geodata
    :members:

scatterplot
+++++++++++
.. automodule:: majatools.majatools.scatterplot
    :members:

single_scatterplot
++++++++++++++++++
.. automodule:: majatools.majatools.single_scatterplot
    :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

