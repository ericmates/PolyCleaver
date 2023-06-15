.. _installation:
.. index:: Installation

Installation
************


Installing via `pip`
====================

:program:`PolyCleaver` can be simply installed via `pip`::

    pip install polycleaver


Cloning the repository
======================

If you want to get the absolutely latest version you can clone the
repo::

    git clone https://github.com/ericmates/PolyCleaver.git

and then install :program:`PolyCleaver` via::

    cd PolyCleaver
    python3 setup.py install --user

in the root directory. This will set up :program:`PolyCleaver` as a Python module.


Requirements
============

:program:`PolyCleaver` requires Python3 and depends on the following libraries

* `pymatgen <https://pymatgen.org>`_ (structure handling)
* `NumPy <http://www.numpy.org/>`_ (numerical linear algebra)
