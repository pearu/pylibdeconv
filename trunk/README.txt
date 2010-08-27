
====================================================
PyLibDeconv - a Python wrapper of the Deconv library
====================================================

:Author: Pearu Peterson
:Created: August 2010

Introduction
============

Deconv is a deconvolution software package for 3-D quantitative
fluorescence microscopy imaging, it is released with LGPL license. See

  http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2818.2009.03205.x/abstract

for more information.

The home page of the PyLibDeconv project is

  http://code.google.com/p/pylibdeconv/

To import SWIG generated wrapper to Python, use

  >>> import libdeconv

that provides number of Deconv classes such as CCube_float,
CCube_double, CSlice_float, CSlice_double, Fluo3DPSF, FluoRZPSF,
LWdeconvolver, CGdeconvolver, and EMdeconvolver.

Installation
============

Prerequisites
-------------

Before building the wrapper, make sure that you have installed the
following libraries that the Deconv software uses (see also
src/deconv/README):

  fftw3, gsl, blas

The Deconv software is included in PyLibDeconv sources.

You may install the C++ libdeconv to your system according to
libdeconv/src/deconv/README but it will not be usable for Python
because Python requires a shared library.  So, you may skip installing
the libdeconv and continue as follows. Installing the C++ libdeconv
may be still useful because it provides various command line tools.

Installing
----------

Run

  python setup.py install

This will install libdeconv Python package to your system.

Enjoy!
Pearu Peterson
