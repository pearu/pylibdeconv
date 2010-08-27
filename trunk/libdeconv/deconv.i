/* -*- c++ -*- */
/*
  SWIG wrapper of Deconv library for Python
  Author: Pearu Peterson
  Created: August 2010
 */

%module deconv
 %{
#define SWIG_FILE_WITH_INIT
#include "CCube.h"
#include "CSlice.h"
#include "Fluo3DPSF.h"
#include "FluoRZPSF.h"
#include "LWdeconvolver.h"
#include "CGdeconvolver.h"
#include "EMdeconvolver.h"
 %}

%include "numpy.i"

%init %{
import_array();
%}

%feature("autodoc", "1");

%exception {
   try {
      $action
   } 
   catch (std::exception &e)
     {
       PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(e.what()));
       return NULL;
     }
}

%include "CCube.h"
%include "CSlice.h"
%include "Fluo3DPSF.h"
%include "FluoRZPSF.h"
%include "LWdeconvolver.h"
%include "CGdeconvolver.h"
%include "EMdeconvolver.h"

#define GETPIXELMTH(T) \
  T getpixel(int x, int y) \
  {\
    T value; \
    self->getpixel(x,y,value); \
    return value; \
  }


#define GETVOXELMTH(T) \
  T getvoxel(int x, int y, int z) \
  {\
    T value; \
    self->getvoxel(x,y,z,value); \
    return value; \
  }

#define SLICEASARRAYMTH(T,NPY_T)			\
  PyObject* asarray(void)\
  {\
    PyObject *arr = NULL;\
    npy_intp dims[2]; \
    dims[1] = self->width();\
    dims[0] = self->height();\
    arr = PyArray_SimpleNew(2, dims, NPY_T);\
    memcpy(PyArray_DATA(arr), self->data(), self->size()*sizeof(T));\
    return arr;\
  }

#define CUBEASARRAYMTH(T,NPY_T)			\
  PyObject* asarray(void)\
  {\
    PyObject *arr = NULL;\
    npy_intp dims[3]; \
    dims[2] = self->length();\
    dims[1] = self->width();\
    dims[0] = self->height();\
    arr = PyArray_SimpleNew(3, dims, NPY_T);\
    memcpy(PyArray_DATA(arr), self->data(), self->size()*sizeof(T));\
    return arr;\
  }

#define SLICEFROMARRAYMTH(T,NPY_T)			\
  PyObject* fromarray(PyObject* obj)\
  {\
    int width = 1, height = 1;\
    PyObject* arr = PyArray_ContiguousFromAny(obj, NPY_T, 1, 2);\
    switch (PyArray_NDIM(arr))\
      {\
      case 2:\
	height = PyArray_DIM(arr, 0);\
	width = PyArray_DIM(arr, 1);\
	break;\
      case 1:\
	height = PyArray_DIM(arr, 0);\
	break;\
      default:\
	PyErr_SetString(PyExc_ValueError, "array argument must have rank 1,2 but got more or 0");\
	return NULL;\
      }\
    self->free();\
    self->init(width, height, (T*)PyArray_DATA(arr));\
    Py_INCREF(Py_None);\
    return Py_None;\
  }

#define CUBEFROMARRAYMTH(T,NPY_T)			\
  PyObject* fromarray(PyObject* obj)\
  {\
    int length = 1, width = 1, height = 1;\
    PyObject* arr = PyArray_ContiguousFromAny(obj, NPY_T, 1, 3);\
    switch (PyArray_NDIM(arr))\
      {\
      case 3:\
	height = PyArray_DIM(arr, 0);\
	width = PyArray_DIM(arr, 1);\
	length = PyArray_DIM(arr, 2);\
	break;\
      case 2:\
	height = PyArray_DIM(arr, 0);\
	width = PyArray_DIM(arr, 1);\
	break;\
      case 1:\
	height = PyArray_DIM(arr, 0);\
	break;\
      default:\
	PyErr_SetString(PyExc_ValueError, "array argument must have rank 1,2,3 but got more or 0");\
	return NULL;\
      }\
    self->free();\
    self->init(length, width, height, (T*)PyArray_DATA(arr));\
    Py_INCREF(Py_None);\
    return Py_None;\
  }

%template(CCube_float) CCube<float>;
%template(CCube_double) CCube<double>;
%template(CSlice_float) CSlice<float>;
%template(CSlice_double) CSlice<double>;

%extend CCube<float>
{
  GETVOXELMTH(float);
  CUBEASARRAYMTH(float, NPY_FLOAT32);
  CUBEFROMARRAYMTH(float, NPY_FLOAT32);
};

%extend CCube<double>
{
  GETVOXELMTH(double);
  CUBEASARRAYMTH(double, NPY_FLOAT64);
  CUBEFROMARRAYMTH(double, NPY_FLOAT64);
};

%extend CSlice<float>
{
  GETPIXELMTH(float);
  SLICEASARRAYMTH(float, NPY_FLOAT32);
  SLICEFROMARRAYMTH(float, NPY_FLOAT32);
};

%extend CSlice<double>
{
  GETPIXELMTH(double);
  SLICEASARRAYMTH(double, NPY_FLOAT64);
  SLICEFROMARRAYMTH(double, NPY_FLOAT64);
};

%extend Fluo3DPSF
{
  PyObject* create(int Nx, int Ny, int Nz, bool Check = true)
  {
    npy_intp dims[] = {Nz, Ny, Nx};
    PyObject *arr = PyArray_EMPTY(3, dims, NPY_FLOAT64, 0);
    self->create(Nx, Ny, Nz, (double*)PyArray_DATA(arr), Check);
    return arr;
  }

  CCube<double>* create_CCube_double(int Nx, int Ny, int Nz, bool Check = true)
  {
    CCube<double>* psf = new CCube<double>(Nx, Ny, Nz); 
    self->create(Nx, Ny, Nz, psf->data(), Check);
    return psf;
  }

  CCube<float>* create_CCube_float(int Nx, int Ny, int Nz, bool Check = true)
  {
    CCube<float>* psf = new CCube<float>(Nx, Ny, Nz); 
    self->create(Nx, Ny, Nz, psf->data(), Check);
    return psf;
  }
};

%extend FluoRZPSF
{
  CCube<double>* get3Dpsf(int Nx, int Ny, int Nz,
			  double dx, double dy, double dz, int Points2Sum = 5)
  {
    CCube<double>* psf = new CCube<double>(Nx, Ny, Nz);
    self->get3Dpsf(Nx, Ny, Nz, dx, dy, dz,
		   psf->data(), NULL, Points2Sum);
    return psf;
  }

  CCube<float>* get3Dpsf(int Nx, int Ny, int Nz,
			 float dx, float dy, float dz, int Points2Sum = 5)
  {
    CCube<float>* psf = new CCube<float>(Nx, Ny, Nz);
    self->get3Dpsf(Nx, Ny, Nz, dx, dy, dz,
		   psf->data(), NULL, Points2Sum);
    return psf;
  }
}

#define RUNMTH(T, WS_T)\
  CCube<T>* run(CCube<T>& image, CCube<T>& psf)\
  {\
    CCube<T>* object = new CCube<T>(image);		\
    WS_T ws;\
    self->run(image.length(), image.width(), image.height(), image.data(), psf.data(), object->data(), ws);\
    return object;\
  }

%extend LWdeconvolver
{
  RUNMTH(float, LWsws);
  RUNMTH(double, LWdws);
}

%extend CGdeconvolver
{
  RUNMTH(float, CGsws);
  RUNMTH(double, CGdws);
}

%extend EMdeconvolver
{
  RUNMTH(float, EMsws);
  RUNMTH(double, EMdws);
}
