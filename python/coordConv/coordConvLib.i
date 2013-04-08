%define coordConvLibDocString
"Python interface to coordConv"
%enddef

%feature("autodoc", "1");
%module(package="coordConv", docstring=coordConvLibDocString) coordConvLib

%{
#include <stdexcept>
#include <new>
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
#include "coordConv/coordConv.h"
%}

%init %{
    import_array();
%}

%include "std_except.i"
%include "std_string.i"
%include "typemaps.i"
%include "boost_shared_ptr.i"
%include "ndarray.i"

// Specify the default C++ to python exception handling interface
%exception {
    try {
        $action
    } catch (std::invalid_argument & e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        SWIG_fail;
    } catch (std::out_of_range & e) {
        PyErr_SetString(PyExc_LookupError, e.what());
        SWIG_fail;
    } catch (std::logic_error & e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        SWIG_fail;
    } catch (std::range_error & e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        SWIG_fail;
    } catch (std::overflow_error & e) {
        PyErr_SetString(PyExc_OverflowError, e.what());
        SWIG_fail;
    } catch (std::runtime_error & e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        SWIG_fail;
    } catch (std::bad_alloc & e) {
        PyErr_SetString(PyExc_MemoryError, e.what());
        SWIG_fail;
    } catch (std::exception & e) {
        PyErr_SetString(PyExc_StandardError, e.what());
        SWIG_fail;
    } catch (...) {
        SWIG_fail;
    }
}

%declareNumPyConverters(Eigen::Vector2d);
%declareNumPyConverters(Eigen::Vector3d);
%declareNumPyConverters(Eigen::Matrix3d);

%include "coordConv/pvt.h"
%include "coordConv/physConst.h"

%apply double &OUTPUT { double & };

%include "coordConv/mathUtils.h"
%include "coordConv/angSideAng.h"
%include "coordConv/rotEqPol.h"
%include "coordConv/rotXY.h"
%include "coordConv/site.h"
%include "coordConv/time.h"
%include "coordConv/coord.h"
%include "coordConv/pvtCoord.h"
%include "coordSys.i"
