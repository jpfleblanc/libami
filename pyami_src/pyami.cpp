#include <pybind11/pybind11.h>
#define PYBIND11_CPP14

namespace py = pybind11;

void init_pyami_wrapper(py::module &);

namespace pai{

PYBIND11_MODULE(pyami, m) {
    // Optional docstring
    m.doc() = "Pyami: a python wrapper for the libami c++ library";
    init_pyami_wrapper(m);
}
}

