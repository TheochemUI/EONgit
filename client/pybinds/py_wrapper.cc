#include "py_wrapper.hpp"

PYBIND11_MODULE(pyeonclient, m) {
    py_parameters(m);
    py_potential(m);
    py_matter(m);
}