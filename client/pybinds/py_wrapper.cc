#include "py_wrapper.hpp"

PYBIND11_MODULE(pyeonclient, m) {
    py_parameters(m);
    py_matter(m);
    py_log(m);
    // Jobs
    py_job(m);
    py_saddlesearchjob(m);
    // Optimizer
    py_optimizer(m);
    py_lbfgs(m);
    // Objective Functions
    py_objectivefunction(m);
    py_matterobjfunc(m);
    // Potentials
    py_potential(m);
    py_morse(m);
}
