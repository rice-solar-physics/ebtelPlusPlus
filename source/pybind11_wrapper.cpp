// Python bindings for the ebtel++ code

#include <pybind11/pybind11.h>
#include <ebtel.h>

PYBIND11_MODULE(ebtelplusplus, m) {
    m.doc() = "Python bindings for ebtel++ model";
    m.def("run", &run, "Run the ebtel++ model for a given config file");
}