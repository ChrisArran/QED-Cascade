#include "RunManager.hh"
#include "ThreeVector.hh"
#include <pybind11/pybind11.h>

namespace py = pybind11;


PYBIND11_MODULE(QEDCascPy, module)
{
    py::class_<RunManager>(module, "RunManager")
        .def(py::init<>())
        .def("setTime", &RunManager::setTime, "Set the end time and time-step")
        .def("setField", &RunManager::setField, "Set the field properties")
        .def("setGenerator", &RunManager::setGenerator, "Set the particle source")
        .def("setPhysics", &RunManager::setPhysics, "Select physics model")
        .def("setSampleFraction", &RunManager::setSampleFraction,
            "Set sampling fraction for continuos radiation")
        .def("usePairProduction", &RunManager::usePairProduction,
            "Turn pair production on.",
             py::arg("useBW"), py::arg("up_scaling") = 1.0)
        .def("useFiniteBeaming", &RunManager::useFiniteBeaming,
            "Turn finite beaming on.")
        .def("beamOn", &RunManager::beamOn, "Run the simulation", 
            py::arg("events"), py::arg("threads") = 1)
        .def("getInput", &RunManager::getInput,
            "Get particle properties before interaction")
        .def("getElectrons", &RunManager::getElectrons,
            "Get electron properties after interaction")
        .def("getPositrons", &RunManager::getPositrons,
            "Get positron properties after interaction")
        .def("getPhotons", &RunManager::getPhotons,
            "Get photon properties after interaction");

    py::class_<ThreeVector>(module, "ThreeVector")
        .def(py::init<double, double, double>());
}
