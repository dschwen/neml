#include "pyhelp.h" // include first to avoid annoying redef warning

#include "state.h"

namespace py = pybind11;

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(statedata, m) {
  py::module::import("neml.objects");

  m.doc() = "State data structure.";

  py::class_<State, NEMLObject, std::shared_ptr<State>>(m, "State")
      ;
}

} // namespace neml
