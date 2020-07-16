#include "pyhelp.h" // include first to avoid annoying redef warning

#include "larsonmiller.h"

#include "nemlerror.h"

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>)

namespace neml {

PYBIND11_MODULE(larsonmiller, m) {
  py::module::import("neml.objects");
  py::module::import("neml.solvers");

  m.doc() = "Class defining a Larson-Miller relation between stress and rupture.";

  py::class_<LMTrialState, TrialState>(m, "LMTrialState")
      ;

  py::class_<LarsonMillerRelation, NEMLObject, Solvable, std::shared_ptr<LarsonMillerRelation>>(m, "LarsonMillerRelation")
      .def(py::init([](py::args args, py::kwargs kwargs)
                    {
                      return create_object_python<LarsonMillerRelation>(args, kwargs, {"function", "C"});
                    }))
      .def("sR",
           [](LarsonMillerRelation & m, double t, double T) -> double
           {
            double s;
            int ier = m.sR(t, T, s);
            py_error(ier);
            return s;
           }, "The rupture stress as a function of time and temperature.")
      .def("tR",
           [](LarsonMillerRelation & m, double s, double T) -> double
           {
            double t;
            int ier = m.tR(s, T, t);
            py_error(ier);
            return t;
           }, "The rupture time as a function of stress and temperature.")
      .def("dtR_ds",
           [](LarsonMillerRelation & m, double s, double T) -> double
           {
            double dt;
            int ier = m.dtR_ds(s, T, dt);
            py_error(ier);
            return dt;
           }, "The derivative of the time to rupture with respect to stress.")
    ;
}

} // namespace neml
