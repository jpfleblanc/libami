#include "../src/ami_base.hpp"
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>

PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::vector<AmiBase::pole_struct>>>); // for mutability of P_t
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::vector<double>>>); // for mutability of S_t
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::vector<AmiBase::g_struct>>>); // for mutability of R_t

namespace py = pybind11;

// trival structs translating typedefs into python ?????
//struct energy_t{
//    std::vector<int> v;
//};

void init_pyami_wrapper(py::module &m) {
	
  py::bind_vector<std::vector<int>>(m, "VectorInt");
  py::bind_vector<std::vector<double>>(m, "VectorDouble");
  py::bind_vector<std::vector<std::vector<std::vector<AmiBase::pole_struct>>>>(m, "P_t");
  py::bind_vector<std::vector<std::vector<std::vector<double>>>>(m, "S_t");
  py::bind_vector<std::vector<std::vector<std::vector<AmiBase::g_struct>>>>(m, "R_t");

  py::class_<AmiBase> AmiBase(m, "AmiBase");
  AmiBase.def(py::init<>());
  AmiBase.def(py::init<AmiBase::ami_parms &>());
    
  py::class_<AmiBase::ami_vars> (AmiBase, "ami_vars")
    .def(py::init<>())
    .def(py::init<AmiBase::energy_t, AmiBase::frequency_t>())
    .def(py::init<AmiBase::energy_t, AmiBase::frequency_t, double>())
    .def(py::init<AmiBase::energy_t, AmiBase::frequency_t, double, double>())
    .def_readwrite("energy_", &AmiBase::ami_vars::energy_)
    .def_readwrite("frequency_", &AmiBase::ami_vars::frequency_)
    .def_readwrite("prefactor", &AmiBase::ami_vars::prefactor)
    .def_readwrite("BETA_", &AmiBase::ami_vars::BETA_)
    .def_readwrite("gamma_", &AmiBase::ami_vars::gamma_);
  
  py::class_<AmiBase::ami_parms> (AmiBase, "ami_parms")
    .def(py::init<>())
    .def(py::init<int, double>())
    .def(py::init<int, double, AmiBase::graph_type>())
    .def(py::init<int, double, AmiBase::graph_type, AmiBase::int_type, AmiBase::disp_type>())
    .def_readwrite("N_INT_", &AmiBase::ami_parms::N_INT_)
    .def_readwrite("N_EXT_", &AmiBase::ami_parms::N_EXT_)
    .def_readwrite("E_REG_", &AmiBase::ami_parms::E_REG_)
    .def_readwrite("tol_", &AmiBase::ami_parms::tol_)
    .def_readwrite("TYPE_", &AmiBase::ami_parms::TYPE_)
    .def_readwrite("int_type_", &AmiBase::ami_parms::int_type_)
    .def_readwrite("dispersion_", &AmiBase::ami_parms::dispersion_);

  py::class_<AmiBase::g_struct> (AmiBase, "g_struct")
    .def(py::init<AmiBase::epsilon_t, AmiBase::alpha_t, AmiBase::stat_type>())
    .def(py::init<AmiBase::epsilon_t, AmiBase::alpha_t>())
    .def(py::init<>())
    .def_readwrite("eps_", &AmiBase::g_struct::eps_)
    .def_readwrite("alpha_", &AmiBase::g_struct::alpha_)
    .def_readwrite("stat_", &AmiBase::g_struct::stat_)
    .def_readwrite("species_", &AmiBase::g_struct::species_)
    .def_readwrite("eff_stat_", &AmiBase::g_struct::eff_stat_)
    .def_readwrite("pp", &AmiBase::g_struct::pp);

  py::class_<AmiBase::pole_struct> (AmiBase, "pole_struct")
    .def(py::init<>())
    .def(py::init<AmiBase::epsilon_t, AmiBase::alpha_t>())
    .def_readwrite("eps_", &AmiBase::pole_struct::eps_)
    .def_readwrite("alpha_", &AmiBase::pole_struct::alpha_)
    .def_readwrite("index_", &AmiBase::pole_struct::index_)
    .def_readwrite("multiplicity_", &AmiBase::pole_struct::multiplicity_)
    .def_readwrite("der_", &AmiBase::pole_struct::der_)
    .def_readwrite("which_g_", &AmiBase::pole_struct::which_g_)
    .def_readwrite("x_alpha_", &AmiBase::pole_struct::x_alpha_);

// Having issues with passing empty S_t, P_t, R_t and filling them, ie push_back is not working

  AmiBase.def("construct", py::overload_cast<AmiBase::ami_parms &, AmiBase::g_prod_t, AmiBase::R_t &, AmiBase::P_t &, AmiBase::S_t &>(&AmiBase::construct), "Construction function for term-by-term construction.");

  AmiBase.def("evaluate", py::overload_cast<AmiBase::ami_parms &, AmiBase::R_t &, AmiBase::P_t &, AmiBase::S_t &, AmiBase::ami_vars &>(&AmiBase::evaluate), "This is the primary evaluation which takes again `ami_parms`, the outputs from `construct` as well as the `ami_vars` external values that enter into the expression");
}
