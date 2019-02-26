#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "libx2z/atom.hh"
#include "libx2z/chem.hh"
#include <sstream>

namespace py = pybind11;

std::string zmatrix_string(const MolecStruct& mol) {
  std::ostringstream s;

  s << mol.zmatrix();
  //
  s << "\n";

  s << std::left;

  for(int i = 1; i < mol.zmat_coval().size(1); ++i) {
    //
    s << MolecStruct::var_name(MolecStruct::DISTANCE) << std::setw(2) << i
	      << " = " << std::setw(9) << mol.zmat_coval()(MolecStruct::DISTANCE, i);
    
    if(i > 1) {
      //
      s << " " << MolecStruct::var_name(MolecStruct::POLAR) <<  std::setw(2) << i
		<< " = " << std::setw(9) << mol.zmat_coval()(MolecStruct::POLAR, i);
    }

    if(i > 2) {
      //
      s << " " <<  MolecStruct::var_name(MolecStruct::DIHEDRAL) << std::setw(2) << i
		<< " = " << std::setw(15) << mol.zmat_coval()(MolecStruct::DIHEDRAL, i);
    }

    s << "\n";
  }

  return s.str();
}


PYBIND11_MODULE(pyx2z, module) {
    py::class_<AtomBase>(module, "AtomBase")
        .def(py::init<const std::string&>())
        .def(py::init<const std::string&, int>())
        .def("mass", &AtomBase::mass);
    py::class_<Atom>(module, "Atom")
        .def(py::init<const std::string&>())
        .def(py::init<const std::string&, int>())
        .def("__getitem__", [](const Atom& a, size_t i) {
            if (i >= 3) throw py::index_error();
            return a[i];
        })
        .def("__setitem__", [](Atom& a, size_t i, float v) {
            if (i >= 3) throw py::index_error();
            a[i] = v;
        });
    py::class_<MolecGeom>(module, "MolecGeom")
        .def(py::init<>())
        .def("size", [](MolecGeom& m) { return m.size(); })
        .def("push_back", [](MolecGeom& m, const Atom& a) { m.push_back(a); });
    py::class_<MolecOrient>(module, "MolecOrient")
        .def(py::init<const MolecGeom&>())
        .def("sym_num", &MolecOrient::sym_num)
        .def("is_enantiomer", &MolecOrient::is_enantiomer)
        .def("is_plane", &MolecOrient::is_plane)
        .def("is_linear", &MolecOrient::is_linear)
        .def("size", &MolecOrient::size);
    py::class_<PrimStruct>(module, "PrimStruct")
        .def(py::init<const MolecGeom&>())
        .def("connected_group", &PrimStruct::connected_group)
        .def("group_stoicheometry", &PrimStruct::group_stoicheometry)
        .def("is_connected",
             (bool (PrimStruct::*)(int, int) const)
             &PrimStruct::is_connected);
    py::class_<MolecStruct>(module, "MolecStruct")
        .def(py::init<const PrimStruct&>())
        .def("size", [](MolecStruct& m) { return m.size(); })
        .def("atom_ordering", &MolecStruct::atom_ordering)
        .def("rotation_bond", &MolecStruct::rotation_bond)
        .def("resonance_averaged_bond_order",
             (double (MolecStruct::*)(int, int) const)
             &MolecStruct::bond_order)
        .def("bond_order",
             (unsigned (MolecStruct::*)(int, int, int) const)
             &MolecStruct::bond_order)
        .def("resonance_count", &MolecStruct::resonance_count)
        .def("is_radical", &MolecStruct::is_radical);
    module.def("zmatrix_string", &zmatrix_string);
}
