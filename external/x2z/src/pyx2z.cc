#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "libx2z/atom.hh"
#include "libx2z/chem.hh"

namespace py = pybind11;


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
        .def("zmatrix", &MolecStruct::zmatrix)
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
}
