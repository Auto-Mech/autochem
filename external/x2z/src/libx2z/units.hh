#ifndef UNITS_HH
#define UNITS_HH

#include <string>
#include "error.hh"

/************************** Physical Constants **************************/

struct Phys_const
{
  static const double bohr; // atomic length in angstroms
  static const double amu;  // amu in atomic units (electron mass = 1)
  static const double kelv; // Kelvin in atomic units
  static const double kcal; // kcal/mol in atomic units
  static const double incm; // cm^-1 energy units

  static double bohr2ang (double x) { return x*bohr; }
  static double ang2bohr (double x) { return x/bohr; }

  static double amu2au (double x) { return x*amu; }
  static double au2amu (double x) { return x/amu; }

  static double kelv2hart (double x) { return x*kelv; }
  static double hart2kelv (double x) { return x/kelv; }

  static double kcal2hart (double x) { return x*kcal; }
  static double hart2kcal (double x) { return x/kcal; }

  static double str2fac (const std::string&) ;
};

#endif
