#include "units.hh"
#include <iostream>

const double Phys_const::bohr = 0.529177249;
const double Phys_const::amu  = 1822.8885;
const double Phys_const::kelv = 3.166829e-6;
const double Phys_const::kcal = 1.59362e-3;
const double Phys_const::incm = 4.55633e-6;

double Phys_const::str2fac(const std::string& unit) throw(Error::General) 
{
  const char funame [] = "Phys_const::str2fac: ";

  if(unit == "Kelvin" || unit == "kelvin" || unit == "kelv")
    return Phys_const::kelv;
  else if(unit == "kcal" || unit == "kcal/mol")
    return Phys_const::kcal;
  else if(unit == "Angstrom" || unit == "angstrom")
    return 1./Phys_const::bohr;
  else if(unit == "incm" || unit == "invcm")
    return incm;
  else if(unit == "amu")
    return amu;
  else if(unit == "Bohr" || unit == "bohr" || 
	  unit == "hartree" || unit == "Hartree" || 
	  unit == "au" || unit == "nu") 
    return 1.;
  else if(unit == "%" || unit == "percent" || unit == "Percent")
    return .01;

  std::cerr << funame << "unknown unit: " << unit << "\n";
  throw Error::Form();
}
