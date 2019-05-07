#include "atom.hh"
#include "units.hh"

#include <sstream>
#include <cstdlib>
#include <vector>

/************************** Atom description ****************************/

const char* AtomBase::name() const 
{
  const char funame [] = "AtomBase::name: ";

  switch(_num) {
  case HYDROGEN: 
    return "H";
  case CARBON: 
    return "C";
  case NITROGEN:
    return "N";
  case OXYGEN:
    return "O";
  case FLUORINE:
    return "F";
  case SODIUM:
    return "Na";
  case SILICON:
    return "Si";
  case PHOSPHORUS:
    return "P";
  case SULFUR:
    return "S";
  case CHLORINE:
    return "Cl";
  case TITANIUM:
    return "Ti";
  case BROMINE:
    return "Br";
  case DUMMY:
    return "X";
  default:
    std::cerr << funame << "unknown atomic number " << _num << "\n";
    throw Error::Range();
  }
}

void AtomBase::_str2num(const std::string& name) 
{
  const char funame [] = "AtomBase::set(const std::string&): ";

  if (name == "X")
    _num = DUMMY;
  else if (name == "H")
    _num = HYDROGEN;
  else if (name == "C")
    _num = CARBON;
  else if (name == "N")
    _num = NITROGEN;
  else if (name == "O")
    _num = OXYGEN;
  else if (name == "F")
    _num = FLUORINE;
  else if (name == "Na")
    _num = SODIUM;
  else if (name == "Si")
    _num = SILICON;
  else if (name == "P")
    _num = PHOSPHORUS;
  else if (name == "S")
    _num = SULFUR;
  else if (name == "Cl")
    _num = CHLORINE;
  else if (name == "Ti")
    _num = TITANIUM;
  else if (name == "Br")
    _num = BROMINE;
  else {
    std::cerr << funame << "unknown atom name "  << name << "\n";
    throw Error::Range();
  }
  _isot = _default_isot();
}

int AtomBase::_default_isot () const 
{
  const char funame [] = "AtomBase::_default_isot: ";

  int res;
  switch (_num) {
  case DUMMY:
    return 0;
  case HYDROGEN:
    return 1;
  case CARBON:
    return 12;
  case NITROGEN:
    return 14;
  case OXYGEN:
    return 16;
  case FLUORINE:
    return 19;
  case SODIUM:
    return 23;
  case SILICON:
    return 28;
  case PHOSPHORUS:
    return 31;
  case SULFUR:
    return 32;
  case CHLORINE:
    return 35;
  case TITANIUM:
    return 49;
  case BROMINE:
    return 79;
  default:
    std::cerr << funame << "unknown atomic number " << _num << "\n";
    throw Error::Range();
  }
}

double AtomBase::mass () const 
{
  const char funame [] = "AtomBase::mass: ";

  switch (_num) {
  case DUMMY:
    return 0.;
  case HYDROGEN:
    switch (_isot) {
    case 1: return Phys_const::amu * 1.007825;
    case 2: return Phys_const::amu * 2.014;
    case 3: return Phys_const::amu * 3.01605;
    default: 
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case CARBON:
    switch (_isot) {
    case 12: return Phys_const::amu * 12.0;
    case 13: return Phys_const::amu * 13.00335;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case NITROGEN:
    switch (_isot) {
    case 14: return Phys_const::amu * 14.00307;
    case 15: return Phys_const::amu * 15.00011;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case OXYGEN:
    switch (_isot) {
    case 16: return Phys_const::amu * 15.99491;
    case 17: return Phys_const::amu * 17.0;
    case 18: return Phys_const::amu * 18.0;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case FLUORINE:
    switch(_isot) {
    case 19: return Phys_const::amu * 18.9984;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case SODIUM:
    switch(_isot) {
    case 23: return Phys_const::amu * 22.9898;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case SILICON:
    switch(_isot) {
    case 28: return Phys_const::amu * 27.97693;
    case 29: return Phys_const::amu * 28.97649;
    case 30: return Phys_const::amu * 29.97376;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case PHOSPHORUS:
    switch(_isot) {
    case 31: return Phys_const::amu * 30.97376;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case SULFUR:
    switch(_isot) {
    case 32: return Phys_const::amu * 31.97207;
    case 33: return Phys_const::amu * 32.97146;
    case 34: return Phys_const::amu * 33.96786;
    case 36: return Phys_const::amu * 35.96709;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case CHLORINE:
    switch(_isot) {
    case 35: return Phys_const::amu * 34.96885;
    case 37: return Phys_const::amu * 37.;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case TITANIUM:
    switch(_isot) {
    case 46: return Phys_const::amu * 45.95263;
    case 47: return Phys_const::amu * 46.951764;
    case 48: return Phys_const::amu * 47.947947;
    case 49: return Phys_const::amu * 48.947871;
    case 50: return Phys_const::amu * 49.944792;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  case BROMINE:
    switch(_isot) {
    case 79: return Phys_const::amu * 78.9183;
    case 81: return Phys_const::amu * 80.9163;
    default:
      std::cerr << funame << "unknown isotope: " << _isot << "\n";
      throw Error::Range();
    }
  default:
    std::cerr << funame << "unknown atomic number " << _num << "\n";
    throw Error::Range();
  }
}

unsigned AtomBase::valence () const 
{
  const char funame [] = "AtomBase::valence: ";

  switch (_num) {
  case DUMMY:
    return 0;
  case HYDROGEN:
    return 1;
  case CARBON:
    return 4;
  case NITROGEN:
    return 3;
  case OXYGEN:
    return 2;
  case FLUORINE:
    return 1;
  case SILICON:
    return 4;
  case PHOSPHORUS:
    return 3;
  case SULFUR:
    return 2;
  case CHLORINE:
    return 1;
  case BROMINE:
    return 1;
  default:
    std::cerr << funame << "unknown atomic number " << _num << "\n";
    throw Error::Range();
  }
}

void Atom::_read (std::istream& from) 
{
    const char funame [] = "Atom::_read: ";

    int itemp;

    std::string line;
    std::getline(from, line);
    if(!from) {
	std::cerr << funame << "cannot read the line\n";
	throw Error::Form();
    }

    std::istringstream iss(line);

    std::string s;
    std::vector<std::string> vec;
    while(iss >> s)
      vec.push_back(s);

    switch(vec.size()) {
    case 1: // atom name only
      set(vec[0]);
      break;

    case 2: // atom name and isotope number
      itemp = std::atoi(vec[1].c_str());
      if(itemp <= 0) {
	std::cerr << funame << "wrong isotope on the line: " << line << "\n";
	throw Error::Form();
      }
      set(vec[0], itemp);
      break;

    case 4: // atom name and coordinates
      set(vec[0]);
      for(int i = 0; i < 3; ++i)
	(*this)[i] = std::atof(vec[i + 1].c_str());
      break;

    case 5: // atom name, isotope number, and coordinates
      itemp = std::atoi(vec[1].c_str());
      if(itemp <= 0) {
	std::cerr << funame << "wrong isotope on the line: " << line << "\n";
	throw Error::Form();
      }
      set(vec[0], itemp);
      for(int i = 0; i < 3; ++i)
	(*this)[i] = std::atof(vec[i + 2].c_str());
      break;

    default:
      std::cerr << funame << "wrong number of items, " << vec.size() << ", on the line: " << line << "\n";
      throw Error::Form();
    }
}




