#include<vector>
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>

#include "libx2z/units.hh"
#include "libx2z/chem.hh"
#include "libx2z/math.hh"

int main(int argc, const char* argv [])
{
  const char funame [] = "x2z: ";

  const std::string geom_key = "Geometry";

  if(argc != 2) {
    std::cout << funame << "usage: x2z input_file\n";
    return 1;
  }
  
  std::ifstream from(argv[1]);
  if(!from) {
    std::cout << funame << "cannot open " << argv[1] << " file\n";
    return 1;
  }

  int itemp;
  double dtemp;
  
  std::string token, line, s;

  MolecGeom start_geom;

  std::getline(from, line);
  std::istringstream iss(line);
  int n;
  if(!(iss >> n)) {
    std::cout << funame << "cannot read number of atoms\n";
    return 1;
  }
  
  std::getline(from, line);
  
  Atom a;
  for(int i = 0; i < n; ++i) {
    from >> a;
    a /= Phys_const::bohr;
    start_geom.push_back(a);
  }
  if(!from) {
    std::cout << funame << " cannot read " << token << " section\n";
    return 1;
  }

  std::string amax_key = "AngleTolerance";


  while(from >> token) {
    if(token == amax_key) {
      //
      if(!(from >> dtemp)) {
	std::cerr << funame << token << ": corrupted\n";
	return 1;
      }
      
      if(dtemp <= 0.) {
	std::cerr << funame << token << ": out of range\n";
	return 1;
      }

      angle_tolerance = dtemp;
    }
    else {
      std::cerr << funame << "unknown token: " << token << "\n";
      return 1;
    }
  }
  
  MolecOrient mo(start_geom);

  std::cout << "molecule is ";
  if(mo.is_linear())
    std::cout << "linear\n";
  else if(mo.is_plane())
    std::cout << "plane\n";
  else {
    std::cout << "nonlinear\n";
    std::cout << "has enantiomer? ";
    if(mo.is_enantiomer())
      std::cout << "yes\n";
    else
      std::cout << "no\n";
  }
  std::cout << "\n";

  std::cout << "rotational symmetry number = " << mo.sym_num() << "\n\n";

  PrimStruct prim(start_geom);

  if(!prim.is_connected()) {
    std::cout << funame << "primary structure is not connected\n";
    return 0;
  }

  MolecStruct mol(prim);

  mol.print(std::cout);

  std::cout << "Z-Matrix:\n";
  //
  std::cout << mol.zmatrix();
  //
  std::cout << "\n";

  std::cout << std::left;

  for(int i = 1; i < mol.zmat_coval().size(1); ++i) {
    //
    std::cout << MolecStruct::var_name(MolecStruct::DISTANCE) << std::setw(2) << i
	      << " = " << std::setw(9) << mol.zmat_coval()(MolecStruct::DISTANCE, i);
    
    if(i > 1) {
      //
      std::cout << " " << MolecStruct::var_name(MolecStruct::POLAR) <<  std::setw(2) << i
		<< " = " << std::setw(9) << mol.zmat_coval()(MolecStruct::POLAR, i);
    }

    if(i > 2) {
      //
      std::cout << " " <<  MolecStruct::var_name(MolecStruct::DIHEDRAL) << std::setw(2) << i
		<< " = " << std::setw(15) << mol.zmat_coval()(MolecStruct::DIHEDRAL, i);
    }

    std::cout << "\n";
  }
  std::cout << "\n";

  // constant parameters
  //
  if(mol.const_var().size()) {
    //
    std::cout << "Constants:";

    for(std::list<int>::const_iterator cit = mol.const_var().begin(); cit != mol.const_var().end(); ++cit) {
      //
      std::cout << "   " << MolecStruct::var_name(*cit % 3) << *cit / 3;
    }

    std::cout << "\n\n";
  }
  else
    std::cout << "no constants\n";

  // rotational dihedral angles
  //
  if(mol.rotation_bond().size()) {
    //
    std::cout << "Rotational bond dihedral angles: ";

    for(std::map<int, std::list<std::list<int> > >::const_iterator bit = mol.rotation_bond().begin();
	bit != mol.rotation_bond().end(); ++bit) {
      //
      if(bit != mol.rotation_bond().begin())
	//
	std::cout << ", ";

      std::cout << MolecStruct::var_name(MolecStruct::DIHEDRAL) << bit->first;
    }

    std::cout << "\n\n";

    std::cout << "Rotational groups:\n";

    for(std::map<int, std::list<std::list<int> > >::const_iterator bit = mol.rotation_bond().begin();
	bit != mol.rotation_bond().end(); ++bit) {
      
      std::cout << MolecStruct::var_name(MolecStruct::DIHEDRAL) << std::setw(2) << bit->first;

      for(std::list<std::list<int> >::const_iterator git = bit->second.begin(); git != bit->second.end(); ++git)
	//
	std::cout << " " << std::setw(10) << mol.group_stoicheometry(*git);

      std::cout << "\n";
    }

    std::cout << "\n";

    for(std::map<int, std::list<std::list<int> > >::const_iterator bit = mol.rotation_bond().begin();
	bit != mol.rotation_bond().end(); ++bit) {
      
      std::cout << MolecStruct::var_name(MolecStruct::DIHEDRAL) << std::setw(2) << bit->first;

      for(std::list<std::list<int> >::const_iterator git = bit->second.begin(); git != bit->second.end(); ++git) {
	//
	std::cout << "   ";
	
	for(std::list<int>::const_iterator it = git->begin(); it != git->end(); ++it) {
	  //
	  if(it != git->begin())
	    //
	    std::cout << ",";
	  
	  std::cout << mol[*it].name() << mol.atom_map(*it) + 1;
	}
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
  else
    std::cout << "no rotational bonds\n";
  
  // beta bonds
  //
  if(mol.beta_bond().size()) {
    //
    std::cout << "Beta-scission bonds: ";

    for(int i = 0; i < mol.beta_bond().size(); ++i) {
      //
      if(i)
	//
	std::cout << ", ";

      std::cout << MolecStruct::var_name(MolecStruct::DISTANCE) << mol.beta_bond()[i];
    }

    std::cout << "\n\n";
  }
  else
    std::cout << "no beta-scision bonds\n";
  
  return 0;

}
