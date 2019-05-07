#include <cmath>
#include "math.hh"

// polar angle
double angle (const D3::Vector& a1, const D3::Vector& a2, const D3::Vector& a3)
{
  double dtemp;
  const D3::Vector v1 = a1 - a2;
  const D3::Vector v2 = a3 - a2;
  dtemp = vdot(v1, v2)/std::sqrt(v1.vdot() * v2.vdot());
  if(dtemp < -1.)
	dtemp = -1.;
  if(dtemp > 1.)
	dtemp = 1.;

  dtemp = std::acos(dtemp) * 180. / M_PI;

  //if(dtemp > 175.)
  //  dtemp = 175.;

  return dtemp;
}

// dihedral angle
double angle (const D3::Vector& a1, const D3::Vector& a2, const D3::Vector& a3, 
	      const D3::Vector& a4) 
{
  double dtemp;

  D3::Vector n  = a3 - a2;
  D3::Vector v1 = a1 - a2;
  D3::Vector v2 = a4 - a3;

  v1.orthogonalize(n);
  v2.orthogonalize(n);
  dtemp = std::sqrt(v1.vdot() * v2.vdot());
  if(dtemp < 1.e-8)
    return 0.;

  dtemp = vdot(v1, v2)/dtemp;
  if(dtemp < -1.)
	dtemp = -1.;
  if(dtemp > 1.)
	dtemp = 1.;

  double res = std::acos(dtemp)
      * 180. / M_PI;

  if(vol(n, v1, v2) > 0.)
    return res;
  else
    return 360. - res;
}


const double MultiIndex::_ceil = 1.e9;

MultiIndex::MultiIndex (const std::vector<int>& l)  
  :  _len(l), _val(l.size()), _end(false)
{
  const char funame [] = "MultiIndex::MultiIndex: "; 

  if(!l.size()) {
    std::cerr << funame << "zeroth order\n";
    throw Error::Range();
  }
  for(int i = 0; i < _len.size(); ++i)
    if(_len[i] < 1) {
      std::cerr << funame << i << "-th dimension is out of range\n";
      throw Error::Range();
    }


  double dtemp = 1.;
  for(int i = 0; i < _len.size(); ++i)
    dtemp *= (double)_len[i];

  if(dtemp > _ceil) {
    std::cerr << funame << "the overall dimension is too big\n";
    throw Error::Range();
  }
}

void MultiIndex::operator++ () 
{
  const char funame [] = "MultiIndex::operator++: ";

  if(end()) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }

  int r;
  for(r = 0; r < _len.size(); ++r)
    if(++_val[r] < _len[r])
      break;
    else
      _val[r] = 0;
    
  if(r == _len.size())
    _end = true;
  else
    _end = false;
}
