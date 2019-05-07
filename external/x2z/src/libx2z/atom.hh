#ifndef ATOM_HH
#define ATOM_HH

#include "d3.hh"
#include <iomanip>

/*************************** Atom description ***************************/

class AtomBase 
{
public:
  enum AN { // atomic number
    DUMMY = 0,
    HYDROGEN = 1,
    CARBON = 6,
    NITROGEN = 7,
    OXYGEN = 8,
    FLUORINE = 9,
    SODIUM = 11,
    SILICON = 14,
    PHOSPHORUS = 15,
    SULFUR = 16,
    CHLORINE = 17,
    TITANIUM = 22,
    BROMINE = 35
  };

private:
  AN  _num;   // atomic number 
  int _isot; // isotope

  int  _default_isot   () const ;
  
  void _str2num (const std::string&) ;
  
public: 
  void set (AN n) { _num = n; _isot = _default_isot(); }
  void set (AN n, int i) { _num = n; _isot = i; }

  void set (const std::string& s)         { _str2num(s); _isot = _default_isot(); }
  void set (const std::string& s, int i)  { _str2num(s); _isot = i; }

  AtomBase () : _num(DUMMY), _isot(0)  {}
  explicit AtomBase (const std::string& s)  { set(s); }
  explicit AtomBase (const std::string& s, int i)  { set(s, i); }
    
  const char* name    () const ;
  double      mass    () const ;
  unsigned    valence () const ;

  operator AN      () const { return _num; }
  AN        number () const { return _num; }
  int      isotope () const { return _isot; }

  bool operator== (const AtomBase& a) const { return _num == a._num && _isot == a._isot; }
  bool operator!= (const AtomBase& a) const { return _num != a._num || _isot != a._isot; }
  bool operator<  (const AtomBase& a) const;
  bool operator>  (const AtomBase& a) const;
};

inline bool AtomBase::operator< (const AtomBase& a) const 
{ 
  if(_num < a._num)
    return true;
  else if(_num == a._num && _isot < a._isot)
    return true;
  return false;
}

inline bool AtomBase::operator> (const AtomBase& a) const 
{ 
  if(_num > a._num)
    return true;
  else if(_num == a._num && _isot > a._isot)
    return true;
  return false;
}

class Atom : public AtomBase, public D3::Vector
{
  void _read(std::istream&) ;

public:
  Atom () {}
  explicit Atom (const std::string& s) : AtomBase(s) {}
  explicit Atom (std::istream& from)  { _read(from); }
  Atom (const std::string& s, int i) : AtomBase(s, i) {}

  friend std::istream& operator>> (std::istream&, Atom&) ;
};

inline std::istream& operator>> (std::istream& from, Atom& a)  
{ 
  a._read(from); 
  return from; 
}

inline std::ostream& operator<< (std::ostream& out , const Atom& a)
{
    out << std::setw(2) << a.name();
    for(int i = 0; i < 3; ++i)
	out << std::setw(15) << a[i];
    return out;
}


#endif
