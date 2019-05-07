#ifndef MATH_HH
#define MATH_HH

#include "error.hh"
#include "linpack.hh"
#include "d3.hh"

#include <iostream>
#include <vector>
#include <list>

double angle (const D3::Vector&, const D3::Vector&, 
	      const D3::Vector&); // polar angle
double angle (const D3::Vector&, const D3::Vector&, const D3::Vector&, 
	      const D3::Vector&) ; // dihedral angle

// connectivity matrix (symmetric matrix with zero diagonal)
template <typename T> 
class ConMat : private std::vector<T> 
{
  int _dim;

public:
  explicit ConMat(int d)  : _dim(d), std::vector<T>(d*(d-1)/2) {}
    
  T  operator() (int, int) const ;
  T& operator() (int, int) ;

  bool operator== (const ConMat<T>& cm) const 
  { return (const std::vector<T>&)(*this) == (const std::vector<T>&)cm; }
    
  bool operator!= (const ConMat<T>& cm) const 
  { return (const std::vector<T>&)(*this) != (const std::vector<T>&)cm; }
    
  int size () const { return _dim; }

  T sum () const;
  T row_sum (int) const;
};

template <typename T>
T ConMat<T>::sum () const
{
  T res = 0;
  for(int i = 0; i < std::vector<T>::size(); ++i)
    res += (*this)[i];
  return res;
}

template <typename T>
T ConMat<T>::row_sum (int i) const
{
  T res = 0;
  for(int j = 0; j < size(); ++j)
    if(j != i)
      res += (*this)(i, j);
  return res;
}

template <typename T>
T ConMat<T>::operator() (int i, int j) const 
{
  const char funame [] = "ConMat<T>::operator() (int, int) const: ";

#ifdef DEBUG
  if(i < 0 || i >= _dim) {
    std::cout << funame << "first index, " << i << ", is out of range\n";
    throw Error::Range();
  }

  if(j < 0 || j >= _dim) {
    std::cout << funame << "second index, " << j << ", is out of range\n";
    throw Error::Range();
  }

#endif

  if(j > i)
    return (*this)[j*(j-1)/2 + i];
  if(i > j)
    return (*this)[i*(i-1)/2 + j];

  std::cout << funame << "indices should be different\n";
  throw Error::Range();
}

template <typename T>
T& ConMat<T>::operator() (int i, int j) 
{
  const char funame [] = "ConMat<T>::operator() (int, int): ";

#ifdef DEBUG
  if(i < 0 || i >= _dim) {
    std::cout << funame << "first index, " << i << ", is out of range\n";
    throw Error::Range();
  }

  if(j < 0 || j >= _dim) {
    std::cout << funame << "second index, " << j << ", is out of range\n";
    throw Error::Range();
  }

#endif

  if(j > i)
    return (*this)[j*(j-1)/2 + i];
  if(i > j)
    return (*this)[i*(i-1)/2 + j];

  std::cout << funame << "indices should be different\n";
  throw Error::Range();
}

// sort
template <class T>
void my_sort (const std::vector<T>& v, std::vector<int>& perm)
{
    perm.resize(v.size());

    if(!v.size())
	return;

    typedef std::list<std::vector<int> >::iterator git_t;

    std::list<std::vector<int> > sort_list(1, std::vector<int>(v.size()));
    for(int at = 0; at < v.size(); ++at)
      sort_list.front()[at] = at;

    std::vector<int> upper, lower;
    upper.reserve(v.size());
    lower.reserve(v.size());
    
    git_t g = sort_list.begin();
    while(sort_list.size() != v.size()) {
      while(g->size() == 1)
	g++;

      int p0 = (*g)[0];
      upper.clear();
      lower.clear();
      for(int i = 1; i < g->size(); ++i) {
	int pi = (*g)[i]; 
	if(v[p0] < v[pi])
	  upper.push_back(pi);
	else
	  lower.push_back(pi);
      }

      g = sort_list.erase(g);
      if(upper.size())
	g = sort_list.insert(g, upper);
      g = sort_list.insert(g, std::vector<int>(1, p0));
      if(lower.size())
	g = sort_list.insert(g, lower);
    }

    int i = 0;
    for(git_t g = sort_list.begin(); g != sort_list.end(); ++g)
      perm[i++] = (*g)[0];
}

// multidimensional index
class MultiIndex {
  std::vector<int> _len; // dimensions array
  std::vector<int> _val;
  bool _end;

  static const double _ceil;
public:

  MultiIndex (const std::vector<int>&) ;

  int operator[] (int i) const { return _val[i]; }
  int rank () const { return _len.size(); }
  bool end () const { return _end; }

  void operator++ () ;

};

#endif
