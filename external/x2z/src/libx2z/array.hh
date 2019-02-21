#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cstdarg>

#include "error.hh"
#include "linpack.hh"

/**************************************************************************
 ******************************* Array ************************************
 **************************************************************************/

template <typename T>
class Array
{// class Array

  int _capacity;
  int _size;
  T* _begin;
  T* _end;

public:
      
  explicit Array (int =0, const T* =0) throw(Error::General);
  Array (int, const T&) throw(Error::General);
  explicit Array(const Slice<T>&);
  explicit Array(const ConstSlice<T>&);
  template <typename S>
  explicit Array(const std::vector<S>&);

  Array (const Array&);
  ~Array () { if(_begin) delete[] _begin; }
    
  T* begin () { return _begin; }
  T* end   () { return _end; }
  const T* begin () const { return _begin; }
  const T* end   () const { return _end; }

  operator       T* ()       { return _begin; }
  operator const T* () const { return _begin; }

  T& front () throw(Error::General);
  T& back  () throw(Error::General);
  const T& front () const throw(Error::General);
  const T& back  () const throw(Error::General);

  int size     () const { return _size; }
  int capacity () const { return _capacity; }
  void resize  (int) throw(Error::General);
  void reserve (int) throw(Error::General);
  void compact ();

  Array& operator= (const Array&);
  Array& operator= (const T*);
  Array& operator=  (const T&);
  template <typename S>
  Array& operator= (const std::vector<S>&);

  T& operator[] (int) throw(Error::General);
  const T&  operator[] (int) const throw(Error::General);

  Array& operator- ();

  Array& operator+= (const Array&) throw(Error::General);
  Array& operator-= (const Array&) throw(Error::General);
  Array& operator+= (const T*);
  Array& operator-= (const T*);

  Array& operator*= (const T&); 
  Array& operator/= (const T&); 
  Array& operator+= (const T&); 
  Array& operator-= (const T&); 
};// class Array

template <typename T>
Array<T>::Array (int s, const T* p1) throw(Error::General)
  : _size(s), _capacity(s)
{
  const char funame [] = "Array<T>::Array(int, T*): ";

  if(_size < 0) {
    std::cerr << funame << "negative size\n";
    throw Error::Range();
  }

  if(!_size) {
    _end = _begin = 0;
    return;
  }

  _begin = new T[_size];
  _end = _begin + _size;

  if(p1)
    for(T* p = _begin; p != _end; ++p, ++p1)
      *p = *p1;
}

template <typename T>
Array<T>::Array (int s, const T& t) throw(Error::General)
  : _size(s), _capacity(s)
{
  const char funame [] = "Array<T>::Array(int, const T&): ";

  if(_size < 0) {
    std::cerr << funame << "negative size\n";
    throw Error::Range();
  }

  if(!_size) {
    _end = _begin = 0;
    return;
  }

  _begin = new T[_size];
  _end = _begin + _size;

  for(T* p = _begin; p != _end; ++p)
    *p = t;
}

template <typename T>
Array<T>::Array (const Slice<T>& s)
  : _size(s.size()), _capacity(s.size())
{
  if(!_size) {
    _end = _begin = 0;
    return;
  }

  _begin = new T[_size];
  _end = _begin + _size;
  
  const T* p1 = s.start();
  for(T* p = _begin; p != _end; ++p, p1 += s.stride())
    *p = *p1; 
}

template <typename T>
template <typename S>
Array<T>::Array (const std::vector<S>& v)
  : _size(v.size()), _capacity(v.size())
{
  if(!_size) {
    _end = _begin = 0;
    return;
  }
  
  _begin = new T[_size];
  _end = _begin + _size;

  //std::vector<S>::const_iterator p1 = v.begin();
  int i = 0;
  for(T* p = _begin; p != _end; ++p, ++i)
      *p = v[i]; 
}

template <typename T>
Array<T>::Array (const ConstSlice<T>& s)
  : _size(s.size()), _capacity(s.size())
{
  if(!_size) {
    _end = _begin = 0;
    return;
  }

  _begin = new T[_size];
  _end = _begin + _size;

  const T* p1 = s.start();
  for(T* p = _begin; p != _end; ++p, p1 += s.stride())
      *p = *p1; 
}

template <typename T>
Array<T>::Array (const Array& ar) : _size(ar._size), _capacity(ar._size)
{
  if(!_size) {
    _end = _begin = 0;
    return;
  }

  _begin = new T[_size];
  _end = _begin + _size;

  const T* it1 = ar.begin();
  for (T* it = begin(); it != end(); ++it)
    *it = *it1++;
}

template <typename T>
void Array<T>::compact ()
{
    if(_capacity == _size)
	return;

    if(!_size) {
	delete[] _begin;
	_end = _begin = 0;
	_capacity = 0;
	return;
    }

    T* _new = new T[_size];
    _end = _new + _size;
    _capacity = _size;

    T* p1 = _begin;
    for(T* p = _new; p != _end; ++p)
	*p = *p1++;

    delete[] _begin;
    _begin = _new;
}

template <typename T>
void Array<T>::resize (int s) throw(Error::General)
{
  const char funame [] = "Array<T>::resize: ";

  if(s < 0) {
    std::cerr << funame << "negative array size\n";
    throw Error::Range();
  }

  if(s <= _capacity) {
    _size = s;
    _end = _begin + _size;
    return;
  }

  T* _new = new T[s];
  T* p1 = _new;
  for(T* p = begin(); p != end(); ++p)
    *p1++ = *p;

  if(_begin)
    delete[] _begin;

  _begin = _new;
  _capacity = _size = s;
  _end = _begin + _size;
}

template <typename T>
void Array<T>::reserve (int s) throw(Error::General)
{
  const char funame [] = "Array<T>::reserve: ";

  if(s < 0) {
    std::cerr << funame << "negative array size\n";
    throw Error::Range();
  }

  if(s <= _capacity)
      return;
  
  T* _new = new T[s];
  T* p1 = _new;
  for(T* p = begin(); p != end(); ++p)
    *p1++ = *p;

  if(_begin)
    delete[] _begin;

  _begin = _new;
  _capacity  = s;
  _end = _begin + _size;
}

template <typename T>
Array<T>& Array<T>::operator= (const Array& a1)
{
  const char funame [] = "Array<T>::operator= (const Array&): ";

  if(!a1.size()) {
    if(_begin)
      delete[] _begin;
    _begin = _end = 0;
    _size = _capacity = 0;
    return *this;
  }

  if(_size != a1.size()) {
    _size = a1.size();
    if(_capacity < a1.size()) {
      if(_begin)
	delete[] _begin;
      _begin = new T[a1.size()];
      _capacity = a1.size();
    }
    _end = _begin + _size;
  }

  T* it1 = begin();
  for(const T* it = a1.begin(); it != a1.end(); ++it, ++it1)
    *it1 = *it;
  return *this;
}

template <typename T>
template <typename S>
Array<T>& Array<T>::operator= (const std::vector<S>& a1)
{
  const char funame [] = "Array<T>::operator=(const std::vector<S>&): ";

  if(!a1.size()) {
    if(_begin)
      delete[] _begin;
    _begin = _end = 0;
    _size = _capacity = 0;
    return *this;
  }

  if(_size != a1.size()) {
    _size = a1.size();
    if(_capacity < a1.size()) {
      if(_begin)
	delete[] _begin;
      _begin = new T[a1.size()];
      _capacity = a1.size();
    }
    _end = _begin + _size;
  }

  T* it = begin();
  for(int i = 0; i < a1.size(); ++i, ++it)
    *it = a1[i];
  
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator= (const T* it1)
{
  for (T* it = begin(); it != end(); ++it)
    {
      *it = *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it = val;
  return *this;
}

template <typename T>
T& Array<T>::operator[] (int i) throw(Error::General)
{
  const char funame [] = "Array<T>::operator[]: ";

#ifdef DEBUG
  if (i >= _size || i < 0) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
#endif

  return _begin[i];
}

template <typename T>
const T& Array<T>::operator[] (int i) const throw(Error::General)
{
  const char funame [] = "Array<T>::operator[]: ";

#ifdef DEBUG
  if (i >= _size || i < 0) {
    std::cerr << funame << "out of range\n";
    throw Error::Range();
  }
#endif

  return _begin[i];
}

template <typename T>
T& Array<T>::front () throw(Error::General)
{
  if(!_size) {
    std::cerr << "Array<T>::front(): array is empty\n";
    throw Error::Range();
  }
  return *_begin;
}

template <typename T>
const T& Array<T>::front () const throw(Error::General)
{
  if(!_size) {
    std::cerr << "Array<T>::front(): array is empty\n";
    throw Error::Range();
  }
  return *_begin;
}

template <typename T>
T& Array<T>::back () throw(Error::General)
{
  if(!_size) {
    std::cerr << "Array<T>::back(): array is empty\n";
    throw Error::Range();
  }
  return _begin[_size-1];
}

template <typename T>
const T& Array<T>::back () const throw(Error::General)
{
  if(!_size) {
    std::cerr << "Array<T>::back(): array is empty\n";
    throw Error::Range();
  }
  return _begin[_size-1];
}

template <typename T>
Array<T>& Array<T>::operator- ()
{
  for (T* it = begin(); it != end(); ++it)
    *it *= -*it;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator+= (const Array& a1) throw(Error::General)
{
  const char funame [] = "Array<T>::operator+=: ";

  if (size() != a1.size()) {
    std::cerr << funame << "different size\n";
    throw Error::Range();
  }

  const T* it1 = a1.begin();
  for(T* it = begin(); it != end(); ++it) {
      *it += *it1;
      ++it1;
  }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const Array& a1) throw(Error::General)
{
  const char funame [] = "Array<T>::operator-=: ";

  if (size() != a1.size()) {
    std::cerr << funame << "different size\n";
    throw Error::Range();
  }

  const T* it1 = a1.begin();
  for (T* it = begin(); it != end(); ++it)
    {
      *it -= *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator+= (const T* it1)
{
  for (T* it = begin(); it != end(); ++it)
    {
      *it += *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const T* it1)
{
  for (T* it = begin(); it != end(); ++it)
    {
      *it -= *it1;
      ++it1;
    }
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator*= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it *= val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator/= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it /= val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator+= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it += val;
  return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const T& val)
{
  for (T* it = begin(); it != end(); ++it)
    *it -= val;
  return *this;
}

/**************************************************************************
 *************************** Array be reference  **************************
 **************************************************************************/

template <typename T>
class RefArr {

  Array<T>* _data;
  int*      _count;  // number of references

  void delete_ref ();
  void create_ref (const RefArr&);

public:

  RefArr () : _data(new Array<T>), _count(new int(1)) {}
  explicit RefArr (int s) : _data(new Array<T>(s)), _count(new int(1)) {}
  RefArr (int s, const T* p) : _data(new Array<T>(s, p)), _count(new int(1)) {}

  RefArr (const RefArr& a) { create_ref(a); }
  ~RefArr () { delete_ref(); }

  RefArr operator= (const RefArr& a) { delete_ref(); create_ref(a); }
  RefArr copy() const;// make a new copy
    
  void resize  (int s) { _data->resize(s);  }
  void reserve (int s) { _data->reserve(s); }
  void compact ()      { _data->compact();  }

  int size ()      const { return _data->size();     }
  int capacity ()  const { return _data->capacity(); }
  int ref_count () const { return *_count; }

  T* begin () { return _data->begin(); }
  T* end   () { return _data->end(); }
  const T* begin () const { return _data->begin(); }
  const T* end   () const { return _data->end(); }

  // conversions
  operator T* () { return *_data; }
  operator const T* () const { return *_data; }

  // indexing
  const T& operator[] (int i) const { return (*_data)[i]; }
  T&       operator[] (int i)       { return (*_data)[i]; }

  // comparisons
  bool operator== (const RefArr& a) const { return _data == a._data; }
  bool operator!= (const RefArr& a) const { return _data != a._data; }

  // arithmetic operations
  RefArr operator+ (const RefArr&) const throw(Error::General);
  RefArr operator- (const RefArr&) const throw(Error::General);

  RefArr operator+= (const RefArr& a) throw(Error::General)
  { *_data += *a._data; return *this; }
  RefArr operator-= (const RefArr& a) throw(Error::General)
  { *_data -= *a._data; return *this; }

  RefArr operator-  ()           { -(*_data)  ; return *this; }
  RefArr operator=  (const T& t) { *_data =  t; return *this; }
  RefArr operator+= (const T& t) { *_data += t; return *this; }
  RefArr operator-= (const T& t) { *_data -= t; return *this; }
  RefArr operator*= (const T& t) { *_data *= t; return *this; }
  RefArr operator/= (const T& t) { *_data /= t; return *this; }

};// class RefArr

template <typename T>
inline void RefArr<T>::create_ref (const RefArr& a)
{
    _data = a._data;
    _count =  a._count;
    ++(*_count);
}

template <typename T>
inline void RefArr<T>::delete_ref ()
{
  if(*_count > 1)
      --(*_count);
  else {
      delete _data;
      delete _count;
  }
}


template <typename T>
RefArr<T> RefArr<T>::copy() const
{
  const char funame [] = "RefArr<T>::copy: ";

  RefArr res(size());
  const T* p1 = *this;
  for(T* p = res.begin(); p != res.end(); ++p, ++p1)
    *p = *p1;

  return res;
}

template <typename T>
RefArr<T> RefArr<T>::operator+ (const RefArr& a) const throw(Error::General)
{
  const char funame [] = "RefArr<T>::operator+: ";

  if(a.size() != size()) {
    std::cerr << funame << "added array has a different size\n";
    throw Error::Range();
  }

  RefArr res(size());
  const T* p1 = *this;
  const T* p2 = a;
  for(T* p = res.begin(); p != res.end(); ++p, ++p1, ++p2)
    *p = *p1 + *p2;

  return res;
}

template <typename T>
RefArr<T> RefArr<T>::operator- (const RefArr& a) const throw(Error::General)
{
  const char funame [] = "RefArr<T>::operator-: ";

  if(a.size() != size()) {
    std::cerr << funame << "subtratcted array has a different size\n";
    throw Error::Range();
  }

  RefArr res(size());
  const T* p1 = *this;
  const T* p2 = a;
  for(T* p = res.begin(); p != res.end(); ++p, ++p1, ++p2)
    *p = *p1 - *p2;

  return res;
}

template<class C>
class Array_2 : public Array<C>
{// class Array_2

      int dd_1;
      int dd_2;

   public:

      Array_2 (int d1, int d2) : Array<C>(d1*d2), dd_1(d1), dd_2(d2) {}

      C operator () (int i1, int i2) const
      {
         return (*this) [i1 + dd_1 * i2];
      }

      C& operator () (int i1, int i2)
      {
         return (*this) [i1 + dd_1 * i2];
      }

      int dim_1 () const { return dd_1; }
      int dim_2 () const { return dd_2; }
};

template <class C>
class Array_3 : public Array<C>
{// class Array_3

      int dd_1;
      int dd_2;
      int dd_3;

   public:

      Array_3 (int d1, int d2, int d3)
        : Array<C>(d1*d2*d3), dd_1(d1), dd_2(d2), dd_3(d3)
      {}

      C operator () (int i1, int i2, int i3) const
      {
         return (*this)[i1 + dd_1 * i2 + dd_1 *dd_2 * i3];
      }

      C& operator () (int i1, int i2, int i3)
      {
         return (*this)[i1 + dd_1 * i2 + dd_1 *dd_2 * i3];
      }

      int dim_1 () const { return dd_1; }
      int dim_2 () const { return dd_2; }
      int dim_3 () const { return dd_3; }

};

/************************************************************************************************************
 ****************************************** MULTI-DIMENSIONAL ARRAY *****************************************
 ************************************************************************************************************/

template <typename T>
class MultiArray : private Array<T> {
  //
  std::vector<int> _size;

  void             _isinit ()         const;

  void             _resize (va_list&);
  int              _index  (va_list&) const;
  std::vector<int> _slice  (va_list&) const;
  
public:
  //
  MultiArray () {}
  
  explicit MultiArray (int ...);

  void resize (int ...);

  const T& operator () (...) const;
  T&       operator () (...);

  ConstSlice<T> slice (...) const;
  Slice<T>      slice (...);

  int  rank ()      const { return _size.size(); }
  int  size (int i) const { return _size[i]; }
  int  size ()      const { return Array<T>::size(); }

  typedef       T*   iterator;
  typedef const T*   const_iterator;
  typedef       T    value_type;
  typedef       int  size_type;

  operator       T* ()       { return *this; }
  operator const T* () const { return *this; }

  T*       begin ()       { return Array<T>::begin(); }
  const T* begin () const { return Array<T>::begin(); }

  T*       end ()       { return Array<T>::end(); }
  const T* end () const { return Array<T>::end(); }
};

template <typename T>
void MultiArray<T>::_isinit () const
{
  Exception::Base funame = "MultiArray::_isinit: ";

  if(!rank())
    //
    throw funame << "not initialized";
}

// resizing
//
template <typename T>
void MultiArray<T>::_resize(va_list& ap) 
{
  Exception::Base funame = "MultiArray::_resize: ";

  double dtemp = 1.;
  
  int    itemp = 1;

  for(typename std::vector<int>::iterator it = _size.begin(); it != _size.end(); itemp *= *it++) {
    //
    *it = va_arg(ap, int);
    
    if(*it < 1)
      //
      throw funame << it - _size.begin() << "-th dimension out of range: " << *it;

    dtemp *= double(*it);
  }

  if(dtemp > 2.e9)
    //
    throw funame << "linear size out of range: " << dtemp;
  
  Array<T>::resize(itemp);
}

template <typename T>
MultiArray<T>::MultiArray (int r ...)
{
  Exception::Base funame = "MultiArray::MultiArray: ";

  if(r < 1)
    //
    throw funame << "rank out of range: " << r;

  _size.resize(r);
  
  va_list ap;
  
  va_start(ap, r);

  _resize(ap);
  
  va_end(ap);
}

template <typename T>
void MultiArray<T>::resize (int r ...)
{
  Exception::Base funame = "MultiArray::resize: ";

  if(r < 1)
    //
    throw funame << "rank out of range: " << r;

  _size.resize(r);
  
  va_list ap;
  
  va_start(ap, r);

  _resize(ap);
  
  va_end(ap);
}

// indexing
//
template <typename T>
int MultiArray<T>::_index (va_list& ap) const
{
  Exception::Base funame = "MultiArray::_index: ";

  _isinit();

  int res = 0;
  
  int step  = 1;

  for(int r = 0; r < rank(); step *= _size[r++]) {
    //
    int i = va_arg(ap, int);

    if(i < 0 || i >= _size[r])
      throw funame << r << "-th index out of range: " << i;

    res += i * step;
  }

  return res;
}

template <typename T>
const T& MultiArray<T>::operator() (...) const
{
  va_list ap;

  va_start(ap, *this);

  int i = _index(ap);

  va_end(ap);
  
  return (*this)[i];
}

template <typename T>
T& MultiArray<T>::operator() (...)
{
  va_list ap;
  
  va_start(ap, *this);

  int i = _index(ap);

  va_end(ap);
  
  return (*this)[i];
}

// Slice
//
template <typename T>
std::vector<int> MultiArray<T>::_slice (va_list& ap) const
{
  Exception::Base funame = "MultiArray::_slice: ";

  _isinit();

  int itemp;
  
  int start  =  0;
  
  int stride = -1;
  
  int step   =  1;

  
  for(int r = 0; r < rank(); step *= _size[r++]) {
    //
    int i = va_arg(ap, int);

    if(i >= _size[r])
      //
      throw funame << r << "-th index out of range: " << i;

    if(i >= 0) {
      //
      start += i * step;
    }
    else {
      //
      if(stride > 0)
	//
	throw funame << "second sliced index: " << itemp << ", " << r;
      
      stride = step;

      itemp = r;
    }
  }

  if(stride < 0)
    //
    throw funame << "no sliced index";

  std::vector<int> res(3);

  res[0] = start;
  res[1] = _size[itemp];
  res[2] = stride;

  return res;
}

template <typename T>
Slice<T> MultiArray<T>::slice (...)
{
  va_list ap;
  
  va_start(ap, *this);

  std::vector<int> s = _slice(ap);
    
  va_end(ap);
  
  return Slice<T>(*this + s[0], s[1], s[2]);
}

template <typename T>
ConstSlice<T> MultiArray<T>::slice (...) const
{
  va_list ap;
  
  va_start(ap, *this);

  std::vector<int> s = _slice(ap);
    
  va_end(ap);
  
  return ConstSlice<T>(*this + s[0], s[1], s[2]);
}

#endif
