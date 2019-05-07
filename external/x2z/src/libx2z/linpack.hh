#ifndef LINPACK_HH
#define LINPACK_HH

#include "error.hh"
#include <iostream>

/********************** General Vector Operations ***************************/

double normalize (double*, int);
double orthogonalize (double*, const double*, int) ;
double vdistance (const double*, const double*, int);
double vlength (const double*, int);
double vdot (const double*, const double*, int);

/*********** PseudoArray, the assignment is by reference **************************/

template <typename T> 
class PseudoArray 
{
    T*  _begin;
    int _size;

public:
    PseudoArray (T*, int) ;
    T& operator[] (int i) const ;
    int size() { return _size; }
};

template <typename T>
PseudoArray<T>::PseudoArray (T* p, int s)  : _begin(p), _size(s)
{
  const char funame [] = "PseudoArray<T>::PseudoArray: ";

#ifdef DEBUG
  if(!p) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(s < 0) {
    std::cerr << funame << "negative size = " << s << std::endl;
    throw Error::Range();
  }
#endif
}

template <typename T>
T& PseudoArray<T>::operator[] (int i) const 
{
    const char funame [] = "PseudoArray<T>::operator[]: ";

#ifdef DEBUG
    if(i < 0 || i >= _size) {
	std::cerr << funame << "index i = " << i
		  << " out of range (size = " << _size
		  << ")\n";
	throw Error::Range();
  }
#endif

    return _begin[i];
}

/*********************************** Stride Pointer **************************/

template <typename T> 
class StridePointer 
{
  T* _pnt;
  int _stride;

public:
  StridePointer (T*, int =1) ;
  int stride () const { return _stride; }
  T& operator*  () const { return *_pnt; }
  T* operator-> () const { return  _pnt; }
  operator T* () const { return _pnt; }
  StridePointer& operator++ ()     { _pnt += _stride; return *this; }
  StridePointer& operator-- ()     { _pnt -= _stride; return *this; }
  StridePointer  operator++ (int)  { StridePointer p = *this; _pnt += _stride; return p; }
  StridePointer  operator-- (int)  { StridePointer p = *this; _pnt -= _stride; return p; }
  StridePointer& operator+= (int d) { _pnt += d * _stride; return *this; }
  StridePointer& operator-= (int d) { _pnt -= d * _stride; return *this; }

  StridePointer operator+   (int d) const { return StridePointer(_pnt + d * _stride, _stride); }
  StridePointer operator-   (int d) const { return StridePointer(_pnt - d * _stride, _stride); }

  int operator-   (const StridePointer&) const ;
};

template <typename T> 
StridePointer<T> operator+ (int d, const StridePointer<T>& p) 
{ 
  return StridePointer<T>(p + d * p.stride(), p.stride());
}

template <typename T> 
int StridePointer<T>::operator- (const StridePointer& sp) const 
{
  const char funame [] = "StridePointer::operator- (const StridePointer&): ";

  if(_stride != sp._stride) {
    std::cerr << funame << "strides are different\n";
    throw Error::Range();
  }

  if(((_pnt - sp._pnt) % _stride)) {
    std::cerr << funame << "different slices\n";
    throw Error::Range();
  }
 
  return (_pnt - sp._pnt) / _stride;

}

template <typename T>
StridePointer<T>::StridePointer (T* p, int s)  : _pnt(p), _stride(s) 
{
  const char funame [] = "StridePointer<T>::StridePointer: ";

#ifdef DEBUG
  if(!p) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(s <= 0) {
    std::cerr << funame << "non-positive stride = " << s << std::endl;
    throw Error::Range();
  }
#endif
}

/****************************** Constant Stride Pointer **************************/

template <typename T> 
class ConstStridePointer 
{
  const T* _pnt;
  int _stride;

public:
  ConstStridePointer (const T*, int =1) ;
  ConstStridePointer (const StridePointer<T>& p) : _pnt(p), _stride(p.stride()) {}
  int stride () const { return _stride; }
  const T& operator* () const { return *_pnt; }
  const T* operator-> () const { return  _pnt; }
  operator const T* () const { return _pnt; }
  ConstStridePointer& operator++ ()      { _pnt += _stride; return *this; }
  ConstStridePointer& operator-- ()      { _pnt -= _stride; return *this; }
  ConstStridePointer  operator++ (int)   
  { ConstStridePointer p = *this; _pnt += _stride; return p; }
  ConstStridePointer  operator-- (int)   
  { ConstStridePointer p = *this; _pnt -= _stride; return p; }
  ConstStridePointer& operator+= (int d) { _pnt += d * _stride; return *this; }
  ConstStridePointer& operator-= (int d) { _pnt -= d * _stride; return *this; }

  ConstStridePointer operator+   (int d) const 
  { return ConstStridePointer(_pnt + d * _stride, _stride); }
  ConstStridePointer operator-   (int d) const 
  { return ConstStridePointer(_pnt - d * _stride, _stride); }

  int operator- (const ConstStridePointer&) const ;
};

template <typename T> 
ConstStridePointer<T> operator+ (int d, const ConstStridePointer<T>& p) 
{ 
  return ConstStridePointer<T>(p + d * p.stride(), p.stride());
}

template <typename T> 
int ConstStridePointer<T>::operator- (const ConstStridePointer& sp) const 
{
  const char funame [] = "ConstStridePointer::operator- (const ConstStridePointer&): ";

  if(_stride != sp._stride) {
    std::cerr << funame << "strides are different\n";
    throw Error::Range();
  }

  if(((_pnt - sp._pnt) % _stride)) {
    std::cerr << funame << "different slices\n";
    throw Error::Range();
  }
 
  return (_pnt - sp._pnt) / _stride;

}

template <typename T>
ConstStridePointer<T>::ConstStridePointer (const T* p, int s)  
    : _pnt(p), _stride(s) 
{
  const char funame [] = "ConstStridePointer<T>::ConstStridePointer: ";

#ifdef DEBUG
  if(!p) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(s <= 0) {
    std::cerr << funame << "non-positive stride = " << s << std::endl;
    throw Error::Range();
  }
#endif
}

template<typename T>
T vdot (StridePointer<T> p1, StridePointer<T> p2, int n)
{
  T res = 0;
  for(int i = 0; i < n; ++i) {
    res += *p1 * *p2; 
    ++p1;
    ++p2;
  }
  return res;
}

template<typename T>
T vdot (ConstStridePointer<T> p1, ConstStridePointer<T> p2, int n)
{
  T res = 0;
  for(int i = 0; i < n; ++i) {
    res += *p1 * *p2; 
    ++p1;
    ++p2;
  }
  return res;
}

/**************************************************************************
 **************************** Constant Slice ******************************
 **************************************************************************/

template <typename T>
class Slice;

template <typename T>
class ConstSlice {
    const T* _start;
    const T* _stop;
    int _size;
    int _stride;

    // no assignments
    ConstSlice& operator= (const ConstSlice&);

public:
    typedef ConstStridePointer<T> iterator;
    iterator begin () { return ConstStridePointer<T>(_start, _stride); }
    iterator end   () { return ConstStridePointer<T>(_stop , _stride); }

    ConstSlice(const T* st, int sz, int sd =1) ;

    const T*  start  () const { return _start; }

    int size   () const { return _size; }
    int stride () const { return _stride; }

    const T&  operator[] (int) const ;
    
    // scalar product
    T operator* (const ConstSlice&) const ;
    T operator* (const Slice<T>&)   const ;
    T operator* (const T*)    const;

    T max (int* =0) const;
    T min (int* =0) const;
    T sum ()        const;
    T product ()    const;
    T vdot ()       const;
};

template <typename T>
ConstSlice<T>::ConstSlice(const T* st, int sz, int sd) 
  : _start(st), _size(sz), _stride(sd), _stop(st + sz * sd)
{
  const char funame [] = "ConstSlice<T>::ConstSlice: ";

#ifdef DEBUG
  if(!st) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(sz < 0) {
    std::cerr << funame << "negative size = " << sz << std::endl;
    throw Error::Range();
  }

  if(sd <= 0) {
    std::cerr << funame << "non-positive stride = " << sd << std::endl;
    throw Error::Range();
  }
#endif

}

template <typename T>
const T& ConstSlice<T>::operator[] (int i) const 
{
  const char funame [] = "ConstSlice<T>::operator[]: ";

#ifdef DEBUG
  if(i < 0 || i >= _size) {
    std::cerr << funame << "out of range: index = " << i << std::endl;
    throw Error::Range();
  }
#endif

  return _start[i * _stride];
}

template <typename T>
T ConstSlice<T>::operator* (const Slice<T>& s) const 
{
  const char funame [] = "ConstSlice<T>::operator*(const Slice<T>&): ";
  
#ifdef DEBUG
  if(_size != s.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s.start();

  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride) {
    res += *p * *p1;
    p1 += s.stride();
  }
  return res;
}

template <typename T>
T ConstSlice<T>::operator* (const ConstSlice& s) const 
{
  const char funame [] = "ConstSlice<T>::operator* (const ConstSlice&): ";

#ifdef DEBUG  
  if(_size != s._size) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s._start;

  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride) {
    res += *p * *p1;
    p1 += s._stride;
  }
  return res;
}

template <typename T>
T ConstSlice<T>::operator* (const T* p1) const
{
    T res = 0;
    for(const T* p = _start; p != _stop; p += _stride) {
	res += *p * *p1++;
    }
    return res;
}

template <typename T>
T ConstSlice<T>::max(int* ip) const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    if(p == _start || *p > res) {
      res = *p;
      if(ip)
	*ip = p - _start;
    }
  return res;
}

template <typename T>
T ConstSlice<T>::min(int* ip) const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    if(p == _start || *p < res) {
      res = *p;
      if(ip)
	*ip = p - _start;
    }
  return res;
}

template <typename T>
T ConstSlice<T>::sum() const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    res += *p;
  return res;
}

template <typename T>
T ConstSlice<T>::product() const
{
  T res = 1;
  for(const T* p = _start; p != _stop; p += _stride)
    res *= *p;
  return res;
}

template <typename T>
T ConstSlice<T>::vdot() const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    res += *p * *p;
  return res;
}

/**************************************************************************
 ********************************** Slice *********************************
 **************************************************************************/

template <typename T>
class Slice 
{
  T* _start;
  T* _stop;
  int _size;
  int _stride;

public:
  typedef StridePointer<T> iterator;
  iterator begin () { return StridePointer<T>(_start, _stride); }
  iterator end   () { return StridePointer<T>(_stop , _stride); }

  Slice(T* st, int sz, int sd =1) ;

  const T*  start  () const { return _start; }
  T*        start  ()       { return _start; }

  int size   () const { return _size; }
  int stride () const { return _stride; }

  const T&  operator[] (int) const ;
  T&        operator[] (int)       ;

  // scalar product
  T operator* (const Slice&)         const ;
  T operator* (const ConstSlice<T>&) const ;
  T operator* (const T*)     const;

  T max (int* =0) const;
  T min (int* =0) const;
  T sum ()        const;
  T product ()    const;
  T vdot ()       const;
  
  // block operations
  void operator=  (const Slice&) ;// copy values
  void operator+= (const Slice&) ;
  void operator-= (const Slice&) ;
  void operator*= (const Slice&) ;
  void operator/= (const Slice&) ;

  void operator=  (const ConstSlice<T>&) ;// copy values
  void operator+= (const ConstSlice<T>&) ;
  void operator-= (const ConstSlice<T>&) ;
  void operator*= (const ConstSlice<T>&) ;
  void operator/= (const ConstSlice<T>&) ;

  void operator=  (const T*);
  void operator+= (const T*);
  void operator-= (const T*);
  void operator*= (const T*);
  void operator/= (const T*);

  void operator-  ();
  void operator=  (const T&);
  void operator+= (const T&);
  void operator-= (const T&);
  void operator*= (const T&);
  void operator/= (const T&);
};

template <typename T>
Slice<T>::Slice(T* st, int sz, int sd) 
  : _start(st), _size(sz), _stride(sd), _stop(st + sz * sd)
{
  const char funame [] = "Slice<T>::Slice: ";

#ifdef DEBUG
  if(!st) {
    std::cerr << funame << "null data pointer\n";
    throw Error::Range();
  }

  if(sz < 0) {
    std::cerr << funame << "negative size = " << sz << std::endl;
    throw Error::Range();
  }

  if(sd <= 0) {
    std::cerr << funame << "negative stride = " << sd << std::endl;
    throw Error::Range();
  }
#endif

}

template <typename T>
const T& Slice<T>::operator[] (int i) const 
{
  const char funame [] = "Slice<T>::operator[]: ";

#ifdef DEBUG
  if(i < 0 || i >= _size) {
    std::cerr << funame << "out of range: index = " << i << std::endl;
    throw Error::Range();
  }
#endif

  return _start[i * _stride];
}

template <typename T>
T& Slice<T>::operator[] (int i) 
{
  const char funame [] = "Slice<T>::operator[]: ";

#ifdef DEBUG
  if(i < 0 || i >= _size) {
    std::cerr << funame << "out of range: index = " << i << std::endl;
    throw Error::Range();
  }
#endif

  return _start[i * _stride];
}

template <typename T>
T Slice<T>::operator* (const Slice& s) const 
{
  const char funame [] = "Slice<T>::operator*(const Slice&): ";

#ifdef DEBUG  
  if(_size != s._size) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s._start;
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride) {
    res += *p * *p1;
    p1 += s._stride;
  }
  return res;
}

template <typename T>
T Slice<T>::operator* (const ConstSlice<T>& s) const 
{
  const char funame [] = "Slice<T>::operator*(const ConstSlice<T>&): ";
  
#ifdef DEBUG
  if(_size != s.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s.size()
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s.start();
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride) {
    res += *p * *p1;
    p1 += s.stride();
  }
  return res;
}

template <typename T>
T Slice<T>::operator* (const T* p1) const
{
    T res = 0;
    for(const T* p = _start; p != _stop; p += _stride) {
	res += *p * *p1++;
    }
    return res;
}

template <typename T>
T Slice<T>::max(int* ip) const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    if(p == _start || *p > res) {
      res = *p;
      if(ip)
	*ip = p - _start;
    }
  return res;
}

template <typename T>
T Slice<T>::min(int* ip) const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    if(p == _start || *p < res) {
      res = *p;
      if(ip)
	*ip = p - _start;
    }
  return res;
}

template <typename T>
T Slice<T>::sum() const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    res += *p;
  return res;
}

template <typename T>
T Slice<T>::product() const
{
  T res = 1;
  for(const T* p = _start; p != _stop; p += _stride)
    res *= *p;
  return res;
}

template <typename T>
T Slice<T>::vdot() const
{
  T res = 0;
  for(const T* p = _start; p != _stop; p += _stride)
    res += *p * *p;
  return res;
}

template <typename T>
void Slice<T>::operator= (const Slice& s) 
{
  const char funame [] = "Slice<T>::operator= (const Slice&): ";
 
#ifdef DEBUG 
  if(_size != s._size) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif 
 
  const T* p1 = s._start;
  for(T* p = _start; p != _stop; p += _stride) {
    *p = *p1;
    p1 += s._stride;
  }
}

template <typename T>
void Slice<T>::operator+= (const Slice& s) 
{
  const char funame [] = "Slice<T>::operator+= (const Slice&): ";
  
#ifdef DEBUG
  if(_size != s._size) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s._start;
  for(T* p = _start; p != _stop; p += _stride) {
    *p += *p1;
    p1 += s._stride;
  }
}

template <typename T>
void Slice<T>::operator-= (const Slice& s) 
{
  const char funame [] = "Slice<T>::operator-= (const Slice&): ";
  
#ifdef DEBUG
  if(_size != s._size) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s._start;
  for(T* p = _start; p != _stop; p += _stride) {
    *p -= *p1;
    p1 += s._stride;
  }
}

template <typename T>
void Slice<T>::operator*= (const Slice& s) 
{
  const char funame [] = "Slice<T>::operator*= (const Slice&): ";
  
#ifdef DEBUG
  if(_size != s._size) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s._start;
  for(T* p = _start; p != _stop; p += _stride) {
    *p *= *p1;
    p1 += s._stride;
  }
}

template <typename T>
void Slice<T>::operator/= (const Slice& s) 
{
  const char funame [] = "Slice<T>::operator/= (const Slice&): ";
  
#ifdef DEBUG
  if(_size != s._size) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s._size
	      << std::endl;
    throw Error::Range();
  }
#endif

  const T* p1 = s._start;
  for(T* p = _start; p != _stop; p += _stride) {
    *p /= *p1;
    p1 += s._stride;
  }
}

template <typename T>
void Slice<T>::operator= (const ConstSlice<T>& s) 
{
  const char funame [] = "Slice<T>::operator= (const ConstSlice<T>&): ";
  
#ifdef DEBUG
  if(_size != s.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s.size()
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s.start();
  for(T* p = _start; p != _stop; p += _stride) {
    *p = *p1;
    p1 += s.stride();
  }
}

template <typename T>
void Slice<T>::operator+= (const ConstSlice<T>& s) 
{
  const char funame [] = "Slice<T>::operator+= (const ConstSlice<T>&): ";
  
#ifdef DEBUG
  if(_size != s.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s.size()
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s.start();
  for(T* p = _start; p != _stop; p += _stride) {
    *p += *p1;
    p1 += s.stride();
  }
}

template <typename T>
void Slice<T>::operator-= (const ConstSlice<T>& s) 
{
  const char funame [] = "Slice<T>::operator-= (const ConstSlice<T>&): ";
  
#ifdef DEBUG
  if(_size != s.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s.size()
	      << std::endl;
    throw Error::Range();
  }
#endif
  
  const T* p1 = s.start();
  for(T* p = _start; p != _stop; p += _stride) {
    *p -= *p1;
    p1 += s.stride();
  }
}

template <typename T>
void Slice<T>::operator*= (const ConstSlice<T>& s) 
{
  const char funame [] = "Slice<T>::operator*= (const ConstSlice<T>&): ";
  
#ifdef DEBUG
  if(_size != s.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s.size()
	      << std::endl;
    throw Error::Range();
  }
#endif

  const T* p1 = s.start();
  for(T* p = _start; p != _stop; p += _stride) {
    *p *= *p1;
    p1 += s.stride();
  }
}

template <typename T>
void Slice<T>::operator/= (const ConstSlice<T>& s) 
{
  const char funame [] = "Slice<T>::operator/= (const ConstSlice<T>&): ";
  
#ifdef DEBUG
  if(_size != s.size()) {
    std::cerr << funame << "dimensions are different:"
	      << " size1 = " << _size
	      << " size2 = " << s.size()
	      << std::endl;
    throw Error::Range();
  }
#endif  

  const T* p1 = s.start();
  for(T* p = _start; p != _stop; p += _stride) {
    *p /= *p1;
    p1 += s.stride();
  }
}

template <typename T>
void Slice<T>::operator= (const T* p1)
{
    for(T* p = _start; p != _stop; p += _stride)
	*p = *p1++;
}

template <typename T>
void Slice<T>::operator+= (const T* p1)
{
    for(T* p = _start; p != _stop; p += _stride)
	*p += *p1++;
}

template <typename T>
void Slice<T>::operator-= (const T* p1)
{
    for(T* p = _start; p != _stop; p += _stride)
	*p -= *p1++;
}

template <typename T>
void Slice<T>::operator*= (const T* p1)
{
    for(T* p = _start; p != _stop; p += _stride)
	*p *= *p1++;
}

template <typename T>
void Slice<T>::operator/= (const T* p1)
{
    for(T* p = _start; p != _stop; p += _stride)
	*p /= *p1++;
}

template <typename T>
void Slice<T>::operator- ()
{
  for(T* p = _start; p != _stop; p += _stride)
    *p = -*p;
}

template <typename T>
void Slice<T>::operator= (const T& t)
{
  for(T* p = _start; p != _stop; p += _stride)
    *p = t;
}

template <typename T>
void Slice<T>::operator+= (const T& t)
{
  for(T* p = _start; p != _stop; p += _stride)
    *p += t;
}

template <typename T>
void Slice<T>::operator-= (const T& t)
{
  for(T* p = _start; p != _stop; p += _stride)
    *p -= t;
}

template <typename T>
void Slice<T>::operator*= (const T& t)
{
  for(T* p = _start; p != _stop; p += _stride)
    *p *= t;
}

template <typename T>
void Slice<T>::operator/= (const T& t)
{
  for(T* p = _start; p != _stop; p += _stride)
    *p /= t;
}

#endif
