#include "d3.hh"
#include <iostream>

/*********************** 3-D reference frame ***********************/

void D3::Frame::fv2lv (const double* fv, double* lv) const // frame vector to lab vector
{
  vprod(fv, orient, lv);
}

void D3::Frame::fp2lp (const double* fp, double* lp) const // frame to lab
{
  vprod(fp, orient, lp);
  for(int i = 0; i < 3; ++i)
    lp[i] += origin[i];
}

void D3::Frame::lv2fv (const double* lv, double* fv) const // lab to frame
{
  vprod(orient, lv, fv);
}

void D3::Frame::lp2fp (const double* lp, double* fp) const // lab to frame
{
  vprod(orient, lp - origin, fp);
}


/*********************** 3-D plane ***********************/

void D3::Plane::set (const Vector& n, double d)
{
  _normal = n;
  _dist = d / _normal.normalize();
  find_orth(_normal, _orth);
}

D3::Plane::Plane () 
{
  _normal = 0.;
  _normal[0] = 1.;
  _dist = 1.;
  find_orth(_normal, _orth);
}

/*********************** 3-D rotational matrix ***********************/

D3::Matrix::Matrix (Vector n1, Vector n2)  
{
  n1.normalize();
  row(0) = n1;

  n2.orthogonalize(n1);
  n2.normalize();  
  row(1) = n2;

  row(2) = vprod(n1, n2);  
}

double& D3::Matrix::operator() (int i, int j) 
{
  const char funame [] = "D3::Matrix::operator() (int, int):";

#ifdef DEBUG
  if(i < 0 || i > 2 || j < 0 || j > 2) {
    std::cerr << funame << "indices out of range\n";
    throw Error::Range();
  }
#endif

  return *(_data + (i * 3 +  j));
}

double  D3::Matrix::operator() (int i, int j) const 
{
  const char funame [] = "D3::Matrix::operator() (int, int) const:";

#ifdef DEBUG
  if(i < 0 || i > 2 || j < 0 || j > 2) {
    std::cerr << funame << "indices out of range\n";
    throw Error::Range();
  }
#endif

  return *(_data + (i * 3 +  j));
}

Slice<double> D3::Matrix::column (int i)  
{
  const char funame [] = "D3::Matrix::column (int):";

#ifdef DEBUG
  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }
#endif

  return Slice<double>(_data + i, 3, 3);
}

ConstSlice<double> D3::Matrix::column (int i) const 
{
  const char funame [] = "D3::Matrix::column (int) const:";

#ifdef DEBUG
  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }
#endif

  return ConstSlice<double>(_data + i, 3, 3);
}

Slice<double> D3::Matrix::row (int i) 
{
  const char funame [] = "D3::Matrix::row (int):";

#ifdef DEBUG
  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }
#endif

  return Slice<double>(_data + i*3, 3);

}
ConstSlice<double> D3::Matrix::row (int i) const 
{
  const char funame [] = "D3::Matrix::row (int) const:";

#ifdef DEBUG
  if(i < 0 || i > 2) {
    std::cerr << funame << "index out of range\n";
    throw Error::Range();
  }
#endif

  return ConstSlice<double>(_data + i*3, 3);
}

D3::Matrix D3::Matrix::operator* (const Matrix& r) const
{
  Matrix res;
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      res(i, j) = row(i) * r.column(j);

  return res;
}

D3::Vector D3::Matrix::operator* (const double* v) const
{
  Vector res;
  for(int i = 0; i < 3; ++i)
    res[i] = row(i) * v;

  return res;
}

void D3::Matrix::orthogonality_check () const 
{
  const char funame [] = "D3::Matrix::orthogonality_check: ";

  static const double tol = 1.e-12;

  double dtemp;
  for(int i = 0; i < 3; ++i) {
    dtemp = row(i) * row(i) - 1.;
    if(dtemp < -tol || dtemp > tol) {
      std::cerr << funame << "failed\n";
      throw Error::Run();
    }
    for(int j = i + 1; j < 3; ++j) {
      dtemp = row(i) * row(j);
      if(dtemp < -tol || dtemp > tol) {
	std::cerr << funame << "failed\n";
	throw Error::Run();
      }
    }
  } 
}

/******************************************
 ************ 3-D Real Vector *************
 ******************************************/

D3::Vector::Vector(const Vector& v)
{
  for(int i = 0; i < 3; ++i)
    _data[i] = v._data[i];
}

D3::Vector::Vector (const double* p)
{
  for(int i = 0; i < 3; ++i)
    _data[i] = *p++;
}

D3::Vector& D3::Vector::operator= (const Vector& v)
{

  for(int i = 0; i < 3; ++i)
    _data[i] = v._data[i];
  
  return *this;
}

D3::Vector& D3::Vector::operator= (const double* p)
{

  for(int i = 0; i < 3; ++i)
    _data[i] = *p++;
  
  return *this;
}

D3::Vector& D3::Vector::operator= (double d)
{

  for(int i = 0; i < 3; ++i)
    _data[i] = d;
  
  return *this;
}

D3::Vector D3::Vector::operator+ (const double* a2) const
{
  Vector res;
  double* p = res._data;
  const double* p1 = _data;

  for(int i = 0; i < 3; ++i)
    *p++ = *p1++ + *a2++;

  return res;
}

D3::Vector D3::Vector::operator- (const double* a2) const
{
  Vector res;
  double* p = res._data;
  const double* p1 = _data;

  for(int i = 0; i < 3; ++i)
    *p++ = *p1++ - *a2++;

  return res;
}

D3::Vector D3::Vector::operator* (double d) const
{
  Vector res;
  double* p = res._data;
  const double* p1 = _data;

  for(int i = 0; i < 3; ++i)
    *p++ = *p1++ * d;

  return res;
}

D3::Vector D3::Vector::operator* (const Matrix& m) const
{
  Vector res;
  vprod(*this, m, res);
  return res;
}

D3::Vector& D3::Vector::operator+= (const double* a)
{
  double* p = _data;

  for(int i = 0; i < 3; ++i)
    *p++ += *a++;

  return *this;
}

D3::Vector& D3::Vector::operator-= (const double* a)
{
  double* p = _data;

  for(int i = 0; i < 3; ++i)
    *p++ -= *a++;

  return *this;
}

D3::Vector& D3::Vector::operator*= (double d)
{
  double* p = _data;

  for(int i = 0; i < 3; ++i)
    *p++ *= d;

  return *this;
}

D3::Vector& D3::Vector::operator/= (double d)
{
  double* p = _data;

  for(int i = 0; i < 3; ++i)
    *p++ /= d;

  return *this;
}

D3::Vector& D3::Vector::operator*= (const Matrix& r) //rotation
{
  return *this = r * (*this);
}

double D3::Vector::vdot () const 
{
  return ::vdot(_data, _data, 3);
}

double D3::Vector::vlength () const 
{
  return ::vlength(_data, 3);
}

double D3::Vector::normalize () 
{
  return ::normalize(_data, 3);
}

double D3::Vector::orthogonalize (const double* n) 
{
  return ::orthogonalize(_data, n, 3);
}

D3::Vector D3::operator* (double d, const Vector& a)
{
  return a * d;
}

D3::Vector D3::operator+ (const double* b, const D3::Vector& a)
{
  return a + b;
}

D3::Vector D3::operator- (const double* p, const D3::Vector& v)
{
  D3::Vector res;
  for(int i = 0; i < 3; ++i)
    res[i] = *p++ - v[i];
  return res;
}

double vdistance (const D3::Vector& a1, const D3::Vector& a2)
{
  return ::vdistance((const double*)a1, (const double*)a2, 3);
}

double vdot (const D3::Vector& a1, const D3::Vector& a2)
{
  return ::vdot((const double*) a1, (const double*) a2, 3);
}

D3::Vector vprod (const double* a1, const double* a2)
{
  D3::Vector res;
  vprod(a1, a2, res);
  return res;
}

void vprod (const double* a1, const double* a2, double* res)
{
  res[0] = a1[1] * a2[2] - a1[2] * a2[1];
  res[1] = a1[2] * a2[0] - a1[0] * a2[2];
  res[2] = a1[0] * a2[1] - a1[1] * a2[0];
}

void vprod (const D3::Matrix& r, const double* v, double* res)
{
  for(int i = 0; i < 3; ++i)
    res[i] = r.row(i) * v;
}

void vprod (const double* v, const D3::Matrix& r, double* res)
{
  for(int i = 0; i < 3; ++i)
    res[i] = r.column(i) * v;
}

void find_orth (const double* n0, D3::Vector ort[2]) // get two vectors orthogonal to the given one
{
  int imin;
  double nmin;
  double dtemp;
  for(int i = 0; i < 3; ++i) {
    dtemp = n0[i] > 0. ? n0[i] : -n0[i];
    if(!i || dtemp < nmin) {
      imin = i;
      nmin = dtemp;
    }
  }
  
  for(int i = 0; i < 3; ++i)
    if(i == imin)
      ort[0][i] = 1.;
    else
      ort[0][i] = 0.;

  ort[0].orthogonalize(n0);
  ort[0].normalize();

  vprod(n0, ort[0], ort[1]);
  ort[1].normalize();
}

double vol (const D3::Vector& a1, const D3::Vector& a2, const D3::Vector& a3)
{
  return vdot(vprod(a1, a2), a3);
}

