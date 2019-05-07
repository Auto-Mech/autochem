#ifndef D3_HH
#define D3_HH

#include "linpack.hh"

/******************************************
 ************ 3-D Space  *****************
 ******************************************/

namespace D3 {
  class Matrix;

  class Vector
  {
    double _data [3];

  public:

    Vector () {}
    Vector (const Vector&);
    explicit Vector (const double*);

    Vector& operator= (const Vector&);
    Vector& operator= (const double*);
    Vector& operator= (double);

    double& operator [] (int i)       { return _data[i]; }
    double  operator [] (int i) const { return _data[i]; }

    operator       double* ()       { return _data; }
    operator const double* () const { return _data; }

    Vector operator+ (const double*)   const;
    Vector operator+ (const Vector& v) const { return *this + (const double*)v; }
    Vector operator- (const double*)   const;
    Vector operator- (const Vector& v) const { return *this - (const double*)v; }
    Vector operator* (double)          const;
    Vector operator* (const Matrix&)   const; // (row)v*M

    Vector& operator+= (const double*);
    Vector& operator-= (const double*);

    Vector& operator*= (double);
    Vector& operator/= (double);

    Vector& operator*= (const Matrix&); // orthogonal transformation M*v

    double vlength () const;
    double vdot    () const;

    double normalize ();
    double orthogonalize (const double*) ;
  };

  Vector operator+ (const double*, const D3::Vector&);
  Vector operator- (const double*, const D3::Vector&);
  Vector operator* (double, const D3::Vector&);

  inline std::ostream& operator<< (std::ostream& to, const Vector& v)
  {
    to << "{";
    for(int i = 0; i < 3; ++i) {
      if(i)
	to << ", ";
      to << v[i];
    }
    to << "}";

    return to;
  }

  /********************** 3-D rotational matrix *************************/

  // C style matrix
  class Matrix 
  {
    double _data [9];

  public:
    Matrix () {}
    Matrix (Vector, Vector) ; // standard orientation

    double& operator() (int, int)       ;   // C style indexing
    double  operator() (int, int) const ;   

    Slice<double>      column (int i)       ;
    ConstSlice<double> column (int i) const ;
    Slice<double>      row    (int i)       ;
    ConstSlice<double> row    (int i) const ;

    Matrix operator* (const Matrix&) const; // matrix multiplication
    Vector operator* (const double*) const; // matrix vector product M*v
    Vector operator* (const Vector& v) const { return *this * (const double*)v; }

    void orthogonality_check () const ;
  };

  class Plane
  {
    Vector _normal;
    double     _dist;
    Vector _orth [2];
  public:
	
    void set (const Vector& n, double d);

    Plane ();
    Plane(const Vector& n, double d) { set(n, d); }

    const Vector& normal ()   const { return _normal; }
    double dist ()            const { return _dist; }
    const Vector& orth(int i) const { return _orth[i]; }
  };

  inline std::ostream& operator<< (std::ostream& to, const Plane& p)
  {
    to << "{normal = " << p.normal() << ", offset = " << p.dist() << "}";
    return to;
  }
   

  // reference frame
  struct Frame
  {
    Vector origin; // displacement
    Matrix orient; // orientation matrix

    void set(const Vector& v1, const Vector& v2, const Vector& v3) 
    { origin = v1; orient = Matrix(v2 - v1, v3 - v1); }

    Frame () {}
    Frame (const Vector& v1, const Vector& v2, const Vector& v3) { set(v1, v2, v3); }

    void fv2lv (const double*, double*) const; // frame vector to lab vector
    void lv2fv (const double*, double*) const; // lab vector to frame vector
    void fp2lp (const double*, double*) const; // frame position to lab position
    void lp2fp (const double*, double*) const; // lab position to frame position
  };

} // D3 namespace

double vdot      (const D3::Vector&, const D3::Vector&);
double vdistance (const D3::Vector&, const D3::Vector&);

D3::Vector   vprod (const double*, const double*);
void         vprod (const double*, const double*, double*);

void vprod (const D3::Matrix&, const double*, double*);
void vprod (const double*, const D3::Matrix&, double*);
void find_orth(const double*, D3::Vector []); // get two vectors orthogonal to the given one

double vol (const D3::Vector&, const D3::Vector&, const D3::Vector&);

#endif
