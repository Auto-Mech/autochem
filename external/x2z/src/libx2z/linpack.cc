#include "linpack.hh"
#include <iostream>
#include <cmath>

/*********************** Vector operations *******************************/

double normalize (double* v, int dim)
{
    const char funame [] = "normalize (double*, int): ";
    double norm = 0.0;
    double dtemp;
    const double* p = v;

    for (int i = 0; i < dim; ++i) {
	dtemp = *p++;
	norm += dtemp * dtemp;
    }

    if (norm == 0.) {
	std::cerr << funame << "WARNING: the vector length is zero\n";
	return norm;
    }
    norm = std::sqrt(norm);

    for (int i = 0; i < dim; ++i)
    {
	*v++ /= norm;
    }

    return norm;
}

double orthogonalize (double* v, const double* n, int dim) 
{
    const char funame [] = "orthogonalize: ";
    static const double tol = 1.e-12;

    double nlen = 0.0, proj = 0.0, vlen = 0.0;
    for (int i = 0; i < dim; ++i) {
	nlen += n[i] * n[i];
	vlen += v[i] * v[i];
	proj += n[i] * v[i];
    }

    if(nlen == 0.) {
	std::cerr << funame << "the ortogonalizing vector length is zero\n";
	throw Error::Range();
    }

    if(vlen == 0.) {
	std::cerr << funame << "WARNING: the ortogonalized vector length is zero\n";
	return 0.;
    }

    proj /= nlen;
    for (int i = 0; i < dim; ++i) {
	v[i] -= proj * n[i];
    }

    double cos2 = proj * proj * nlen / vlen;

#ifdef DEBUG
    if(cos2 > 1. - tol) {
      std::cerr << funame << "vectors are colinear, sin**2 = " << 1. - cos2 << "\n";
	throw Error::Range();
    }
#endif

    return cos2;
}

double vdistance (const double* v1, const double* v2, int dim)
{
  double dist2 = 0.0;
  double dtemp;
  for (int i = 0; i < dim; ++i)
    {
      dtemp = *v1++ - *v2++;
      dist2 += dtemp * dtemp;
    }
  return std::sqrt(dist2);
}

double vlength (const double* v, int dim)
{
  double norm = 0.0;
  double dtemp;

  for (int i = 0; i < dim; ++i)
    {
      dtemp = *v++;
      norm += dtemp * dtemp;
    }

  return std::sqrt(norm);
}

double vdot (const double* v1, const double* v2, int n)
{
  double res = 0.;
  for (int i = 0; i < n; ++i)
    res += *v1++ * *v2++;
  return res;
}

