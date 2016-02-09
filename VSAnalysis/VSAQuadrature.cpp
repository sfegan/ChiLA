//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAQuadrature.cpp

  Classes for numerical integration.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       12/01/2008
*/

#include <VSAQuadrature.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAMath;



const double QuadratureRuleGauss7::s_abscissas[] = 
  {
    -0.949107912342759,
    -0.741531185599394,
    -0.405845151377397,
    0,
    0.405845151377397,
    0.741531185599394,
    0.949107912342759
  };

const double QuadratureRuleGauss7::s_weights[] = 
  {
    0.129484966168870,
    0.279705391489277,
    0.381830050505119,
    0.417959183673469,
    0.381830050505119,
    0.279705391489277,
    0.129484966168870
  };

const double QuadratureRuleKronrod15::s_abscissas[] = 
  {
    -0.991455371120813,
    -0.949107912342759,
    -0.864864423359769,
    -0.741531185599394,
    -0.586087235467691,
    -0.405845151377397,
    -0.207784955007898,
    0,
    0.207784955007898,
    0.405845151377397,
    0.586087235467691,
    0.741531185599394,
    0.864864423359769,
    0.949107912342759,
    0.991455371120813 
  };

const double QuadratureRuleKronrod15::s_weights[] = 
  { 
    0.022935322010529,
    0.063092092629979,
    0.104790010322250,
    0.140653259715525,
    0.169004726639267,
    0.190350578064785,
    0.204432940075298,
    0.209482141084728,
    0.204432940075298,
    0.190350578064785,
    0.169004726639267,
    0.140653259715525,
    0.104790010322250,
    0.063092092629979,
    0.022935322010529
  };

#ifdef TEST_MAIN

// g++ -O3 -g -Wall -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DUSEALLOCA -fPIC -I/usr/include/mysql -I../VSUtility -I../VSSimDB -I../Physics -I../VSShower -I../VSOptics -I../VSAnalysis -I../VSNSpace -I../VSDataReduction -I../SEphem -I. -I../VSCommon -I/home/sfegan/include -I/usr/local/veritas/include -DTEST_MAIN -o test VSAQuadrature.cpp ../Physics/RandomNumbers_TOS.o VSAMath.o VSAAlgebra.o -lpthread

#include <iostream>
#include <vector>

#include <VSACoord.hpp>

class Poly
{
public:
  Poly() {}

  unsigned ndim() const { return 1; }

  double val(const VSACoord::Coord1D x)
  {
    return 3.47*x[0]*x[0]-4*x[0];
  }

};

class Poly2
{
public:
  Poly2() {}

  unsigned ndim() const { return 1; }

  double val(const VSACoord::Coord2D x) const
  {
    return x[0]*x[1];
  }

};

class Volume : public QuadratureVolume<VSACoord::Coord2D>
{
public:
  Volume(VSACoord::Coord2D lo, VSACoord::Coord2D hi):
    QuadratureVolume<VSACoord::Coord2D>(lo,hi)
  {

  }

  
  bool isInVolume(const VSACoord::Coord2D& c) const
  {
    //    std::cout << c.y() << " " << c.x() << std::endl;

    if(c.y() > c.x()) return false;
    else return true;
  }
};

int main(int argc, char**argv)
{
  Poly poly;
  Poly2 poly_fn;

  VSACoord::Coord2D lo(0.1,0.1);
  VSACoord::Coord2D hi(10,10);

  Quadrature2D< Poly2, VSACoord::Coord2D, QuadratureRuleGauss7 > qd2d(poly_fn);

  double v2 = qd2d.integrate(lo,hi);
  std::cout << v2 << std::endl;


  QuadratureND<Poly2, VSACoord::Coord2D, QuadratureRuleGauss7> quad7(poly_fn);
  QuadratureND<Poly2, VSACoord::Coord2D, QuadratureRuleKronrod15> 
    quad15(poly_fn);

  QuadratureGaussKronrod<Poly2,VSACoord::Coord2D> gk_quad(poly_fn);



  double a1 = quad7.integrate(lo,hi);
  double a2 = quad15.integrate(lo,hi);
  double a3 = gk_quad.integrate(lo,hi);

  double a4;
  
  a4 = Quadrature::integrate(poly_fn,lo,hi);

  Volume v(lo,hi);

  double a5 = Quadrature::integrate<Poly2,VSACoord::Coord2D,Volume>(poly_fn,v);


  std::cout << "----------------------------------" << std::endl;
  std::cout << a1 << " " << a2 << " " << a3 << " " << a4 << " " << a5
	    << " " << exp(-0.01) - exp(-100) << std::endl;
}

#endif
