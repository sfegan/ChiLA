//-*-mode:c++; mode:font-lock;-*-

/*! \file VSALinearLeastSquares.cpp

  Generalized Linear Least Squares Fitting

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       12/01/2007
*/

#include <vector>
#include <VSAAlgebra.hpp>
#include <VSALinearLeastSquares.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAMath;
using namespace VERITAS::VSAAlgebra;

// ----------------------------------------------------------------------------
// PolyFit
// ----------------------------------------------------------------------------
double PolyFit::fit(unsigned n, const Data<double>& data, VecND& param,
		      MatrixND* cov, const PreferredFitter::Options& opt)
{
  Fn fn(n);
  PreferredFitter fit(data, fn, opt);
  fit.fit();
  param = fit.param();
  if(cov)*cov = fit.cov();
  return fit.chi2();
}

double PolyFit::val(const VecND& a, double x)
{
#if 0
  Fn fn(a.ndim());
  VecND poly(a.ndim());
  fn(x, poly);
  return a * poly;
#else
  const unsigned na = a.ndim();  
  int ia = na; --ia;
  double y = a[ia];
  while(--ia >= 0) y = x*y+a[ia];
  return y;
#endif
}

double PolyFit::var(const VecND& a, const MatrixND& cov, double x)
{
  VSAAlgebra::VecND dyda;
  VSAMath::PolyFit::dyda(a,x,dyda);

  return dyda*(cov*dyda);
}

void PolyFit::val(const VecND& a, const std::vector<double> x,
			std::vector<double>& y)
{
  unsigned nx = x.size();
  y.resize(nx);
#if 0
  Fn fn(a.ndim());
  VecND poly(a.ndim());
  for(unsigned ix=0;ix<nx;ix++) { fn(x[ix], poly); y[ix] = a * poly; }
#else
  const unsigned na = a.ndim();  
  for(unsigned ix=0;ix<nx;ix++) 
    {
      const double xx = x[ix];
      int ia = na; --ia;
      double yy = a[ia];
      while(--ia >= 0) yy = xx*yy+a[ia];
      y[ix] = yy;
    }
#endif
}

void PolyFit::differentiate(const VecND& a, VecND& dydx_a)
{
  unsigned nd = a.ndim();
  if(nd==0)dydx_a = a;
  else if(nd==1)dydx_a = VecND(1,0.0);
  else
    {
      dydx_a.resize(nd-1);
      for(unsigned id=1;id<nd;id++)dydx_a[id-1] = double(id)*a[id];
    }
}

void PolyFit::dyda(const VecND& a, double x, VecND& dyda)
{
  const unsigned nd = a.ndim();
  dyda.resize(nd);

  double s = dyda[0] = 1.0;
  for(unsigned i=1;i<nd;i++)dyda[i] = s *= x;
}

#ifdef TEST_MAIN

// g++ -march=core2 -O3 -g -mssse3 -mfpmath=sse -Wall -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DUSEALLOCA -fPIC -I/usr/include/mysql -I../VSUtility -I../VSSimDB -I../Physics -I../VSShower -I../VSOptics -I../VSAnalysis -I../VSNSpace -I../VSDataReduction -I../SEphem -I. -I/home/sfegan/include -I/usr/local/veritas/include -DTEST_MAIN -o test VSALinearLeastSquares.cpp VSASVD.cpp ../Physics/RandomNumbers_TOS.o VSAMath.o VSAAlgebra.o -lpthread

#include<iostream>
#include<iomanip>
#include<VSAMath.hpp>
#include<RandomNumbers.hpp>

int main(int argc, char**argv)
{
  RandomNumbers rng("test.dat");

  unsigned ndata = 50;
  std::vector<DataPoint> data(ndata);
  for(unsigned idata=0;idata<ndata;idata++)
    {
      double x = double(idata)/10.0;
      double s = 0.4;
      double y = (((0.5*x)-2)*x+1)*x-7 + rng.Normal()*s;
      data[idata] = DataPoint(x,y,s);
    }

  try
    {
      unsigned nfit = 3;
      VecND param;
      MatrixND cov;
      double chi2;
     
      // Either choose one of the fitters (SVD or Matrix inverse)
      // explicitly below or use default as defined by static "fit"
      // function in PolyFit
#if 0
      PolyFit::Fn fn(nfit);
      Fitlin<PolyFit::Fn> fit(data, fn);
      //Fitsvd<PolyFit::Fn> fit(data, fn);
      fit.fit();
      param = fit.param();
      cov = fit.cov();
      chi2 = fit.chi2();
#else
      chi2 = PolyFit::fit(nfit,data,param,/* optional */ &cov);
#endif

      unsigned np = param.ndim();

      for(unsigned idata=0;idata<ndata;idata++)
	std::cout << data[idata].x << ' ' << data[idata].y << ' ' 
		  << data[idata].sigma << ' '
		  << PolyFit::val(param,data[idata].x) << '\n';
      std::cout << '\n';


      std::cout << chi2 << ' ' << chi2/double(100-np) << "\n\n";
      for(unsigned ip=0; ip<np; ip++)
	std::cout << param[ip] << " +/- " << sqrt(cov(ip,ip)) << '\n';
      
    }
  catch(const std::string& s)
    {
      std::cerr << s << '\n';
    }
}

#endif
