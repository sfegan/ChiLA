//-*-mode:c++; mode:font-lock;-*-

/*! \file VSANewtonRaphsonRootFinder.hpp
  Find a root of a function of using the Newton-Raphson method + Bisection.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    1.0
  \date       03/08/2010
*/

#ifndef VSANEWTONRAPHSONROOTFINDER_HPP
#define VSANEWTONRAPHSONROOTFINDER_HPP

#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <stdexcept>

#ifdef VSA_ASSERT
#include<cassert>
#endif

#include <VSAMinimizationWithDerivative.hpp>

namespace VERITAS
{

  namespace VSAMath
  {
    class RootNotBracketed
    {
    public:
    };

    // template<typename Fn>
    // inline double findNewtonRaphsonRoot(Fn& fn, double lo_x, double hi_x, 
    // 					double eps = 1E-6)
    // {
    //   return findRoot(fn,lo_x,hi_x,0.0,eps);
    // }

    template<typename Fn>
    inline double findNewtonRaphsonRoot(Fn& fn, 
					double lo_x, double hi_x, 
					double fn_val, double eps = 1E-6)
    {
      const unsigned int MAXIT = 100;
      const double xacc = eps*(hi_x+lo_x)*0.5;
	
      double flo = fn.val(lo_x) - fn_val;
      double fhi = fn.val(hi_x) - fn_val;
      
      if((flo > 0.0 && fhi > 0.0) || (flo < 0.0 && fhi < 0.0))
	{
	  std::ostringstream os;
	  os << "Root not bracketed. " << std::endl
	     << std::setprecision(16)
	     << "lo_x   = " << std::setw(25) << lo_x << std::endl
	     << "hi_x   = " << std::setw(25) << hi_x << std::endl
	     << "flo    = " << std::setw(25) << flo << std::endl
	     << "fhi    = " << std::setw(25) << fhi << std::endl
	     << "fn_val = " << std::setw(25) << fn_val << std::endl;

	  throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ ": " +
				   os.str());

	}
      //	throw RootNotBracketed();

      if(flo == 0.0) return lo_x;
      if(fhi == 0.0) return hi_x;

      double xl, xh;

      if(flo < 0.0)
	{
	  xl = lo_x;
	  xh = hi_x;
	}
      else
	{
	  xh = lo_x;
	  xl = hi_x;
	}

      double rts = 0.5*(lo_x + hi_x);
      double dxold = fabs(hi_x - lo_x);
      double dx = dxold;

      double f, df, temp;

      f = fn.val(rts) - fn_val;
      df = fn.dfdx(rts);
      for(unsigned j = 0; j < MAXIT; j++)
	{
	  // Bisection
	  if(!std::isfinite(df) ||  
	     (((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) ||
	     (fabs(2.0*f) > fabs(dxold*df)) )
	    {
	      dxold=dx;
	      dx = 0.5*(xh-xl);
	      rts=xl+dx;
	      if(xl == rts) return rts;
	    }
	  // Newton step
	  else
	    {
	      dxold = dx;
	      dx = f/df;
	      temp = rts;
	      rts -= dx;
	      if(temp == rts) return rts;
	    }

	  if(fabs(dx) < xacc) return rts;

	  f = fn.val(rts) - fn_val;
	  df = fn.dfdx(rts);

	  if(f < 0.0) xl = rts;
	  else xh = rts;
	}

      std::ostringstream os;
      os << ": No convergence after " << MAXIT << " iterations." << std::endl;
      throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ ": " +
			       os.str());

      return 0.0;
    }


    template<typename Fn>
    inline double findBrentRoot(Fn& fn, 
				double lo_x, double hi_x, 
				double fn_val, double eps = 1E-6)
    {
      const int MAXIT = 100;
      const double tol = eps*fabs(hi_x+lo_x)*0.5;

      const double EPS = std::numeric_limits<double>::epsilon();

      int iter;
      double a = lo_x;
      double b = hi_x;
      double c = hi_x;
      double d = 0, e = 0, min1, min2;
      double fa = fn.val(a) - fn_val;
      double fb = fn.val(b) - fn_val;
      double fc,p,q,r,s,tol1,xm;

      if((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
	throw RootNotBracketed();
      fc=fb;

      for(iter=0;iter<MAXIT;iter++)
	{
	  if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
	    {
	      c=a;
	      fc=fa;
	      e=d=b-a;
	    }

	  if(fabs(fc) < fabs(fb))
	    {
	      a=b;
	      b=c;
	      c=a;
	      fa=fb;
	      fb=fc;
	      fc=fa;
	    }
	 
	  tol1=2.0*EPS*fabs(b)+0.5*tol;
	  xm=0.5*(c-b);
	  if(fabs(xm) <= tol1 || fb == 0.0) return b;
	  if(fabs(e) >= tol1 && fabs(fa) > fabs(fb))
	    {
	      s=fb/fa;
	      if(a == c)
		{
		  p=2.0*xm*s;
		  q=1.0-s;
		}
	      else
		{
		  q=fa/fc;
		  r=fb/fc;
		  p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
		  q=(q-1.0)*(r-1.0)*(s-1.0);	     
		}
	      if(p > 0.0) q = -q;
	      p=fabs(p);
	      min1=3.0*xm*q-fabs(tol1*q);
	      min2=fabs(e*q);
	      if(2.0*p < (min1 < min2 ? min1 : min2))
		{
		  e=d;
		  d=p/q;
		}
	      else
		{
		  d=xm;
		  e=d;
		}
	    }
	  else
	    {
	      d=xm;
	      e=d;
	    }
	  a=b;
	  fa=fb;
	  if(fabs(d) > tol1) b += d;
	  else b += SIGN(tol1,xm);
	  fb=fn.val(b) - fn_val;
	}

      std::ostringstream os;
      os << ": No convergence after " << MAXIT << " iterations." << std::endl;
      throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ ": " +
			       os.str());

      return 0;
    }

  }
}  

#endif // VSANEWTONRAPHSONROOTFINDER
