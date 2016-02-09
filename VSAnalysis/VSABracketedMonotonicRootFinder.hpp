//-*-mode:c++; mode:font-lock;-*-

/*! \file VSABracketedMonotonicRootFinder.hpp
  Find a root of a monotonic function which has been bracketed by two values

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#ifndef VSABRACKETEDMONOTONICROOTFINDER_HPP
#define VSABRACKETEDMONOTONICROOTFINDER_HPP

#include<cmath>

#ifdef VSA_ASSERT
#include<cassert>
#endif

namespace VERITAS
{

  namespace VSAMath
  {

    class FunctionNotMonotonic
    {
    public:
    };

    template<typename Function>
    inline double findMonotonicHiBracket(Function& function, double x, 
					 double dx = 1.0)
    {
      double f0 = function(x);
      double f;
      if(f0>0)
	do { x += dx; dx *= 2.0; f = function(x); 
	  if(f>f0)throw FunctionNotMonotonic(); f0=f; } while(f0>0);
      else
	do { x += dx; dx *= 2.0; f = function(x); 
	  if(f<f0)throw FunctionNotMonotonic(); f0=f; } while(f0>0);
      return x;
    }

    template<typename Function>
    inline double findMonotonicLoBracket(Function& function, double x, 
					 double dx = 1.0)
    {
      double f0 = function(x);
      double f;
      if(f0>0)
	do { x -= dx; dx *= 2.0; f = function(x); 
	  if(f>f0)throw FunctionNotMonotonic(); f0=f; } while(f0>0);
      else
	do { x -= dx; dx *= 2.0; f = function(x); 
	  if(f<f0)throw FunctionNotMonotonic(); f0=f; } while(f0>0);
      return x;
    }

    template<typename Function, typename T>
    inline double findBracketedMonotonicRoot(Function& function,
					     T lo_x, T hi_x,
					     unsigned idim,
					     double val,
					     const double& tolerence)
    {
      double lo_f = function(lo_x)-val;
      double hi_f = function(hi_x)-val;
      while(fabs(lo_f-hi_f)>tolerence)
	{
	  T mid_x = lo_x;
	  mid_x[idim] = 0.5*(lo_x[idim]+hi_x[idim]);
	  const double mid_f = function(mid_x)-val;
	  if(lo_f<hi_f)
	    {
	      if((mid_f<=lo_f)||(mid_f>=hi_f))throw FunctionNotMonotonic();
	      if(mid_f<0)lo_x = mid_x, lo_f = mid_f;
	      else hi_x = mid_x, hi_f = mid_f;
	    }
	  else
	    {
	      if((mid_f>=lo_f)||(mid_f<=hi_f))throw FunctionNotMonotonic();
	      if(mid_f<0)hi_x = mid_x, hi_f = mid_f;
	      else lo_x = mid_x, lo_f = mid_f;
	    }
	}
      return 0.5*(lo_x[idim]+hi_x[idim]);
    }

    template<typename Function>
    inline double findBracketedMonotonicRoot(Function& function,
					     double lo_x, double hi_x,
					     const double& tolerence)
    {
      double lo_f = function(lo_x);
      double hi_f = function(hi_x);
      while(fabs(lo_f-hi_f)>tolerence)
	{
	  const double mid_x = 0.5*(lo_x+hi_x);
	  const double mid_f = function(mid_x);
	  if(lo_f<hi_f)
	    {
	      if((mid_f<=lo_f)||(mid_f>=hi_f))throw FunctionNotMonotonic();
	      if(mid_f<0)lo_x = mid_x, lo_f = mid_f;
	      else hi_x = mid_x, hi_f = mid_f;
	    }
	  else
	    {
	      if((mid_f>=lo_f)||(mid_f<=hi_f))throw FunctionNotMonotonic();
	      if(mid_f<0)hi_x = mid_x, hi_f = mid_f;
	      else lo_x = mid_x, lo_f = mid_f;
	    }
	}
      return 0.5*(lo_x+hi_x);
    }

    template<typename Function>
    inline double findBracketedMonotonicRootXTOL(Function& function,
						 double lo_x, double hi_x,
						 const double& tolerence)
    {
      double lo_f = function(lo_x);
      double hi_f = function(hi_x);
      while(fabs(lo_x-hi_x)>tolerence)
	{
	  const double mid_x = 0.5*(lo_x+hi_x);
	  const double mid_f = function(mid_x);
	  if(lo_f<hi_f)
	    {
	      if((mid_f<=lo_f)||(mid_f>=hi_f))throw FunctionNotMonotonic();
	      if(mid_f<0)lo_x = mid_x, lo_f = mid_f;
	      else hi_x = mid_x, hi_f = mid_f;
	    }
	  else
	    {
	      if((mid_f>=lo_f)||(mid_f<=hi_f))throw FunctionNotMonotonic();
	      if(mid_f<0)hi_x = mid_x, hi_f = mid_f;
	      else lo_x = mid_x, lo_f = mid_f;
	    }
	}
      return 0.5*(lo_x+hi_x);
    }
    
  }

}

#endif // VSABRACKETEDMONOTONICROOTFINDER
