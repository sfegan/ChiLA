//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAMinimizationWithDerivative.hpp
  Minimization of function with derivative

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/06/2007
*/

#ifndef VSAMINIMIZATIONWITHDERIVATIVE_HPP
#define VSAMINIMIZATIONWITHDERIVATIVE_HPP

#include <cassert>

#include <VSAMinimization.hpp>

namespace VERITAS 
{
  namespace VSAMath
  {
    inline void MOV3(double& a, double& b, double& c,
		     const double& d, const double& e, const double& f)
    {
      a=d; b=e; c=f;
    }

    template<typename FunctionAndDerivative>
    double dbrent(double ax, double bx, double cx, 
		  double& xmin, const double& tol,
		  FunctionAndDerivative& function, 
		  bool find_minimum = true)
    {
      static const unsigned ITERMAX = 100;
      static const double   ZEPS    = 1.0e-10;

      //int iter,ok1,ok2; Will be used as flags for whether profloat
      //a,b,d,d1,d2,du,dv,dw,dx,e=0.0; posed steps are acceptable or not.
      //float fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
  
      double a = (ax < cx ? ax : cx);
      double b = (ax > cx ? ax : cx);

      double x(bx);
      double w(bx);
      double v(bx);

      double fx(0.0);
      double dx(0.0); 

      function.val(x,fx,dx);
      if(!find_minimum) { fx *= -1; dx *= -1; }

      double fw(fx);
      double fv(fx);

      double dw(dx);
      double dv(dx);

      double e(0.0);
      double d(0.0);
      double u(0.0);
      double fu(0.0);
      double du(0.0);

#if 0
      std::cerr << "DBRNT " << ax << ' ' << bx << ' ' << cx << ' ';
#endif

      for(unsigned iter=0;iter<=ITERMAX;iter++) 
	{
	  double xm = 0.5*(a+b);
	  double tol1 = tol*fabs(x)+ZEPS;
	  double tol2 = 2.0*tol1;
	  if(fabs(x-xm) <= (tol2-0.5*(b-a))) 
	    {
	      xmin=x;
#if 0
	      std::cerr << fx << std::endl;
#endif
	      return fx;
	    }

	  if(fabs(e) > tol1) 
	    {
	      double d1=2.0*(b-a);
	      double d2=d1;
	      if (dw != dx) d1=(w-x)*dx/(dx-dw);
	      if (dv != dx) d2=(v-x)*dx/(dx-dv);

	      double u1 = x+d1;
	      double u2 = x+d2;

	      bool ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
	      bool ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;

	      double olde = e;
	      e = d;

	      if (ok1 || ok2) 
		{
		  if (ok1 && ok2)
		    d = (fabs(d1) < fabs(d2) ? d1 : d2);
		  else if (ok1)
		    d = d1;
		  else
		    d = d2;
		  if (fabs(d) <= fabs(0.5*olde))
		    {
		      u = x+d;
		      if (u-a < tol2 || b-u < tol2)
			d = SIGN(tol1,xm-x);
		    } 
		  else 
		    {
		      e = (dx >= 0.0 ? a-x : b-x);
		      d = 0.5*e;
		    }
		} 
	      else
		{
		  e = (dx >= 0.0 ? a-x : b-x);
		  d = 0.5*e;
		}
	    } 
	  else
	    {
	      e = (dx >= 0.0 ? a-x : b-x);
	      d = 0.5*e;
	    }
	  if (fabs(d) >= tol1)
	    {
	      u = x+d;
	      fu = function.val(u);
	      if(!find_minimum) { fu *= -1; }
	    } 
	  else
	    {
	      u = x+SIGN(tol1,d);
	      fu = function.val(u);
	      if(!find_minimum) { fu *= -1; }
	      if (fu > fx) 
		{
		  xmin=x;
#if 0
		  std::cerr << fx << std::endl;
#endif
		  return fx;
		}
	    }

	  du = function.dfdx(u);
	  if(!find_minimum) { du *= -1; }
	  if (fu <= fx) 
	    {
	      if (u >= x)a = x; 
	      else b = x;
	      MOV3(v,fv,dv, w,fw,dw);
	      MOV3(w,fw,dw, x,fx,dx);
	      MOV3(x,fx,dx, u,fu,du);
	    } 
	  else 
	    {
	      if (u < x)a = u; 
	      else b = u;

	      if (fu <= fw || w == x) 
		{
		  MOV3(v,fv,dv, w,fw,dw);
		  MOV3(w,fw,dw, u,fu,du);
		} 
	      else if (fu < fv || v == x || v == w)
		{
		  MOV3(v,fv,dv, u,fu,du);
		}
	    }
	}
      throw MinimizationFailed();
    }

    template<typename FunctionAndDerivative2>
    class FunctionAndDerivative2Projection
    {
    public:
      FunctionAndDerivative2Projection(const double& _px, const double& _py,
				       const double& _xix, const double& _xiy,
				       FunctionAndDerivative2& _function):
	px(_px), py(_py), xix(_xix), xiy(_xiy), function(_function) { }

      double val(const double& x)
      {
	setX(x);
	double f;
	eval(f);
	return f;
      }

      double dfdx(const double& x)
      {
	setX(x);
	double _dfdx;
	deriv(_dfdx);
	return _dfdx;
      }

      void val(const double& x, double& _val, double& _dfdx)
      {
	setX(x);
	eval(_val);
	deriv(_dfdx);
      }

      void setX(const double& x) { function.setXY(px + x*xix, py + x*xiy); }
      void eval(double& f) { function.eval(f); }
      void deriv(double& df) 
      { double dfx,dfy; function.deriv(dfx,dfy); df = dfx*xix + dfy*xiy; }
    private:
      const double& px;
      const double& py;
      const double& xix;
      const double& xiy;
      FunctionAndDerivative2& function;
    };

    template<typename FuncionAndDerivative2>
    void linmin(double& px, double& py, double& xix, double& xiy,
		double& fret, FuncionAndDerivative2& function,
		const double DBRENT_TOL = 2e-4)
    {
      double ax=0.0; 
      double xx=1.0;

      double bx;
      double fa;
      double fx;
      double fb;

      FunctionAndDerivative2Projection<FuncionAndDerivative2>
	f1(px, py, xix, xiy, function);

      mnbrak(ax,xx,bx,fa,fx,fb,f1);

      double xmin;
      fret = dbrent(ax,xx,bx,xmin,DBRENT_TOL,f1);

      xix *= xmin;
      xiy *= xmin;
      px += xix;
      py += xiy;
    }

    template<typename FuncionAndDerivative2> 
    void minimizeWithDerivative2(double& _px, double& _py, double& _fmin,
				 const double& ftol,
				 FuncionAndDerivative2& function,
				 const double DBRENT_TOL = 2e-4)
    {
      static const unsigned MAXITER = 200;

      double px = _px;
      double py = _py;
      double gx;
      double gy;
      double hx;
      double hy;
      double xix;
      double xiy;
      double fp;
    
      function.setXY(px,py);
      function.eval(fp);
      function.deriv(xix,xiy);

      gx = -xix;
      gy = -xiy;
      xix = hx = gx;
      xiy = hy = gy;

      for(unsigned iter=0;iter<MAXITER;iter++)
	{
#if 0
	  std::cerr << px << ' ' << py << ' ' << xix << ' ' << xiy << ' '
		    << fp << std::endl;
#endif

	  double fmin;
	  linmin(px, py, xix, xiy, fmin, function, DBRENT_TOL);

	  if(2.0*fabs(fmin-fp) <= 
	     ftol*(fabs(fmin)+fabs(fp)+std::numeric_limits<double>::epsilon()))
	    {
	      _px = px;
	      _py = py;
	      _fmin = fmin;
	      return;
	    }

	  fp = fmin;
	  function.setXY(px,py);
	  function.deriv(xix,xiy);

	  const double gg  = gx*gx + gy*gy;
	  if(gg == 0.0)
	    {
	      _px = px;
	      _py = py;
	      _fmin = fmin;
	      return;
	    }

	  const double dgg = (xix+gx)*xix + (xiy+gy)*xiy;
	  const double gamma = dgg/gg;

	  gx = -xix;
	  gy = -xiy;
	  xix = hx = gx + gamma*hx;
	  xiy = hy = gy + gamma*hy;
	}
      throw MinimizationFailed();
    }

  } // namespace VSAMath
} // namespace VERITAS

#endif // VSAMINIMIZATIONWITHDERIVATIVE_HPP
