//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAMinimization.hpp
  Minimization general

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/06/2007
*/

#ifndef VSAMINIMIZATION_HPP
#define VSAMINIMIZATION_HPP

#include <VSAAlgebra.hpp>

namespace VERITAS 
{
  namespace VSAMath
  {
    class MinimizationFailed
    {
    public:
      MinimizationFailed() { }
    };

    inline void SHFT(double& a, double& b, double& c, const double& d)
    {
      a = b; b = c; c = d;
    }

    inline double SIGN(const double& a, const double& b) 
    {
      return ((b) >= 0.0 ? fabs(a) : -fabs(a));
    }

    template<typename Fn>
    void mnbrak(double& ax, double& bx, double& cx, 
		double& fa, double& fb, double& fc,
		Fn& fn)
    {
      static const double GOLD    = 1.618034;
      static const double GLIMIT  = 100.0;
      static const double TINY    = 1.0e-20;
    
      fa = fn.val(ax);
      fb = fn.val(bx);

      if (fb > fa) 
	{
	  double dummy;
	  dummy = ax; ax = bx; bx = dummy;
	  dummy = fa; fa = fb; fb = dummy;
	}

      cx = bx + GOLD*(bx-ax); 

      fc = fn.val(cx);

#if 0
      std::cerr << "MNBRAK_START: " << ax << ' ' << bx << ' ' << cx
		<< ' ' << fa+1 << ' ' << fb+1 << ' ' << fc+1 << std::endl;
#endif

      double fu;
      while (fb > fc) 
	{ 
	  double r = (bx-ax)*(fb-fc);
	  double q = (bx-cx)*(fb-fa);
	  double u = (bx-((bx-cx)*q-(bx-ax)*r)
		      /(2.0*SIGN(std::max(fabs(q-r),TINY),q-r)));
	  double ulim = bx+GLIMIT*(cx-bx);

#if 0
	  std::cerr << "MNBRAK_ITER: " << ax << ' ' << bx << ' ' << cx
		    << ' ' << u << std::endl;
#endif

	  if ((bx-u)*(u-cx) > 0.0) 
	    { 
#if 0
	      std::cerr << "MNBRAK_A " << std::endl;
#endif
	      fu = fn.val(u);

	      if (fu < fc) 
		{
		  ax = bx;
		  bx = u;
		  fa = fb;
		  fb = fu;
		  return;
		} 
	      else if (fu > fb) 
		{
		  cx = u;
		  fc = fu;
		  return;
		}

	      u = cx+GOLD*(cx-bx); 
	      fu = fn.val(u);
	    }
	  else if ((cx-u)*(u-ulim) > 0.0) 
	    { 
	      fu = fn.val(u);
#if 0
	      std::cerr << "MNBRAK_B " << fu+1 << std::endl;
#endif
	      if (fu < fc) 
		{
		  bx = cx;
		  cx = u;
		  u  = cx+GOLD*(cx-bx);
		  fb = fc;
		  fc = fu;
		  fu = fn.val(u);
#if 0
		  std::cerr << "MNBRAK_B2 " << u << ' ' << fu+1 << std::endl;
#endif
		}
	    } 
	  else if ((u-ulim)*(ulim-cx) >= 0.0) 
	    {
#if 0
	      std::cerr << "MNBRAK_C " << std::endl;
#endif
	      u=ulim; 
	      fu = fn.val(u);
	    } 
	  else 
	    {
#if 0
	      std::cerr << "MNBRAK_D " << std::endl;
#endif
	      u=cx+GOLD*(cx-bx); 
	      fu = fn.val(u);
	    }
	  SHFT(ax,bx,cx,u);
	  SHFT(fa,fb,fc,fu);
	}
    }

    template<typename Fn>
    double brent(double ax, double bx, double cx, 
		 double& xmin, const double& tol,
		 Fn& fn)
    {
      static const unsigned ITMAX=100;
      static const double CGOLD = 0.3819660;
      static const double ZEPS = std::numeric_limits<double>::epsilon()*1.0E-3;
      
      double a, b, d = 0.0, etemp;
      double fu,fw,fx,fv;
      double p,q,r,tol1,tol2,u,v,w,x,xm;
      double e = 0.0;

      a = (ax < cx ? ax : cx);
      b = (ax > cx ? ax : cx);
      x = w = v = bx;
      fw = fv = fx = fn.val(x);
      
      for(int iter = 0; iter < (int)ITMAX; iter++)
	{
	  xm=0.5*(a+b);
	  tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
	  if(fabs(x-xm) <= (tol2-0.5*(b-a)))
	    {
	      xmin = x;
	      return fx;
	    }

	  if(fabs(e) > tol1)
	    {
	      r=(x-w)*(fx-fv);
	      q=(x-v)*(fx-fw);
	      p=(x-v)*q-(x-w)*r;
	      q=2.0*(q-r);
	      if(q > 0.0) p = -p;
	      q=fabs(q);
	      etemp=e;
	      e=d;

	      if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
		d=CGOLD*(e=(x >= xm ? a-x : b-x));
	      else
		{
		  d=p/q;
		  u=x+d;
		  if(u-a < tol2 || b-u < tol2)
		    d = SIGN(tol1,xm-x);
		}
	    }
	  else
	    d=CGOLD*(e=(x >= xm ? a-x : b-x));

	  u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
	  fu=fn.val(u);
	  if(fu <= fx)
	    {
	      if(u >= x) a = x; else b = x;
	      SHFT(v,w,x,u);
	      SHFT(fv,fw,fx,fu);
	    }
	  else
	    {
	      if(u < x) a = u; else b = u;
	      if(fu <= fw || w == x)
		{
		  v=w;
		  w=u;
		  fv=fw;
		  fw=fu;
		}
	      else if(fu <= fv || v == x || v == w)
		{
		  v=u;
		  fv=fu;
		}
	    }	      
	}
    

      std::ostringstream os;
      os << ": No convergence after " << ITMAX << " iterations." << std::endl;
      throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ ": " +
			       os.str());

      return 0;
    }

    template<typename Fn>
    class Fn1Dim
    {
    public:
      
      Fn1Dim(const Fn& fn,
	     const VSAAlgebra::VecND& p,
	     const VSAAlgebra::VecND& xi): 
	m_fn(fn), m_p(p), m_xi(xi) { }
      

      double val(double x)
      {
	VSAAlgebra::VecND xt = m_p+x*m_xi;
	return m_fn.val(xt);
      }

    private:

      const Fn&         m_fn;
      VSAAlgebra::VecND m_p;
      VSAAlgebra::VecND m_xi;

    };

    template<typename Fn>
    void linmin(VSAAlgebra::VecND& p, VSAAlgebra::VecND& xi, double& fret,
		Fn& fn)
    {
      static const double TOL=1.0E-8;

      Fn1Dim<Fn> fn1d(fn,p,xi);

      //       double xmin, fx, fb, fa, bx;

      //       int n = p.ndim();

      double ax=0.0;
      double xx=1.0;

      double bx;
      double fa;
      double fx;
      double fb;


      mnbrak(ax,xx,bx,fa,fx,fb,fn1d);
      double xmin;
      fret=brent(ax,xx,bx,xmin,TOL,fn1d);
      xi *= xmin;
      p += xi;       
    }

    template<typename Fn>
    void powell(VSAAlgebra::VecND& p, VSAAlgebra::MatrixND& xi,
		const double ftol, int& iter, double& fret,
		Fn& fn)
    {
      static const int ITMAX = 200;
      static const double TINY = 1.0e-20;
      
      int ibig;
      double del, fp, fptt, t;

      int n = p.ndim();

      VSAAlgebra::VecND ptt(n), xit(n);
      fret = fn.val(p);
      VSAAlgebra::VecND pt = p;
      
      for(iter = 0;;++iter)
	{
	  fp=fret;
	  ibig=0;
	  del=0.0; // Will be the biggest function decrease

	  for(int i=0;i<n;i++)
	    {
	      for(int j=0;j<n;j++) xit[j]=xi[j][i];
	      fptt=fret;
	      linmin(p,xit,fret,fn);
	      if(fptt-fret > del)
		{
		  del=fptt-fret;
		  ibig=i+1;
		}
	    }

	  if(2.0*(fp-fret) <= ftol*(fabs(fp)+fabs(fret))+TINY) return;


	  if(iter == ITMAX)
	    {
	      std::ostringstream os;
	      os << ": No convergence after " << ITMAX 
		 << " iterations." << std::endl;
	      throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ ": " +
				       os.str());

	    }

	  for(int j=0;j<n;j++)
	    {
	      ptt[j]=2.0*p[j]-pt[j];
	      xit[j]=p[j]-pt[j];
	      pt[j]=p[j];
	    }

	  fptt=fn.val(pt);
	  if(fptt < fp)
	    {
	      t = 2.0*(fp-2.0*fret+fptt)*std::pow(fp-fret-del,2) -
		del*std::pow(fp-fptt,2);
	      if(t < 0.0)
		{
		  linmin(p,xit,fret,fn);
		  for(int j=0;j<n;j++)
		    {
		      xi[j][ibig-1]=xi[j][n-1];
		      xi[j][n-1]=xit[j];
		    }
		}

	    }	  
	}
    }
  }
}

#endif // define VSAMINIMIZATION_HPP
