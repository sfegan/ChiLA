//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAConfidenceInterval.hpp

  Routines for calculating confidence intervals.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       01/20/2009
*/

#include <VSANewtonRaphsonRootFinder.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSAMath.hpp>

namespace VERITAS
{
  namespace VSAMath
  {

    template< typename Fn >
    class ConfidenceInterval
    {
    public:

      class FnA
      {
      public:
	FnA(NLFitter<Fn>* fitter, const VSAAlgebra::VecND& a, unsigned ip): 
	  m_fitter(fitter->clone()), m_a(a), m_ip(ip) 
	{ 
	  m_fitter->free();
	  m_fitter->hold(ip);
	}

	double val(double x) 
	{
	  m_a[m_ip] = x;
	  m_fitter->fit(m_a);

	  double chi2 = m_fitter->chi2();

	  // std::cout << std::setw(25) << x 
	  //		    << std::setw(25) << chi2 << std::endl;

	  return chi2;
	}


	NLFitter<Fn>* m_fitter;
	VSAAlgebra::VecND m_a;
	unsigned          m_ip;
      };

      class FnB
      {
      public:
	FnB(NLFitter<Fn>* fitter, const VSAAlgebra::VecND& a, 
	    unsigned ip1, unsigned ip2): 
	  m_fitter(fitter), m_a(a), m_ip1(ip1), m_ip2(ip2)
	{ 
	  m_fitter->free();
	  m_fitter->hold(ip1);
	  m_fitter->hold(ip2);
	}

	double val(double x) 
	{
	  m_a[m_ip2] = x;
	  m_fitter->fit(m_a);

	  double chi2 = m_fitter->chi2();

	  // std::cout << std::setw(25) << x 
	  //		    << std::setw(25) << chi2 << std::endl;

	  return chi2;
	}

      private:

	NLFitter<Fn>*     m_fitter;
	VSAAlgebra::VecND m_a;
	unsigned          m_ip1;
	unsigned          m_ip2;
      };

      ConfidenceInterval(NLFitter<Fn>* fitter, const Fn& fn): 
	m_fitter(fitter->clone()), m_fn(fn) 
      { 
	//	m_fitter = new NLFitterPowell<Fn>(fn);
      }

      ConfidenceInterval(const Fn& fn): 
	m_fn(fn) 
      { 
	m_fitter = new NLFitterPowell<Fn>(fn);
      }

      void limits(double prob, 
		  const VSAAlgebra::VecND& param,
		  unsigned ip, double dx,
		  double& lo, double &hi);
      
      // void pdf(const VSAAlgebra::VecND& param, unsigned ip)
      // {
	
      // }

      void limits(double prob, 
      		  const VSAAlgebra::VecND& param,
      		  unsigned ip1, double dx1,
      		  unsigned ip2, double dx2,
      		  std::vector< std::pair<double,double> >& xy);


    private:

      NLFitter<Fn>*  m_fitter;
      Fn             m_fn;

    };

    template< typename Fn >
    void ConfidenceInterval<Fn>::limits(double prob, 
					const VSAAlgebra::VecND& param,
					unsigned ip, double dx,
					double& lo, double &hi)
    {
      FnA fn(m_fitter,param,ip);
      
      double fval = m_fn.val(param);
      double flo = 0, fhi = 0;

      double xm = param(ip);

      double dxhi = dx;
      double hix = xm+dxhi;

      double df = VSAMath::chi2q(prob,1);


      while(fn.val(hix)-fval < df)
	{
	  dxhi *= 2;
	  hix = xm+dxhi;
	}


      hi = findBrentRoot(fn,xm,hix,fval+df,1E-3);

      //      std::cout << "LO " << lo << std::endl;

      //      hix = lox - dx;
      double dxlo = dx;
      double lox = xm-dxlo;

      // std::cout << fn.val(lox) << " " << lox << std::endl;
      // std::cout << fn.val(hix) << " " << hix << std::endl;

      while(fn.val(lox)-fval < df)
	{
	  dxlo *= 2;
	  lox = xm - dxlo;
	}
      // std::cout << "FNLO " << dx << " " << fn.val(lox) << " " << lox << " "
      // 		<< fval+df << std::endl;

      lo = findBrentRoot(fn,lox,xm,fval+df,1E-3);
    }

    template< typename Fn >
    void ConfidenceInterval<Fn>::
    limits(double prob, 
	   const VSAAlgebra::VecND& param,
	   unsigned ip1, double dx1,
	   unsigned ip2, double dx2,
	   std::vector< std::pair<double,double> >& xy)
    {
      FnA fna1(m_fitter,param,ip1);
      FnA fna2(m_fitter,param,ip2);
      
      double fval = m_fn.val(param);
      double df = VSAMath::chi2q(prob,2);

      double x1m = param(ip1);
      double dx1lo = dx1;
      double x1lo = x1m-dx1lo;

      while(fna1.val(x1lo)-fval < df)
	{
	  dx1lo *=1.25;
	  x1lo = x1m-dx1lo;
	}

      double dx1hi = dx1;
      double x1hi = x1m + dx1hi;
      while(fna1.val(x1hi)-fval < df)
	{
	  dx1hi *=1.25;
	  x1hi = x1m+dx1hi;
	}
      
      double x2m = param(ip2);
      double dx2lo = dx2;
      double x2lo = x2m-dx2lo;

      //     std::cout << "X2M " << x2m << std::endl;
      
      while(fna2.val(x2lo)-fval < df)
	{
	  dx2lo *=1.25;
	  x2lo = x2m-dx2lo;
	}

      fna2.m_a = param;

      double dx2hi = dx2;
      double x2hi = x2m + dx2hi;
      while(fna2.val(x2hi)-fval < df)
	{
	  dx2hi *=1.25;
	  x2hi = x2m+dx2hi;
	}

      //      x2hi = x2m + 1.1*dx2;

      // std::cout << x1lo << " " << x1hi << std::endl;
      // std::cout << x2lo << " " << x2hi << std::endl;

      std::vector< std::pair<double,double> > xy1;
      std::vector< std::pair<double,double> > xy2;

      for(double x1 = x1lo; x1 < x1hi; x1 += (x1hi-x1lo)/200.)
      	{
      	  m_fitter->free();
      	  m_fitter->set(param);
      	  m_fitter->hold(ip1,x1);
      	  m_fitter->fit();

      	  double x2m = m_fitter->param(ip2);

      	  if(m_fitter->chi2() - fval > df) continue;

      	  FnB fn2(m_fitter,m_fitter->param(),ip1,ip2);
	  
	  //std::cout << x1 << " " << x2m << " ";

	  double x2 = findBrentRoot(fn2,x2m,x2hi,fval+df,1E-4);
	  xy1.push_back(std::make_pair(x1,x2));

	  x2 = findBrentRoot(fn2,x2lo,x2m,fval+df,1E-4);
	  xy2.push_back(std::make_pair(x1,x2));
	  //std::cout << x2 << std::endl;
      	}

      std::reverse(xy2.begin(),xy2.end());
      
      xy.insert(xy.end(),xy1.begin(),xy1.end());
      xy.insert(xy.end(),xy2.begin(),xy2.end());
      xy.insert(xy.end(),xy1.begin(),xy1.begin()+1);

    }

  }

};
