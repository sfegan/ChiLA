//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAQuadrature.hpp

  Classes for numerical integration.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       12/01/2008
*/

#ifndef VSAQUADRATURE_HPP
#define VSAQUADRATURE_HPP

#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iostream>

#include <VSACoord.hpp>

namespace VERITAS
{
  namespace VSAMath
  {
    class QuadratureVar1D
    {
    public:
      QuadratureVar1D(): m_x() { }
      QuadratureVar1D(double x): m_x(x) { }
      
      void set(const std::vector<double>& x) { m_x = x[0]; }
      void set(unsigned idim, double x) { m_x = x; }

      double& operator[] (unsigned idim) { return m_x; }
      
      const double& operator[] (unsigned idim) const { return m_x; }


      static unsigned ndim() { return 1; }
      double jacobian() const { return 1; }
      double jacobian(unsigned idim) const { return 1; }    

    protected:

      double m_x;
    };


    // ========================================================================
    // QuadratureVolume
    // ========================================================================
    template< typename T >
    class QuadratureVolume
    {
    public:
      QuadratureVolume(): 
	m_ndim(T::ndim()), m_idim(), m_lo(), m_hi(), s() 
      { 
	for(unsigned idim = 0; idim < m_ndim; idim++) m_idim.push_back(idim);
      }

      QuadratureVolume(T lo, T hi): 
	m_ndim(T::ndim()), m_idim(), m_lo(lo), m_hi(hi), s() 
      { 
	for(unsigned idim = 0; idim < m_ndim; idim++) m_idim.push_back(idim);
      }
      
      QuadratureVolume(T lo, T hi, unsigned idim): 
	m_ndim(1), m_idim(), m_lo(lo), m_hi(hi), s() 
      { 
	m_idim.push_back(idim);
      }

      virtual ~QuadratureVolume() { }

      bool operator <(const QuadratureVolume& rhs) const
      {
	return s < rhs.s;
      }

      // Accessors ------------------------------------------------------------
      const T& lo() const { return m_lo; }
      const T& hi() const { return m_hi; }
      unsigned ndim() const { return m_ndim; }
      const std::vector<unsigned>& idim() const { return m_idim; }

      double lo(unsigned idim) const { return m_lo[idim]; }
      double hi(unsigned idim) const { return m_hi[idim]; }

      double lo(const T& x, unsigned idim) const
      {
	return m_lo[idim];
      }

      double hi(const T& x, unsigned idim) const
      {
	return m_hi[idim];
      }

      virtual bool isInVolume(const T& x) const { return true; }

      // Setters --------------------------------------------------------------
      T& lo() { return m_lo; }
      T& hi() { return m_hi; }
      void set(const std::vector<double>& lo,
	       const std::vector<double>& hi)
      {
	m_lo.set(lo);
	m_hi.set(hi);
      }

      unsigned m_ndim;
      std::vector<unsigned> m_idim;
      T        m_lo;
      T        m_hi;
      double   s;
    };

    // ========================================================================
    // QuadratureVolume<double>
    // ========================================================================
    template< >
    class QuadratureVolume<double>
    {
    public:
      QuadratureVolume(): m_ndim(1), m_idim(), m_lo(), m_hi(), s() 
      { 
	m_idim.push_back(0);
      }
      QuadratureVolume(double lo, double hi): 
	m_ndim(1), m_idim(), m_lo(lo), m_hi(hi), s() 
      { 
	m_idim.push_back(0);
      }
      
      virtual ~QuadratureVolume() { }

      bool operator <(const QuadratureVolume& rhs) const
      {
	return s < rhs.s;
      }

      // Accessors ------------------------------------------------------------
      double lo() const { return m_lo; }
      double hi() const { return m_hi; }
      unsigned ndim() const { return m_ndim; }
      const std::vector<unsigned>& idim() const { return m_idim; }

      double lo(unsigned idim) const { return m_lo; }
      double hi(unsigned idim) const { return m_hi; }

      void set(const std::vector<double>& lo,
	       const std::vector<double>& hi)
      {
	m_lo = lo[0];
	m_hi = hi[0];
      }

      virtual bool isInVolume(const double& x) const { return true; }

      unsigned m_ndim;
      std::vector<unsigned> m_idim;
      double   m_lo;
      double   m_hi;
      double   s;
    };

    template<typename Fn, typename T, typename Method>
    class Quadrature1D 
    {
    public:
      Quadrature1D(): m_fn(), m_vol(), m_idim() { }
       
      Quadrature1D(const Fn& fn, unsigned idim,
		   const QuadratureVolume<T>& vol):
	m_fn(fn), m_vol(vol), m_idim(idim), m_fnval()
      {

      }

      double val() const
      {
	T x = m_vol.lo();
	return val(x);
      }

      double val(T x) const
      {
	m_vol.lo() = x;
	m_vol.hi() = x;
	QuadratureVolume<T> vol(m_vol);
	double fsum = 0, fdiff = 0;
	std::vector< std::pair<bool,double> > fnval;
	return val(m_vol,fsum,fdiff,fnval);
      }

      double val(const QuadratureVolume<T>& vol, double &fsum, double& fdiff,
		 const std::vector< std::pair<bool,double> >& fnval) const
      {
	const unsigned np = Method::npoints();
	
	T x = vol.lo();

	double lo = vol.lo(m_idim);
	double hi = vol.hi(m_idim);

	const double half_length = 0.5*(hi-lo);
	const double abs_half_length = fabs(half_length);
	const double xcenter = 0.5*(hi+lo);
	
	fsum = 0;
	for(unsigned ip = 0; ip < np; ip++)
	  {	 
	    x.set(m_idim,xcenter + half_length*Method::getAbscissa(ip));

	    if(fnval.empty()) m_fnval[ip] = m_fn.val(x);
	    else if(fnval[ip].first) m_fnval[ip] = fnval[ip].second;
	    else m_fnval[ip] = m_fn.val(x);

	    fsum += m_fn.val(x)*Method::getWeight(ip)*x.jacobian(m_idim);
	  }

	fdiff = 0;
	double favg = 0.5*fsum;
	for(unsigned ip = 0; ip < np; ip++)
	  fdiff += Method::weight(ip)*fabs(m_fnval[ip]-favg);
	
	if(!std::isfinite(fsum)) 
	  throw std::domain_error(std::string(__PRETTY_FUNCTION__)+
				  ": integral is not finite");
	
	fsum *= half_length;
	fdiff *= abs_half_length;

	return fsum;
      }      

    private:      
      Fn                      m_fn;
      QuadratureVolume<T>     m_vol;
      unsigned                m_idim;
      std::vector< double >   m_fnval;
    };

    template<typename Fn, typename T, typename Method>
    class Quadrature2D
    {
    public:
      Quadrature2D(const Fn& fn)
      {
	
      }
      
      void setFn(const Fn& fn) { m_fn = fn; }

      double integrate(T lo, T hi)
      {
	QuadratureVolume<T> bnd_fn(lo,hi);
	Quadrature1D<Fn,T,Method> fn1(m_fn,0,bnd_fn);	
	Quadrature1D<Quadrature1D<Fn,T,Method>, T, Method> fn2(fn1,1,bnd_fn);

	return fn2.val();
      }

    private:      
      Fn                         m_fn;
    };

    // ========================================================================
    // QuadratureND
    //
    // Quadrature class which performs integration over arbitrary N
    // dimensional volume.
    //
    // ========================================================================
    template<typename Fn, typename T, typename Method, 
	     typename VolFn = QuadratureVolume<T> >
    class QuadratureND
    {
    public:
      QuadratureND(const Fn& fn); 
      
      double integrate(const VolFn& vol)
      {
	double fsum, fdiff;
	std::vector< std::pair<bool,double> > fnval;
	return integrate(vol,fsum,fdiff,fnval);
      }

      double integrate(const VolFn& vol,
		       const std::vector< std::pair<bool,double> >& fnval)
      {
	double fsum, fdiff;
	return integrate(vol,fsum,fdiff,fnval);
      }

      double integrate(const VolFn& vol,
		       double & fsum, double& fdiff,
		       const std::vector< std::pair<bool,double> >& fnval);

      double integrate(T lo, T hi)
      {
	return integrate(VolFn(lo,hi));
      }

      void setFn(const Fn& fn) { m_fn = fn; }
       
      unsigned npoints(unsigned ndim) const 
      { 
	return lround(std::pow((double)Method::npoints(),(int)ndim));
      }
      double fnval(unsigned ip) const { return m_fnval[ip]; }

    private:

      //void computeWeights();

      std::vector< double >       m_weights;
      std::vector< double >       m_fnval;
      Fn                          m_fn;
      Method                      m_method;
    };

    template<typename Fn, typename T, typename Method, typename VolFn>
    QuadratureND<Fn,T,Method,VolFn>::QuadratureND(const Fn& fn): 
      m_fnval(), m_fn(fn), m_method()
    { 
      
    }

    // template<typename Fn, typename T, typename Method, typename VolFn>
    // void QuadratureND<Fn,T,Method,VolFn>::computeWeights()
    // {
    //   std::vector< unsigned > n(m_ndim);
    //   m_weights.clear();
    //   m_weights.resize(m_np,1);
    //   for(unsigned ip = 0; ip < m_np; ip++)
    // 	{
    // 	  for(unsigned idim = 0; idim < m_ndim; idim++)
    // 	    m_weights[ip] *= Method::getWeight(n[idim]);

    // 	  for(unsigned idim = 0; idim < m_ndim; idim++)
    // 	    {
    // 	      n[idim]++;
    // 	      n[idim] = n[idim] % Method::npoints();
    // 	      if(n[idim]) break;
    // 	    }
    // 	}
    // }
				   
    // template<typename Fn, typename T, typename Method, typename VolFn>
    // double QuadratureND<Fn,T,Method,VolFn>::integrate(const VolFn& vol)
    // {
    //   std::vector< unsigned > n(m_ndim);
    //   T x = vol.lo();
    //   double integral = 0;

    //   for(unsigned ip = 0; ip < m_np; ip++)
    // 	{	 
    // 	  std::vector<double> xc(m_ndim);
    // 	  for(unsigned idim = 0; idim < m_ndim; idim++)
    // 	    xc[idim] = Method::abscissa(n[idim],vol.lo(idim),vol.hi(idim));
    // 	  x.set(xc);

    // 	  double w = 1;
    // 	  for(unsigned idim = 0; idim < m_ndim; idim++)
    // 	    w *= Method::weight(n[idim],vol.lo(idim),vol.hi(idim));

    // 	  for(unsigned idim = 0; idim < m_ndim; idim++)
    // 	    {
    // 	      n[idim]++;
    // 	      n[idim] = n[idim] % Method::npoints();
    // 	      if(n[idim]) break;
    // 	    }

    // 	  if(!vol.isInVolume(x)) continue;

    // 	  m_fnval[ip] = m_fn.val(x);

    // 	  if(!std::isfinite(m_fnval[ip])) 
    // 	    {	      
    // 	      std::ostringstream os;
    // 	      os << ": Function is not finite. " << std::endl
    // 		 << std::setprecision(12)
    // 		 << "x = "      << std::setw(20) << x 
    // 		 << " w = "     << std::setw(20) << w 
    // 		 << " fnval = " << std::setw(20) << m_fnval[ip] 
    // 		 << std::endl;
    // 	      throw std::domain_error(std::string(__PRETTY_FUNCTION__)+	
    // 				      os.str());
    // 	    }

    // 	  integral += m_fnval[ip]*w*x.jacobian();
    // 	}

    //   return integral;
    // }

    template<typename Fn, typename T, typename Method, typename VolFn>
    double QuadratureND<Fn,T,Method,VolFn>::
    integrate(const VolFn& vol, double& fsum, double& fdiff,
	      const std::vector<std::pair<bool,double> >& fnval)
    {
      const unsigned ndim = vol.ndim();
      const unsigned np = 
	lround(std::pow((double)Method::npoints(),(int)ndim));
      m_fnval.resize(np,0.0);

      std::vector< unsigned > n(ndim,0);
      std::vector< double > half_length(ndim);
      std::vector< double > abs_half_length(ndim);
      std::vector< double > xcenter(ndim);

      unsigned idim = 0;
      for(typename std::vector<unsigned>::const_iterator itr = 
	    vol.idim().begin(); itr != vol.idim().end(); ++itr)
	{
	  half_length[idim] = 0.5*(vol.hi(*itr) - vol.lo(*itr));
	  abs_half_length[idim] = fabs(half_length[*itr]);
	  xcenter[idim] = 0.5*(vol.hi(*itr) + vol.lo(*itr));
	  idim++;
	}

      T x = vol.lo();
      fsum = 0;

      for(unsigned ip = 0; ip < np; ip++)
	{	 
	  //std::vector<double> xc(m_ndim);
	  //for(unsigned idim = 0; idim < m_ndim; idim++)
	  // for(typename std::vector<unsigned>::const_iterator itr = 
	  // 	vol.idim().begin(); itr != vol.idim().end(); ++itr)
	  for(unsigned idim = 0; idim < ndim; idim++)
	    x.set(vol.idim()[idim],xcenter[idim] + 
		  half_length[idim]*Method::abscissa(n[idim]));

	  //	    xc[idim] = xcenter[idim] + 
	  //	      half_length[idim]*Method::abscissa(n[idim]);
	  //Method::abscissa(n[idim],vol.lo(idim),vol.hi(idim));
	  //	  x.set(xc);

	  double w = 1;
	  for(unsigned idim = 0; idim < ndim; idim++)
	    w *= Method::weight(n[idim])*x.jacobian(vol.idim()[idim]);

	  for(unsigned idim = 0; idim < ndim; idim++)
	    {
	      n[idim]++;
	      n[idim] = n[idim] % Method::npoints();
	      if(n[idim]) break;
	    }

	  if(!vol.isInVolume(x)) continue;

	  if(fnval.empty()) m_fnval[ip] = m_fn.val(x);
	  else if(fnval[ip].first) m_fnval[ip] = fnval[ip].second;
	  else m_fnval[ip] = m_fn.val(x);

	  if(!std::isfinite(m_fnval[ip])) 
	    {	      
	      std::ostringstream os;
	      os << ": Function is not finite. " << std::endl
		 << std::setprecision(12)
		 << "x = "      << std::setw(20) << x 
		 << " w = "     << std::setw(20) << w 
		 << " fnval = " << std::setw(20) << m_fnval[ip] 
		 << std::endl;
	      throw std::domain_error(std::string(__PRETTY_FUNCTION__)+	
				      os.str());
	    }

	  fsum += m_fnval[ip]*w;
	}

      fdiff = 0;
      if(ndim == 1)
	{
	  double favg = std::pow(0.5,(int)ndim)*fsum;
	  for(unsigned ip = 0; ip < np; ip++)
	    fdiff += Method::weight(ip)*fabs(m_fnval[ip]-favg);
	}

      for(unsigned idim = 0; idim < ndim; idim++)
	{
	  fsum *= half_length[idim];
	  fdiff *= abs_half_length[idim];
	}

      return fsum;
    }

    // ========================================================================
    // QuadratureND -- double specialization
    // ========================================================================
    template<typename Fn, typename Method>
    class QuadratureND<Fn, double, Method, QuadratureVolume<double> >
    {
    public:
      QuadratureND(const Fn& fn); 
      
      //      double integrate(const QuadratureVolume<double>& vol);
      double integrate(const QuadratureVolume<double>& vol)
      {
	double fsum, fdiff;
	std::vector< std::pair<bool,double> > fnval;
	return integrate(vol,fsum,fdiff,fnval);
      }

      double integrate(const QuadratureVolume<double>& vol,
		       const std::vector< std::pair<bool,double> >& fnval)
      {
	double fsum, fdiff;
	return integrate(vol,fsum,fdiff,fnval);
      }

      double integrate(const QuadratureVolume<double>& vol,
		       double & fsum, double& fdiff,
		       const std::vector< std::pair<bool,double> >& fnval);

      double integrate(double lo, double hi)
      {
	return integrate(QuadratureVolume<double>(lo,hi));
      }

      unsigned npoints(unsigned ndim) const 
      { 
	return Method::npoints();
      }

      void setFn(const Fn& fn) { m_fn = fn; }
      double fnval(unsigned ip) const { return m_fnval[ip]; }

    private:
      //      void computeWeights();

      std::vector< double >       m_weights;
      std::vector< double >       m_fnval;
      Fn                          m_fn;
      Method                      m_method;
    };

    template<typename Fn, typename Method>
    QuadratureND<Fn,double,Method,QuadratureVolume<double> >::
    QuadratureND(const Fn& fn): 
      m_weights(), m_fnval(), m_fn(fn), m_method()
    { 

    }

    // template<typename Fn, typename Method>
    // void QuadratureND<Fn,double,Method,QuadratureVolume<double> >::
    // computeWeights()
    // {
    //   m_weights.clear();
    //   m_weights.resize(m_np,1);
    //   for(unsigned ip = 0; ip < m_np; ip++)
    // 	m_weights[ip] = Method::getWeight(ip);
    // }
		
    // template<typename Fn, typename Method>
    // double QuadratureND<Fn,double,Method,QuadratureVolume<double> >::
    // integrate(const QuadratureVolume<double>& vol)
    // {
    //   double fsum = 0;
    //   const double half_length = 0.5*(vol.hi() - vol.lo());

    //   for(unsigned ip = 0; ip < m_np; ip++)
    // 	{	 
    // 	  double x = Method::abscissa(ip,vol.lo(),vol.hi());
    // 	  double w = Method::weight(ip);

    // 	  if(!vol.isInVolume(x)) continue;

    // 	  m_fnval[ip] = m_fn.val(x);

    // 	  if(!std::isfinite(m_fnval[ip])) 
    // 	    {	      
    // 	      std::ostringstream os;
    // 	      os << ": Function is not finite. " << std::endl
    // 		 << std::setprecision(12)
    // 		 << "x = "      << std::setw(20) << x 
    // 		 << " w = "     << std::setw(20) << w 
    // 		 << " fnval = " << std::setw(20) << m_fnval[ip] 
    // 		 << std::endl;
    // 	      throw std::domain_error(std::string(__PRETTY_FUNCTION__)+	
    // 				      os.str());
    // 	    }

    // 	  fsum += m_fnval[ip]*w;
    // 	}

    //   double favg = 0.5*fsum;

    //   return fsum*half_length;
    // }

    template<typename Fn, typename Method>
    double QuadratureND<Fn,double,Method,QuadratureVolume<double> >::
    integrate(const QuadratureVolume<double>& vol, double& fsum, double& fdiff,
	      const std::vector<std::pair<bool,double> >& fnval )
    {
      const unsigned np = Method::npoints();
      m_fnval.resize(np,0.0);
      fsum = 0;

      const double half_length = 0.5*(vol.hi() - vol.lo());
      const double abs_half_length = fabs(half_length);
      const double xcenter = 0.5*(vol.hi() + vol.lo());

      for(unsigned ip = 0; ip < np; ip++)
	{	 
	  double x = xcenter + half_length*Method::abscissa(ip);
	  if(!vol.isInVolume(x)) continue;
	 
	  if(fnval.empty()) m_fnval[ip] = m_fn.val(x);
	  else if(fnval[ip].first) m_fnval[ip] = fnval[ip].second;
	  else m_fnval[ip] = m_fn.val(x);

	  if(!std::isfinite(m_fnval[ip])) 
	    {	      
	      std::ostringstream os;
	      os << ": Function is not finite. " << std::endl
		 << std::setprecision(12)
		 << "x = "      << std::setw(20) << x 
		 << " w = "     << std::setw(20) << Method::weight(ip)
		 << " fnval = " << std::setw(20) << m_fnval[ip] 
		 << std::endl;
	      throw std::domain_error(std::string(__PRETTY_FUNCTION__)+	
				      os.str());
	    }

	  fsum += m_fnval[ip]*Method::weight(ip);
	}

      double favg = 0.5*fsum;

      fdiff = 0;
      for(unsigned ip = 0; ip < np; ip++)
	fdiff += Method::weight(ip)*fabs(m_fnval[ip]-favg);

      fsum *= half_length;
      fdiff *= abs_half_length;

      return fsum;
    }

    // ========================================================================
    // QuadratureRuleGauss7
    //
    // Defines weights and abscissas for 7-point gaussian quadrature.
    //
    // ========================================================================
    class QuadratureRuleGauss7 
    {
    public:
      QuadratureRuleGauss7(): m_abscissa(s_npoints,0) { }

      static double abscissa(unsigned i, double lo, double hi)
      {
	return ((hi-lo)*s_abscissas[i] + (hi+lo))/2.;
      }

      static double weight(unsigned i, double lo, double hi)
      {
	return s_weights[i]*(hi-lo)/2.;
      }

      static double abscissa(unsigned i)
      {
	return s_abscissas[i];
      }

      static double weight(unsigned i)
      {
	return s_weights[i];
      }

      static unsigned npoints() { return s_npoints; }
      static double getAbscissa(unsigned i) { return s_abscissas[i]; }
      static double getWeight(unsigned i) { return s_weights[i]; }

    private:

      std::vector<double>   m_abscissa;

      static const unsigned s_npoints = 7;
      static const double   s_abscissas[s_npoints];
      static const double   s_weights[s_npoints];
    };

    // ========================================================================
    // QuadratureRuleKronrod15
    //
    // Defines weights and abscissas for 15-point kronrod quadrature.
    //
    // ========================================================================
    class QuadratureRuleKronrod15 
    {
    public:
      QuadratureRuleKronrod15() { }
     
      static double abscissa(unsigned i, double lo, double hi)
      {
	return ((hi-lo)*s_abscissas[i] + (hi+lo))/2.;
      }

      static double weight(unsigned i, double lo, double hi)
      {
	return s_weights[i]*(hi-lo)/2.;
      }

      static double abscissa(unsigned i)
      {
	return s_abscissas[i];
      }

      static double weight(unsigned i)
      {
	return s_weights[i];
      }

      static unsigned npoints() { return s_npoints; }
      static double getAbscissa(unsigned i) { return s_abscissas[i]; }
      static double getWeight(unsigned i) { return s_weights[i]; }

    private:
      static const unsigned s_npoints = 15;
      static const double   s_abscissas[s_npoints];
      static const double   s_weights[s_npoints];
    };

    // ========================================================================
    // Gauss-Kronrod Adaptive Quadrature
    //
    // Adaptive quadrature routine that uses 7-point Gauss and
    // 15-point Kronrod rules to calculate an estimate of the integral
    // and the error on that estimate.  On each iteration if the
    // fractional error is greater than the tolerance, the integration
    // volume is split into 2^N sub volumes.  Integration volume can be
    // defined using a user-defined functor class (VolFn).
    //
    // ========================================================================
    template<typename Fn, typename T, typename VolFn = QuadratureVolume<T> >
    class QuadratureGaussKronrod
    {
    public:
      QuadratureGaussKronrod(const Fn& fn = Fn(), 
			     double tol = s_default_tolerance,
			     unsigned max_niter = s_default_max_niter); 
      
      double integrate(T lo, T hi)
      {
	VolFn vol(lo,hi);
	m_ndim = vol.ndim();
	m_vol.push_back(vol);
	return integrate();
      }

      double integrate(const VolFn& v)
      {
	m_ndim = v.ndim();
	m_vol.push_back(v);
	return integrate();
      }

      double integrate();

      void setFn(const Fn& fn) 
      { 
	m_quad_gauss7.setFn(fn); 
	m_quad_kronrod15.setFn(fn); 
      }

      double sum() const { return m_sum; }
      double err() const { return m_err; }
      unsigned niter() const { return m_niter; }

      void setMaxIterations(unsigned niter) { m_max_niter = niter; }

    private:
      void split(VolFn& vol);

      QuadratureND<Fn,T,QuadratureRuleGauss7>    m_quad_gauss7;
      QuadratureND<Fn,T,QuadratureRuleKronrod15> m_quad_kronrod15;

      std::vector< VolFn >                   m_vol;
      unsigned                               m_ndim;

      double                                 m_sum;
      double                                 m_err;
      unsigned                               m_niter;

      double                                 m_tol;
      unsigned                               m_max_niter;

      static double                          s_default_tolerance;
      static unsigned                        s_default_max_niter;
    };

    template<typename Fn, typename T, typename VolFn>
    QuadratureGaussKronrod<Fn,T,VolFn>::
    QuadratureGaussKronrod(const Fn& fn,
			   double tol, unsigned max_niter):
      m_quad_gauss7(fn), m_quad_kronrod15(fn), m_vol(), m_ndim(),
      m_sum(), m_err(), m_niter(),
      m_tol(tol), m_max_niter(max_niter)
    {
      
    }
    
    template<typename Fn, typename T, typename VolFn>
    double QuadratureGaussKronrod<Fn,T,VolFn>::integrate() 
    {
      m_sum = 0;
      m_err = 0;
      m_niter = 0;

      std::vector< std::pair<bool,double> > 
	fnval(m_quad_kronrod15.npoints(m_ndim));

      while(!m_vol.empty())
	{
	  VolFn v = m_vol.back();
	  m_vol.pop_back();

	  double g7_sum = m_quad_gauss7.integrate(v);

	  // ------------------------------------------------------------------
	  // Reuse function evaluations from 7-point gaussian quadrature
	  // ------------------------------------------------------------------
	  for(unsigned i = 0; i < m_quad_gauss7.npoints(m_ndim); i++)
	    {
	      unsigned k = 0;
	      for(unsigned idim = 0, j = i, n = 1; idim < m_ndim; idim++)
	  	{
	  	  k += (2*(j%QuadratureRuleGauss7::npoints())+1)*n;
	  	  j /= QuadratureRuleGauss7::npoints();
	  	  n *= QuadratureRuleKronrod15::npoints();
	  	}

	      fnval[k].first = true;
	      fnval[k].second = m_quad_gauss7.fnval(i);
	    }

	  double k15_sum, k15_diff;
	  m_quad_kronrod15.integrate(v,k15_sum,k15_diff,fnval);
	  //	  double k = m_quad_kronrod15.integrate(v);

	  v.s = k15_sum;

	  double gk_diff = fabs(k15_sum - g7_sum);
	  double gk_err = 0;

	  if(m_ndim == 1)
	    {
	      double scale = std::pow((200*gk_diff/k15_diff),1.5);
	      if(scale < 1) gk_err = k15_diff*scale;
	      else gk_err = k15_diff;
	    }
	  else gk_err = gk_diff;
	  
	  double sum = m_sum+k15_sum;
	  for(typename std::vector<VolFn>::iterator itr = m_vol.begin(); itr !=
		m_vol.end(); ++itr)
	    sum += itr->s;

 	  // std::cout 
	  //   << std::setw(22) << g7_sum 	    
	  //   << std::setw(22) << k15_sum 
	  //   << std::setw(22) << gk_err 
	  //   << std::setw(22) << gk_diff
	  //   << std::setw(22) << m_sum
	  //   << std::setw(22) << sum
	  //   << std::setw(22) << m_vol.size()
	  //   << std::endl;
// 		    << std::setw(22) << m_tol*sum 
// 		    << std::setw(22) << k15_sum 
// 		    << std::setw(22) << gk_err 
// 		    << std::setw(22) << gk_err/k15_sum
// 		    << std::setw(22) << v.lo()
// 		    << std::setw(22) << v.hi()
// 		    << std::endl;

	  if(gk_err <= m_tol*fabs(sum))
	    {
	      m_sum += k15_sum;
	      m_err += gk_err;
	    }
	  else split(v);

	  // ------------------------------------------------------------------
	  // Sort volumes by integral estimate
	  // ------------------------------------------------------------------
	  std::sort(m_vol.begin(),m_vol.end());

	  m_niter++;
	  if(m_niter == m_max_niter) 
	    {
	      std::ostringstream os;
	      os << ": No convergence after " << m_max_niter 
		 << " iterations." << std::endl
		 << std::setprecision(14)
		 << "g7       = " << std::setw(22) << g7_sum << std::endl
		 << "k15      = " << std::setw(22) << k15_sum << std::endl
		 << "gk_err   = " << std::setw(22) << gk_err << std::endl
		 << "gk_diff  = " << std::setw(22) << gk_diff << std::endl
		 << "sum      = "  << std::setw(22) << sum << std::endl
		 << "m_vol.size() = " << m_vol.size() << std::endl
		 << "m_vol.back().lo() = " 
		 << std::setw(22) << m_vol.back().lo() << std::endl
		 << "m_vol.back().hi() = " 
		 << std::setw(22) << m_vol.back().hi() << std::endl;

	      throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
				       os.str());
	    }
	  
	}

      return m_sum;
    }
    


    // ========================================================================
    // Monte Carlo Quadrature
    //
    // Description...
    //
    // ========================================================================
    template<typename Fn, typename T, typename VolFn = QuadratureVolume<T> >
    class QuadratureMC
    {
    public:
      QuadratureMC() { }
      
      // For ring bkgd model integration over non-trivial volume, 
      // bc of star exclusion regions.
      double integrate(const VolFn& v);
      
    private:
      
  };
    
    
    
    template<typename Fn, typename T, typename VolFn>
    void QuadratureGaussKronrod<Fn,T,VolFn>::split(VolFn& vol) 
    {
      const unsigned ndim = vol.ndim();
      const unsigned nvol = lround(std::pow(2.0,(int)ndim));
      for(unsigned ivol = 0; ivol < nvol; ivol++)
	{
	  VolFn v = vol;
	  std::vector<double> lo(ndim);
	  std::vector<double> hi(ndim);
	  for(typename std::vector<unsigned>::const_iterator itr = 
		vol.idim().begin(); itr != vol.idim().end(); ++itr)
	    {
	      unsigned idim = *itr;
	    //	  for(unsigned idim = 0; idim < ndim; idim++)
	    //	    {
	      double dx = (vol.hi(idim) - vol.lo(idim))/2.;
	      lo[idim] = vol.lo(idim) + dx*((ivol>>idim)&1);
	      hi[idim] = vol.lo(idim) + dx*(((ivol>>idim)&1)+1);
	    }
	  v.set(lo,hi);
	  v.s = vol.s/(double)nvol;
	  m_vol.push_back(v);
	}
    }
  
    template<typename Fn, typename T, typename VolFn>
    double QuadratureGaussKronrod<Fn,T,VolFn>::s_default_tolerance = 1E-4;
    template<typename Fn, typename T, typename VolFn>
    unsigned QuadratureGaussKronrod<Fn,T,VolFn>::s_default_max_niter = 1000;

    template< typename Fn >
    class QuadFn
    {
    public:
      QuadFn(const Fn &fn): m_fn(fn) { }

      double val(const VSACoord::Coord1D& x) const { return m_fn.val(x.x()); }

    private:
      Fn m_fn;
    };

    // ========================================================================
    // Quadrature Helper Class
    //
    // The static methods of this class can be used in lieu of
    // manually instantiating one of the quadrature templates.  The
    // function to be integrated is provided as a functor defining a
    // val(T) method which returns the value of the function at the
    // given coordinate.
    // ========================================================================
    class Quadrature
    {
    public:

      template< typename Fn, typename T >
      static double integrate1D(const Fn& fn, T lo, T hi, unsigned idim,
				double tol = 1E-4)
      {
	QuadratureVolume<T> vol(lo,hi,idim);
	QuadratureGaussKronrod<Fn,T> quad(fn,tol);
	return quad.integrate(vol);
      }

      template< typename Fn, typename T >
      static double integrate(const Fn& fn, T lo, T hi, double tol = 1E-4)
      {
	QuadratureGaussKronrod<Fn,T> quad(fn,tol);
	return quad.integrate(lo,hi);
      }

      template< typename Fn, typename T >
      static void integrate(const Fn& fn, T lo, T hi, 
			    double& sum, double& err,
			    double tol = 1E-4)
      {
	QuadratureGaussKronrod<Fn,T> quad(fn,tol);
	quad.integrate(lo,hi);
	sum = quad.sum();
	err = quad.err();
      }

      // template< typename Fn >
      // static double integrate(const Fn& fn, double lo, double hi, 
      // 			      double tol = 1E-4)
      // {	
      // 	QuadFn<Fn> quad_fn(fn);
      // 	QuadratureGaussKronrod<QuadFn<Fn>,VSACoord::Coord1D> 
      // 	  quad(quad_fn,tol);
      // 	return quad.integrate(lo,hi);
      // }

      template< typename Fn, typename T, typename VolFn >
      static double integrate(const Fn& fn, const VolFn& v)
      {
	QuadratureGaussKronrod<Fn,T,VolFn> quad(fn);
	return quad.integrate(v);
      }

    };

    class QuadratureFixed
    {
    public:

      template< typename Fn, typename T >
      static double integrate(const Fn& fn, T lo, T hi)
      {
	QuadratureND<Fn,T,QuadratureRuleKronrod15> qd(fn);
	return qd.integrate(lo,hi);
      }

    };
  }
}

#endif // VSAQUADRATURE_HPP
