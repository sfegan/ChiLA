//-*-mode:c++; mode:font-lock;-*-

/*! \file VSANonlinearFitting.hpp

  Routines for nonlinear function minimization.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       01/20/2009
*/

#ifndef VSANONLINEARFITTING_HPP
#define VSANONLINEARFITTING_HPP

#include <vector>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include <VSAAlgebra.hpp>
#include <VSASVD.hpp>
#include <VSAMath.hpp>
#include <VSAData.hpp>
#include <VSAssert.hpp>
#include <VSAFunction.hpp>
#include <VSAMinimization.hpp>
#include <VSABracketedMonotonicRootFinder.hpp>

namespace VERITAS
{
  namespace VSAMath
  {
    // ------------------------------------------------------------------------
    // Documentation 
    //
    // The nonlinear fitting routines provided here are designed to
    // perform function minimization with respect to one or more
    // parameters.  
    // ------------------------------------------------------------------------

    class ConvergenceFailed : public std::runtime_error
    {
    public:
      ConvergenceFailed(const std::string& e): runtime_error(e) { }
    };

    //! 
    //! NLFitter: Abstract base class for non-linear fitter.  Objective
    //! function (e.g. Log-Likelihood or Chi-Squared) to be minimized is 
    //! provided as a template parameter.
    //! 
    template< typename Fn >
    class NLFitter
    {
    public:
      struct Options 
      { 
	Options(): 
	  min_tolerance(0.01), 
	  max_tolerance(0.01), 
	  min_iterations(2) { }

	double    min_tolerance;
	double    max_tolerance;
	unsigned  min_iterations;
      };


      NLFitter(const Options& opt = defaultOptions()):
	m_options(opt),
	m_scale(),
	m_has_lo_bound(), m_lo_bound(),
	m_has_hi_bound(), m_hi_bound()
      { }

      virtual ~NLFitter() { }

      //! Set the parameter starting value.
      //! @param ip Parameter index.
      //! @param val Parameter value.
      virtual void set(unsigned ip, double val) = 0;

      //! Set the parameter starting values.
      //! @param a Vector of parameter values.
      virtual void set(const VSAAlgebra::VecND& a) = 0;

      //! Set the objective function.
      //! @param fn Reference to objective function.
      virtual void setFn(const Fn& fn)
      {
	m_scale.resize(fn.nparm(),1.0);
	m_has_lo_bound.resize(fn.nparm(),false); 
	m_lo_bound.resize(fn.nparm(),0.0); 
	m_has_hi_bound.resize(fn.nparm(),false); 
	m_hi_bound.resize(fn.nparm(),0.0); 
      }

      //! Fix a parameter during fitting.
      //! @param ip Parameter index.
      virtual void hold(unsigned ip) = 0;

      //! Fix a parameter at the specified value during fitting.
      //! @param ip Parameter index.
      //! @param val Parameter value.
      virtual void hold(unsigned ip, double val) = 0;

      //! Free all parameters.
      virtual void free() = 0;

      //! Free a held parameter.
      //! @param ip Parameter index.
      virtual void free(unsigned ip) = 0;

      //! Set rescaling factor.
      void setScaling(unsigned ip, double scale)
      {
	m_scale[ip] = scale;
      }

      //! Set a lower bound on a parameter.
      void setLoBound(unsigned ip, double val)
      {
	m_has_lo_bound[ip] = true;
	m_lo_bound[ip] = val;
      }

      //! Set an upper bound on a parameter.
      void setHiBound(unsigned ip, double val)
      {
	m_has_hi_bound[ip] = true;
	m_hi_bound[ip] = val;
      }

      //! Initialize fitter using current parameter values.  Call this
      //! function before calling NLFitter::fit().
      virtual void initialize() = 0;

      //! Initialize fitter using given parameter values.  Call this
      //! function before calling NLFitter::fit().
      virtual void initialize(const VSAAlgebra::VecND& a) = 0;

      //! Initialize fitter using given parameter values and array specifying
      //! whether each parameter should be fixed.  Call this
      //! function before calling NLFitter::fit().
      virtual void initialize(const VSAAlgebra::VecND& a, 
			      const std::vector< bool >& fixed) = 0;

      //! Minimize the objective function with respect to all free 
      //! parameters.
      virtual void fit() = 0;

      virtual void fit(const VSAAlgebra::VecND& a) = 0;

      //! Scan for the smallest value of the the objective function
      //! over the given parameter range.
      //! @param ip Parameter index.
      //! @param lo Lower search bound.
      //! @param hi Upper search bound.
      //! @param step Search step size.
      virtual void scan(unsigned ip, double lo, double hi, double step)
      {
	vsassert(ip < nparam());
	for(double p = lo; p <= hi; p += step)
	  {
	    VSAAlgebra::VecND a = param();
	    a[ip] = p;
	    if(chi2(a) < chi2(param())) set(a);
	  }
	initialize();
      }

      //! Scan for the smallest value of the the objective function
      //! over the given range of two parameters.
      //! @param ip1 Parameter index.
      //! @param lo1 Lower search bound.
      //! @param hi1 Upper search bound.
      //! @param step1 Search step size.
      virtual void scan(unsigned ip1, double lo1, double hi1, double step1,
			unsigned ip2, double lo2, double hi2, double step2)
      {
	vsassert(ip1 < nparam() && ip2 < nparam());
	for(double p1 = lo1; p1 <= hi1; p1 += step1)
	  {
	    for(double p2 = lo2; p2 <= hi2; p2 += step2)
	      {
		VSAAlgebra::VecND a = param();
		a[ip1] = p1;
		a[ip2] = p2;

		if(chi2(a) < chi2(param())) set(a);
	      }
	  }
	initialize();
      }

      //! Return objective function.
      virtual Fn& fn() = 0;

      void setTolerance(double tol) 
      { 
	m_options.min_tolerance = tol; 
	m_options.max_tolerance = tol; 
      }

      void setMinTolerance(double tol) 
      { 
	vsassert(tol <= m_options.max_tolerance);
	m_options.min_tolerance = tol; 
      }

      void setMaxTolerance(double tol) 
      { 
	vsassert(tol >= m_options.min_tolerance);
	m_options.max_tolerance = tol; 
      }

      void setMinIterations(unsigned niter) 
      { m_options.min_iterations = niter; }

      // ----------------------------------------------------------------------
      // Accessors 
      // ----------------------------------------------------------------------

      //! Return number of parameters.  
      virtual unsigned nparam() const = 0;

      //! Return vector of parameters.  
      virtual const VSAAlgebra::VecND& param() const = 0;

      //! Return parameter value.  
      virtual double param(unsigned ip) const = 0;

      //! Return vector of parameter errors.  
      virtual const VSAAlgebra::VecND& err() const = 0;

      //! Return parameter error.  
      virtual double err(unsigned ip) const = 0;

      //! Return parameter covariance matrix.  
      virtual const VSAAlgebra::MatrixND& cov() const = 0;

      //! Return value of objective function.      
      virtual double chi2() const = 0;

      //! Return value of objective function for given parameter values.
      virtual double chi2(const VSAAlgebra::VecND& a) const = 0;
      virtual void chi2(const VSAAlgebra::VecND& a, 
			VSAAlgebra::VecND& dyda) const = 0;

      //! Return number of degrees of freedom.
      virtual unsigned ndf() const = 0;

      unsigned minIterations() const { return m_options.min_iterations; }
      double minTolerance() const { return m_options.min_tolerance; }
      double maxTolerance() const { return m_options.max_tolerance; }

      double scale(unsigned ip) const { return m_scale[ip]; }

      double hiBound(unsigned ip) const { return m_hi_bound[ip]; }
      bool hasHiBound(unsigned ip) const { return m_has_hi_bound[ip]; }
      double loBound(unsigned ip) const { return m_lo_bound[ip]; }
      bool hasLoBound(unsigned ip) const { return m_has_lo_bound[ip]; }

      static Options defaultOptions() { return Options(); }

      //! Return objective function.
      virtual const Fn& fn() const = 0;

      //! Virtual constructor
      virtual NLFitter<Fn>* clone() const = 0;

    private:

      Options                m_options;

      std::vector<double>    m_scale;
      std::vector<bool>      m_has_lo_bound;
      std::vector<double>    m_lo_bound;
      std::vector<bool>      m_has_hi_bound;
      std::vector<double>    m_hi_bound;
    };

    template<typename Fn> 
    class NLFitterPowell : public NLFitter<Fn>
    {
    public:
      // struct Options 
      // { 
      // 	Options(): tolerance(1E-6), min_iterations(2) { }

      // 	double    tolerance;
      // 	unsigned  min_iterations;
      // };

      NLFitterPowell(const Fn& fn, 
		     const typename NLFitter<Fn>::Options& opt = 
		     NLFitter<Fn>::defaultOptions()):
	NLFitter<Fn>(opt),
	m_fn(), m_nparam(), m_nfit(), m_fit(), 
	m_a(), m_err(), m_cov(), m_chi2(), m_tol(1E-6)
      {
	setFn(fn);	
      }
      
      void set(const VSAAlgebra::VecND& a) { m_a = a; }
      void set(unsigned i, double val) { m_a[i]=val; }
      void setFn(const Fn& fn) 
      {   
	NLFitter<Fn>::setFn(fn);
	m_fn = fn;
	m_nparam = fn.nparm();
	m_a      = fn.param();
	m_fit   .resize(m_nparam,true); 
	m_err   .resize(m_nparam,0.0);
	m_cov   .resize(m_nparam,m_nparam,0.0); 

	for(unsigned ip = 0; ip < m_nparam; ip++)
	  if(m_fn.fixed(ip)) hold(ip,m_fn.param(ip));
      }

      void hold(unsigned i) { m_fit[i] = false; }
      void hold(unsigned i, double val) { m_fit[i] = false; m_a[i]=val; }
      void free() { std::fill_n(m_fit.begin(),m_fit.size(),true); }
      void free(unsigned i) { m_fit[i] = true; }

      void initialize() { initialize(m_a); }      
      void initialize(const VSAAlgebra::VecND& a)
      {
	m_a = a;
      }

      void initialize(const VSAAlgebra::VecND& a, 
		      const std::vector< bool >& fixed)
      {
	for(unsigned ip=0;ip<fixed.size();ip++) m_fit[ip] = !fixed[ip];
	initialize(a);
      }

      void fit()
      {
	VSAAlgebra::MatrixND xi(nparam());
	unsigned nfit = 0;
	for(unsigned ip = 0; ip < nparam(); ip++)
	  if(m_fit[ip])
	    {
	      xi(ip,ip) = 1.0;
	      nfit++;
	    }

	if(nfit == 0)
	  {
	    m_chi2 = m_fn.val(m_a);
	    return;
	  }

	int iter = 0;
	VSAMath::powell<Fn>(m_a,xi,m_tol,iter,m_chi2,m_fn);
 

	xi.set(0.0);

	for(unsigned ip = 0; ip < nparam(); ip++)
	  if(m_fit[ip]) xi(ip,ip) = 1.0;

	VSAMath::powell<Fn>(m_a,xi,m_tol,iter,m_chi2,m_fn);
     }

      void fit(const VSAAlgebra::VecND& a)
      {
	m_a = a;
	fit();
      }

      Fn& fn() { return m_fn; }

      // Accessors ------------------------------------------------------------
      unsigned nparam() const { return m_a.size(); }
      const VSAAlgebra::VecND& param() const { return m_a; }
      double param(unsigned ip) const { return m_a[ip]; }
      const VSAAlgebra::VecND& err() const { return m_err; }
      double err(unsigned ip) const { return m_err[ip]; }
      const VSAAlgebra::MatrixND& cov() const { return m_cov; }
      double chi2() const { return m_chi2; }
      double chi2(const VSAAlgebra::VecND& a) const 
      { 
	return m_fn.val(a);
      }

      void chi2(const VSAAlgebra::VecND& a,
		VSAAlgebra::VecND& dyda) const 
      { 
	
      }

      unsigned ndf() const
      {
	return 0;
      }

      const Fn& fn() const { return m_fn; }
      
      //! Virtual constructor
      virtual NLFitterPowell<Fn>* clone() const
      {
	return new NLFitterPowell<Fn>(*this);
      }

    private:

      static const unsigned MAX_ITERATIONS = 1000;

      Fn                     m_fn;

      unsigned               m_nparam;
      unsigned               m_nfit;
      std::vector<bool>      m_fit;
      VSAAlgebra::VecND      m_a;
      VSAAlgebra::VecND      m_err;
      VSAAlgebra::MatrixND   m_cov;
      double                 m_chi2;
      double                 m_tol;
    };

    ///////////////////////////////////////////////////////////////////////////
    //! NLFitterLM: Nonlinear fitter class that implements Levenberg-Marquardt
    //! algorithm from Numerical Recipes Section 15.5
    template<typename Fn> 
    class NLFitterLM : public NLFitter<Fn>
    {
    public:


      //      static Options defaultOptions() { return Options(); }

      NLFitterLM(const typename NLFitter<Fn>::Options& opt = 
		 NLFitter<Fn>::defaultOptions()):
	NLFitter<Fn>(opt),
	m_fn(), m_nparam(), m_nfit(), m_fit(), 
	m_a(), m_err(), m_da(), m_cov(), m_lambda(), m_chi2()
      {
	setFn(m_fn);
      }

      NLFitterLM(const Fn& fn, 
		 const typename NLFitter<Fn>::Options& opt = 
		 NLFitter<Fn>::defaultOptions()):
	NLFitter<Fn>(opt),
	m_fn(), m_nparam(), m_nfit(), m_fit(), 
	m_a(), m_err(), m_da(), m_cov(), m_lambda(), m_chi2()
      {
	setFn(fn);	
      }

      ~NLFitterLM()
      {

      }

      void set(const VSAAlgebra::VecND& a) { m_a = a; }
      void set(unsigned i, double val) { m_a[i]=val; }
      void setFn(const Fn& fn) 
      {   
	NLFitter<Fn>::setFn(fn);
	m_fn = fn;
	m_nparam = fn.nparm();
	m_a      = fn.param();
	m_da    .resize(m_nparam,0.0);
	m_fit   .resize(m_nparam,true); 
	m_err   .resize(m_nparam,0.0);
	m_cov   .resize(m_nparam,m_nparam,0.0); 
	m_alpha .resize(m_nparam,m_nparam,0.0);
	m_beta  .resize(m_nparam,0.0); 

	for(unsigned ip = 0; ip < m_nparam; ip++)
	  if(m_fn.fixed(ip)) hold(ip,m_fn.param(ip));
      }

      void hold(unsigned i) { m_fit[i] = false; }
      void hold(unsigned i, double val) { m_fit[i] = false; m_a[i]=val; }
      void free() { std::fill_n(m_fit.begin(),m_fit.size(),true); }
      void free(unsigned i) { m_fit[i] = true; }

      void initialize() { initialize(m_a); }      
      void initialize(const VSAAlgebra::VecND& a); 
      void initialize(const VSAAlgebra::VecND& a, 
		      const std::vector< bool >& fixed)
      {
	for(unsigned ip=0;ip<fixed.size();ip++) m_fit[ip] = !fixed[ip];
	initialize(a);
      }

      void fit();

      void fit(const VSAAlgebra::VecND& a)
      {
	initialize(a);
	fit();
      }

      void iterate();

      // Setters --------------------------------------------------------------
      Fn& fn() { return m_fn; }

      // Accessors ------------------------------------------------------------
      unsigned nparam() const { return m_a.size(); }
      const VSAAlgebra::VecND& param() const { return m_a; }
      double param(unsigned ip) const { return m_a[ip]; }
      const VSAAlgebra::VecND& err() const { return m_err; }
      double err(unsigned ip) const { return m_err[ip]; }
      const VSAAlgebra::MatrixND& cov() const { return m_cov; }
      double chi2() const { return m_chi2; }
      double chi2(const VSAAlgebra::VecND& a) const 
      { 
	return m_fn.val(a);
      }

      void chi2(const VSAAlgebra::VecND& a,
		VSAAlgebra::VecND& dyda) const 
      { 
	m_fn.dyda(a,dyda);
      }

      unsigned ndf() const
      {
	if(m_nfit > m_fn.ndata()) return 0;
	else return m_fn.ndata() - m_nfit;
      }

      const Fn& fn() const { return m_fn; }

      //! Virtual constructor
      virtual NLFitterLM<Fn>* clone() const
      {
	return new NLFitterLM<Fn>(*this);
      }

      void generateContour(double dchi2, 
			   unsigned ip1, double ip1lo, double ip1hi,
			   unsigned ip2, double ip2lo, double ip2hi,
			   std::vector< std::pair<double,double> >& xy);

      void findRoot(VSAAlgebra::VecND param,
		    unsigned ip1,
		    unsigned ip2,
		    double fnval,
		    std::vector< std::pair<double,double> >& xy1,
		    std::vector< std::pair<double,double> >& xy2);


    private:
      bool mrqcof(VSAAlgebra::MatrixND& alpha, VSAAlgebra::VecND& beta,
		  VSAAlgebra::VecND& a, double& chi2);

      static const unsigned MAX_ITERATIONS = 1000;

      Fn                     m_fn;

      unsigned               m_nparam;
      unsigned               m_nfit;
      std::vector<bool>      m_fit;
      VSAAlgebra::VecND      m_a;
      VSAAlgebra::VecND      m_err;
      VSAAlgebra::VecND      m_da;
      VSAAlgebra::MatrixND   m_cov;
      VSAAlgebra::MatrixND   m_alpha;
      VSAAlgebra::VecND      m_beta;
      double                 m_lambda;
      double                 m_chi2;
    };

    // ------------------------------------------------------------------------
    // NLFitterLM::initialize() - Call this before fit() with the starting
    // values for the fitted parameters.
    // ------------------------------------------------------------------------
    template<typename Fn> 
    void NLFitterLM<Fn>::initialize(const VSAAlgebra::VecND& a)
    {
      vsassert(m_fn.ndata() > 0);

      m_nfit = 0;
      for(unsigned iparm=0;iparm<m_nparam;iparm++)if(m_fit[iparm])m_nfit++;

      vsassert(a.ndim() <= m_a.ndim());
      for(unsigned ip = 0; ip < a.ndim(); ip++) 
	m_a[ip] = a[ip];

      m_lambda = 0.001;
      m_chi2 = 0;

      m_alpha = VSAAlgebra::MatrixND(m_nfit,m_nfit,0.0);
      m_beta = VSAAlgebra::VecND(m_nfit,0.0);
      m_da = VSAAlgebra::VecND(m_nfit,0.0);

      try
	{
	  mrqcof(m_alpha,m_beta,m_a,m_chi2);
	}
      catch(const std::domain_error& e)
	{
	  m_chi2 = 0;
	}

      vsassert(std::isfinite(m_chi2));
    }

    // ------------------------------------------------------------------------
    // NLFitterLM::fit() - Main fitting routine which minimizes the
    // objective function to the specified fit tolerance.
    // ------------------------------------------------------------------------
    template<typename Fn> 
    void NLFitterLM<Fn>::fit()
    {
      vsassert(m_fn.ndata() > 0);

      const double min_tol = NLFitter<Fn>::minTolerance();
      const double max_tol = NLFitter<Fn>::maxTolerance();

      unsigned niter = 0;
      unsigned nzero = 0;
      double dchisq = 0;
      double chi2 = m_chi2;
      double tol = min_tol;
      double dtol = (max_tol - min_tol)/10.;
      while(1)
	{
	  if(niter >= MAX_ITERATIONS && tol < max_tol)
	    {
	      tol += dtol;
	      std::cerr << std::string(__PRETTY_FUNCTION__)
			<< ": Fit failed to converge after " << MAX_ITERATIONS
			<< " iterations. Tolerance = " << tol << std::endl;
	      niter = 0;	      
	    }
	  else if(niter >= MAX_ITERATIONS && tol >= max_tol)
	    {
	      std::ostringstream os;
	      os << std::string(__PRETTY_FUNCTION__)
		 << ": Fit failed to converge after " << MAX_ITERATIONS
		 << " iterations. Tolerance = " << tol << std::endl;
	      throw ConvergenceFailed(os.str());
	    }

	  chi2 = m_chi2;

	  iterate();
	  
	  dchisq = (m_chi2-chi2);	

	  if(dchisq == 0) nzero++;
	  else nzero = 0;

	  // 	  	  std::cout << std::setw(10) << niter
	  // 	  		    << std::setprecision(10)
	  // 	  		    << std::setw(20) << chi2
	  // 	  		    << std::setw(20) << m_chi2
	  // 	  		    << std::setw(20) << dchisq
	  // 	  		    << std::setw(20) << fabs(dchisq/m_chi2)
	  // 	  		    << std::setw(20) << m_lambda
	  // 	  		    << std::endl;


	  if(niter >= NLFitter<Fn>::minIterations() && fabs(dchisq) < tol &&
	     dchisq < 0)
	    break;
	  else if(nzero > 10 && fabs(dchisq) < tol)
	    break;

	  niter++;
	}

      m_lambda = 0;

      iterate();
    }
    
    // ------------------------------------------------------------------------
    // NLFitterLM::iterate() - Called by fit() for each fit iteration.
    // Attempts to improve chi2 by choosing a step size for each free
    // parameter using the objective function derivatives.
    // ------------------------------------------------------------------------
    template<typename Fn> 
    void NLFitterLM<Fn>::iterate()
    {
      VSAAlgebra::MatrixND alpha(m_nfit,m_nfit,0.0);
      VSAAlgebra::MatrixND scale(m_nfit,m_nfit,0.0);

      for(unsigned ip = 0, ifit = 0; ip < m_nparam; ip++)
	{
	  if(m_fit[ip])
	    {
	      scale(ifit,ifit) = 1./NLFitter<Fn>::scale(ip);
	      ifit++;
	    }
	}

      for(unsigned ifit=0; ifit<m_nfit; ifit++)
	{
	  for(unsigned jfit=0; jfit<m_nfit; jfit++)
	    alpha(ifit, jfit) = m_alpha(ifit, jfit);

	  alpha(ifit,ifit)=m_alpha(ifit,ifit)*(1+m_lambda);
	}

      VSAAlgebra::VecND b = m_beta;

      b = scale*b;
      alpha = scale*alpha*scale;

      // std::cout << "-----------------------------------" << std::endl;
      // std::cout << "b = " << std::endl << b << std::endl;
      // std::cout << "alpha = " << std::endl << alpha << std::endl;
      // std::cout << "scale = " << std::endl << scale << std::endl;

      try
	{	  
	  // Solve for inverse of alpha and solution to normal equations ------
	  VSAAlgebra::SVD svd(alpha);
	  svd.solve(b,m_da);
	  
	  // // Generate inverse by SVD if matrix has less than full rank
	  // if(svd.rank() < alpha.nrow()) 
	  //   {
	  //     for(unsigned ifit = 0; ifit < m_nfit; ifit++)
	  // 	beta(ifit,0) = b2(ifit);
	  //   }
	  // // Generate inverse by Gauss-Jordan otherwise
	  // else
	  //   VSAAlgebra::MatrixND::gaussJordan(alpha2, beta);
	}
      catch(const std::exception& e)
	{
	  std::cout << e.what() << std::endl;
	  std::cout << alpha << std::endl;
	}  

      m_da = scale*m_da;

      // Evaluate covariance matrix -------------------------------------------
      if(m_lambda == 0)
	{
	  VSAAlgebra::SVD svd(alpha);
	  VSAAlgebra::MatrixND cov = scale*svd.inverse()*scale;

	  for(unsigned ifit=0; ifit<m_nfit; ifit++)
	    for(unsigned jfit=0; jfit<m_nfit; jfit++)
	      m_cov(ifit, jfit) = cov(ifit, jfit);

	  for(unsigned iparm=m_nfit; iparm<m_nparam; iparm++)
	    for(unsigned jparm=0; jparm<iparm+1; jparm++)
	      m_cov(iparm,jparm) = m_cov(jparm,iparm) = 0.0;
      
	  unsigned kparm = m_nfit-1;
	  for(int jparm = m_nparam-1; jparm>=0; jparm--)
	    if(m_fit[jparm])
	      {
		for(unsigned iparm=0; iparm<m_nparam; iparm++)
		  std::swap(m_cov(iparm,kparm), m_cov(iparm, jparm));
		for(unsigned iparm=0; iparm<m_nparam; iparm++)
		  std::swap(m_cov(kparm,iparm), m_cov(jparm, iparm));
		kparm--;
	      }

	  for(unsigned iparm=0; iparm<m_nparam; iparm++)
	    m_err(iparm) = sqrt(m_cov(iparm,iparm));

	  return;
	}      

      double chi2 = 0;
      VSAAlgebra::VecND atry(m_nparam);
      for(unsigned iparm=0, jparm=0; iparm<m_nparam; iparm++)
	{
	  if(m_fit[iparm]) 
	    {
	      atry[iparm] = m_a[iparm]+m_da(jparm++);

	      if(NLFitter<Fn>::hasLoBound(iparm) && 
		 atry[iparm] <= NLFitter<Fn>::loBound(iparm))
		atry[iparm] = NLFitter<Fn>::loBound(iparm);
	      else if(NLFitter<Fn>::hasHiBound(iparm) && 
		      atry[iparm] >= NLFitter<Fn>::hiBound(iparm))
		atry[iparm] = NLFitter<Fn>::hiBound(iparm);
	    }
	  else atry[iparm] = m_a[iparm];
	}

      bool ret = false;

      try
	{
	  ret = mrqcof(alpha,m_da,atry,chi2);
	}
      catch(const std::domain_error& e)
	{
	  ret = false;
	}
      
      if(chi2 < m_chi2 && ret)
	{
	  m_lambda *= 0.1;
	  m_chi2 = chi2;	  

	  m_beta = m_da;
	  m_alpha = alpha;
	  m_a = atry;
	}
      else m_lambda *= 10.;
    }

    template<typename Fn> 
    void NLFitterLM<Fn>::
    generateContour(double dchi2, 
		    unsigned ip1, double ip1lo, double ip1hi,
		    unsigned ip2, double ip2lo, double ip2hi,
		    std::vector< std::pair<double,double> >& xy)
    {
      vsassert(ip1 != ip2 && ip1 < nparam() && ip2 < nparam());

      VSAAlgebra::VecND param = m_a;
      VSAAlgebra::VecND err = m_err;

      const double min_chi2 = m_chi2;

      double chi2 = m_chi2;
      double x1 = param(ip1);
      VSAAlgebra::VecND lo = param;
      VSAAlgebra::VecND hi = param;

      std::vector< std::pair<double,double> > xy1;
      std::vector< std::pair<double,double> > xy2;
      std::vector< std::pair<double,double> > xy3;
      std::vector< std::pair<double,double> > xy4;

      double dx1 = 0.01*err(ip1);
      if(ip1lo == ip1hi && dx1 == 0) return;
      else if(dx1 == 0) dx1 = (ip1hi-ip1lo)/100.;
      else dx1 = std::min(dx1,(ip1hi-ip1lo)/100.);

      free();

      while(1)
	{
	  if((x1 < ip1lo || x1 > ip1hi) && ip1lo != ip1hi) break;

	  set(param);
	  hold(ip1,x1);	  
	  initialize();

	  try
	    {
	      fit();
	    }
	  catch(const std::exception& e)
	    {
	      std::cout << e.what() << std::endl;
	    }

	  chi2 = m_chi2;
	  
	  if(chi2 > min_chi2 + dchi2) break;
// 	  else if((x1 <= m_lo_bound[ip1] && m_has_lo_bound[ip1]) || 
// 		  (x1 >= m_hi_bound[ip1] && m_has_hi_bound[ip1])) break;

	  findRoot(m_a,ip1,ip2,min_chi2 + dchi2,xy1,xy2);	  
	  x1 += dx1;
	}

      std::reverse(xy2.begin(),xy2.end());

      x1 = param(ip1) - dx1;

      while(1)
	{
	  if((x1 < ip1lo || x1 > ip1hi) && ip1lo != ip1hi) break;

	  set(param);
	  hold(ip1,x1);
	  initialize();

	  try
	    {
	      fit();
	    }
	  catch(const std::exception& e)
	    {
	      std::cout << e.what() << std::endl;
	    }

	  chi2 = m_chi2;
	  
	  if(chi2 > min_chi2 + dchi2) break;
// 	  else if((x1 <= m_lo_bound[ip1] && m_has_lo_bound[ip1]) || 
// 		  (x1 >= m_hi_bound[ip1] && m_has_hi_bound[ip1])) break;

	  findRoot(m_a,ip1,ip2,min_chi2 + dchi2,xy4,xy3);
	  x1 -= dx1;
	}

      std::reverse(xy4.begin(),xy4.end());

      xy.insert(xy.end(),xy1.begin(),xy1.end());
      xy.insert(xy.end(),xy2.begin(),xy2.end());
      xy.insert(xy.end(),xy3.begin(),xy3.end());
      xy.insert(xy.end(),xy4.begin(),xy4.end());

      m_a = param;
      m_err = err;
      m_chi2 = chi2;
    }

    template<typename Fn> 
    void NLFitterLM<Fn>::findRoot(VSAAlgebra::VecND param,
				  unsigned ip1,
				  unsigned ip2,
				  double fnval,
				  std::vector< std::pair<double,double> >& xy1,
				  std::vector< std::pair<double,double> >& xy2)
    {
      VSAAlgebra::VecND lo = param;
      VSAAlgebra::VecND hi = param;

      hi(ip2) = lo(ip2) - 5*m_err(ip2);

      if(NLFitter<Fn>::hasLoBound(ip2)) 
	hi(ip2) = std::max(hi(ip2),NLFitter<Fn>::loBound(ip2));

      if(m_fn.val(hi) > fnval)
	{
	  try
	    {
	      double x2lo = 
		VSAMath::findBracketedMonotonicRoot(m_fn,lo,hi,ip2,
						    fnval,1E-4);
	      xy1.push_back(std::make_pair(lo(ip1),x2lo));
	    }
	  catch(const VSAMath::FunctionNotMonotonic& x)
	    {
	      std::cout << "Function not monotonic " << std::endl;
	      std::cout << lo << std::endl;
	      std::cout << hi << std::endl;
	    }
	}


      hi(ip2) = lo(ip2) + 5*m_err(ip2);
      
      if(NLFitter<Fn>::hasHiBound(ip2)) 
	hi(ip2) = std::min(hi(ip2),NLFitter<Fn>::hiBound(ip2));
      
      if(m_fn.val(hi) > fnval)
	{
	  try
	    {
	      double x2hi = 
		VSAMath::findBracketedMonotonicRoot(m_fn,lo,hi,ip2,
						    fnval, 1E-4);
	      
	      xy2.push_back(std::make_pair(lo(ip1),x2hi));
	    }
	  catch(const VSAMath::FunctionNotMonotonic& x)
	    {
	      std::cout << "Function not monotonic " << std::endl;
	      std::cout << lo << std::endl;
	      std::cout << hi << std::endl;
	    }
	}
    }

    template<typename Fn> 
    bool NLFitterLM<Fn>::mrqcof(VSAAlgebra::MatrixND& alpha, 
				VSAAlgebra::VecND& beta,
				VSAAlgebra::VecND& a,
				double& chi2)
    {
      alpha.set(0.0);
      beta.set(0.0);

      VSAAlgebra::MatrixND alpha2;
      VSAAlgebra::MatrixND beta2;

      m_fn.val(a,chi2,beta2,alpha2);

      for(unsigned iparm=0, ifit=0; iparm<m_nparam; iparm++)
	{
	  if(m_fit[iparm])
	    {
	      for(unsigned jparm=0, jfit=0; jparm<m_nparam; jparm++)
		if(m_fit[jparm])
		  alpha(ifit,jfit++) = alpha2(iparm,jparm);
	      beta(ifit++) = beta2(iparm,0);
	    }	      
	}

      // Fill in symmetric side -----------------------------------------------
      for(unsigned ifit=1; ifit<m_nfit; ifit++)
	for(unsigned jfit=0; jfit<ifit; jfit++)
	  alpha(jfit,ifit) = alpha(ifit,jfit);

      return true;
    }

    class NLFitterFactory
    {
    public:
      
      template< typename Fn, typename T > 
      static NLFitter< VSAFunction::LnPoissonLikelihood<Fn,T> >* 
      createLnL(const Fn* fn, const Data<T>& data = Data<T>())
      {
	return new NLFitterLM< VSAFunction::LnPoissonLikelihood<Fn,T> >
	  (VSAFunction::LnPoissonLikelihood<Fn,T>(data,fn));
      }

      template< typename Fn, typename T > 
      static NLFitter< VSAFunction::ChiSquared<Fn,T> >* 
      createChiSquared(const Data<T>& data, const Fn* fn)
      {
	return new NLFitterLM< VSAFunction::ChiSquared<Fn,T> >
	  (VSAFunction::ChiSquared<Fn,T>(data,fn));
      }

    };


//     template<typename Fn> 
//     class LikelihoodContourGenerator 
//     {
//     public:
      
//       LikelihoodContourGenerator(NLFitter<Fn>* fitter):
// 	m_fitter(fitter->clone())
//       {

//       }

//       void generateContour(double dchi2, 
// 			   unsigned ip1, unsigned ip2,
// 			   std::vector< std::pair<double,double> >& xy);


//     private:

//       NLFitter<Fn>*     m_fitter;

//       std::vector<double> m_lo_bound;

//     };

//     template<typename Fn> 
//     void LikelihoodContourGenerator<Fn>::
//     generateContour(double dchi2, 
// 		    unsigned ip1, 
// 		    unsigned ip2,
// 		    std::vector< std::pair<double,double> >& xy)
//     {
//       vsassert(ip1 != ip2 && ip1 < nparam() && ip2 < nparam());

//       VSAAlgebra::VecND param = m_a;
//       VSAAlgebra::VecND err = m_err;

//       const double min_chi2 = m_chi2;

//       double chi2 = m_chi2;
//       double x1 = param(ip1);
//       VSAAlgebra::VecND lo = param;
//       VSAAlgebra::VecND hi = param;

//       std::vector< std::pair<double,double> > xy1;
//       std::vector< std::pair<double,double> > xy2;
//       std::vector< std::pair<double,double> > xy3;
//       std::vector< std::pair<double,double> > xy4;

//       double dx1 = 0.01*err(ip1);
//       if(ip1lo != ip1hi) dx1 = std::min(dx1,(ip1hi-ip1lo)/100.);

//       free();

//       while(1)
// 	{
// 	  if((x1 < ip1lo || x1 > ip1hi) && ip1lo != ip1hi) break;

// 	  set(param);
// 	  hold(ip1,x1);	  
// 	  initialize();

// 	  try
// 	    {
// 	      fit();
// 	    }
// 	  catch(const std::exception& e)
// 	    {
// 	      std::cout << e.what() << std::endl;
// 	    }

// 	  chi2 = m_chi2;
	  

// 	  if(chi2 > min_chi2 + dchi2) break;
// 	  else if((x1 <= m_lo_bound[ip1] && m_has_lo_bound[ip1]) || 
// 		  (x1 >= m_hi_bound[ip1] && m_has_hi_bound[ip1])) break;

// 	  findRoot(m_a,ip1,ip2,min_chi2 + dchi2,xy1,xy2);	  
// 	  x1 += dx1;
// 	}

//       std::reverse(xy2.begin(),xy2.end());

//       x1 = param(ip1) - dx1;

//       while(1)
// 	{
// 	  if((x1 < ip1lo || x1 > ip1hi) && ip1lo != ip1hi) break;

// 	  set(param);
// 	  hold(ip1,x1);
// 	  initialize();

// 	  try
// 	    {
// 	      fit();
// 	    }
// 	  catch(const std::exception& e)
// 	    {
// 	      std::cout << e.what() << std::endl;
// 	    }

// 	  chi2 = m_chi2;
	  
// 	  if(chi2 > min_chi2 + dchi2) break;
// 	  else if((x1 <= m_lo_bound[ip1] && m_has_lo_bound[ip1]) || 
// 		  (x1 >= m_hi_bound[ip1] && m_has_hi_bound[ip1])) break;

// 	  findRoot(m_a,ip1,ip2,min_chi2 + dchi2,xy4,xy3);
// 	  x1 -= dx1;
// 	}

//       std::reverse(xy4.begin(),xy4.end());

//       xy.insert(xy.end(),xy1.begin(),xy1.end());
//       xy.insert(xy.end(),xy2.begin(),xy2.end());
//       xy.insert(xy.end(),xy3.begin(),xy3.end());
//       xy.insert(xy.end(),xy4.begin(),xy4.end());
//     }
    

  }
}

#endif // #ifndef VSANONLINEARFITTING_HPP
