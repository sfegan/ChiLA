//-*-mode:c++; mode:font-lock;-*-

/*! \file VSALinearLeastSquares.hpp

  Generalized Linear Least Squares Fitting

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       12/01/2007
*/

#ifndef VSALINEARLEASTSQUARES_HPP
#define VSALINEARLEASTSQUARES_HPP

#include <vector>
#include <stdexcept>
#include <iostream>

#include <VSAAlgebra.hpp>
#include <VSASVD.hpp>
#include <VSAMath.hpp>
#include <VSAData.hpp>
#include <VSAssert.hpp>
#include <VSAFunction.hpp>

#define VSA_ASSERT

namespace VERITAS
{
  namespace VSAMath
  {
    //! Linear least squares fitting by solving normal equations.  Finds
    //! the set of parameter values that minimize chi-squared by solving 
    //! the matrix equation:
    //!
    //! a = (A'*A)^(-1)*A'*y
    //! 
    //! where A is the m x n design matrix and y is the vector of n
    //! data points.  This method can fail if A'*A is singular and
    //! generally the SVD fitter should be preferred for this reason.
    //! However the normal equations method is marginally faster and
    //! uses less memory than the SVD method.
    template<typename Fn, typename T = double> class Fitlin
    {
    public:
      struct Options { /* nothing to see here */ };
      static Options defaultOptions() { return Options(); }

      Fitlin(const Data<T>& data, Fn& fn, 
	     const Options& opt = defaultOptions()):
	m_ndata(data.size()), m_data(data), m_fn(fn), 
	m_nparm(), m_fit(), m_a(), m_cov(), m_chi2()
      {
	VSAAlgebra::VecND atemp;
	m_fn(m_data[0].x, atemp);

	m_nparm = atemp.ndim();
	m_a     .resize(m_nparm); 
	m_fit   .resize(m_nparm,true); 
	m_cov   .resize(m_nparm,m_nparm); 
      }

      void hold(unsigned i, double val) { m_fit[i] = false; m_a[i]=val; }
      void free(unsigned i) { m_fit[i] = true; }
      
      void fit();

      const VSAAlgebra::VecND& param() const { return m_a; }
      const VSAAlgebra::MatrixND& cov() const { return m_cov; }
      double chi2() const { return m_chi2; }

    private:
      unsigned               m_ndata;
      Data<T>                m_data;
      Fn                     m_fn;

      unsigned               m_nparm;
      std::vector<bool>      m_fit;
      VSAAlgebra::VecND      m_a;
      VSAAlgebra::MatrixND   m_cov;
      double                 m_chi2;
    };

    template<typename Fn,typename T> void Fitlin<Fn,T>::fit()
    {
      // Based on NR3 algorithm - Section 15.4
      
      unsigned nfit = 0;
      for(unsigned iparm=0;iparm<m_nparm;iparm++)if(m_fit[iparm])nfit++;
      if(nfit==0)throw
	std::domain_error(std::string(__PRETTY_FUNCTION__)+
			  ": no parameters to be fitted");
      
      VSAAlgebra::MatrixND temp(nfit,nfit,0.0);
      VSAAlgebra::MatrixND beta(nfit,1,0.0);
      
      VSAAlgebra::VecND yfunc(m_nparm);      
      for(unsigned idata=0;idata<m_ndata;idata++)
	{
	  m_fn(m_data[idata].x, yfunc);
	  double ym = m_data[idata].y;
	  double sm = m_data[idata].sigma;
	  if(sm<=0)throw 
	    std::domain_error(std::string(__PRETTY_FUNCTION__) +
			      ": RMS of data point is zero or negative");

	  if(nfit < m_nparm)
	    for(unsigned iparm=0;iparm<m_nparm;iparm++)
	      if(!m_fit[iparm])ym -= m_a[iparm]*yfunc[iparm];

	  double sig2inv = 1.0/(sm*sm);
	  for(unsigned iparm=0, ifit=0; iparm<m_nparm; iparm++)
	    if(m_fit[iparm])
	      {
		const double wt = yfunc[iparm]*sig2inv;
		for(unsigned jparm=0, jfit=0; jparm<m_nparm; jparm++)
		  if(m_fit[jparm])temp(ifit,jfit++) += wt*yfunc[jparm];
		beta(ifit++,0) += wt*ym;
	      }
	}

      for(unsigned ifit=1; ifit<nfit; ifit++)
	for(unsigned jfit=0; jfit<ifit; jfit++)
	  temp(jfit,ifit) = temp(ifit,jfit);
      
      VSAAlgebra::MatrixND::gaussJordan(temp, beta);

      for(unsigned iparm=0, ifit=0; iparm<m_nparm; iparm++)
	if(m_fit[iparm])m_a[iparm] = beta(ifit++,0);

      m_chi2 = 0.0;
      for(unsigned idata=0; idata<m_ndata; idata++)
	{
	  m_fn(m_data[idata].x, yfunc);
	  double sum = 0.0;
	  for(unsigned iparm=0; iparm<m_nparm; iparm++)
	    sum += m_a[iparm]*yfunc[iparm];
	  sum = (m_data[idata].y - sum)/m_data[idata].sigma;
	  m_chi2 += sum*sum;
	}

      for(unsigned ifit=0; ifit<nfit; ifit++)
	for(unsigned jfit=0; jfit<nfit; jfit++)
	  m_cov(ifit, jfit) = temp(ifit, jfit);

      for(unsigned iparm=nfit; iparm<m_nparm; iparm++)
	for(unsigned jparm=0; jparm<iparm+1; jparm++)
	  m_cov(iparm,jparm) = m_cov(jparm,iparm) = 0.0;
      
      unsigned kparm = nfit-1;
      for(int jparm = m_nparm-1; jparm>=0; jparm--)
	if(m_fit[jparm])
	  {
	    for(unsigned iparm=0; iparm<m_nparm; iparm++)
	      std::swap(m_cov(iparm,kparm), m_cov(iparm, jparm));
	    for(unsigned iparm=0; iparm<m_nparm; iparm++)
	      std::swap(m_cov(kparm,iparm), m_cov(jparm, iparm));
	    kparm--;
	  }
    }


    //! Linear least squares fitting using SVD.  When the least squares
    //! solution is not unique the SVD method finds the solution with
    //! minimum norm in the fit parameters.  The solution vector is given by
    //!
    //! a = V*W^(-1)*U'*y
    //!
    //! where U*W*V' = A is the SVD decomposition of the design matrix.
    //! The reciprocal of singular and/or small values in W are set to
    //! 0 in this procedure.
    //!
    template<typename Fn, typename T = double> class Fitsvd
    {
    public:
      typedef double Options;
      static Options defaultOptions() { return 1e-12; }

      Fitsvd(const Data<T>& data, Fn& fn, 
	     const double& tol = defaultOptions()):
	m_ndata(data.size()), m_data(data), m_fn(fn), m_tol(tol),
	m_nparm(), m_a(), m_cov(), m_chi2()
      {
	VSAAlgebra::VecND atemp;
	m_fn(m_data[0].x, atemp);

	m_nparm = atemp.ndim();
	m_a     .resize(m_nparm); 
	m_cov   .resize(m_nparm,m_nparm); 
      }

      double tol() const { return m_tol; }
      void setTol(const double& tol) { m_tol = tol; }

      void fit()
      {
	VSAAlgebra::MatrixND aa;
	VSAAlgebra::VecND b;
	initialize(aa,b);
	fit(aa,b,m_a,m_cov,m_tol);
	m_chi2 = chi2(aa,b,m_a);
      }

      static void fit(const VSAAlgebra::MatrixND& aa,
		      const VSAAlgebra::VecND& b,
		      VSAAlgebra::VecND& a,
		      VSAAlgebra::MatrixND& cov,
		      double tol = 1E-12);

      const VSAAlgebra::VecND& param() const { return m_a; }
      const VSAAlgebra::MatrixND& cov() const { return m_cov; }
      double chi2() const { return m_chi2; }

      void initialize(VSAAlgebra::MatrixND& aa, VSAAlgebra::VecND& b)
      {
	aa = VSAAlgebra::MatrixND(m_ndata,m_nparm);
	b = VSAAlgebra::VecND(m_ndata);
	VSAAlgebra::VecND yfunc(m_nparm);      
	for(unsigned idata=0;idata<m_ndata;idata++)
	  {
	    m_fn(m_data[idata].x, yfunc);
	    const double ym = m_data[idata].y;
	    const double sm = m_data[idata].sigma;
	    if(sm<=0)
	      throw 
		std::domain_error(std::string(__PRETTY_FUNCTION__) +
				  ": RMS of data point is zero or negative");
	    
	    const double siginv = 1.0/sm;
	    for(unsigned iparm=0; iparm<m_nparm; iparm++)
	      aa(idata,iparm) = yfunc(iparm)*siginv;
	    b[idata] = ym*siginv;
	  }
      }

      static double chi2(const VSAAlgebra::MatrixND& aa, 
			 const VSAAlgebra::VecND& b,
			 const VSAAlgebra::VecND& a)
      {
	const unsigned ndata = aa.nrow();
	const unsigned nparm = aa.ncol();
	double chi2 = 0.0;
	for(unsigned idata=0; idata<ndata; idata++)
	{
	  double sum = 0.0;
	  for(unsigned iparm=0; iparm<nparm; iparm++)
	    sum += aa(idata,iparm)*a[iparm];
	  sum -= b[idata];
	  chi2 += sum*sum;
	}

	return chi2;
      }

    private:

      unsigned               m_ndata;
      Data<T>                m_data;
      Fn                     m_fn;
      double                 m_tol;

      unsigned               m_nparm;
      VSAAlgebra::VecND      m_a;
      VSAAlgebra::MatrixND   m_cov;
      double                 m_chi2;
    };

    template<typename Fn, typename T> 
    void Fitsvd<Fn,T>::fit(const VSAAlgebra::MatrixND& aa,
			   const VSAAlgebra::VecND& b,
			   VSAAlgebra::VecND& a,
			   VSAAlgebra::MatrixND& cov,
			   double tol)
    {
      // Based on NR3 algorithm - Section 15.4
      const unsigned nparm = aa.ncol();

      VSAAlgebra::SVD svd(aa);
      double thresh = (tol > 0 ? tol*svd.w()[0] : -1);
      svd.solve(b,a,thresh);

      for(unsigned iparm=0; iparm<nparm; iparm++)
	for(unsigned jparm=0; jparm<iparm+1; jparm++)
	  {
	    double sum = 0;
	    for(unsigned kparm=0; kparm<nparm; kparm++)
	      {
		const double wk = svd.w()[kparm];
		if(wk > svd.thresh())
		  sum += svd.v()(iparm,kparm)*svd.v()(jparm,kparm)/(wk*wk);
	      }
	    cov(jparm,iparm) = cov(iparm,jparm) = sum;
	  }
    }
    
    //! Constrained linear least squares with SVD.  Solves linear
    //! least squares problem given a system of constraint equations
    //! (equality or inequality).
    template<typename Fn, typename T = double> class Fitcsvd
    {
    public:
      typedef double Options;
      static Options defaultOptions() { return 1e-12; }
      
      Fitcsvd(const Data<T>& data, const Fn& fn, 
	      const double& tol = defaultOptions()):
    	m_ndata(data.size()), m_data(data), m_fn(fn), m_tol(tol),
    	m_nparm(), m_a(), m_cov(), m_chi2(), m_nconst()
      {
    	VSAAlgebra::VecND atemp;
    	m_fn(m_data[0].x, atemp);
	
    	m_nparm = atemp.ndim();
    	m_a     .resize(m_nparm); 
    	m_cov   .resize(m_nparm,m_nparm); 
      }
      
      const VSAAlgebra::VecND& param() const { return m_a; }
      const VSAAlgebra::MatrixND& cov() const { return m_cov; }
      double chi2() const { return m_chi2; }

      //! Return number of active constraints
      unsigned nconst() const { return m_nconst; }

      void initialize(VSAAlgebra::MatrixND& aa, VSAAlgebra::VecND& b)
      {
	aa = VSAAlgebra::MatrixND(m_ndata,m_nparm);
	b = VSAAlgebra::VecND(m_ndata);
	VSAAlgebra::VecND yfunc(m_nparm);      
	for(unsigned idata=0;idata<m_ndata;idata++)
	  {
	    m_fn(m_data[idata].x, yfunc);
	    const double ym = m_data[idata].y;
	    const double sm = m_data[idata].sigma;
	    if(sm<=0)
	      throw 
		std::domain_error(std::string(__PRETTY_FUNCTION__) +
				  ": RMS of data point is zero or negative");
	    
	    const double siginv = 1.0/sm;
	    for(unsigned iparm=0; iparm<m_nparm; iparm++)
	      aa(idata,iparm) = yfunc(iparm)*siginv;
	    b[idata] = ym*siginv;
	  }
      }

      //! Add a sequential set of constraints on the derivative which will
      //! constrain the function to be monotonically positive/negative.
      void addMonotonicConstraints(T lo, T step, unsigned n)
      {
	for(unsigned i = 0; i < n; i++)
	  addDerivativeConstraint(lo+i*step);
      }

      //! Add a constraint on the derivative to be positive/negative at x.
      void addDerivativeConstraint(T x, bool gtr = true)
      {
	VSAAlgebra::VecND dydx;
	m_fn.dydx(x,dydx);

	const unsigned nrow = m_g.nrow();

	VSAAlgebra::MatrixND g(nrow+1,m_nparm);
	VSAAlgebra::VecND h(nrow+1);
	g.setSubMatrix(0,0,m_g);
	h.setSubVector(0,m_h);

	for(unsigned icol = 0; icol < g.ncol(); icol++)
	  {
	    if(gtr) g(nrow,icol) = dydx(icol);
	    else g(nrow,icol) = -dydx(icol);
	  }

	h(nrow) = 0;

	m_g = g;
	m_h = h;
      }

      void setInequalityConstraint(const VSAAlgebra::MatrixND& g,
				   const VSAAlgebra::VecND& h)
      {
	vsassert(g.nrow() == h.ndim());
	m_g = g;
	m_h = h;
      }

      void fit()
      {
	VSAAlgebra::MatrixND aa;
	VSAAlgebra::VecND b;
	initialize(aa,b);

	m_nconst = 0;

	// If the set of constraint equations is empty then compute
	// lls using svd fitter
	if(m_g.nrow() == 0) Fitsvd<Fn,T>::fit(aa,b,m_a,m_cov);
	else
	  {
	    if(!icls(aa,b,m_g,m_h,m_a,m_cov,m_nconst))
	      throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ 
				       ": Error in icls");
	  }

	m_chi2 = Fitsvd<Fn,T>::chi2(aa,b,m_a);
      }

      //! Non-negative least squares algorithm taken from Lawson and
      //! Hanson, Solving Least Squares Problems, Chap. 23, pg. 160-165
      //! 
      //! Solves for min ||A*x-b|| with constraint x >= 0
      //!
      static void nnls(const VSAAlgebra::MatrixND& aa,
		       const VSAAlgebra::VecND& b,
		       VSAAlgebra::VecND& a, 
		       VSAAlgebra::VecND& lambda, 
		       double tol = 1E-12);

      //! Least Distance algorithm taken from Lawson and Hanson,
      //! Solving Least Squares Problems, Chap. 23, pg. 165-167
      //!
      //! Solves for min ||x|| with constraint G*x >= h
      //!
      static bool lsd(const VSAAlgebra::MatrixND& g,
		      const VSAAlgebra::VecND& h,
		      VSAAlgebra::VecND& a);

      //! Inequality Constrained least squares algorithm 
      //!
      //! Solves for min ||A*x-b|| with constraint G*x >= h
      //!
      //! @param aa     Design matrix
      //! @param b      Data vector
      //! @param g      Constraint matrix
      //! @param h      Constraint vector
      //! @param a      LS Solution vector
      //! @param cov    Covariance matrix of LS solution
      //! @param nconst Number of active constraints at the LS solution
      static bool icls(const VSAAlgebra::MatrixND& aa,
		       const VSAAlgebra::VecND& b,
		       const VSAAlgebra::MatrixND& g,
		       const VSAAlgebra::VecND& h,
		       VSAAlgebra::VecND& a,
		       VSAAlgebra::MatrixND& cov,
		       unsigned& nconst);

    private:
      unsigned               m_ndata;
      Data<T>                m_data;
      Fn                     m_fn;
      double                 m_tol;

      VSAAlgebra::MatrixND   m_g;
      VSAAlgebra::VecND      m_h;

      unsigned               m_nparm;
      VSAAlgebra::VecND      m_a;
      VSAAlgebra::MatrixND   m_cov;
      double                 m_chi2;
      unsigned               m_nconst;
    };
    
    template<typename Fn, typename T> 
    void Fitcsvd<Fn,T>::nnls(const VSAAlgebra::MatrixND& aa,
			     const VSAAlgebra::VecND& b,
			     VSAAlgebra::VecND& a, 
			     VSAAlgebra::VecND& lambda, 
			     double tol)
    {
      const unsigned nparm = aa.ncol();
      tol = 10*tol*VSAAlgebra::SVD::norm(aa);
      
      std::set<unsigned> zset; // set of indices for negative coefficients
      std::set<unsigned> pset; // set of indices for positive coefficients
      for(unsigned ip = 0; ip < nparm; ip++) zset.insert(ip);
      VSAAlgebra::VecND x(nparm);

      // dual vector (lagrange multipliers)
      VSAAlgebra::VecND w = aa.getTranspose()*(b-aa*x); 

      const unsigned nmax = 3*nparm;
      unsigned n = 0;

      while(!zset.empty())
	{
#if 0	 
	  std::cout << "---------------------------" << std::endl;
	  std::cout << "N " << n << std::endl << "x = "
		    << std::endl << x << std::endl;

	  std::cout << "z = ";
	  for(typename std::set<unsigned>::iterator itr = zset.begin(); 
		itr != zset.end(); ++itr)
	    std::cout << *itr << " ";
	  std::cout << std::endl;

	  std::cout << "p = ";
	  for(typename std::set<unsigned>::iterator itr = pset.begin(); 
		itr != pset.end(); ++itr)
	    std::cout << *itr << " ";
	  std::cout << std::endl;
#endif

	  // Exit if the number of iterations exceeds the maximum
	  if(n > nmax)
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": " 
			<< "Exceeded max iterations." << std::endl;
	      break;
	    }
	  n++;
      
	  // Find the largest positive value in w
	  int tmax = 0;
	  double wmax = 0;
	  for(unsigned ip = 0; ip < nparm; ip++)
	    for(typename std::set<unsigned>::iterator itr = zset.begin(); 
		itr != zset.end(); ++itr)
	    {
	      if(wmax < w(ip))
		{
		  tmax = ip;
		  wmax = w(ip);
		}
	    }

#if 0	
	  std::cout << "w = " << std::endl << w << std::endl;	  
	  std::cout << "tmax " << tmax << " " << w(tmax) << std::endl;
	  std::cout << "tol " << tol << std::endl;
#endif

	  // ------------------------------------------------------------------
	  // Exit when all values in w are negative or equal to zero within
	  // the tolerance value
	  if(wmax < tol) break;

	  // Move basis vector with largest w from negative to positive set
	  zset.erase(tmax);
	  pset.insert(tmax);
	  	 
#if 0	  
	  std::cout << "z = ";
	  for(typename std::set<unsigned>::iterator itr = zset.begin(); 
		itr != zset.end(); ++itr)
	    std::cout << *itr << " ";
	  std::cout << std::endl;

	  std::cout << "p = ";
	  for(typename std::set<unsigned>::iterator itr = pset.begin(); 
		itr != pset.end(); ++itr)
	    std::cout << *itr << " ";
	  std::cout << std::endl;
#endif	  

	  unsigned it = 0;
	  while(1)
	    {
	      // --------------------------------------------------------------
	      // Compute the LLS solution in the subspace of basis vectors
	      // in the positive set
	      VSAAlgebra::MatrixND ep(aa.nrow(),nparm);
	  
	      for(typename std::set<unsigned>::iterator itr = pset.begin(); 
		  itr != pset.end(); ++itr)
		for(unsigned irow = 0; irow < ep.nrow(); irow++)
		  ep(irow,*itr) = aa(irow,*itr);	       
	     
	      VSAAlgebra::VecND z(nparm);

	      VSAAlgebra::SVD svd(ep);
	      svd.solve(b,z);

	      for(typename std::set<unsigned>::iterator itr = zset.begin(); 
		  itr != zset.end(); ++itr)
		z(*itr) = 0;

#if 0
	      std::cout << "zt = " << z(tmax) << std::endl;
	      std::cout << "z = " << std::endl << z << std::endl;
#endif

	      bool zp = true;
	      // Check that all elements in pset are positive
	      for(typename std::set<unsigned>::iterator itr = pset.begin(); 
		  itr != pset.end(); ++itr)
		{
		  if(z(*itr) < 0)
		    {
		      zp = false;
		      break;
		    }
		}

	      //	      std::cout << "zp " << zp << std::endl;

	      if(zp)
		{
		  x = z;
		  break;
		}	      

	      int q = -1;
	      double qmin = 0;
	      
	      for(typename std::set<unsigned>::iterator itr = pset.begin(); 
		  itr != pset.end(); ++itr)
		{
		  if(z(*itr) > 0) continue;
		  
		  double qtmp = x(*itr)/(x(*itr)-z(*itr));
		  if(q == -1 || qtmp < qmin)
		    {
		      q = *itr;
		      qmin = qtmp;
		    }
		}
	      

	      double alpha = x(q)/(x(q)-z(q));
	      x += alpha*(z-x);

	      // std::cout << "q = " << q << " " << qmin << std::endl;
	      // std::cout << n << " alpha = " << alpha << std::endl;
	      // std::cout << "x = " << std::endl << x << std::endl;

	      // Move any negative or zero elements from the pset to the
	      // zset
	      for(typename std::set<unsigned>::iterator itr = pset.begin(); 
		  itr != pset.end(); )
		{
		  if(x(*itr) <= std::numeric_limits<double>::epsilon())
		    {
		      zset.insert(*itr);
		      pset.erase(itr++);
		    }
		  else ++itr;
		}
	      it++;
	    }
	  
	  w = aa.getTranspose()*(b-aa*x);
	}

      a = x;
      lambda = w;
    }

    template<typename Fn, typename T> 
    bool Fitcsvd<Fn,T>::lsd(const VSAAlgebra::MatrixND& g,
			    const VSAAlgebra::VecND& h,
			    VSAAlgebra::VecND& a)
    {
      const unsigned nparm = g.ncol();

      a = VSAAlgebra::VecND(nparm);

      VSAAlgebra::MatrixND gt = g.getTranspose();
      VSAAlgebra::MatrixND e(gt.nrow()+1,gt.ncol());
      VSAAlgebra::VecND f(e.nrow());

      for(unsigned irow = 0; irow < e.nrow(); irow++)
	{
	  if(irow + 1 == e.nrow()) f(irow) = 1;
	  for(unsigned icol = 0; icol < e.ncol(); icol++)
	    {
	      if(irow + 1 == e.nrow()) e(irow,icol) = h(icol);
	      else e(irow,icol) = gt(irow,icol);
	    }
	}

      // ----------------------------------------------------------------------
      // Solve lsd problem by solving the equivalent nnls problem
      VSAAlgebra::VecND u,lambda;
      nnls(e,f,u,lambda);

      VSAAlgebra::VecND r = e*u-f;
      double rnorm = r*r;

      // std::cout << "e = " << std::endl << e << std::endl;
      // std::cout << "f = " << std::endl << f << std::endl;
      // std::cout << "u = " << std::endl << u << std::endl;
      // std::cout << "rnorm = " << rnorm << std::endl;
      // std::cout << "r = " << std::endl << r << std::endl;
      // std::cout << "lambda = " << std::endl << lambda << std::endl;

      // ----------------------------------------------------------------------
      // If ||r|| <= 0 then the constraint equations are incompatible
      // and no solution is returned.
      // ----------------------------------------------------------------------
      if(rnorm <= 1E-12) return false;
      for(unsigned ip = 0; ip < a.ndim(); ip++) a(ip) = -r(ip)/r(a.ndim());
      
      return true;
    }

    template<typename Fn, typename T> 
    bool Fitcsvd<Fn,T>::icls(const VSAAlgebra::MatrixND& aa,
			     const VSAAlgebra::VecND& b,
			     const VSAAlgebra::MatrixND& g,
			     const VSAAlgebra::VecND& h,
			     VSAAlgebra::VecND& a,
			     VSAAlgebra::MatrixND& cov,
			     unsigned& nconst)
    {
      vsassert(aa.nrow() == b.ndim() && g.nrow() == h.ndim());
      const unsigned nparm = aa.ncol();

      // std::cout << "ICLS" << std::endl;
      // std::cout << aa << std::endl;
      // std::cout << b << std::endl;
      // std::cout << g << std::endl;
      // std::cout << h << std::endl;


      a = VSAAlgebra::VecND(nparm);
      cov = VSAAlgebra::MatrixND(nparm,nparm);

      // Perform SVD decomposition of the design matrix
      VSAAlgebra::SVD svd(aa);

      VSAAlgebra::MatrixND s = VSAAlgebra::MatrixND::makeDiagonal(svd.w());
      VSAAlgebra::MatrixND u = svd.u();
      VSAAlgebra::MatrixND v = svd.v();

#if 0
      std::cout << "s = " << std::endl << s << std::endl;
      std::cout << "v = " << std::endl << v << std::endl;
      std::cout << "u = " << std::endl << u << std::endl;
#endif

      // Convert to lsd problem by change of variables
      VSAAlgebra::MatrixND si = s.inverse();
      VSAAlgebra::VecND f = u.getTranspose()*b;
      VSAAlgebra::VecND f1 = f.subVector(0,nparm);
      VSAAlgebra::MatrixND gg = g*v*si;
      VSAAlgebra::VecND hh = h - gg*f1;

#if 0
      std::cout << "si = " << std::endl << si << std::endl;
      std::cout << "f = " << std::endl << f << std::endl;
      std::cout << "gg = " << std::endl << gg << std::endl;
      std::cout << "hh = " << std::endl << hh << std::endl;
#endif

      VSAAlgebra::VecND z;

      // If lsd returns false then one or more of the constraint equations
      // are incompatible
      if(!lsd(gg,hh,z)) return false;

      // std::cout << "z = " << z << std::endl;

      a = v*si*(z+f1);

      // Covariance matrix of unconstrained solution
      cov = v*si*si*v.getTranspose();

      VSAAlgebra::VecND gr = g*a-h;
      VSAAlgebra::MatrixND d;

      // std::cout << "a " << std::endl << a << std::endl;
      // std::cout << "g " << std::endl << g << std::endl;
      // std::cout << "gr " << std::endl << gr << std::endl;

      // Determine the set of constraints that are in effect at the LS
      // solution point and populate the matrix d with this subset
      nconst = 0;
      for(unsigned i = 0; i < gr.ndim(); i++)
	if(gr(i) <= 1E-12)
	  {
	    nconst++;
	    d.concatenate(g.subMatrix(i,0,1,g.ncol()));
	  }      

      if(d.nrow() == 0) return true;

      VSAAlgebra::MatrixND vd = d*cov*d.getTranspose();
      VSAAlgebra::SVD svd2(vd);
      vd = svd2.inverse();

      VSAAlgebra::MatrixND cov2 = d.getTranspose()*vd*d*cov;
      //      VSAAlgebra::MatrixND cov2 = cov*d.getTranspose()*vd*d*cov;
      VSAAlgebra::MatrixND i = 
	VSAAlgebra::MatrixND::makeDiagonal(VSAAlgebra::VecND(nparm,1));
      cov = cov*(i-cov2);
      return true;
    }

    class PolyFit
    {
    public:

      // ----------------------------------------------------------------------
      // Functor used by fitter to calculate value of basis functions at data
      // x-values. Returns value of th n+1 polynomial functions (x^0 .. x^n)
      // in a VSAAlgebra::VecND
      // ----------------------------------------------------------------------

      class Fn
      {
      public:
	Fn(unsigned n): m_n(n+1) { /* nothing to see here */ }
	void operator() (double x, VSAAlgebra::VecND& v)
	{
	  v.resize(m_n);
	  double s = v[0] = 1.0;
	  for(unsigned i=1;i<m_n;i++)v[i] = s *= x;
	}
      private:
	unsigned m_n;
      };

      // ----------------------------------------------------------------------
      // Static interface class to free user from need to manually construct 
      // the fitter, call it and extract the results. Simply call
      // "PolyFit::fit" and let it do all the work
      // ----------------------------------------------------------------------

      typedef Fitsvd<Fn> PreferredFitter;

      static double fit(unsigned n, const Data<double>& data, 
			VSAAlgebra::VecND& param,
			VSAAlgebra::MatrixND* cov = 0, 
			const PreferredFitter::Options& opt =
			PreferredFitter::defaultOptions());

      static VSAAlgebra::VecND fit(unsigned n, const Data<double>& data)
      {
	VSAAlgebra::VecND param;
	fit(n, data, param);
	return param;
      }

      static double val(const VSAAlgebra::VecND& a, double x);
      static double var(const VSAAlgebra::VecND& a, 
			const VSAAlgebra::MatrixND& cov,
			double x);

      static void val(const VSAAlgebra::VecND& a, const std::vector<double> x,
		      std::vector<double>& y);

      static std::vector<double> 
      val(const VSAAlgebra::VecND& a, const std::vector<double> x)
      {
	std::vector<double> y;
	val(a,x,y);
	return y;
      }

      static void differentiate(const VSAAlgebra::VecND& a, 
				VSAAlgebra::VecND& dydx_a);

      static VSAAlgebra::VecND differentiate(const VSAAlgebra::VecND& a)
      {
	VSAAlgebra::VecND dydx_a;
	differentiate(a, dydx_a);
	return dydx_a;
      }

      static void dyda(const VSAAlgebra::VecND& a, 
		       double x, VSAAlgebra::VecND& dyda);

    };
  }
}

#endif // #ifndef VSALINEARLEASTSQUARES_HPP
