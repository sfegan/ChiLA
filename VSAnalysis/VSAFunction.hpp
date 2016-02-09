//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAFunction.hpp
  Various predefined functions.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#ifndef VSAFUNCTION_HPP
#define VSAFUNCTION_HPP

#include <typeinfo>
#include <sstream>

#include <VSAAlgebra.hpp>
#include <VSAData.hpp>
#include <VSACoord.hpp>
#include <VSSimple2DHist.hpp>
#include <VSNSpace.hpp>

namespace VERITAS
{
  namespace VSAFunction
  {
    template< typename T >
    class BaseFn
    {
    public:
      BaseFn() { }
      virtual ~BaseFn() { }
 
      // Accessors ------------------------------------------------------------
      virtual unsigned nparm() const = 0;
      virtual VSAAlgebra::VecND param() const  = 0;
      virtual double param(unsigned ip) const = 0;
      virtual std::vector<bool> fixed() const  = 0;
      virtual bool fixed(unsigned ip) const = 0;

      // Evaluation -----------------------------------------------------------
      virtual double val(const T& x) const = 0;
      virtual double val(const T& x, const VSAAlgebra::VecND& a) const = 0;
      virtual void dyda(const T& x, VSAAlgebra::VecND& dyda) const = 0;
      virtual void dyda(const T& x, const VSAAlgebra::VecND& a,
			VSAAlgebra::VecND& dyda) const = 0;

      // Setters --------------------------------------------------------------
      virtual void setParam(const VSAAlgebra::VecND& a) = 0; 
      virtual void setParam(const VSAAlgebra::VecND& a, 
			    const std::vector< bool >& fixed) = 0;
      virtual void setParam(unsigned ip, double a) = 0; 
      virtual void fixParam(unsigned ip, bool fixed = true) = 0; 
      virtual void fixParam(bool fixed = true) = 0;

      // Virtual Constructor --------------------------------------------------
      virtual BaseFn<T>* clone() const = 0;
    };

    template< typename T >
    class ParamFn : public BaseFn<T>
    {
    public:
      ParamFn(): BaseFn<T>(), m_a(), m_fixed() { }
      ParamFn(unsigned nparm): m_a(nparm), m_fixed(nparm,false) { }
      virtual ~ParamFn() { }

      // Accessors ------------------------------------------------------------
      virtual unsigned nparm() const { return m_a.ndim(); }
      virtual VSAAlgebra::VecND param() const { return m_a; }
      virtual double param(unsigned ip) const { return m_a[ip]; }
      virtual std::vector<bool> fixed() const { return m_fixed; }
      virtual bool fixed(unsigned ip) const { return m_fixed[ip]; }

      // Evaluation -----------------------------------------------------------
      virtual double val(const T& x) const { return 0; }
      virtual double val(const T& x, const VSAAlgebra::VecND& a) const
      { return 0; }
      virtual void dyda(const T& x, VSAAlgebra::VecND& dyda) const { }
      virtual void dyda(const T& x, const VSAAlgebra::VecND& a,
			VSAAlgebra::VecND& dyda) const { }

      // Setters --------------------------------------------------------------
      virtual void setNParam(unsigned n)
      { m_a = VSAAlgebra::VecND(n); m_fixed = std::vector<bool>(n,false); }
      virtual void setParam(const VSAAlgebra::VecND& a) 
      { m_a = a; m_fixed.resize(a.ndim(),false); }
      virtual void setParam(const VSAAlgebra::VecND& a, 
			    const std::vector< bool >& fixed) 
      { m_a = a; m_fixed = fixed; }
      virtual void setParam(unsigned ip, double a) { m_a[ip] = a; }
      virtual void fixParam(unsigned ip, bool fixed = true)
      { m_fixed[ip] = fixed; }
      virtual void fixParam(bool fixed = true) 
      { m_fixed.clear(); m_fixed.assign(nparm(),fixed); }

      // Virtual Constructor --------------------------------------------------
      virtual ParamFn<T>* clone() const { return new ParamFn<T>(*this); }

    private:
      VSAAlgebra::VecND      m_a;
      std::vector<bool>      m_fixed;
    };

    // ========================================================================
    // CompositeFn
    // ========================================================================

    template< typename T >
    class CompositeFn : public BaseFn<T>
    {
    public:
      CompositeFn(): m_fn1(), m_fn2() { }
      CompositeFn(BaseFn<T>* fn1, BaseFn<T>* fn2): 
	m_fn1(fn1->clone()), m_fn2(fn2->clone()) { }

      ~CompositeFn() 
      { 
	delete m_fn1;
	delete m_fn2;
      }
      
      // Accessors ------------------------------------------------------------
      unsigned nparm() const 
      { 
	return m_fn1->nparm() + m_fn2->nparm();
      }

      VSAAlgebra::VecND param() const
      {
	VSAAlgebra::VecND p(nparm());
	p.set(m_fn1->param(),0);
	p.set(m_fn2->param(),m_fn1->nparm());
	return p;
      }

      double param(unsigned ip) const 
      { 
	if(ip < m_fn1->nparm()) return m_fn1->param(ip);
	else return m_fn2->param(ip-m_fn1->nparm());
      }

      std::vector<bool> fixed() const 
      { 
	std::vector<bool> fixed;
	
	const unsigned nparm1 = m_fn1->nparm();
	for(unsigned i = 0; i < nparm1; i++)
	  fixed.push_back(m_fn1->fixed(i));

	const unsigned nparm2 = m_fn2->nparm();
	for(unsigned i = 0; i < nparm2; i++)
	  fixed.push_back(m_fn2->fixed(i));

	return fixed;
      }

      bool fixed(unsigned ip) const 
      { 
	if(ip < m_fn1->nparm()) return m_fn1->fixed(ip);
	else return m_fn2->fixed(ip-m_fn1->nparm());
      }

      BaseFn<T>* getFn1() { return m_fn1; }
      BaseFn<T>* getFn2() { return m_fn2; }

      // Evaluation -----------------------------------------------------------
      virtual double val(const T& x) const = 0;
      virtual double val(const T& x, const VSAAlgebra::VecND& a) const = 0;
      virtual void dyda(const T& x, VSAAlgebra::VecND& _dyda) const = 0;
      virtual void dyda(const T& x, const VSAAlgebra::VecND& a, 
			VSAAlgebra::VecND& _dyda) const = 0;

      // Setters --------------------------------------------------------------
      void setParam(const VSAAlgebra::VecND& a)
      {
	vsassert(a.ndim() == nparm());
	m_fn1->setParam(a.subVector(0,m_fn1->nparm()));
	m_fn2->setParam(a.subVector(m_fn1->nparm(),nparm()-m_fn1->nparm()));
      }
      
      void setParam(const VSAAlgebra::VecND& a, 
		    const std::vector< bool >& fixed)
      {
	vsassert(a.ndim() == nparm());

	std::vector< bool > tmp1(fixed.begin(),fixed.begin()+m_fn1->nparm());
	std::vector< bool > tmp2(fixed.begin()+m_fn1->nparm(),fixed.end());
	
	VSAAlgebra::VecND a1 = a.subVector(0,m_fn1->nparm());
	VSAAlgebra::VecND a2 = a.subVector(m_fn1->nparm(),
					   nparm()-m_fn1->nparm());

	m_fn1->setParam(a1,tmp1);
	m_fn1->setParam(a2,tmp2);
      }
    
      void setParam(unsigned ip, double a)
      {
	vsassert(ip < nparm());

	if(ip < m_fn1->nparm()) m_fn1->setParam(ip,a);
	else m_fn2->setParam(ip-m_fn1->nparm(),a);
      }

      void fixParam(unsigned ip, bool fixed = true)
      {
	vsassert(ip < nparm());

	if(ip < m_fn1->nparm()) m_fn1->fixParam(ip,fixed);
	else m_fn2->fixParam(ip-m_fn1->nparm(),fixed);
      }

      void fixParam(bool fixed = true)
      {
	m_fn1->fixParam(fixed);
	m_fn2->fixParam(fixed);
      }

      void setFn1(BaseFn<T>* fn) 
      {       
	delete m_fn1; 
	m_fn1 = fn;
      }
    
      void setFn2(BaseFn<T>* fn) 
      {       
	delete m_fn2; 
	m_fn2 = fn;
      }

      // Virtual Constructor --------------------------------------------------
//       virtual CompositeFn<T>* clone() const
//       {
// 	return new CompositeFn<T>(*this);
//       }

      // Copy-Constructor -----------------------------------------------------
      CompositeFn(const CompositeFn& o)
      {
	m_fn1 = o.m_fn1->clone();
	m_fn2 = o.m_fn2->clone();
      }

      // Assignment Operator --------------------------------------------------
      CompositeFn& operator= (const CompositeFn& o)
      {
	if(this == &o) return *this;

	delete m_fn1;
	delete m_fn2;
	m_fn1 = o.m_fn1->clone();
	m_fn2 = o.m_fn2->clone();

	return *this;
      }

    protected:
      BaseFn<T>* m_fn1;
      BaseFn<T>* m_fn2;
    };

    template< typename T >
    class CompositeSumFn : public CompositeFn<T>
    {
    public:
      CompositeSumFn(): CompositeFn<T>() { }
      CompositeSumFn(BaseFn<T>* fn1, BaseFn<T>* fn2):
	CompositeFn<T>(fn1,fn2) { }
      ~CompositeSumFn() {}

      // Evaluation -----------------------------------------------------------
      double val(const T& x) const 
      { 
	return CompositeFn<T>::m_fn1->val(x) + CompositeFn<T>::m_fn2->val(x);
      }

      double val(const T& x, const VSAAlgebra::VecND& a) const
      { 
	return 
	  CompositeFn<T>::
	  m_fn1->val(x,a.subVector(0,CompositeFn<T>::m_fn1->nparm())) + 
	  CompositeFn<T>::
	  m_fn2->val(x,a.subVector(CompositeFn<T>::m_fn1->nparm(),
				   CompositeFn<T>::nparm()-
				   CompositeFn<T>::m_fn1->nparm()));
      }
	
      void dyda(const T& x, VSAAlgebra::VecND& _dyda) const 
      {  
	_dyda.resize(CompositeFn<T>::nparm());

	VSAAlgebra::VecND tmp1;
	VSAAlgebra::VecND tmp2;

	CompositeFn<T>::m_fn1->dyda(x,tmp1);
	CompositeFn<T>::m_fn2->dyda(x,tmp2);

	_dyda.set(tmp1,0);
	_dyda.set(tmp2,CompositeFn<T>::m_fn1->nparm());
      }

      void dyda(const T& x, const VSAAlgebra::VecND& a, 
		VSAAlgebra::VecND& _dyda) const 
      {  
	_dyda.resize(CompositeFn<T>::nparm());

	VSAAlgebra::VecND dyda1;
	VSAAlgebra::VecND dyda2;

	VSAAlgebra::VecND a1 = a.subVector(0,CompositeFn<T>::m_fn1->nparm());
	VSAAlgebra::VecND a2 = a.subVector(CompositeFn<T>::m_fn1->nparm(),
					   CompositeFn<T>::nparm()-
					   CompositeFn<T>::m_fn1->nparm());
	  
	CompositeFn<T>::m_fn1->dyda(x,a1,dyda1);
	CompositeFn<T>::m_fn2->dyda(x,a2,dyda2);
	_dyda.set(dyda1,0);
	_dyda.set(dyda2,CompositeFn<T>::m_fn1->nparm());
      }

      // Virtual Constructor --------------------------------------------------
      virtual CompositeSumFn<T>* clone() const
      {
	return new CompositeSumFn<T>(*this);
      }

    };

    // ========================================================================
    // MemberFn
    // ========================================================================
    template< class FN, class T >
    class MemberFn
    {
    public:

      typedef double (FN::* FnPtr) (const T&) const;
      
      MemberFn(const FN* fn, FnPtr pmf):
	m_fn(fn->clone()),m_pmf(pmf)  
      {}
      
      virtual ~MemberFn() { delete m_fn; }
      
      virtual double val(const T& x) const
      {
	return (m_fn->*m_pmf)(x);
      }
      
      virtual double operator() (const T& x) const { return val(x); }

      // Copy-Constructor -----------------------------------------------------
      MemberFn(const MemberFn& o)
      {
	m_fn = o.m_fn->clone();
	m_pmf = o.m_pmf;
      }

    private:
      
      FN*           m_fn;
      FnPtr         m_pmf;            
    };

    template< class FN, class T >
    class ParamMemberFn
    {
    public:

      typedef double (FN::* FnPtr) 
	(const T&, const VSAAlgebra::VecND& ) const;
      
      ParamMemberFn(const FN* fn, FnPtr pmf, const VSAAlgebra::VecND& a):
	m_fn(fn->clone()),m_pmf(pmf), m_a(a)
      {}
      
      virtual ~ParamMemberFn() { delete m_fn; }
      
      virtual double val(const T& x) const
      {
	return (m_fn->*m_pmf)(x,m_a);
      }
      
      // Copy-Constructor -----------------------------------------------------
      ParamMemberFn(const ParamMemberFn& o)
      {
	m_fn = o.m_fn->clone();
	m_pmf = o.m_pmf;
	m_a = o.m_a;
      }

    private:
      
      FN*           m_fn;
      FnPtr         m_pmf;      
      VSAAlgebra::VecND m_a;
    };

    template< class FN, class T >
    class MemberGradFn
    {
    public:

      typedef double (FN::* FnPtr) (const T&, unsigned ip) const;
      
      MemberGradFn(const FN* fn, FnPtr pmf, unsigned ip):
	m_ip(ip), m_fn(fn->clone()),m_pmf(pmf)  
      {}
      
      virtual ~MemberGradFn() { delete m_fn; }
      
      virtual double val(const T& x) 
      {
	return (m_fn->*m_pmf)(x,m_ip);
      }
      
      // Copy-Constructor -----------------------------------------------------
      MemberGradFn(const MemberGradFn& o)
      {
	m_ip = o.m_ip;
	m_fn = o.m_fn->clone();
	m_pmf = o.m_pmf;
      }

    private:      
      unsigned      m_ip;
      FN*           m_fn;
      FnPtr         m_pmf;            
    };

    // ------------------------------------------------------------------------
    // HistFn
    // ------------------------------------------------------------------------
    template< typename HIST, typename T >
    class HistFn : public ParamFn<T>
    {
    public:
      HistFn(): ParamFn<T>(0), m_hist() { }
      HistFn(const HIST& h): ParamFn<T>(0), m_hist(h) { }

      double val(const T& x) const
      {
	return m_hist.countForVal(x[0],x[1]);
      }
      
      double val(const T& x, const VSAAlgebra::VecND& a) const
      {
	return m_hist.countForVal(x[0],x[1]);
      }

      void dyda(const T& x, VSAAlgebra::VecND& dyda) const { }
      void dyda(const T& x, const VSAAlgebra::VecND& a, 
		VSAAlgebra::VecND& dyda) const
      { }

      // Virtual Constructor --------------------------------------------------
      virtual HistFn<HIST,T>* clone() const
      {
	return new HistFn<HIST,T>(*this);
      }

    private:
      HIST   m_hist;
    };

    // ------------------------------------------------------------------------
    // Bicubic Interpolation
    //
    // Implementation taken from NR Ch. 3.6.
    //
    // ------------------------------------------------------------------------
    class BicubicInterpolation
    {
    public:
      BicubicInterpolation();
      BicubicInterpolation(double x1lo, double x1hi, unsigned nx1,
			   double x2lo, double x2hi, unsigned nx2,
			   const std::vector< std::vector<double> >& y,
			   const std::vector< std::vector<double> >& y1,
			   const std::vector< std::vector<double> >& y2,
			   const std::vector< std::vector<double> >& y12);
			   

      double val(double x1, double x2) const;

    private:

      void bcucof(unsigned ix, unsigned iy,
		  std::vector< std::vector<double> >& c) const;

      
      double                             m_x1lo;
      double                             m_x1hi;
      double                             m_nx1;
      double                             m_dx1;
      double                             m_x2lo;
      double                             m_x2hi;
      double                             m_nx2;
      double                             m_dx2;
      std::vector< std::vector<double> > m_y;
      std::vector< std::vector<double> > m_y1;
      std::vector< std::vector<double> > m_y2;
      std::vector< std::vector<double> > m_y12;
    };


    // ------------------------------------------------------------------------
    // Cubic Spline
    //
    // Function which returns the value of a piecewise cubic spline
    // interpolated from a set of user-specified data points.  Note
    // that after adding additional data points with setPoint() the
    // spline must be reevaluated with the spline() method.
    // ------------------------------------------------------------------------
    class Spline 
    {
    public:
      // ----------------------------------------------------------------------
      // Boundary Conditions:
      //
      // BC_NATURAL:   Second Derivative is zero at endpoints.
      // BC_LAGRANGE:  
      // BC_CLAMPED:   Set first derivatives at endpoints to user-defined
      //               values.
      // ----------------------------------------------------------------------
      enum BoundaryCondition
	{
	  BC_NATURAL,
	  BC_LAGRANGE,
	  BC_CLAMPED
	};

      Spline(BoundaryCondition bc = BC_NATURAL);
      Spline(const std::vector< std::pair<double,double> >& xy,
	     BoundaryCondition bc = BC_NATURAL);
      double val(const double& x) const;
      double val(const VSACoord::Coord1D& x) const;

      void setPoint(double x, double y);

      void setBoundaryCondition(BoundaryCondition bc, double yp1 = 0,
				double ypn = 0)
      {
	m_bc = bc;
	m_yp1 = yp1;
	m_ypn = ypn;
	spline();
      }

      void spline();
      void clear()
      {
	m_xy.clear();
	m_y2.clear();
      }

      const std::vector< std::pair<double,double> >& getPoints() const
      {
	return m_xy;
      }

    private:

      
      std::vector< std::pair<double,double> >  m_xy;
      std::vector<double>                      m_y2;
      double                                   m_yp1;
      double                                   m_ypn;
      BoundaryCondition                        m_bc;
    };

    // ------------------------------------------------------------------------
    // ChiSquared
    // ------------------------------------------------------------------------
    template< typename Fn, typename T >
    class ChiSquared
    {
    public:
      ChiSquared(): m_data(), m_fn()
      {

      }

      ChiSquared(const VSAMath::Data<T>& data, const Fn* fn):
	m_data(), m_fn(fn->clone())
      {
	setData(data);
      }

      ~ChiSquared()
      {
	delete m_fn;
      }

      // Accessors ------------------------------------------------------------
      unsigned nparm() const { return m_fn->nparm(); }
      VSAAlgebra::VecND param() const { return m_fn->param(); }
      double param(unsigned ip) const { return m_fn->param(ip); }
      std::vector<bool> fixed() const { return m_fn->fixed(); }
      bool fixed(unsigned ip) const { return m_fn->fixed(ip); }

      // Evaluation -----------------------------------------------------------
      double val() const
      {
	return val(m_fn->param());
      }

      double val(const VSAAlgebra::VecND& a) const
      {
	m_fn->setParam(a);
	double chi2 = 0;

	const unsigned ndata = m_data.size();
	for(unsigned idata=0;idata<ndata;idata++)
	  {
	    double ym = m_fn->val(m_data[idata].x);
	    double y = m_data[idata].y;
	    chi2 += std::pow((y-ym)/m_data[idata].sigma,2);
	  }
	return chi2;
      }

      void val(const VSAAlgebra::VecND& a, 
	       double& chi2,
	       VSAAlgebra::MatrixND& beta,
	       VSAAlgebra::MatrixND& alpha) const
      {
	m_fn->setParam(a);
	const unsigned nparm = a.ndim();
	beta.resize(nparm,1);
	alpha.resize(nparm,nparm);

	alpha.set(0.0);
	beta.set(0.0);
	chi2 = 0;

	VSAAlgebra::VecND dyda(nparm);
	const unsigned ndata = m_data.size();
	for(unsigned idata=0;idata<ndata;idata++)
	  {
	    m_fn->dyda(m_data[idata].x,dyda);
	    double ym = m_fn->val(m_data[idata].x);
	    double y = m_data[idata].y;
	    double sig2i = 1./std::pow(m_data[idata].sigma,2);
	    double dy = y-ym;
	    
	    chi2 += std::pow(dy,2)*sig2i;

	    for(unsigned iparm = 0; iparm < nparm; iparm++)
	      {
		const double wt = dyda[iparm]*sig2i;
		for(unsigned jparm = 0; jparm < nparm; jparm++)
		  alpha(iparm,jparm) += wt*dyda(jparm);

		beta(iparm,0) += dy*wt;
	      }
	  }
      }

      void dyda(const VSAAlgebra::VecND& a, VSAAlgebra::VecND& dyda) const
      {
	m_fn->setParam(a);
	const unsigned nparm = a.ndim();
	dyda.resize(nparm);
	dyda.clear();

	VSAAlgebra::VecND dfda(nparm);
	const unsigned ndata = m_data.size();
	for(unsigned idata=0;idata<ndata;idata++)
	  {
	    m_fn->dyda(m_data[idata].x,dfda);
	    double ym = m_fn->val(m_data[idata].x);
	    double y = m_data[idata].y;
	    double dy = y-ym;
	    double sig2i = 1./std::pow(m_data[idata].sigma,2);

	    dyda += dfda*sig2i*dy;
	  }
      }
		
      unsigned ndata() const { return m_data.size(); }

      // Setters --------------------------------------------------------------
      void setData(const VSAMath::Data<T>& data)
      {
	m_data = data;
      }

      void setFn(const Fn* fn) { delete m_fn; m_fn = fn->clone(); }

      void setFnParam(const VSAAlgebra::VecND& a) { m_fn->setParam(a); }

      void set(const VSAMath::Data<T>& data, const Fn& fn)
      {
	m_fn = fn;
	setData(data);
      }

      // Copy-Constructor -----------------------------------------------------
      ChiSquared(const ChiSquared& o)
      {
	m_data = o.m_data;
	m_fn = o.m_fn->clone();
      }
      
      // Assignment Operator --------------------------------------------------
      ChiSquared& operator=(const ChiSquared& o)
      {
	m_data = o.m_data;
	delete m_fn;
	m_fn = o.m_fn->clone();

	return *this;
      }

    private:
      VSAMath::Data<T>       m_data;
      Fn*                    m_fn;
    };

    // ------------------------------------------------------------------------
    //! LnPoissonLikelihood: Function class that calculates -2lnL for a given 
    //! dataset and a parameterized function specified as a template 
    //! parameter.  Assumes that each datapoint 
    //! follows a poisson distribution and thus is appropriate to use when the 
    //! data consists of a binned distribution of counts.
    // ------------------------------------------------------------------------
    template< typename Fn, typename T >
    class LnPoissonLikelihood
    {
    public:
      LnPoissonLikelihood(): m_data(), m_lny(), m_fn() { }

      LnPoissonLikelihood(const VSAMath::Data<T>& data, const Fn* fn):
	m_data(), m_lny(), m_fn(fn->clone())
      {
	setData(data);
      }

      ~LnPoissonLikelihood()
      {
	delete m_fn;
      }

      // Accessors ------------------------------------------------------------
      unsigned nparm() const { return m_fn->nparm(); }
      VSAAlgebra::VecND param() const { return m_fn->param(); }
      double param(unsigned ip) const { return m_fn->param(ip); }
      std::vector<bool> fixed() const { return m_fn->fixed(); }
      bool fixed(unsigned ip) const { return m_fn->fixed(ip); }

      // Evaluation -----------------------------------------------------------
      double val() const
      {
	return val(m_fn->param());
      }

      //! Return -2lnL for the given parameters values.
      double val(const VSAAlgebra::VecND& a) const
      {
	m_fn->setParam(a);
	double lnl = 0;

	const unsigned ndata = m_data.size();
	for(unsigned idata=0;idata<ndata;idata++)
	  {
	    double ym = m_fn->val(m_data[idata].x);
	    double y = m_data[idata].y;

	    if(ym <= 0)
	      {
		std::ostringstream os;

		os << std::endl
		   << "y = " << std::setw(15) << y
		   << " ym = " << std::setw(15) << ym
		   << " x = " << std::setw(15) << m_data[idata].x
		   << std::endl;
		
		throw 
		  std::domain_error(std::string(__PRETTY_FUNCTION__)+
				    ": Model value less than or equal to 0."+
				    os.str());
	      }

	    lnl += -2*(y*log(ym)-ym);
	  }

	lnl += m_lny;

	return lnl;
      }

      //! Return -2lnL, first derivates, and curvature matrix for 
      //! the given parameters values.
      void val(const VSAAlgebra::VecND& a, 
	       double& lnl,
	       VSAAlgebra::MatrixND& beta,
	       VSAAlgebra::MatrixND& alpha) const
      {
	m_fn->setParam(a);
	const unsigned nparm = a.ndim();
	beta.resize(nparm,1);
	alpha.resize(nparm,nparm);

	alpha.set(0.0);
	beta.set(0.0);
	lnl = 0;

	VSAAlgebra::VecND dyda(nparm);
	const unsigned ndata = m_data.size();
	for(unsigned idata=0;idata<ndata;idata++)
	  {
	    m_fn->dyda(m_data[idata].x,dyda);
	    double ym = m_fn->val(m_data[idata].x);
	    double y = m_data[idata].y;

	    if(ym <= 0)
	      {
		std::ostringstream os;

		os << std::endl
		   << "y = " << std::setw(15) << y
		   << " ym = " << std::setw(15) << ym
		   << " x = " << std::setw(15) << m_data[idata].x
		   << std::endl;
		
		throw 
		  std::domain_error(std::string(__PRETTY_FUNCTION__)+
				    ": Model value less than or equal to 0."+
				    os.str());
	      }

	    lnl += -2*(y*log(ym)-ym);

	    const double wt = y/(ym*ym);
	    for(unsigned iparm = 0; iparm < nparm; iparm++)
	      {
		for(unsigned jparm = 0; jparm < nparm; jparm++)
		  alpha(iparm,jparm) += wt*dyda(jparm)*dyda(iparm);

		beta(iparm,0) += dyda(iparm)*(y/ym-1);
	      }
	  }

	lnl += m_lny;
      }

      void dyda(const VSAAlgebra::VecND& a, VSAAlgebra::VecND& dyda) const
      {
	m_fn->setParam(a);
	const unsigned nparm = a.ndim();
	dyda.resize(nparm);
	dyda.clear();

	VSAAlgebra::VecND dfda(nparm);
	const unsigned ndata = m_data.size();
	for(unsigned idata=0;idata<ndata;idata++)
	  {
	    m_fn->dyda(m_data[idata].x,dfda);
	    double ym = m_fn->val(m_data[idata].x);
	    double y = m_data[idata].y;
	    dyda += 2*dfda*(y/ym-1);
	  }
      }
		
      //      const VSAMath::Data<T>& getData() const { return m_data; }
      unsigned ndata() const { return m_data.size(); }

      // Setters --------------------------------------------------------------
      void setData(const VSAMath::Data<T>& data)
      {
	m_data = data;
	m_lny = 0;
	const unsigned ndata = m_data.size();
	for(unsigned idata=0;idata<ndata;idata++)
	  m_lny += 2*VSAMath::lnGammaLanczos(m_data[idata].y+1);	
      }

      void setFn(const Fn* fn) { delete m_fn; m_fn = fn->clone(); }

      void setFnParam(const VSAAlgebra::VecND& a) { m_fn->setParam(a); }

      void set(const VSAMath::Data<T>& data, const Fn& fn)
      {
	m_fn = fn;
	setData(data);
      }

      // Copy-Constructor -----------------------------------------------------
      LnPoissonLikelihood(const LnPoissonLikelihood& o)
      {
	m_data = o.m_data;
	m_lny  = o.m_lny;
	m_fn   = o.m_fn->clone();
      }

      // Assignment Operator --------------------------------------------------
      LnPoissonLikelihood& operator=(const LnPoissonLikelihood& o)
      {
	m_data = o.m_data;
	m_lny = o.m_lny;
	delete m_fn;
	m_fn = o.m_fn->clone();
	
	return *this;
      }

    private:
      VSAMath::Data<T>       m_data;
      double                 m_lny;
      Fn*                    m_fn;
    };

    // ------------------------------------------------------------------------
    // Constant
    // ------------------------------------------------------------------------
    template< typename T >
    class Constant : public ParamFn<T>
    {
    public:
      Constant(double a = 1.0): ParamFn<T>(1) { ParamFn<T>::setParam(0,a); }

      // Function Evalulation -------------------------------------------------
      double val(const T& x) const { return ParamFn<T>::param(0); }
      double val(const T& x, const VSAAlgebra::VecND& a) const { return a[0]; }

      void dyda(const T& x, VSAAlgebra::VecND& dyda) const 
      {
	dyda.resize(1), dyda(0) = 1.0;
      }

      void dyda(const T& x, const VSAAlgebra::VecND& a, 
		VSAAlgebra::VecND& dyda) const 
      {
	dyda.resize(1), dyda(0) = 1.0;
      }

      // Virtual Constructor --------------------------------------------------
      Constant* clone() const 
      { return new Constant(*this); }
    };

    // ------------------------------------------------------------------------
    // Poly - 1D
    // ------------------------------------------------------------------------
    class Poly : public ParamFn<double>
    {
    public:
      Poly(unsigned n = 0): ParamFn<double>(n+1), m_n(n+1) { }

      // Accessors ------------------------------------------------------------
      void operator() (const double& x, VSAAlgebra::VecND& v) const;
      double val(const double& x) const;
      double val(const double& x, const VSAAlgebra::VecND& a) const;
      void dyda(const double& x, VSAAlgebra::VecND& dyda) const;
      void dyda(const double& x, const VSAAlgebra::VecND& a,
		VSAAlgebra::VecND& dyda) const;
      
      void dydx(const double& x, VSAAlgebra::VecND& dydx) const;

      // Setters --------------------------------------------------------------
      void set(unsigned n) { m_n = n+1; ParamFn<double>::setNParam(n+1); }

      // Virtual Constructor --------------------------------------------------
      Poly* clone() const 
      { return new Poly(*this); }

    private:
      unsigned               m_n;
    };

    // ------------------------------------------------------------------------
    // PolyR2
    // ------------------------------------------------------------------------
    template< typename T >
    class PolyR2 : public ParamFn<T>
    {
    public:
      PolyR2(unsigned n = 0): ParamFn<T>(n+1), m_n(n+1) { }

      // Evaluation -----------------------------------------------------------
      void operator() (const T& x, VSAAlgebra::VecND& v) const;
      double val(const T& x) const;
      double val(const T& x, const VSAAlgebra::VecND& a) const;
      void dyda(const T& x, VSAAlgebra::VecND& dyda) const;
      void dyda(const T& x, const VSAAlgebra::VecND& a,
		VSAAlgebra::VecND& dyda) const;

      // Setters --------------------------------------------------------------
      void set(unsigned n) { m_n = n+1; ParamFn<T>::setNParam(n+1); }

      // Virtual Constructor --------------------------------------------------
      PolyR2* clone() const 
      { return new PolyR2(*this); }

    private:
      unsigned               m_n;
    };

    template< typename T >
    void PolyR2<T>::operator() (const T& x, VSAAlgebra::VecND& v) const
    {
      v.resize(m_n);
      double s = v[0] = 1.0;
      for(unsigned i=1;i<m_n;i++) v[i] = s *= std::pow(x.r(),2);
    }

    template< typename T >
    double PolyR2<T>::val(const T& x) const
    {
      VSAAlgebra::VecND v;
      (*this)(x,v);
      return ParamFn<T>::param()*v;
    }

    template< typename T >
    double PolyR2<T>::val(const T& x, const VSAAlgebra::VecND& a) const
    {
      VSAAlgebra::VecND v;
      (*this)(x,v);
      return a * v;
    }

    template< typename T >
    void PolyR2<T>::dyda(const T& x, VSAAlgebra::VecND& dyda) const
    {
      (*this)(x,dyda);
    }

    template< typename T >
    void PolyR2<T>::dyda(const T& x, const VSAAlgebra::VecND& a,
			 VSAAlgebra::VecND& dyda) const
    {
      (*this)(x,dyda);
    }

    // ------------------------------------------------------------------------
    // Poly2D
    // ------------------------------------------------------------------------
    template< typename T >
    class Poly2D : public ParamFn<T>
    {
    public:
      Poly2D(unsigned n): 
	ParamFn<T>(0), m_n()
      { 
	initialize(n,n);
      }

      Poly2D(unsigned nx, unsigned ny): 
	ParamFn<T>(), m_n()
      { 
	initialize(nx,ny);
      }

      void initialize(unsigned nx, unsigned ny)
      {
	m_nxy.clear();
	for(unsigned ix = 0; ix <= nx; ix++)
	  for(unsigned iy = 0; iy <= ny; iy++)
	    m_nxy.push_back(std::make_pair(ix,iy));

	m_n = m_nxy.size();
	ParamFn<T>::setNParam(m_n);
      }

      // Evaluation -----------------------------------------------------------
      void operator() (const T& x, VSAAlgebra::VecND& v) const;
      double val(const T& x) const;
      double val(const T& x, const VSAAlgebra::VecND& a) const;
      void dyda(const T& x, VSAAlgebra::VecND& dyda) const;
      void dyda(const T& x, const VSAAlgebra::VecND& a, 
		VSAAlgebra::VecND& dyda) const;

      // Virtual Constructor --------------------------------------------------
      Poly2D* clone() const 
      { return new Poly2D(*this); }

    private:
      unsigned               m_n;

      std::vector< std::pair<unsigned,unsigned> > m_nxy;
    };

    template< typename T >
    void Poly2D<T>::operator() (const T& x, VSAAlgebra::VecND& v) const
    {
      v.resize(m_n);

      unsigned n = 0;
      for(std::vector< std::pair<unsigned,unsigned> >::const_iterator itr =
	    m_nxy.begin(); itr != m_nxy.end(); ++itr)
	{
	  v[n] = std::pow(x[0],(int)itr->first)*std::pow(x[1],(int)itr->second);
	  n++;
	}
    }
    
    template< typename T >
    double Poly2D<T>::val(const T& x) const
    {
      return val(x,ParamFn<T>::param());
    }

    template< typename T >
    double Poly2D<T>::val(const T& x, const VSAAlgebra::VecND& a) const
    {
      VSAAlgebra::VecND v;
      (*this)(x,v);
      return a*v;
    }

    template< typename T >
    void Poly2D<T>::dyda(const T& x, VSAAlgebra::VecND& dyda) const
    {
      (*this)(x,dyda);
    }

    template< typename T >
    void Poly2D<T>::dyda(const T& x, const VSAAlgebra::VecND& a, 
			 VSAAlgebra::VecND& dyda) const
    {
      (*this)(x,dyda);
    }

    // ------------------------------------------------------------------------
    // Gauss1D
    // ------------------------------------------------------------------------
    template< typename T = VSACoord::Coord1D >
    class Gauss1D : public ParamFn<T>
    {
    public:
      Gauss1D(): ParamFn<T>(3) { }

      // Evaluation -----------------------------------------------------------
      double val(const T& x) const;
      double val(const T& x, double n, double s, double mu) const;
      double val(const T& x, const VSAAlgebra::VecND& a) const;
      void dyda(const T& x, VSAAlgebra::VecND& dyda) const;
      void dyda(const T& x, const VSAAlgebra::VecND& a,
		VSAAlgebra::VecND& dyda) const;

      static double val(const double& x,
			double n, double s, double mu) 
      {
	double ex = std::pow(x-mu,2)/(2*std::pow(s,2));
	return n/(sqrt(2*M_PI)*fabs(s))*exp(-ex);
      }
      
      // Virtual Constructor --------------------------------------------------
      Gauss1D* clone() const 
      { return new Gauss1D(*this); }
    };

    template< typename T >
    double Gauss1D<T>::val(const T& x) const
    {
      return val(x,ParamFn<T>::param());
    }

    template< typename T >
    double Gauss1D<T>::val(const T& x, double n, double s, double mu) const
    {
      VSAAlgebra::VecND a(3);
      a(0) = n;
      a(1) = s;
      a(2) = mu;
      return val(x,a);
    }

    template< typename T >
    double Gauss1D<T>::val(const T& x, const VSAAlgebra::VecND& a) const
    {
      double ex = std::pow(x.x()-a(2),2)/(2*std::pow(a(1),2));
      return a(0)/(sqrt(2*M_PI)*fabs(a(1)))*exp(-ex);
    }

    template< typename T >
    void Gauss1D<T>::dyda(const T& x, VSAAlgebra::VecND& dyda) const
    {
      Gauss1D<T>::dyda(x,ParamFn<T>::param(),dyda);
    }      

    template< typename T >
    void Gauss1D<T>::dyda(const T& x, const VSAAlgebra::VecND& a,
			  VSAAlgebra::VecND& dyda) const
    {
      dyda.resize(3);
      double dx = x.x()-a(2);
      double ex = dx*dx/(2*std::pow(a(1),2));
      double c = 1./sqrt(2*M_PI);

      dyda(0) = exp(-ex)*c/fabs(a(1));
      dyda(1) = -a(0)*c*exp(-ex)*(a(1)*a(1)-dx*dx)*std::pow(a(1),-4);
      dyda(2) = a(0)*c*exp(-ex)*std::pow(fabs(a(1)),-3)*dx;
    }      

    // ------------------------------------------------------------------------
    // Gauss2D
    // ------------------------------------------------------------------------
    template< typename T = VSACoord::Coord2D >
    class Gauss2D : public ParamFn<T>
    {
    public:
      Gauss2D(): ParamFn<T>(4) { }

      // Evaluation -----------------------------------------------------------
      double val(const T& x) const;
      double val(const T& x, const VSAAlgebra::VecND& a) const;
      void dyda(const T& x, VSAAlgebra::VecND& dyda) const;
      void dyda(const T& x, const VSAAlgebra::VecND& a,
		VSAAlgebra::VecND& dyda) const;

      // Virtual Constructor --------------------------------------------------
      Gauss2D* clone() const 
      { return new Gauss2D(*this); }
    };

    template< typename T >
    double Gauss2D<T>::val(const T& x) const
    {
      return val(x,ParamFn<T>::param());
    }

    template< typename T >
    double Gauss2D<T>::val(const T& x, const VSAAlgebra::VecND& a) const
    {
      double s2i = std::pow(a(1),-2);
      double ex = (std::pow(x.x()-a(2),2) + std::pow(x.y()-a(3),2))*s2i/2.;
      return a(0)*s2i/(2*M_PI)*exp(-ex);
    }

    template< typename T >
    void Gauss2D<T>::dyda(const T& x, VSAAlgebra::VecND& dyda) const
    {
      Gauss2D<T>::dyda(x,ParamFn<T>::param(),dyda);
    }      

    template< typename T >
    void Gauss2D<T>::dyda(const T& x, const VSAAlgebra::VecND& a,
			  VSAAlgebra::VecND& dyda) const
    {
      dyda.resize(4);
      dyda.clear();
      double s2i = std::pow(a(1),-2);
      double ex = (std::pow(x.x()-a(2),2) + std::pow(x.y()-a(3),2))*s2i/2.;

      if(!ParamFn<T>::fixed(0)) dyda(0) = exp(-ex)*s2i/(2*M_PI);
      if(!ParamFn<T>::fixed(1)) 
	dyda(1) = 
	  exp(-ex)*(a(0)/(2*M_PI*std::pow(a(1),5))*
		    (std::pow(x.x()-a(2),2) + std::pow(x.y()-a(3),2)) -
		    a(0)/(M_PI*std::pow(a(1),3)));
      if(!ParamFn<T>::fixed(2)) 
	dyda(2) = a(0)/(2*M_PI)*s2i*s2i*(x.x()-a(2))*exp(-ex);
      if(!ParamFn<T>::fixed(3)) 
	dyda(3) = a(0)/(2*M_PI)*s2i*s2i*(x.y()-a(3))*exp(-ex);
    }      

    // ------------------------------------------------------------------------
    // FourierBesselSeries
    // ------------------------------------------------------------------------
    class FourierBesselSeries : public ParamFn<VSACoord::Coord2D>
    {
    public:
      FourierBesselSeries(unsigned m, unsigned n, double R);
      FourierBesselSeries(const std::vector<std::pair<int,int> >& mn,
			  double R);

      // Accessors ------------------------------------------------------------
      const VSAAlgebra::VecND& norm() { return m_norm; };

      // Function Evalulation -------------------------------------------------
      void operator() (const VSACoord::Coord2D& x, VSAAlgebra::VecND& v) const;
      double val(const VSACoord::Coord2D& x, const VSAAlgebra::VecND& a) const;
      double val(const VSACoord::Coord2D& x) const;
      void dyda(const VSACoord::Coord2D& x, VSAAlgebra::VecND& dyda) const;
      void dyda(const VSACoord::Coord2D& x, const VSAAlgebra::VecND& a, 
		VSAAlgebra::VecND& dyda) const; 

      // Setters --------------------------------------------------------------
      void set(const std::vector<std::pair<int,int> >& mn, double R);

      // Virtual Constructor --------------------------------------------------
      FourierBesselSeries* clone() const 
      { return new FourierBesselSeries(*this); }

    private:
      std::vector<std::pair<int,int> > m_mn;
      double                           m_R;
      VSAAlgebra::VecND                m_norm;
    };

    // ------------------------------------------------------------------------
    // FourierBesselSeries2
    // ------------------------------------------------------------------------
    class FourierBesselSeries2 : public ParamFn<VSACoord::Coord2D>
    {
    public:

      struct Param
      {
	Param(int ip, int m, int n): ip(ip), m(m), n(n) { }

	int ip;
	int m;
	int n;
      };

      typedef std::vector<Param>::const_iterator const_iterator;

      FourierBesselSeries2(double R);
      FourierBesselSeries2(unsigned m, unsigned n, double R);
      FourierBesselSeries2(const std::vector<std::pair<unsigned,
			   unsigned> >& mn, double R);

      // Accessors ------------------------------------------------------------
      const_iterator begin() const
      {
	return m_mn.begin();
      }

      const_iterator end() const
      {
	return m_mn.end();
      }

      double getC() const { return m_c; }

      // Function Evalulation -------------------------------------------------
      double val(const VSACoord::Coord2D& x, const VSAAlgebra::VecND& a) const;
      double val(const VSACoord::Coord2D& x) const;
      
      void dyda(const VSACoord::Coord2D& x, VSAAlgebra::VecND& dyda) const; 
      void dyda(const VSACoord::Coord2D& x, const VSAAlgebra::VecND& a, 
		VSAAlgebra::VecND& dyda) const; 

      // Setters --------------------------------------------------------------
      virtual void setNParam(unsigned n);

      virtual void setParam(const VSAAlgebra::VecND& a) 
      { 
	ParamFn<VSACoord::Coord2D>::setParam(a);
	calc(param(),m_c,m_dc);
      }

      virtual void setParam(const VSAAlgebra::VecND& a, 
			    const std::vector< bool >& fixed)
      {
	ParamFn<VSACoord::Coord2D>::setParam(a,fixed);
	calc(param(),m_c,m_dc);
      }

      virtual void setParam(unsigned ip, double a)
      {
	ParamFn<VSACoord::Coord2D>::setParam(ip,a);
	calc(param(),m_c,m_dc);
      } 

      void set(unsigned m, unsigned n);
      void set(const std::vector<std::pair<unsigned,unsigned> >& mn);

      void set(double c) { m_c = c; }

      // Virtual Constructor --------------------------------------------------
      FourierBesselSeries2* clone() const 
      { return new FourierBesselSeries2(*this); }

    private:

      void calc(const VSAAlgebra::VecND& a,double& c,
		VSAAlgebra::VecND& dc) const;

      void getJ(const VSACoord::Coord2D& x, std::vector<double>& J) const;

      void computeSpline();

      std::vector<Param>                         m_mn;
      std::vector<Spline>                        m_bessel_spline;
      std::vector<VSNSpace>                      m_bessel_hist;
      double                                     m_R;
      double                                     m_norm;
      double                                     m_c;
      VSAAlgebra::VecND                          m_dc;
    };
  } 
}

#endif // VSAFUNCTION_HPP
