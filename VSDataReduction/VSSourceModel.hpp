//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSourceModel.hpp

  Class for fitting a source model to data.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.13 $
  \date       07/29/2007

  $Id: VSSourceModel.hpp,v 3.13 2010/10/25 00:42:03 matthew Exp $

*/

#ifndef VSSOURCEMODEL_HPP
#define VSSOURCEMODEL_HPP

#include <vector>
#include <VSAFunction.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSOctaveH5Writer.hpp>
#include <VSModelCoord.hpp>

namespace VERITAS
{

  //! Base class for parametric source models used by the MLM
  //! analysis.
  class VSSourceModel : public VSAFunction::ParamFn<VSModelCoord>
  {
  public:
    VSSourceModel(double bin_size_deg);
    virtual ~VSSourceModel() { }

    // Evaluation -------------------------------------------------------------
    double val(const VSAAlgebra::Vec2D& xy) const;

    //! Evaluate the amplitude of the source model in counts at the
    //! given coordinate.
    virtual double val(const VSModelCoord& x) const = 0;

    //! Evaluate the amplitude of the source model in counts at the
    //! given coordinate.
    virtual double val(const VSModelCoord& x, 
		       const VSAAlgebra::VecND& a) const = 0;

    //! Evaluate the derivative of the source model amplitude in
    //! counts with respect to the model parameters at the given
    //! coordinate.
    //! @param dyda Vector of derivatives of model amplitude with
    //! respect to model parameters.
    virtual void dyda(const VSModelCoord& x, 
		      VSAAlgebra::VecND& dyda) const = 0;

    //! Evaluate the derivative of the source model amplitude in
    //! counts with respect to the model parameters at the given
    //! coordinate.
    //! @param dyda Vector of derivatives of model amplitude with
    //! respect to model parameters.
    virtual void dyda(const VSModelCoord& x, const VSAAlgebra::VecND& a,
		      VSAAlgebra::VecND& dyda) const = 0;
    
    virtual double dyda(const VSModelCoord& x, unsigned ip) const = 0;


    virtual VSAAlgebra::Vec2D getXY() const = 0;

    // Setters ----------------------------------------------------------------
    virtual void setSourceParam(const std::vector<std::string>& param) = 0;
    virtual void setXY(const VSAAlgebra::Vec2D& xy) = 0;
    virtual void setAmplitude(double x) = 0;
    
    virtual void setParam(const VSAAlgebra::VecND& a) 
    {
      VSAFunction::ParamFn<VSModelCoord>::setParam(a);
    }

    virtual void setParam(const VSAAlgebra::VecND& a,
			  const std::vector< bool >& fixed) 
    { 
      VSAFunction::ParamFn<VSModelCoord>::setParam(a,fixed);
    }

    virtual void setParam(unsigned ip, double a) 
    { 
      VSAFunction::ParamFn<VSModelCoord>::setParam(ip,a);
    }

    //! Integrate total counts in source model.
    //! @param a        Model parameters.
    //! @param cov      Model parameters covariance matrix.
    //! @param integral Value of integral.
    //! @param err      Statistical error on integral value.
    virtual void integrate(const VSAAlgebra::VecND& a,
			   const VSAAlgebra::MatrixND& cov,
			   double& integral, double& err) const;

    //! Integrate total counts in source model inside the annular
    //! region defined by xy, R1, R2.
    //! @param xy      Coordinate of annulus center.
    //! @param R1      Inner radius of integration region.
    //! @param R2      Outer radius of integration region.
    virtual double integrate(const VSAAlgebra::Vec2D& xy, double R1,
			     double R2) const;

    //! Integrate total counts in source model inside the annular
    //! region defined by R1 and R2.  The center of the annulus is
    //! assumed to be centered on the source position.  This method
    //! explicitly assumes that the source function is azimuthally
    //! symmetric.  
    //! @param R1 Inner radius of integration region.
    //! @param R2 Outer radius of integration region.
    virtual double integrate(double R1, double R2) const;

    virtual void clear();

    //! Define an observation coordinate.
    void setObs(const VSAAlgebra::Vec2D& xy, 
		const VSNSpace& effarea,
		const VSNSpace& psf);

    virtual void initialize(const VSNSpace& psf) { }


    virtual VSSourceModel* clone() const = 0;

  protected:
    double                                       m_bin_size;
    double                                       m_domega;
    double                                       m_rmax;
    std::vector< VSAAlgebra::Vec2D >             m_obs_xy;
    std::vector< VSNSpace >                      m_effarea;
    std::vector< VSNSpace >                      m_psf;
  };

  //! Source model using a two-dimensional symmetric gaussian.
  class VSSourceModelGauss : public VSSourceModel
  {
  public:
    VSSourceModelGauss(double bin_size_deg);
    virtual  ~VSSourceModelGauss() { }

    // Evaluation -------------------------------------------------------------
    double val(const VSModelCoord& x) const;
    double val(const VSModelCoord& x, const VSAAlgebra::VecND& a) const;
    void dyda(const VSModelCoord& x, VSAAlgebra::VecND& dyda) const;
    void dyda(const VSModelCoord& x, const VSAAlgebra::VecND& a,
	      VSAAlgebra::VecND& dyda) const;
    
    double dyda(const VSModelCoord& x, unsigned ip) const;
    
    double getAcceptance(const VSModelCoord& x, 
			 const VSAAlgebra::VecND& a) const;

    VSAAlgebra::Vec2D getXY() const
    {
      VSAAlgebra::Vec2D xy(param(2),param(3));
      return xy;
    }

    // Setters ----------------------------------------------------------------
    void setParam(const VSAAlgebra::VecND& a) 
    {
      VSAFunction::ParamFn<VSModelCoord>::setParam(a);
      m_fn.setParam(a); 
    }

    void setParam(const VSAAlgebra::VecND& a,
		  const std::vector< bool >& fixed) 
    { 
      VSAFunction::ParamFn<VSModelCoord>::setParam(a,fixed);
      m_fn.setParam(a,fixed); 
    }

    void setParam(unsigned ip, double a) 
    { 
      VSAFunction::ParamFn<VSModelCoord>::setParam(ip,a);
      m_fn.setParam(ip,a); 
    }

    void setSourceParam(const std::vector<std::string>& param);

    void setXY(const VSAAlgebra::Vec2D& xy) 
    {
      setParam(2,xy.x());
      setParam(3,xy.y());
    }

    void setAmplitude(double x) { setParam(0,x); }

    // Virtual Constructors ---------------------------------------------------
    virtual VSSourceModelGauss* clone() const 
    { return new VSSourceModelGauss(*this); }

  private:
    VSAFunction::Gauss2D<VSModelCoord>       m_fn;
  };

  //! Point source model which uses the simulated PSF to generate a
  //! model for the distribution of signal counts.
  class VSSourceModelPointSource : public VSSourceModel
  {
  public:
    VSSourceModelPointSource(double bin_size_deg);
    virtual  ~VSSourceModelPointSource() { }

    // Evaluation -------------------------------------------------------------
    double val(const VSModelCoord& x) const;
    double val(const VSModelCoord& x, const VSAAlgebra::VecND& a) const;
    void dyda(const VSModelCoord& x, VSAAlgebra::VecND& dyda) const;
    void dyda(const VSModelCoord& x, const VSAAlgebra::VecND& a,
	      VSAAlgebra::VecND& dyda) const;
    
    double dyda(const VSModelCoord& x, unsigned ip) const;

    double getAcceptance(const VSModelCoord& x, 
			 const VSAAlgebra::VecND& a) const;

    VSAAlgebra::Vec2D getXY() const
    {
      VSAAlgebra::Vec2D xy(param(1),param(2));
      return xy;
    }

    // Setters ----------------------------------------------------------------
    void setSourceParam(const std::vector<std::string>& param);

    void setXY(const VSAAlgebra::Vec2D& xy) 
    {
      setParam(1,xy.x());
      setParam(2,xy.y());
    }

    void setAmplitude(double x) { setParam(0,x); }

    // Virtual Constructors ---------------------------------------------------
    virtual VSSourceModelPointSource* clone() const 
    { return new VSSourceModelPointSource(*this); }

  };


  //! Azimuthally symmetric extended source model.
  class VSSourceModelExtSource : public VSSourceModel
  {
  public:
    struct ConvolveFn
    {
      ConvolveFn(const VSNSpace& src_model, 
		 const VSNSpace& effarea_domega, 
		 double theta): 
	m_src_model(src_model), m_effarea_domega(effarea_domega),
	m_theta(theta) 
      { }

      double val(const VSACoord::Coord2D &a);

      const VSNSpace& m_src_model;
      const VSNSpace& m_effarea_domega;
      double m_theta;
    };


    VSSourceModelExtSource(double bin_size_deg, const std::string& ext_model);
    virtual  ~VSSourceModelExtSource() { }

    // Evaluation -------------------------------------------------------------
    double val(const VSModelCoord& x) const;
    double val(const VSModelCoord& x, const VSAAlgebra::VecND& a) const;
    void dyda(const VSModelCoord& x, VSAAlgebra::VecND& dyda) const;
    void dyda(const VSModelCoord& x, const VSAAlgebra::VecND& a,
	      VSAAlgebra::VecND& dyda) const;
    
    double dyda(const VSModelCoord& x, unsigned ip) const;

    double getAcceptance(const VSModelCoord& x, 
			 const VSAAlgebra::VecND& a) const;

    VSAAlgebra::Vec2D getXY() const
    {
      VSAAlgebra::Vec2D xy(param(1),param(2));
      return xy;
    }

    // Setters ----------------------------------------------------------------
    void setSourceParam(const std::vector<std::string>& param);

    void setXY(const VSAAlgebra::Vec2D& xy) 
    {
      setParam(1,xy.x());
      setParam(2,xy.y());
    }

    void setAmplitude(double x) { setParam(0,x); }

    void initialize(const VSNSpace& psf);

    virtual void clear();

    // Virtual Constructors ---------------------------------------------------
    virtual VSSourceModelExtSource* clone() const 
    { return new VSSourceModelExtSource(*this); }

  private:
    VSNSpace                                 m_model;
    std::vector< VSNSpace >                  m_effarea_domega;
  };
}

#endif // VSSOURCEMODEL_HPP
