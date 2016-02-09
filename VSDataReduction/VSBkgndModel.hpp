//-*-mode:c++; mode:font-lock;-*-

/*! \file VSBkgndModel.hpp

  Class for generating a CR acceptance model.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.15 $
  \date       07/29/2007

  $Id: VSBkgndModel.hpp,v 3.15 2010/06/18 18:10:42 matthew Exp $

*/

#ifndef VSBKGNDMODEL_HPP
#define VSBKGNDMODEL_HPP

#include <vector>
#include <VSAFunction.hpp>
#include <VSModelCoord.hpp>
#include <VSAnalysisStage3Data.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSAcceptanceModel
  // ==========================================================================
  class VSAcceptanceModel : public VSAFunction::ParamFn<VSModelCoord>
  {
  public:
    VSAcceptanceModel(double bin_size_deg, double offset_max);
    virtual ~VSAcceptanceModel() { }

    // Accessors --------------------------------------------------------------

    //! Return the number of sky pointings in the acceptance model
    unsigned nobs() const { return m_obs_xy.size(); }
    const VSAAlgebra::Vec2D& getObsXY(unsigned iobs) { return m_obs_xy[iobs]; }

    // Evaluation -------------------------------------------------------------

    //! Return cumulative occupation for all bins at the given sky
    //! coordinate
    virtual double val(const VSAAlgebra::Vec2D& xy) const = 0;
    virtual void val(const VSAAlgebra::Vec2D& xy, const VSAAlgebra::VecND& a,
		     const VSAAlgebra::MatrixND& cov,
		     double& val, double& err) const = 0;

    //! Return mean occuption for sky bin with the given model coordinate
    virtual double val(const VSModelCoord& x) const = 0;
    virtual double val(const VSModelCoord& x, const VSAAlgebra::VecND& a) 
      const = 0;

    //! Return the value of the counts density function (dN/dOmega) at the
    //! given 2D coordinate.
    virtual double getFoVCounts(const VSAAlgebra::Vec2D& xy) const = 0;

    //! Return the value of the acceptance function (dP/dOmega) at the
    //! given 2D coordinate.
    virtual double getFoVAcceptance(const VSAAlgebra::Vec2D& xy) const = 0;

    //! Return the value of the acceptance function (dP/dR^2) averaged
    //! over azimuthal angle.
    virtual double getFoVAcceptance(double r2) const = 0;

    
    virtual void integrate(const VSAAlgebra::VecND& a,
			   const VSAAlgebra::MatrixND& cov,
			   double& integral, double& integral_err) const = 0;

    //! Integrate background model in sky coordinate system.
    virtual double integrate(const VSAAlgebra::Vec2D& xy, double R1,
			     double R2) const = 0;
    
    //! Integrate background model in FoV coordinate system.
    virtual double integrateFoV(double R1, double R2) const = 0;

    // Setters ----------------------------------------------------------------
    virtual void setAcceptanceParam(const VSAAlgebra::VecND& a) = 0;

    // Utility Functions ------------------------------------------------------
    virtual VSAcceptanceData* fit(const VSAMath::Data<VSModelCoord>& data) = 0;
    virtual void clear() 
    {
      m_obs_xy.clear();
    }

    void setObs(const VSAAlgebra::Vec2D& xy);

    // Virtual Constructor ----------------------------------------------------
    virtual VSAcceptanceModel* clone() const = 0;

  protected:

    double                                       m_bin_size;
    double                                       m_domega;
    double                                       m_offset_max;
    std::vector< VSAAlgebra::Vec2D >             m_obs_xy;
  };

  // ==========================================================================
  // VSAcceptanceModelFixed
  // ==========================================================================
  class VSAcceptanceModelFixed : public VSAcceptanceModel
  {
  public:
    VSAcceptanceModelFixed(VSAcceptanceModel* model,
			   const VSSimple2DHist<double,double>& counts_hist);
    virtual ~VSAcceptanceModelFixed() { }

    // Evaluation -------------------------------------------------------------
    double val(const VSAAlgebra::Vec2D& xy) const;
    void val(const VSAAlgebra::Vec2D& xy, const VSAAlgebra::VecND& a,
	     const VSAAlgebra::MatrixND& cov, double& val, double& err) const
    { }

    double val(const VSModelCoord& x) const;
    double val(const VSModelCoord& x, const VSAAlgebra::VecND& a) const;

    void dyda(const VSModelCoord& x, VSAAlgebra::VecND& dyda) const;
    void dyda(const VSModelCoord& x, const VSAAlgebra::VecND& a,
	      VSAAlgebra::VecND& dyda) const;

    // Setters ----------------------------------------------------------------
    double getFoVCounts(const VSAAlgebra::Vec2D& xy) const;
    double getFoVAcceptance(const VSAAlgebra::Vec2D& xy) const;
    double getFoVAcceptance(double r2) const { return 0; }

    void integrate(const VSAAlgebra::VecND& a,
		   const VSAAlgebra::MatrixND& cov,
		   double& integral, double& integral_err) const;

    double integrate(const VSAAlgebra::Vec2D& xy, double R1, double R2) const;

    double integrateFoV(double R1, double R2) const { return 0; }

    virtual void setAcceptanceParam(const VSAAlgebra::VecND& a)
    {

    }

    // Utility Functions ------------------------------------------------------
    VSAcceptanceData* fit(const VSAMath::Data<VSModelCoord>& data)
    {
      return NULL;
    }

    // Virtual Constructor ----------------------------------------------------
    virtual VSAcceptanceModelFixed* clone() const 
    { return new VSAcceptanceModelFixed(*this); }

  protected:

    std::vector< VSSimple2DHist<double,double> > m_fn_val;
  };

  // ==========================================================================
  // VSAcceptanceModelPoly
  // ==========================================================================
  class VSAcceptanceModelPoly : public VSAcceptanceModel
  {
  public:
    
    VSAcceptanceModelPoly(const std::string& param,
			  double bin_size_deg, double offset_max);    
    ~VSAcceptanceModelPoly() { }

    // Evaluation -------------------------------------------------------------
    double val(const VSAAlgebra::Vec2D& xy) const { return 0; }
    void val(const VSAAlgebra::Vec2D& xy, const VSAAlgebra::VecND& a,
	     const VSAAlgebra::MatrixND& cov, double& val, double& err) const
    { }

    double val(const VSModelCoord& x) const;
    double val(const VSModelCoord& x, const VSAAlgebra::VecND& a) const;
    
    void dyda(const VSModelCoord& x, VSAAlgebra::VecND& dyda) const;
    void dyda(const VSModelCoord& x, const VSAAlgebra::VecND& a,
	      VSAAlgebra::VecND& dyda) const;
    
    // Setters ----------------------------------------------------------------
    double getFoVCounts(const VSAAlgebra::Vec2D& xy) const { return 0; }
    double getFoVAcceptance(const VSAAlgebra::Vec2D& xy) const;
    double getFoVAcceptance(double r2) const { return 0; }

    void integrate(const VSAAlgebra::VecND& a,
		   const VSAAlgebra::MatrixND& cov,
		   double& integral, double& integral_err) const;
    double integrate(const VSAAlgebra::Vec2D& xy, double R1, double R2) const;
    double integrate(const VSAAlgebra::Vec2D& xy, double R) const;

    double integrateFoV(double R1, double R2) const { return 0; }

    virtual void setAcceptanceParam(const VSAAlgebra::VecND& a)
    {
      vsassert(a.ndim() == m_fn.nparm());
      for(unsigned ip = 0; ip < a.ndim(); ip++) setParam(ip,a(ip));
    }

    // Utility Functions -----------------------------------------------------
    virtual VSAcceptanceData* fit(const VSAMath::Data<VSModelCoord>& data)
    {
      return NULL;
    }

    virtual void clear() 
    {
      m_obs_xy.clear();
      setNParam(m_fn.nparm());
    }

    // Virtual Constructors ---------------------------------------------------
    virtual VSAcceptanceModelPoly* clone() const 
    { return new VSAcceptanceModelPoly(*this); }

  private:
    VSAFunction::PolyR2<VSACoord::Coord2D> m_fn;
  };

  // ==========================================================================
  // VSAcceptanceModelBessel2
  // ==========================================================================
  class VSAcceptanceDataBessel2 : public VSAcceptanceData
  {
  public:    

    struct Param
    {
      Param(): n(), c(), c_err() { }

      unsigned n;
      double   c;
      double   c_err;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,Param,n);
	H5_ADDMEMBER(c,Param,c);
	H5_ADDMEMBER(c,Param,c_err);
      }
    };

    struct AzimuthalParam
    {
      AzimuthalParam(): 
	n(), 
	cx(), cx_err(), cy(), cy_err(),
	cr(), cr_err(), cphi(), cphi_err()
      { }

      unsigned n;
      double   cx;
      double   cx_err;
      double   cy;
      double   cy_err;
      double   cr;
      double   cr_err;
      double   cphi;
      double   cphi_err;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,AzimuthalParam,n);
	H5_ADDMEMBER(c,AzimuthalParam,cx);
	H5_ADDMEMBER(c,AzimuthalParam,cx_err);
	H5_ADDMEMBER(c,AzimuthalParam,cy);
	H5_ADDMEMBER(c,AzimuthalParam,cy_err);
	H5_ADDMEMBER(c,AzimuthalParam,cr);
	H5_ADDMEMBER(c,AzimuthalParam,cr_err);
	H5_ADDMEMBER(c,AzimuthalParam,cphi);
	H5_ADDMEMBER(c,AzimuthalParam,cphi_err);
      }

    };

    VSAcceptanceDataBessel2(double offset_max,
			    const VSAAlgebra::VecND& p = VSAAlgebra::VecND(),
			    const VSAAlgebra::MatrixND& cov =
			    VSAAlgebra::MatrixND());

    virtual VSAcceptanceDataBessel2* clone()
    {
      return new VSAcceptanceDataBessel2(*this);
    }

    // Load/Save --------------------------------------------------------------
    virtual bool load(VSOctaveH5ReaderStruct* reader);
    virtual void save(VSOctaveH5WriterStruct* writer) const;

    double                              c_param;
    unsigned                            nparam_shape;

    std::vector< Param >                         m0_param;
    std::vector< std::vector< AzimuthalParam > > mn_param;

    VSLimitedErrorsHist<double, double>               m0_param_hist;
    std::vector<VSLimitedErrorsHist<double, double> > mx_param_hist;
    std::vector<VSLimitedErrorsHist<double, double> > my_param_hist;
    std::vector<VSLimitedErrorsHist<double, double> > mr_param_hist;
    std::vector<VSLimitedErrorsHist<double, double> > mphi_param_hist;
  };

  class VSAcceptanceModelBessel2 : public VSAcceptanceModel
  {
  public:

    typedef VSAFunction::LnPoissonLikelihood<VSAcceptanceModelBessel2,
					     VSModelCoord> LnLFn;

    VSAcceptanceModelBessel2(const std::string& param,
			     double bin_size_deg, double offset_max);

    // Evaluation -------------------------------------------------------------
    double val(const VSAAlgebra::Vec2D& xy) const;
    void val(const VSAAlgebra::Vec2D& xy, const VSAAlgebra::VecND& a,
	     const VSAAlgebra::MatrixND& cov, double& val, double& err) const;

    double val(const VSModelCoord& x) const;
    double val(const VSModelCoord& x, const VSAAlgebra::VecND& a) const;
    
    void dyda(const VSModelCoord& x, VSAAlgebra::VecND& dyda) const;
    void dyda(const VSModelCoord& x, const VSAAlgebra::VecND& a,
	      VSAAlgebra::VecND& dyda) const;

    double getFoVCounts(const VSAAlgebra::Vec2D& xy) const; 
    double getFoVAcceptance(const VSAAlgebra::Vec2D& xy) const; 
    double getFoVAcceptance(double r2) const; 

    void integrate(const VSAAlgebra::VecND& a,
		   const VSAAlgebra::MatrixND& cov,
		   double& integral, double& integral_err) const;
    double integrate(const VSAAlgebra::Vec2D& xy, double R1, double R2) const;
    
    double integrateFoV(double R1, double R2) const;

    // Setters ----------------------------------------------------------------
    virtual void setParam(const VSAAlgebra::VecND& a) 
    { 
      VSAFunction::ParamFn<VSModelCoord>::setParam(a);
      m_fn.setParam(param().subVector(0,m_fn.nparm()));
    }

    virtual void setParam(const VSAAlgebra::VecND& a, 
			  const std::vector< bool >& fixed)
    {
      VSAFunction::ParamFn<VSModelCoord>::setParam(a,fixed);
      m_fn.setParam(param().subVector(0,m_fn.nparm()));
    }

    virtual void setParam(unsigned ip, double a)
    {
      VSAFunction::ParamFn<VSModelCoord>::setParam(ip,a);
      if(ip < m_fn.nparm()) m_fn.setParam(ip,a);
    } 

    void set(unsigned m, unsigned n);

    virtual void setAcceptanceParam(const VSAAlgebra::VecND& a)
    {
      vsassert(a.ndim() == m_fn.nparm());
      for(unsigned ip = 0; ip < a.ndim(); ip++) setParam(ip,a(ip));
    }
    
    // Utility Functions -----------------------------------------------------
    virtual VSAcceptanceData* fit(const VSAMath::Data<VSModelCoord>& data);

    virtual void clear() 
    {
      m_obs_xy.clear();
      setNParam(m_fn.nparm());
    }

    // Virtual Constructors ---------------------------------------------------
    virtual VSAcceptanceModelBessel2* clone() const 
    { return new VSAcceptanceModelBessel2(*this); }

  private:
    void fillData(VSAcceptanceDataBessel2* data);
    void calc(const VSAMath::Data<VSModelCoord>& data,
	      VSAAlgebra::VecND& p);
    VSAFunction::FourierBesselSeries2 m_fn;    
  };

}

#endif // VSBKGNDMODEL_HPP
