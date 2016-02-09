//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMLMAnalysis.hpp

  Maximum likelihood method.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.6 $
  \date       07/29/2007

  $Id: VSMLMAnalysis.hpp,v 3.6 2010/10/20 03:33:16 matthew Exp $

*/

#ifndef VSMLMANALYSIS_HPP
#define VSMLMANALYSIS_HPP

#include <vector>

#include <VSIntegralAnalysis.hpp>
#include <VSPointing.hpp>
#include <VSSimple2DHist.hpp>
#include <VSAAlgebra.hpp>
#include <VSSourceModel.hpp>
#include <VSANonlinearFitting.hpp>

namespace VERITAS
{
  //! Integral stage3 analysis class responsible for maximum
  //! likelihood method (MLM).
  class VSMLMAnalysis : public VSIntegralAnalysis
  {
  public:
    struct FitResults
    {
      FitResults(): chi2(), param(), cov(), 
		    bkgnd_param(), bkgnd_cov(),
		    source_param(), source_cov()
      { }

      double                chi2;
      VSAAlgebra::VecND     param;
      VSAAlgebra::MatrixND  cov;

      VSAAlgebra::VecND     bkgnd_param;
      VSAAlgebra::MatrixND  bkgnd_cov;

      VSAAlgebra::VecND     source_param;
      VSAAlgebra::MatrixND  source_cov;
    };

    typedef VSAFunction::LnPoissonLikelihood<VSDataModel,VSModelCoord> LnLFn;

    VSMLMAnalysis(const std::string& bkgnd_model,	
		  double bin_size_deg,
		  double theta_cut, 
		  const std::pair< double, double >& ring_cut,
		  double theta_max,
		  const std::string& source_spectrum,
		  bool fit_skymap,
		  const std::string& source_model,
		  const std::string& ext_model,
		  double skymap_fit_radius);
    virtual ~VSMLMAnalysis();

    //! Generate 2D sky maps.
    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisData& o);

    //! Generate results for a single source position.
    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisDatum& o);
    
  private:
    
    bool loadFitData(const VSAnalysisStage3Data& data,
		     const VSAcceptanceData& acceptance);

    void fitSource(const VSAnalysisStage3Data& data,
		   VSIntegralAnalysisDatum& o);

    void fit(const VSAMath::Data<VSModelCoord>& data,
	     const VSAAlgebra::Vec2D& xy, double R, bool fit_bkgnd,
	     FitResults& bkgnd_fit_data, FitResults& source_fit_data);

    void excludeData(const VSAMath::Data<VSModelCoord>& data,
		     VSAMath::Data<VSModelCoord>& d);
    void excludeData(const VSAMath::Data<VSModelCoord>& data,
		     const VSAAlgebra::Vec2D& xy, double R,
		     VSAMath::Data<VSModelCoord>& d);

    double                                           m_domega;
    bool                                             m_fit_skymap;
    double                                           m_fit_radius;

    std::vector<double>                              m_livetime;
    VSAMath::Data<VSModelCoord>                      m_fit_data;
    FitResults                                       m_bkgnd_fit_data;
    FitResults                                       m_source_fit_data;
    std::vector< std::pair<double,bool> >            m_source_model_param;

    VSSourceModel*                                   m_source_model;
    VSDataModel*                                     m_data_model;
    VSAMath::NLFitter< LnLFn >*                      m_fitter; 
  };
}

#endif // VSMLMANALYSIS_HPP
