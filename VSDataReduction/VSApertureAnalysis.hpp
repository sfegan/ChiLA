//-*-mode:c++; mode:font-lock;-*-

/*! \file VSApertureAnalysis.hpp

  Base call for aperture-based integral analysis calcultors (ring
  background and reflected region).

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.7 $
  \date       07/29/2007

  $Id: VSApertureAnalysis.hpp,v 3.7 2010/10/20 17:51:57 matthew Exp $

*/

#ifndef VSAPERTUREANALYSIS_HPP
#define VSAPERTUREANALYSIS_HPP

#include <VSIntegralAnalysis.hpp>

namespace VERITAS
{
  //! Base class for integral analysis calculators using a circular
  //! signal aperture (reflected region and ring background).
  class VSApertureAnalysis : public VSIntegralAnalysis
  {
  public:
    VSApertureAnalysis(const std::string bkgnd_model,
		       double bin_size_deg,
		       double theta_cut,
		       const std::pair< double, double >& ring_cut,
		       double max_offset_cut,
		       unsigned max_nregion,
		       const std::string& spmodel);  

    virtual ~VSApertureAnalysis();

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

  protected:    
    
    std::vector< VSSimple2DHist<double,double> > m_sky_on_counts_hist;
    std::vector< VSSimple2DHist<double,double> > m_sky_off_counts_hist;
    std::vector< VSSimple2DHist<double,double> > m_sky_domega_hist;
    std::vector< VSSimple2DHist<double,double> > m_sky_livetime_hist;
    std::vector< VSSimple2DHist<double,double> > m_sky_exposure_hist;
    std::vector< VSSimple2DHist<double,double> > m_sky_alpha_hist;
    std::vector< VSSimple2DHist<double,double> > m_sky_effarea_hist;
    std::vector< double >                        m_livetime;
    std::vector<VSAAlgebra::Vec2D>               m_ptg_xy;
  };
};


#endif // VSAPERTUREANALYSIS_HPP
