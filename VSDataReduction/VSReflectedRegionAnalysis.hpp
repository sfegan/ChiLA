//-*-mode:c++; mode:font-lock;-*-

/*! \file VSReflectedRegionAnalysis.hpp

  Integral analysis with reflected region algorithm.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.4 $
  \date       07/29/2007

  $Id: VSReflectedRegionAnalysis.hpp,v 3.4 2010/10/20 03:10:26 matthew Exp $

*/

#ifndef VSREFLECTEDREGIONANALYSIS_HPP
#define VSREFLECTEDREGIONANALYSIS_HPP

#include <VSApertureAnalysis.hpp>
#include <VSPointing.hpp>
#include <VSSimple2DHist.hpp>
#include <VSIntegralAnalysisData.hpp>
#include <VSAAlgebra.hpp>

namespace VERITAS
{
  class VSReflectedRegionAnalysis : public VSApertureAnalysis
  {
  public:
    VSReflectedRegionAnalysis(const std::string bkgnd_model,
			      double bin_size_deg,
			      double theta_cut,
			      const std::pair< double, double >& ring_cut,
			      double theta_max,
			      unsigned max_nregion,
			      const std::string& spmodel);  

    virtual ~VSReflectedRegionAnalysis();

    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisData& o);

    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisDatum& o);

  private:
    
    void accumulate(unsigned iptg, const VSAnalysisStage3Data::RunData& data);

    void getBkgndBins(const VSSimple2DHist<double,double>& h,
		      const VSAAlgebra::Vec2D& coord,
		      const VSAAlgebra::Vec2D& obs_xy,
		      std::vector< unsigned > &bins);

    void getSignalBins(const VSSimple2DHist<double,double>& h,
		       const VSAAlgebra::Vec2D& coord,
		       const VSAAlgebra::Vec2D& obs_xy,
		       std::vector< unsigned > &bins);
  };
};

#endif // VSREFLECTEDREGIONANALYSIS_HPP
