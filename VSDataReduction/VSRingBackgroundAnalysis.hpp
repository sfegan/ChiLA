//-*-mode:c++; mode:font-lock;-*-

/*! \file VSRingBackgroundAnalysis.hpp

  Integral analysis with ring background algorithm.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.5 $
  \date       07/29/2007

  $Id: VSRingBackgroundAnalysis.hpp,v 3.5 2010/10/20 03:10:06 matthew Exp $

*/

#ifndef VSRINGBACKGROUNDANALYSIS_HPP
#define VSRINGBACKGROUNDANALYSIS_HPP

#include <VSIntegralAnalysis.hpp>
#include <VSApertureAnalysis.hpp>
#include <VSPointing.hpp>
#include <VSSimple2DHist.hpp>
#include <VSAAlgebra.hpp>

namespace VERITAS
{

  //! Integral stage3 analysis class responsible for ring background
  //! method (RBM).
  class VSRingBackgroundAnalysis : public VSApertureAnalysis
  {
  public:
    VSRingBackgroundAnalysis(const std::string bkgnd_model,
			     double bin_size_deg,
			     double theta_cut,
			     const std::pair< double, double >& ring_cut,
			     double max_offset_cut,
			     unsigned max_nregion,
			     const std::string& spmodel);  

    virtual ~VSRingBackgroundAnalysis();

    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisData& o);

    virtual void analyze(const std::vector<VSIntegralAnalysis::Data>& d,
			 const VSAnalysisStage3Data& data,
			 const VSAcceptanceData& acceptance,
			 VSIntegralAnalysisDatum& o);

    void accumulate(unsigned iptg, const VSAnalysisStage3Data::RunData& data);
    
  private:    
    double integrateAcceptance(const VSSimple2DHist<double,double>& h,
			       const std::vector<unsigned>& bins);
    
    std::vector< std::pair< int, int > > m_signal_offsets;
    std::vector< std::pair< int, int > > m_bkgnd_offsets;
    VSSimple2DHist<double,double>        m_acceptance_hist;


    VSBinCalcLinear<double>              m_xbinner;
    VSBinCalcLinear<double>              m_ybinner;
  };
};

#endif // VSRESULTSCALCRINGBACKGROUND_HPP
