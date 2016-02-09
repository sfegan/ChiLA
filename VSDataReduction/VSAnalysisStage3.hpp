//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage3.hpp

  Stage 3 analysis

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.37 $
  \date       07/29/2007

  $Id: VSAnalysisStage3.hpp,v 3.37 2010/10/20 03:57:55 matthew Exp $

*/

#ifndef VSANALYSISSTAGE3_HPP
#define VSANALYSISSTAGE3_HPP

#include <string>
#include <VSOptions.hpp>
#include <VSCutsCalc.hpp>
#include <VSPointing.hpp>
#include <VSResultsSimData.hpp>
#include <VSCutsEvaluator.hpp>
#include <VSDatumElementExtractor.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSHistAccumulator.hpp>
#include <VSScaledParameterCalc.hpp>
#include <VSEnergyCalcLT.hpp>
#include <VSEventDataVisitor.hpp>
#include <VSAnalysisStage3Visitor.hpp>

namespace VERITAS
{
  class VSAnalysisStage3
  {
  public:
    VSAnalysisStage3();
    ~VSAnalysisStage3();

    void runStage3(const std::list<std::string>& stg2_files,
		   VSOctaveH5Writer* writer);

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="s3_");

  private:
    VSAnalysisStage3(const VSAnalysisStage3&);
    VSAnalysisStage3& operator= (const VSAnalysisStage3&);

    void addHist(const std::string& hist_def);    
 
    VSEventDataReader::MemberSubset  m_array_subset;
    VSEventDataReader::MemberSubset  m_scope_subset;

    VSCutsEvaluator*                 m_cuts_evaluator;

    VSEventDataVisitor*              m_visitor;

    class Options
    {
    public:
      Options();
      std::string                 hist_file;
      double                      max_offset_cut;
      double                      min_energy_gev;
      double                      max_energy_gev;
      double                      sky_bin_width_deg;
      double                      source_exclusion_radius;
      double                      star_exclusion_vmag_limit;
      triple< std::string, double,double > source_pos;
      std::string                 source_name;
      bool                        no_sim;
      bool                        no_run_results;
      bool                        write_run_hists;
      unsigned                    nevents;
      std::string                 coord_system;
      unsigned                    rng_seed;
    };

    SEphem::SphericalCoords                 m_src_radec_J2000;    
    std::vector< SEphem::SphericalCoords >  m_ptg_radec_J2000;

    static Options s_default_options;
    Options m_options;
  };

}  // namespace VERITAS


#endif // VSANALYSISSTAGE3_HPP
