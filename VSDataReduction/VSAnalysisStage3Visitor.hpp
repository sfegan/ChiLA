//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage3Visitor.hpp

  Visitor class for stage3.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       12/17/2009

  $Id: VSAnalysisStage3Visitor.hpp,v 3.7 2010/10/20 03:39:36 matthew Exp $

*/

#ifndef VSANALYSISSTAGE3VISITOR_HPP
#define VSANALYSISSTAGE3VISITOR_HPP

#include <VSEventDataVisitor.hpp>
#include <VSIntegralAnalysis.hpp>

namespace VERITAS
{
  /////////////////////////////////////////////////////////////////////////////
  //! Visitor class for stage3 analysis.
  /////////////////////////////////////////////////////////////////////////////
  class VSAnalysisStage3Visitor : public VSEventDataVisitor
  {
  public:

    struct Run
    {
      Run(): iptg(), ptg_radec() { }

      unsigned                iptg;
      SEphem::SphericalCoords ptg_radec;
    };

    VSAnalysisStage3Visitor(VSCutsEvaluator* cuts_evaluator,
			    VSOctaveH5WriterStruct* writer,
			    double sky_bin_width_deg,
			    double theta_max,
			    const std::string& coord_system,
			    const std::string& src_name,
			    const SEphem::SphericalCoords& origin_radec,
			    const SEphem::SphericalCoords& src_radec,
			    const std::vector< SEphem::SphericalCoords >&
			    ptg_radec,
			    const VSExclusionRegion& exclusion_region,
			    const std::map< unsigned, Run >& run_list,
			    bool no_run_results,
			    bool write_run_hists,
			    unsigned rng_seed);
    virtual ~VSAnalysisStage3Visitor();

    virtual void visitRun(const VSAnalysisStage1Data& stage1,
			  const VSTargetTable::Observation& obs,
			  const VSArrayMergedCalibrationData& cal);

    virtual void leaveRun();

    virtual void visitEvent(const VSEventArrayDatum& event);
    virtual void leaveEvent();

    virtual void visitSimEvent(const VSArraySimulationDatum& sim);
    virtual void leaveSimEvent();

    virtual void visitSimHeader(const VSHeaderSimulationDatum& header);
    virtual void leaveSimHeader();

  private:

    void finalizeSim();
    double fitThreshold(const VSSimpleGraph<double, double>& gr);
    void calcRate(double& rate, double& rate_err,
		  const VSSimpleGraph<double, double>& effarea,
		  VSSpectrumFn* spectrum);

    void calcParallacticAngle(const VSTime& lo_time, const VSTime& hi_time,
			      const SEphem::SphericalCoords& radec,
			      double& pangle_mean, double& pangle_rms);


    VSIntegralAnalysis*               m_integral_analysis;
    VSCutsEvaluator*                  m_cuts_evaluator;
    VSOctaveH5WriterStruct*           m_writer;
    VSScaledParameterCalc*            m_sp_calc;
    VSEnergyCalcLT*                   m_egy_calc;
    VSSimEnergyWeightCalc*            m_egywt_calc;
    VSSourceInjector*                 m_source_injector;
    VSSpectrumCalc*                   m_spectrum_calc;
    RandomNumbers*                    m_rng;

    VSAnalysisStage3Data::RunData*    m_data;
    std::vector< VSAnalysisStage3Data::RunData* > m_run_data;

    VSStage3SimDatum*                 m_sim_data;
    std::vector< VSStage3SimDatum* >  m_sim_run_data;
    std::vector< VSStage3SimArrayTableDatum* > m_sim_table_data;

    VSIntegralAnalysis::Data          m_intanl_data;
    std::vector< VSIntegralAnalysis::Data > m_intanl_run_data;
    VSSpectrumData                    m_spectrum_data;

    VSEventArrayDatum                 m_event_data;
    VSArraySimulationDatum            m_sim;

    std::vector<unsigned>             m_nchan;
    double                            m_sky_bin_width_deg;
    double                            m_offset_max;
    double                            m_theta_cut;    
    double                            m_spectrum_theta_cut;
    std::pair< double, double >       m_ring_cut;
    unsigned                          m_max_nregions;
    std::string                       m_source_spectrum;
    std::string                       m_coord_system;
    std::string                       m_src_name;
    VSAAlgebra::Vec2D                 m_src_xy;
    SEphem::SphericalCoords           m_origin_radec;
    SEphem::SphericalCoords           m_src_radec;
    std::vector< SEphem::SphericalCoords > m_ptg_radec;
    std::map< unsigned, Run >         m_run_list;
    bool                              m_no_run_results;
    bool                              m_write_run_hists;
    unsigned                          m_rng_seed;

    double                            m_log10_egy_min;
    double                            m_log10_egy_max;
    double                            m_egy_bin_width;

    bool                              m_has_sim;
    unsigned                          m_event_code;
  };
}

#endif // VSANALYSISSTAGE3VISITOR_HPP
