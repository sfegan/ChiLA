//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage1.hpp

  Stage 1 analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/12/2006

  $Id: VSAnalysisStage1.hpp,v 3.7 2008/12/15 21:19:18 matthew Exp $

*/

#ifndef VSANALYSISSTAGE1_HPP
#define VSANALYSISSTAGE1_HPP

#include<SphericalCoords.h>
#include<VSOptions.hpp>
#include<VSOctaveIO.hpp>

#include<VSRunInfoData.hpp>
#include<VSTimeDepSupData.hpp>
#include<VSTimeDepPedData.hpp>
#include<VBFTimeDepNSB.hpp>
#include<VSPointing.hpp>
#include<VSLaserData.hpp>
#include<VSAnalysisData.hpp>
#include<VSMiscellaneousDBData.hpp>
#include<VSHiLoData.hpp>

namespace VERITAS
{

  class VSAnalysisStage1
  {
  public:
    VSAnalysisStage1();
    ~VSAnalysisStage1();

    class Data
    {
    public:
      Data()
	:analyze(), run_info(), nsb(), suppress(), pedestals(), 
	 suppress_event_slice(), pedestal_time_slice(), 
	 l3_pointing(), db_pointing(), laser(), target_table(),
	 misc_db(), sim_info(), hilo()
      { /* nothing to see here */ }
      ~Data();

      VSAnalysisData                  analyze;
      VSRunInfoData                   run_info;
      VBFTimeDepNSB::Data             nsb;
      VSTimeDepSupData                suppress;
      VSTimeDepPedData                pedestals;
      std::vector<VBFRunInfo::Slice>  suppress_event_slice;
      std::vector<VBFRunInfo::Slice>  pedestal_time_slice;
      VSL3PointingData*               l3_pointing;
      VSDBPointingData*               db_pointing;
      VSLaserData*                    laser;
      VSTargetTableData*              target_table;
      VSMiscellaneousDBData*          misc_db;
      VSSimInfoData*                  sim_info;
      VSHiLoData*                     hilo;      

      void clear();
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;

    private:
      Data(const Data&);
      Data& operator= (const Data&);
    };

    void runStage1(const std::string& vbf_filename, Data& data,
		   const SEphem::SphericalCoords& earth_position,
		   bool no_db = false, bool no_l3 = false,
		   bool no_threads = false, bool no_verbose = false,
		   bool do_use_overflow = false);

    static void configure(VSOptions& options, 
			  const std::string& profile = "",
			  const std::string& opt_prefix="s1_");

  private:
    VSAnalysisStage1(const VSAnalysisStage1&);
    VSAnalysisStage1& operator= (const VSAnalysisStage1&);

    // CONFIGURATION PARAMETERS -----------------------------------------------
    unsigned              m_print_frequency;
    unsigned              m_npackets;
    bool                  m_no_pedestals_in_core;
    unsigned              m_nsb_window_start;
    unsigned              m_nsb_window_width;
    unsigned              m_events_per_slice;
    std::vector<double>   m_nsb_suppress_dev;
    unsigned              m_nsb_min_events_per_slice;
    bool                  m_no_nsb_suppress;
    double                m_ped_slice_time;
    unsigned              m_ped_min_window;
    unsigned              m_ped_sample_separation;
    bool                  m_ped_no_suppress;

    // STATIC DEFAULT PARAMETERS ----------------------------------------------
    static unsigned              s_default_print_frequency;
    static unsigned              s_default_npackets;
    static bool                  s_default_no_pedestals_in_core;
    static unsigned              s_default_nsb_window_start;
    static unsigned              s_default_nsb_window_width;
    static unsigned              s_default_events_per_slice;
    static std::vector<double>   s_default_nsb_suppress_dev;
    static unsigned              s_default_nsb_min_events_per_slice;
    static bool                  s_default_no_nsb_suppress;
    static double                s_default_ped_slice_time;
    static unsigned              s_default_ped_min_window;
    static unsigned              s_default_ped_sample_separation;
    static bool                  s_default_ped_no_suppress;
  };

  typedef VSAnalysisStage1::Data VSAnalysisStage1Data;

} // namespace VERITAS

#endif // VSANALYSISSTAGE1_HPP
