//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage2.hpp

  Driver for stage 2 analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/10/2006

  $Id: VSAnalysisStage2.hpp,v 3.42 2010/03/23 00:47:47 matthew Exp $

*/

#ifndef VSANALYSISSTAGE2_HPP
#define VSANALYSISSTAGE2_HPP

#include <WhippleCams.h>
#include <VSOptions.hpp>
#include <VSDataReductionTypes.hpp>

#include <VBFAnalysisStage2.hpp>

namespace VERITAS
{

  // --------------------------------------------------------------------------
  // __   _____   _             _         _    ___ _                 ___
  // \ \ / / __| /_\  _ _  __ _| |_  _ __(_)__/ __| |_ __ _ __ _ ___|_  )
  //  \ V /\__ \/ _ \| ' \/ _` | | || (_-< (_-<__ \  _/ _` / _` / -_)/ /
  //   \_/ |___/_/ \_\_||_\__,_|_|\_, /__/_/__/___/\__\__,_\__, \___/___|
  //                              |__/                     |___/
  // --------------------------------------------------------------------------

  class VSAnalysisStage2
  {
  public:
    VSAnalysisStage2();
    ~VSAnalysisStage2();

    void runStage2(const std::string& vbf_filename, 
		   VSAnalysisStage1Data* stage1, 
		   VSOctaveH5WriterStruct* writer,
		   const SEphem::SphericalCoords& earth_position,
		   double earth_elevation,
		   VSAnalysisStage1Data* pad_stage1 = 0,
		   bool no_db = false, bool no_l3 = false,
		   bool no_threads = false, bool no_verbose = false,
		   bool do_use_overflow = false);

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="s2_");

  private:
    VSAnalysisStage2(const VSAnalysisStage2&);
    VSAnalysisStage2& operator= (const VSAnalysisStage2&);

    std::vector<VBFAnalysisStage2::Coords>
    getTargets(std::vector<target_type> targets, 
	       const VSTargetTable::Observation& observation,
	       const VSTargetTable& target_table, const VSTime& approx_time,
	       double wobble_region_radius);

    class Options
    {
    public:
      typedef quad<unsigned,double,double,unsigned> QualityCuts;
      typedef VBFAnalysisStage2::LinearCoefficients LinearCoefficients;

      Options();
      std::vector<unsigned>    integration_lo_zero_sample;
      std::vector<unsigned>    integration_hi_zero_sample;
      double                   integration_threshold_frac;
      double                   integration_window_start;
      unsigned                 integration_window_width;
      bool                     integration_apply_increase;
      double                   integration_threshold_charge;
      double                   integration_window_start_increase;
      double                   integration_window_width_increase;
      std::string              ped_suppress_mode;
      double                   ped_suppress_hi;
      double                   ped_suppress_lo;
      double                   ped_suppress_fraction;
      double                   ped_suppress_scale;
      unsigned                 ped_event_count_min;
      double                   gain_suppress_hi;
      double                   gain_suppress_lo;
      bool                     no_laser;
      bool                     permissive_laser;
      bool                     unity_gains;
      std::vector<std::string> primary_cleaning_args;
      std::vector<std::string> secondary_cleaning_args;
      std::string              cleaning_scale;
      bool                     all_events;
      unsigned                 nimage_cut;
      unsigned                 nscope_cut;
      unsigned                 l2_trigger_threshold;
      unsigned                 l3_trigger_threshold;
      std::vector<double>      scope_gain;
      std::vector<LinearCoefficients> scope_dev_scaling;
      std::vector<double>      scope_hi_lo_gain_ratio;
      std::vector<unsigned>    scope_suppress;
      bool                     auto_suppress_l2_chan;
      bool                     pad_zero_suppressed_chan;
      std::vector<std::pair<unsigned,unsigned> >  chan_no_pmt;
      std::vector<std::pair<unsigned,unsigned> >  chan_suppress;
      std::vector<pos_type>    scope_pos;
      bool                     no_set_scope_pos_from_simulations;
      std::string              camera;
      std::vector<double>      cam_rotation;
      std::string              atmosphere;
      unsigned                 method;
      std::string              weighting;
      double                   theta_cut;
      bool                     no_diagnostics;
      bool                     no_slow_diagnostics;
      unsigned                 print_frequency;
      std::string              pointing_string;
      std::vector<std::pair<unsigned,std::string> > 
                               tracking_recorrections_file;
      std::vector<std::pair<std::string,std::string> >
                               tracking_recorrections_date;
      std::string              tracking_target;
      bool                     no_commanded_target;
      double                   wobble_region_radius;
      bool                     enable_l2_corrections;
      std::vector<target_type> theta_targets;
      bool                     no_reorder;
#ifndef NOTHREADS
      unsigned                 nthreads;
      bool                     no_vbf_reader_thread;
#endif
      unsigned                 npackets;
      bool                     no_muon_analysis;
      double                   muon_raw_ring_radius;
      unsigned                 muon_nimage_cut;
      double                   muon_radius_min_cut;
      double                   muon_radius_max_cut;
      double                   muon_rms_max_cut;
      double                   muon_ring_edge_dist_max_cut;
      double                   muon_centroid_radius_ratio_max_cut;
#ifdef MUON_TEST_NON_UNIFORMITY
      double                   muon_non_uniformity_beta;
#endif
      QualityCuts              primary_qc;
      QualityCuts              secondary_qc;
      std::string              nspace_cuts;
      std::vector<std::string> software_trigger_masks;
      std::vector<std::string> qc_masks;
      std::pair<unsigned,std::string> simple_trigger_masks;
      std::pair<double,double> limited_dt_range;
      std::string              sc_parameter_lookup_file;
      double                   msc_weight_power;
      std::vector<double>      psf_poly_radial;
      std::vector<double>      psf_poly_tangential;
     };
    
    static Options s_default_options;
    Options m_options;
  };

} // namespace VERITAS

#endif // not defined VSANALYSISSTAGE2_HPP
