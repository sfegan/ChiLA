//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage2.hpp

  Stage 2 analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/10/2006

  $Id: VBFAnalysisStage2.hpp,v 3.42 2010/03/23 00:49:12 matthew Exp $

*/

#ifndef VBFANALYSISSTAGE2_HPP
#define VBFANALYSISSTAGE2_HPP

#include <stdint.h>
#include <vector>
#include <map>

#include <SphericalCoords.h>
#include <RandomNumbers.hpp>
#include <VSAReconstruction.hpp>
#include <VSSimpleVBF.hpp>
#include <VSCleaner.hpp>
#include <VSTimingCalc.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSOctaveIO.hpp>
#include <VSSimpleStat.hpp>
#include <VSTime.hpp>
#include <VSEventData.hpp>
#include <VSMuonAnalysis.hpp>
#include <VSDiagnosticsData.hpp>
#include <VSMergedCalibrationData.hpp>
#include <VSNSpaceCutsCalc.hpp>
#include <VSSimulationData.hpp>
#include <VSChannelMap.hpp>
#include <VSLaserData.hpp>
#include <VSHiLoData.hpp>
#include <VSEnergyCalcLT.hpp>
#include <VSScaledParameterCalc.hpp>
#include <VSSimCoordTransform.hpp>

#define vstream std::cout

#define __CPC(x) x(o.x)
#define __CPA(x) x = o.x;

namespace VERITAS
{

  // --------------------------------------------------------------------------
  // __   _____ ___ _             _         _
  // \ \ / / _ ) __/_\  _ _  __ _| |_  _ __(_)___
  //  \ V /| _ \ _/ _ \| ' \/ _` | | || (_-< (_-<
  //   \_/ |___/_/_/ \_\_||_\__,_|_|\_, /__/_/__/
  //                                |__/
  // --------------------------------------------------------------------------

  class VBFAnalysisStage2: public VSSimpleVBFVisitor
  {
  public:
    
    typedef VSTargetTable::Target Coords;
    typedef std::vector<std::pair<unsigned,unsigned> > ChanList;
    enum CleaningScale { CS_UNITY, CS_GAIN, CS_PEDRMS };
    typedef std::pair<VSNSpace*,VSNSpace*> SCParameterPair;
    typedef std::pair<VSCleaner*,VSAReconstruction::ArrayQualityCuts*> SecondaryCleaning;
    typedef std::pair<double,double> LinearCoefficients;

    struct Settings
    {
      Settings():
        earth_position(), process_all_events(), enable_l2_corrections(),
        method(), weighting(),
	theta_targets(), theta_cut(), nimage_cut(), nscope_cut(),
        do_not_write_diagnostics(), no_slow_diagnostics(),
	l2_trigger_threshold(), 
	l3_trigger_threshold(), print_frequency(), cleaning_scale(),
        ped_event_count_min(), permissive_laser(), unity_gain(),
	auto_suppress_l2_chan(true), pad_zero_suppressed_chan(true),
	chan_suppress(), chan_has_no_pmt(), 
	scope_suppress(), scope_gain(), scope_dev_scaling(), 
	scope_hi_lo_gain_ratio(),
        integration_lo_zero_sample(), integration_hi_zero_sample(),
        integration_threshold_frac(), integration_window_start(),
        integration_window_width(), integration_apply_increase(),
        integration_threshold_charge(), integration_window_start_increase(),
        integration_window_width_increase(), neighbor_nchan(), neighbors(),
        verbose(), software_trigger_masks(), qc_masks(),
	limited_dt_range(), 
	psf_poly_radial(), psf_poly_tangential()
      { /* nothing to see here */ }
	
      SEphem::SphericalCoords              earth_position;
      bool                                 process_all_events;
      bool                                 enable_l2_corrections;
      VSAReconstruction::Method            method;
      VSAReconstruction::ScopeWeighting    weighting;
      std::vector<Coords>                  theta_targets;
      double                               theta_cut;
      unsigned                             nimage_cut;
      unsigned                             nscope_cut;
      bool                                 do_not_write_diagnostics;
      bool                                 no_slow_diagnostics;
      unsigned                             l2_trigger_threshold;
      unsigned                             l3_trigger_threshold;
      unsigned                             print_frequency;
      CleaningScale                        cleaning_scale;
      unsigned                             ped_event_count_min;
      bool                                 permissive_laser;
      bool                                 unity_gain;
      bool 	                           auto_suppress_l2_chan;
      bool                                 pad_zero_suppressed_chan;
      ChanList                             chan_suppress;
      ChanList                             chan_has_no_pmt;
      std::vector<bool>                    scope_suppress;
      std::vector<double>                  scope_gain;
      std::vector<LinearCoefficients>      scope_dev_scaling;
      std::vector<double>                  scope_hi_lo_gain_ratio;
      std::vector<unsigned>                integration_lo_zero_sample;
      std::vector<unsigned>                integration_hi_zero_sample;
      double                               integration_threshold_frac;
      double                               integration_window_start;
      unsigned                             integration_window_width;
      bool                                 integration_apply_increase;
      double                               integration_threshold_charge;
      double                               integration_window_start_increase;
      double                               integration_window_width_increase;
      unsigned                             neighbor_nchan;
      const int                          (*neighbors)[NUM_NEIGHBORS];
      bool                                 verbose;
      std::vector<unsigned>                software_trigger_masks;
      std::vector<unsigned>                qc_masks;
      std::pair<double,double>             limited_dt_range;
      VSAAlgebra::VecND                    psf_poly_radial;
      VSAAlgebra::VecND                    psf_poly_tangential;
    };

    // ------------------------------------------------------------------------
    // PRIMARY VBF VISITOR API
    // ------------------------------------------------------------------------

    VBFAnalysisStage2(VSCleaner* clean,
		      const SecondaryCleaning& secondary_clean,
		      VSPointing* pointing,
		      VSAReconstruction* reconstruction, 
		      std::vector<VSMuonAnalysis*> muon_analysis,
		      VSScaledParameterCalc* sp_calc,
		      VSEnergyCalc* energy_calc,
		      VSCutsCalc* cuts,
		      const VSAnalysisStage1Data& stage1, 
		      const VSChannelMap& channel_map,
		      const Settings& settings,
		      VSOctaveH5Writer* io,
		      const VSAnalysisStage1Data* pad_stage1);
    virtual ~VBFAnalysisStage2();

#ifndef NOTHREADS
    virtual void usingThreads(unsigned nthreads);
#endif

    virtual void visitPacket(bool& veto_packet, void*& user_data,
			     uint32_t                   seq_packet_number,
			     uint32_t                   vbf_packet_number,
			     bool                       has_array_event,
			     bool                       has_sim_header,
			     bool                       has_sim_event,
			     uint32_t                   num_overflow_datum,
			     const VBFPacket*           packet);

    virtual void leavePacket(bool veto_packet, void* user_data);

    virtual void visitArrayEvent(bool& veto_array_event, void* user_data,
				 uint32_t               num_scope_events,
				 bool                   has_array_trigger,
				 bool                   has_event_number,
				 uint32_t               event_num,
 				 bool                   has_good_event_time,
				 const VSTime&          best_event_time,
				 EventType              event_type,
				 uint32_t               l2_trigger_mask,
				 const VArrayEvent*     array_event);

    virtual void leaveArrayEvent(bool veto_array_event, void* user_data);
  
    virtual void visitArrayTrigger(bool& veto_array_event, void* user_data,
				   uint32_t             event_num,
				   const VEventType&    event_type,
				   uint32_t             trigger_mask,
				   uint32_t             flags,
				   const VSTime&        raw_time,
				   uint32_t             at_flags,
				   uint32_t             config_mask,
				   uint32_t             num_telescopes,
				   uint32_t             num_trigger_telescopes,
				   uint32_t             run_number,
				   const uint32_t*      ten_mhz_clocks,
				   const uint32_t*      cal_count,
				   const uint32_t*      ped_count,
				   const VArrayTrigger* trigger);
    
    virtual void visitArrayTriggerScope(bool& veto_array_event, 
					void* user_data,
					uint32_t        telescope_num,
					bool            triggered,
					uint32_t        event_type,
					float           altitude,
					float           azimuth,
					uint32_t        tdc,
					uint32_t        shower_delay,
					uint32_t        comp_delay,
					const uint32_t* l2_counts,
					const uint32_t* cal_counts);

    virtual void visitScopeEvent(bool& veto_scope_event, void* user_data,
				 uint32_t               event_num,
				 uint32_t               telescope_num, 
				 const VEventType&      event_type,
				 uint32_t               trigger_mask,
				 uint32_t               flags,
				 const VSTime&          raw_time,
				 uint32_t               num_samples,
				 uint32_t               num_channels_saved,
				 uint32_t               num_channels_total,
				 uint32_t               num_clock_trigger,
				 const VEvent*          event);

    virtual void leaveScopeEvent(bool veto_scope_event, void* user_data);

    virtual void visitChannel(bool& veto_channel, void* user_data,
			      uint32_t                  channel_num, 
			      bool                      hit, 
			      bool                      trigger);

    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);

    virtual void visitSimulationHeader(void* user_data,
				       uint32_t         run_number,
				       const VSTime&    date_of_sims,
				       uint32_t         simulation_package,
				       const std::string& simulator_name,
				       const VSTime&    date_of_array,
				       uint32_t         corsika_atm_model,
				       float            obs_altitude_m,
				       const std::vector<VArrayConfiguration>&
				                        array,
				       const std::string& stringified_config);

    virtual void visitSimulationEvent(bool& veto_packet, void* user_data,
				      uint32_t          run_number,
				      uint32_t          event_num,
				      uint32_t          corsika_particle_id,
				      float             energy_gev,
				      float             obs_zenith_deg,
				      float             obs_azimuth_deg,
				      float             primary_zenith_deg,
				      float             primary_azimuth_deg,
				      float             ref_zenith_deg,
				      float             ref_azimuth_deg,
				      float             ref_position_angle_deg,
				      float             core_east_m,
				      float             core_south_m,
				      float             core_elevation_asl_m);

    virtual void visitChiLAHeader(void* user_data,
				  const std::string& database_name,
				  const std::vector<VChiLASimParamTable>&
				                   sim_param_tables,
				  const std::vector<VChiLAOpticsConfiguration>&
				                   optics_configurations,
				  const std::vector<VChiLAElectronicsConfiguration>& 
				                   electronics_configurations);

    virtual void visitChiLAEvent(bool& veto_packet, void* user_data,
				 uint32_t               table_index,
				 uint32_t               electronics_id,
				 uint32_t               table_event_index,
				 float                  a_omega,
				 float                  a_omega_var);

    virtual void visitKascadeHeader(void* user_data,
				    uint32_t            corsika_particle_id,
				    float               energy_gev,
				    uint32_t            shower_id,
				    float               x_area_width_m,
				    float               y_area_width_m,
				    uint32_t            north_south_grid);
 
    virtual void visitKascadeEvent(bool& veto_packet, void* user_data,
				   int32_t      nx,
				   int32_t      ny,
				   int32_t      direction_index,
				   float        emission_altitude_m,
				   float        emission_altitude_sigma,
				   float        muon_ratio,
				   float        a_omega,
				   float        rel_tel_trigger_time_ns,
				   float        differential_rate_per_event_hz,
				   float        integral_rate_per_event_hz);

  private:
    VBFAnalysisStage2(const VBFAnalysisStage2&);
    VBFAnalysisStage2& operator=(const VBFAnalysisStage2&);

  protected:

    // ------------------------------------------------------------------------
    //  __  __                     _  ___      _
    // |  \/  |___ _ _ __ _ ___ __| |/ __|__ _| |
    // | |\/| / -_) '_/ _` / -_) _` | (__/ _` | |
    // |_|  |_\___|_| \__, \___\__,_|\___\__,_|_|
    //                |___/
    // ------------------------------------------------------------------------

    class ChanMergedCal: public VSChannelMergedCalibrationData
    {
    public:
      ChanMergedCal(unsigned nchan=0):
	VSChannelMergedCalibrationData(nchan),
	clean_multiplier(0), suppress(true), dopad(false), dev(0), pad(0)
      { /* nothing to see here */ }

      // These variables change from timeslice to timeslice -------------------
      
      double                                   clean_multiplier;
      bool                                     suppress;
      bool                                     dopad;
      double                                   dev;
      double                                   pad;
    };

    class ScopeMergedCal: public VSScopeMergedCalibrationBase
    {
    public:
      ScopeMergedCal(unsigned nchan=0):
	VSScopeMergedCalibrationBase(nchan), channel(nchan,nchan),
	ch_has_pmt(nchan), ch_ped_hi(nchan), ch_ped_lo(nchan), 
	ch_hi_lo_gain_ratio(nchan),
	ch_l2chan(nchan), ch_l2time(nchan), ch_chantime(nchan), 
	ch_suppress(nchan), ch_suppress_all_events(nchan), 
	ch_suppress_lo_gain(nchan),
	ch_dopad(nchan), ch_pad(nchan),	ch_clean_multiplier(nchan), 
	ch_gain(nchan)
      { /* nothing to see here */ }

      std::vector<ChanMergedCal>               channel;

      std::vector<bool>                        ch_has_pmt;
      std::vector<double>                      ch_ped_hi;
      std::vector<double>                      ch_ped_lo;
      std::vector<double>                      ch_hi_lo_gain_ratio;
      std::vector<uint16_t>                    ch_l2chan;      
      std::vector<double>                      ch_l2time;	      
      std::vector<double>                      ch_chantime;	      
      std::vector<bool>                        ch_suppress;
      std::vector<bool>                        ch_suppress_all_events;
      std::vector<bool>                        ch_suppress_lo_gain;
      std::vector<bool>                        ch_dopad;
      std::vector<double>                      ch_pad;
      std::vector<double>                      ch_clean_multiplier;
      std::vector<double>                      ch_gain;
    };

    class ArrayMergedCal: public VSArrayMergedCalibrationBase
    {
    public:
      ArrayMergedCal(const VSRunInfoData& run_info,
		     const VSTimeDepSupData* sup,
		     const VSTimeDepPedData* ped,
		     const std::vector<unsigned> hi_sample_zero,
		     const std::vector<unsigned> lo_sample_zero,
		     CleaningScale cs, unsigned window_size,
		     const VSChannelMap* channel_map,
		     const VSAReconstruction::ArrayInfo& array_info,
		     const VSLaserData* laser = 0,
		     bool permissive_laser = false,
		     bool unity_gain = false,
		     const VSHiLoData* hilo = 0,
		     const std::vector<double>& scope_hi_lo_gain_ratio =
		     std::vector<double>(),
		     unsigned ped_event_count_min = 0,
		     const std::vector<bool>& sup_scope = std::vector<bool>(),
		     const std::vector<double>& gain = std::vector<double>(),
		     const bool auto_suppress_l2_chan = true,
		     const ChanList& sup_chan = ChanList(),
		     const ChanList& no_pmt_chan = ChanList(),
		     const VSTimeDepSupData* pad_sup = 0,
		     const VSTimeDepPedData* pad_ped = 0,
		     const VSMiscellaneousDBData* db_data = 0,
		     const std::vector<LinearCoefficients>& dev_scaling
		     = std::vector<LinearCoefficients>());

      ArrayMergedCal(const ArrayMergedCal& o):
	VSArrayMergedCalibrationBase(o),
	__CPC(scope),
	__CPC(m_sup),
	__CPC(m_ped),
	__CPC(m_cs),
	__CPC(m_window_size),
	__CPC(m_pad_sup),
	__CPC(m_pad_ped),
	__CPC(m_sup_islice),
	__CPC(m_ped_islice),
	__CPC(m_ped_event_count_min)
      { /* nothing to see here */ }

      ArrayMergedCal& operator=(const ArrayMergedCal& o)
      {
	*static_cast<VSArrayMergedCalibrationBase*>(this) = o;
	__CPA(scope);
	__CPA(m_sup);
	__CPA(m_ped);
	__CPA(m_cs);
	__CPA(m_window_size);
	__CPA(m_pad_sup);
	__CPA(m_pad_ped);
	__CPA(m_sup_islice);
	__CPA(m_ped_islice);
	__CPA(m_ped_event_count_min);
	return *this;
      }
      
      inline void setEventNumber(unsigned event_num);
      void print(const std::ostream& stream);
      void save(VSOctaveH5WriterStruct* s) const;      

      std::vector<ScopeMergedCal>              scope;

      void setFixedCacheVariables();

    private:
      void mergeSupSlice();
      void mergePedSlice();

      const VSTimeDepSupData*                  m_sup;
      const VSTimeDepPedData*                  m_ped;
      CleaningScale                            m_cs;
      unsigned                                 m_window_size;
      const VSTimeDepSupData*                  m_pad_sup;
      const VSTimeDepPedData*                  m_pad_ped;
      unsigned                                 m_sup_islice;
      unsigned                                 m_ped_islice;
      unsigned                                 m_ped_event_count_min;
    };

    // ------------------------------------------------------------------------
    //  ___  _                       _   _
    // |   \(_)__ _ __ _ _ _  ___ __| |_(_)__ ___
    // | |) | / _` / _` | ' \/ _ (_-<  _| / _(_-<
    // |___/|_\__,_\__, |_||_\___/__/\__|_\__/__/
    //             |___/
    // ------------------------------------------------------------------------

    class IED;
    class ScopeIED;

    class ArrayDiagnostics;

    class PartialScopeDiagnostics: public VSPartialScopeDiagnosticsBase
    {
    public:
      PartialScopeDiagnostics(unsigned nchan, unsigned nsample);
      virtual ~PartialScopeDiagnostics();

      inline void recordTriggeredChannel(unsigned ichan,
					 bool l2_trigger_received);

      inline void recordHitChannel(bool trigger_l2,
				   unsigned ichan, bool lo_gain, bool trigger,
				   unsigned peak, double sij, double tij, 
				   unsigned nsample,
				   const uint32_t* samples,
				   bool slow_diagnostics);

      inline void recordImageChannel(unsigned ichan, bool lo_gain,
				     double sij, double tij, 
				     bool slow_diagnostics);

      inline void recordPedestal(unsigned ichan,
				 unsigned nsample,
				 const uint32_t* samples);

      inline void recordChannelIsTriggeredIsolated(unsigned ichan)
      {
	channel_is_triggered_isolated[ichan]++;
      }

      inline void recordChannelIsTriggeredIsolatedNoL2(unsigned ichan)
      {
	channel_is_triggered_isolated_no_l2[ichan]++;
      }

      inline void recordChannelIsTriggeredInLargestRegion(unsigned ichan)
      {
	channel_is_triggered_in_largest_region[ichan]++;
      }

      inline void recordChannelIsTriggeredInLargestRegionNoL2(unsigned ichan)
      {
	channel_is_triggered_in_largest_region_no_l2[ichan]++;
      }

      inline void recordChannelIsTriggeredInSubThresholdEvent(unsigned ichan)
      {
	channel_is_triggered_in_sub_threshold_event[ichan]++;
      }

      inline void 
      recordChannelIsTriggeredInSubThresholdEventNoL2(unsigned ichan)
      {
	channel_is_triggered_in_sub_threshold_event_no_l2[ichan]++;
      }

      PartialScopeDiagnostics& operator+=(const PartialScopeDiagnostics& o);

    protected:
      PartialScopeDiagnostics(const PartialScopeDiagnostics&);
      PartialScopeDiagnostics& operator= (const PartialScopeDiagnostics&);

      friend class ArrayDiagnostics;

      // DATA -----------------------------------------------------------------

      uint32_t* channel_mean_lo_trace;
      uint32_t* channel_mean_hi_trace;
      uint32_t* channel_mean_lo_trace_no_l2;
      uint32_t* channel_mean_hi_trace_no_l2;
      uint32_t* channel_pedestal_covariance;
      uint32_t* channel_pedestal_mean;
      unsigned* npedestal;

      // NEIGHBORS FOR DO_TRIGGER_DIAGNOSTICS ---------------------------------

      unsigned m_nsample;
    };

    class PartialArrayDiagnostics
    {
    public:
      PartialArrayDiagnostics(const std::vector<bool>& config_mask,
			      const std::vector<unsigned>& nchan, 
			      const std::vector<unsigned>& nsample);
      virtual ~PartialArrayDiagnostics();

      PartialScopeDiagnostics* partialScopeDiagnostics(unsigned iscope)
      {
	return scopes[iscope];
      }

    private:
      PartialArrayDiagnostics(const PartialArrayDiagnostics&);
      PartialArrayDiagnostics& operator= (const PartialArrayDiagnostics&);

      friend class ArrayDiagnostics;

      std::vector<PartialScopeDiagnostics*> scopes;
    };

    class ScopeDiagnostics: 
      public PartialScopeDiagnostics, VSScopeDiagnosticsBase
    {
    public:
      ScopeDiagnostics(unsigned nchan, unsigned nsample, unsigned nscope);
      virtual ~ScopeDiagnostics();
      
      void integrateEventIED(IED* ied, ScopeIED* sied);
      void integratePartial(const std::list<PartialScopeDiagnostics*>& sdiags);
      
      void finalize(const ScopeMergedCal* scal,
		    const VSRunInfoData& run_info,
		    const VSMiscellaneousDBData::Scope* db_scope_data,
		    const std::vector<VATime>* nsb_time,
		    const std::vector<double>* nsb_median,
		    const SEphem::SphericalCoords& earth_position);

      void save(VSOctaveH5WriterStruct* s, bool slow_diagnostics = true) const;

      // DATA -----------------------------------------------------------------

      bool scope_zero_set;
      double scope_az_zero;
      double scope_ra_zero;
      double scope_l_zero;
      std::vector<VSSimpleStat2<double> > scope_zn;
      std::vector<VSSimpleStat2<double> > scope_az;
      std::vector<VSSimpleStat2<double> > scope_ra;
      std::vector<VSSimpleStat2<double> > scope_dec;
      std::vector<VSSimpleStat2<double> > scope_l;
      std::vector<VSSimpleStat2<double> > scope_b;
      std::vector<VSSimpleStat2<double> > scope_gps_diff_hist;

      friend class ArrayDiagnostics;
    };

    class ArrayDiagnostics: public VSArrayDiagnosticsBase
    {
    public:
      ArrayDiagnostics(const std::vector<bool>& config_mask,
		       const std::vector<unsigned>& nchan, 
		       const std::vector<unsigned>& nsample,
		       const std::pair<double,double>& limited_dt_range);
      virtual ~ArrayDiagnostics();

      void integrateEventIED(IED* ied);
      void integratePartial(const std::list<PartialArrayDiagnostics*>& diags);

      void finalize(const ArrayMergedCal* cal,
		    const VSRunInfoData& run_info,
		    const VSMiscellaneousDBData* db_data,
		    const VBFTimeDepNSB::Data* nsb,
		    const SEphem::SphericalCoords& earth_position);

      void save(VSOctaveH5WriterStruct* s, bool slow_diagnostics = true) const;

      SEphem::SphericalCoords getMeanRaDec() const;
      void packetFound() { packets_found++; }

    private:
      ArrayDiagnostics(const ArrayDiagnostics&);
      ArrayDiagnostics& operator= (const ArrayDiagnostics&);

      // SETTINGS -------------------------------------------------------------

      std::pair<double,double>       m_limited_dt_range;

      // DATA -----------------------------------------------------------------

      VSSimpleStat1<double>          m_limited_dt;
      std::vector<VSSimpleStat2<double> > l3_gps_diff_hist;

      std::vector<ScopeDiagnostics*> scope;
    };

    // ------------------------------------------------------------------------
    //  ___     _           _       ___             _   ___       _
    // |_ _|_ _| |_ ___ _ _(_)_ __ | __|_ _____ _ _| |_|   \ __ _| |_ __ _
    //  | || ' \  _/ -_) '_| | '  \| _|\ V / -_) ' \  _| |) / _` |  _/ _` |
    // |___|_||_\__\___|_| |_|_|_|_|___|\_/\___|_||_\__|___/\__,_|\__\__,_|
    // 
    // ------------------------------------------------------------------------

    class ThreadSpecificData
    {
    public:
      ThreadSpecificData(const VSAnalysisStage1Data* stage1,
			 const Settings& settings, const ArrayMergedCal& _cal,
			 const ArrayDiagnostics* diagnostics,
			 bool construct_rng);
      ~ThreadSpecificData();
      static void destructor(void* ptr);

      // DIAGNOSTICS
      PartialArrayDiagnostics*   partial_diagnostics;

      // HELPERS
      std::vector<unsigned>      trace_nsample;
      ArrayMergedCal             cal;
      RandomNumbers*             rng;

    private:
      ThreadSpecificData(const ThreadSpecificData&);
      ThreadSpecificData& operator= (const ThreadSpecificData&);

      const ArrayDiagnostics*    m_diagnostics;
    };

    class ScopeIED
    {
    public:
      ScopeIED(const IED* _ied, unsigned _iscope, unsigned _nchan);
      ~ScopeIED();
      
      // L3 DATA
      bool                       has_l3;
      bool                       l3_sent;
      bool                       l3_l2_received;
      double                     l3_az;
      double                     l3_el;
      uint32_t                   l3_counts_l2;
      uint32_t                   l3_tdc;

      // L2 COUNTS
      int64_t                    counts_l2;

      // VDAQ DATA
      bool                       has_vdaq;
      bool                       processed_scope_event;
      int64_t                    gps_dt;
      unsigned                   nchan_hit;
      unsigned                   sij_max1_raw_j;
      unsigned                   sij_max2_raw_j;
      unsigned                   sij_max3_raw_j;
      unsigned                   trigger_j;
      double                     total_signal;

      // CLEANING
      unsigned                   nchan_logain;
      unsigned                   nchan_trigger;
      unsigned                   nchan_image;
      bool                       has_image;
      unsigned                   sij_max1_sig_j;
      unsigned                   sij_max2_sig_j;
      unsigned                   sij_max3_sig_j;
      double                     sij_max1_sig;
      double                     sij_max2_sig;
      double                     sij_max3_sig;
      unsigned                   largest_region_nchan;

      // TRACKING
      bool                       has_azzn;
      double                     az_rad;
      double                     zn_rad;
      double                     ra_rad;
      double                     dec_rad;
      double                     l_rad;
      double                     b_rad;

      // RECONSTRUCTION
      bool                       has_moments;
      double                     fp_Ni;
      double                     fp_xc;
      double                     fp_yc;
      double                     fp_dist;
      double                     fp_length;
      double                     fp_width;
      double                     fp_psi;
      double                     fp_ex;
      double                     fp_ey;
      double                     fp_disp;
      double                     intrinsic_length;
      double                     intrinsic_width;

      bool                       used_in_reconstruction;
      double                     sc_width;
      double                     sc_length;
      double                     sc_disp;
      double                     lt_log10_energy;
      double                     lt_log10_energy_err;

      // MUON IMAGE
      bool                       has_muon;
      VSMuonAnalysisDatum        muon_data;

      // CONSTANT ITEMS WHICH MAY BE USEFUL
      const IED*const            ied;
      const unsigned             iscope;
      const unsigned             nchan;

      void reset()
      {
	has_l3                 = false;
	l3_sent                = false;
	l3_l2_received         = false;

	has_vdaq               = false;
	processed_scope_event  = false;

	has_image              = false;

	has_azzn               = false;
	has_moments            = false;
	used_in_reconstruction = false;

	has_muon               = false;
      };

      void transferDataToESD(VSEventScopeDatum* esd,
			     const VSAReconstruction::Reconstruction& recon,
			     const VSAReconstruction::ScopeParameters& sp,
			     const VSAReconstruction::ScopeMoments& sm);

    private:
      ScopeIED(const ScopeIED&);
      ScopeIED& operator= (const ScopeIED&);
    };

    class IED
    {
    public:
      typedef VSSimpleVBFVisitor::EventType EventType;

      IED(const VSAnalysisStage1Data* stage1, unsigned ntheta, 
	  VSCutsCalc *cuts_calc);
      ~IED();

      // EVENT STATE
      EventType                  event_type;
      bool                       event_failed_software_trigger;
      uint32_t                   event_l2_trigger_mask;
      uint32_t                   event_sent_trigger_mask;
      uint32_t                   event_has_vdaq_mask;
      unsigned                   scope_num;
      ScopeIED*                  scope_ied;
      const ScopeMergedCal*      scope_cal;
      PartialScopeDiagnostics*   scope_diag;
      unsigned                   scope_nchan;

      // SCOPE IMAGE DATA
      unsigned                   scope_nchan_hit;
      unsigned                   scope_nchan_logain;
      unsigned                   scope_nchan_trigger;
      double                     scope_total_signal;
      bool                       this_chan_trigger;
      std::vector<bool>          chan_has_signal;
      std::vector<bool>          chan_hit;
//      bool*                      chan_has_signal;
//      bool*                      chan_hit;
      double*                    chan_sij;
      double*                    chan_sij_total;
      double*                    chan_tij;
//      bool*                      chan_trigger;
      std::vector<bool>          chan_trigger;
      std::vector<bool>          chan_logain;
      double                     chan_sij_max1_raw;
      double                     chan_sij_max2_raw;
      double                     chan_sij_max3_raw;
      unsigned                   chan_sij_max1_raw_j;
      unsigned                   chan_sij_max2_raw_j;
      unsigned                   chan_sij_max3_raw_j;

      // SCOPE MEAN POSITIONS
      bool                       got_one_position;
      double                     zn_zero;
      double                     az_zero;
      VSSimpleStat1<double>      mean_position_x;
      VSSimpleStat1<double>      mean_position_y;
     
      // SCOPE DATA
      std::vector<ScopeIED*>     scope;

      // THREAD CONTROL
      bool                       ied_processing;
      ThreadSpecificData*        tsd;

      // VBF DATA
      unsigned                   vbf_packet_number;
      unsigned                   seq_packet_number;

      // EVENT DATA
      bool                       has_array_event;
      bool                       processed_array_event;
      unsigned                   event_num;
      bool                       has_event_time;
      VSTime                     event_time;
      double                     event_time_sec;
      double                     event_time_hist;
      int64_t                    ticks_elap;
      int64_t                    ticks_both;
      int64_t                    ticks_vdaq;
      int64_t                    ticks_lev3;

      // L3 DATA
      bool                       has_l3;
      uint32_t                   ten_mhz_elapsed;
      uint32_t                   ten_mhz_veto_both;
      uint32_t                   ten_mhz_veto_vdaq;
      uint32_t                   ten_mhz_veto_l3;
      unsigned                   nscope_trigger;
      unsigned                   nscope_sent_l3;
      int64_t                    gps_dt;

      // VDAQ DATA
      unsigned                   nscope_has_event;

      // CLEANING
      unsigned                   event_has_image_mask;
      unsigned                   nscope_image;
      unsigned                   nscope_quality;

      // RECONSTRUCTION
      unsigned                   event_used_in_reconstruction_mask;
      bool                       reconstruction_attempted;
      bool                       reconstruction_successful;
      std::vector<VSAReconstruction::ScopeImage>     images;
      VSAReconstruction::Reconstruction              recon;
      VSAReconstruction::ArrayParameters             param;
      std::vector<VSAReconstruction::ScopeImage>     secondary_images;
      VSAReconstruction::ArrayMoments                secondary_moments;  

      std::vector<double>        theta;
      double                     mjd;
      SEphem::Angle              lmst;
      double                     mean_az;
      double                     mean_zn;
      double                     recon_az;
      double                     recon_zn;
      double                     recon_ra;
      double                     recon_dec;
      double                     recon_ra_j2000;
      double                     recon_dec_j2000;
      double                     recon_fov_x;
      double                     recon_fov_y;
      double                     recon_derotated_fov_x;
      double                     recon_derotated_fov_y;
      double                     recon_r0_x;
      double                     recon_r0_y;

      double                     N2;
      double                     msc_width;
      double                     msc_length;
      double                     msc_disp;
      double                     mlt_log10_energy;
      double                     mlt_log10_energy_chi2;

      bool                       write_event;

      // ARRAY DATA
      VSEventArrayDatum          ead;

      // CUTS
      VSArrayCutsDatum*          acd;

      // SIMULATIONS
      bool                       has_sim;
      VSArraySimulationDatum*    sim;
      VSHeaderSimulationDatum*   sim_header;
      unsigned                   sim_num;
      double                     sim_aomega;

      // LAST STATE (FOR DIAGNOSTICS)
      double                     last_event_time_sec;
      double                     last_event_time_hist;
      unsigned                   last_event_num;
      uint32_t                   last_ten_mhz_elapsed;
      uint32_t                   last_ten_mhz_veto_both;
      uint32_t                   last_ten_mhz_veto_vdaq;
      uint32_t                   last_ten_mhz_veto_l3;

      void reset()
      {
	ied_processing                = false;
	has_array_event               = false;
	processed_array_event         = false;
	has_l3                        = false;
	nscope_trigger                = 0;
	nscope_sent_l3                = 0;
	nscope_has_event              = 0;
	nscope_image                  = 0;
	nscope_quality                = 0;
	reconstruction_attempted      = false;
	reconstruction_successful     = false;
	write_event                   = false;
	for(std::vector<ScopeIED*>::iterator iscope=scope.begin();
	    iscope!=scope.end(); iscope++)
	  (*iscope)                   ->reset();
	got_one_position              = false;
	mean_position_x               .clear();
	mean_position_y               .clear();
 	has_sim                       = false;
	event_type                    = ET_UNKNOWN;
	event_failed_software_trigger = false;
	event_sent_trigger_mask       = 0;
	event_has_vdaq_mask           = 0;
	event_has_image_mask          = 0;
	event_used_in_reconstruction_mask = 0;

	scope_ied                     = 0;
      }

      void resetForScopeImage(const unsigned iscope, const unsigned nchan)
      {
	scope_num                     = iscope;	
	scope_ied                     = scope[iscope];
	scope_cal                     = &tsd->cal.scope[iscope];
	scope_diag                    = 0;
	if(tsd->partial_diagnostics)
	  scope_diag                  =
	    tsd->partial_diagnostics->partialScopeDiagnostics(iscope);
	scope_nchan                   = nchan;

	scope_nchan_hit               = 0;
	scope_nchan_logain            = 0;
	scope_nchan_trigger           = 0;
	scope_total_signal            = 0.0;

	chan_sij_max1_raw             = 0.0;
	chan_sij_max2_raw             = 0.0;
	chan_sij_max3_raw             = 0.0;
	chan_sij_max1_raw_j           = nchan;
	chan_sij_max2_raw_j           = nchan;
	chan_sij_max3_raw_j           = nchan;
      }

      void finalTransferDataToEAD(uint64_t event_time_ten_mhz_elap,
				  uint64_t event_time_ten_mhz_live);
      void transferDataToEAD();

    private:
      IED(const IED&);
      IED& operator= (const IED&);
    };

    // ------------------------------------------------------------------------
    // IED Raw Data Fetcher
    // ------------------------------------------------------------------------

    class IEDRawDataFetcher: public VSMuonAnalysis::RawDataFetcher
    {
    public:
      IEDRawDataFetcher(const IED* ied, const ScopeMergedCal* scal):
	RawDataFetcher(), m_ied(ied), m_scal(scal)
      { /* nothing to see here */ }
      virtual ~IEDRawDataFetcher();
      virtual bool fetchRawData(unsigned npixel, 
				double* raw_data, double* raw_total_data);
    private:
      const IED*            m_ied;
      const ScopeMergedCal* m_scal;
    };

    // ------------------------------------------------------------------------
    // Simulation energy binner
    // ------------------------------------------------------------------------

    class EBinner
    {
    public:
      EBinner(): m_binner(0), m_dloge(), m_ztable(), m_ntable() { };
      EBinner(double logemin, double logemax, double _dloge,
	      unsigned _ztable = 0):
	m_binner(_dloge,logemin), m_dloge(_dloge), m_ztable(_ztable), 
	m_ntable(unsigned(m_binner.valToBin(logemax, 0))+1) { };
      unsigned table(double loge) const
      { int ibin = m_binner.valToBin(loge, 0);
	return((ibin>=0)&&(ibin<int(m_ntable)))?
	  unsigned(ibin+m_ztable):0xFFFFFFFFU; }
      double energy(unsigned itable) const
      { return m_binner.binToVal(itable, 0)+0.5*m_dloge; }
      unsigned ntable() const { return m_ntable; }
      unsigned ztable() const { return m_ztable; }
    private:
      VSBinCalcLinear<double>         m_binner;
      double                          m_dloge;
      unsigned                        m_ztable;
      unsigned                        m_ntable;      
    };

    // ------------------------------------------------------------------------
    //  ___     _                     _   ___             _   _
    // |_ _|_ _| |_ ___ _ _ _ _  __ _| | | __|  _ _ _  __| |_(_)___ _ _  ___
    //  | || ' \  _/ -_) '_| ' \/ _` | | | _| || | ' \/ _|  _| / _ \ ' \(_-<
    // |___|_||_\__\___|_| |_||_\__,_|_| |_| \_,_|_||_\__|\__|_\___/_||_/__/
    // 
    // ------------------------------------------------------------------------
    
    // ------------------------------------------------------------------------
    // Diagnostics
    // ------------------------------------------------------------------------

    void doTriggerDiagnostics(//const bool* trigger, const double* signal,
			      const std::vector<bool>& trigger, const double* signal,
			      bool l2_trigger_received,
			      PartialScopeDiagnostics* sdiag,
			      unsigned& largest_region_nchan,
			      unsigned& trigger_ichan);

    // ------------------------------------------------------------------------
    // IED
    // ------------------------------------------------------------------------

    inline unsigned allocateIED(unsigned nalloc=1);
    IED* getEmptyIED(unsigned seq_packet_number, unsigned vbf_packet_number);
    void releaseIED(IED* ied);

    void processConstructedIEDs();
    void processConstructedIED(IED* ied);
    void finalizeIED(IED* ied);
    bool writeIED(const IED* ied);

    // ------------------------------------------------------------------------
    // LOCKING
    // ------------------------------------------------------------------------

    inline void lockIEDList();
    inline void unlockIEDList();
    inline void lockIEDListAndWait(unsigned seq_number);
    inline void unlockIEDListAndBroadcast();
    inline void lockPointing();
    inline void unlockPointing();
    inline bool trylockWriter();
    inline void unlockWriter();

#ifndef NOTHREADS
    void doIEDListWaitForSeqNumber(unsigned seq_number);
#endif

    // ------------------------------------------------------------------------
    //  ___       _          __  __           _
    // |   \ __ _| |_ __ _  |  \/  |___ _ __ | |__  ___ _ _ ___
    // | |) / _` |  _/ _` | | |\/| / -_) '  \| '_ \/ -_) '_(_-<
    // |___/\__,_|\__\__,_| |_|  |_\___|_|_|_|_.__/\___|_| /__/
    // 
    // ------------------------------------------------------------------------

    // HELPERS - PASSED TO CLASS
    VSCleaner*                         m_clean;
    SecondaryCleaning                  m_secondary_clean;
    VSPointing*                        m_pointing;
    VSAReconstruction*                 m_reconstruction;
    std::vector<VSMuonAnalysis*>       m_muon_analysis;
    VSScaledParameterCalc*             m_sp_calc;
    VSEnergyCalc*                      m_energy_calc;
    VSCutsCalc*                        m_cuts_calc;
    VSOctaveH5Writer*                  m_io;

    // CALIBRATION - PASSED TO CLASS
    const VSAnalysisStage1Data*        m_stage1;
    const VSAnalysisStage1Data*        m_pad_stage1;
    const VSChannelMap*                m_channel_map;

    // SETTINGS - PASSED TO CLASS
    Settings                           m_settings;
    bool                               m_has_some_scp_lookups;
    bool                               m_slow_diagnostics;

    // HELPERS - CREATED BY CLASS
    VSConstantFractionCalc*            m_qcalc;
    ArrayDiagnostics*                  m_diagnostics;
    VSEventDataWriter*                 m_edw;
    VSArrayCutsWriter*                 m_acw;
    VSArraySimulationWriter*           m_asw;
    VSOctaveH5WriterVector<unsigned>*  m_s2ew;
    VSOctaveH5WriterVector<unsigned>*  m_e2sw;
    ArrayMergedCal*                    m_cal;
    VSMuonAnalysisData*                m_muon_data;
    VSSimCoordTransform*               m_sim_transform;

    // SETTINGS - CALCULATED BY CLASS
    SEphem::Angle                      m_j2000_phi;
    SEphem::Angle                      m_j2000_theta;
    SEphem::Angle                      m_j2000_psi;
    unsigned                           m_sim_ntable;
    std::map<unsigned,EBinner>         m_sim_ebinner;
    VBFRunInfo::SimData::SimPackage    m_sim_package;

    // INTERIM EVENT DATA - CREATED BY CLASS
    unsigned                           m_ied_maximum;
    unsigned                           m_ied_created;
    std::list<IED*>                    m_ied_free;
    std::list<IED*>                    m_ied_processing;
    bool                               m_ied_writer_active;
    unsigned                           m_ied_writer_next_seq_packet_number;
    std::list<ThreadSpecificData*>     m_tsd_list;

    // THREAD PROCESSING - CREATED BY CLASS
#ifndef NOTHREADS
    bool                               m_threaded;
    pthread_mutex_t                    m_ied_list_mutex;
    pthread_cond_t                     m_ied_list_wait_cond;
    pthread_mutex_t                    m_pointing_mutex;
    pthread_key_t                      m_tsd_key;
    unsigned                           m_nthreads;
    unsigned                           m_ied_wait_priority_next;
    std::list<unsigned>                m_ied_wait_priority_list;
    VSSimpleStat2<uint64_t>            m_ied_free_stat;
    unsigned                           m_ied_wait_count;
#endif
    
    // LINEARIZED PROCESSING - CREATED BY CLASS
    VSTime                             m_last_event_time;
    double                             m_last_event_time_sec;
    double                             m_last_event_time_hist;
    uint64_t                           m_last_event_time_ten_mhz_elap;
    uint64_t                           m_last_event_time_ten_mhz_live;
    unsigned                           m_last_event_num;
    unsigned                           m_last_l3_event_num;
    uint32_t                           m_last_ten_mhz_elapsed;
    uint32_t                           m_last_ten_mhz_veto_both;
    uint32_t                           m_last_ten_mhz_veto_vdaq;
    uint32_t                           m_last_ten_mhz_veto_l3;
    std::vector<uint32_t>              m_last_counts_l2;
    bool                               m_mean_have_ra_zero;
    double                             m_mean_ra_zero;
    VSSimpleStat1<double>              m_mean_ra;
    VSSimpleStat1<double>              m_mean_dec;
    
    bool                               m_have_one_array_event;
    bool                               m_have_one_array_event_with_l3;
    VSTime                             m_first_time;
    std::vector<unsigned>              m_first_scope_has_event;

    VSHeaderSimulationDatum*           m_sim_header;
    unsigned                           m_nevent_written;
    unsigned                           m_nsim_written;
    unsigned                           m_sim_ntel;
  };

  // --------------------------------------------------------------------------
  //  ___      _ _            ___             _   _
  // |_ _|_ _ | (_)_ _  ___  | __|  _ _ _  __| |_(_)___ _ _  ___
  //  | || ' \| | | ' \/ -_) | _| || | ' \/ _|  _| / _ \ ' \(_-<
  // |___|_||_|_|_|_||_\___| |_| \_,_|_||_\__|\__|_\___/_||_/__/
  // 
  // --------------------------------------------------------------------------
  
  // ==========================================================================
  // ARRAY MERGED CAL
  // ==========================================================================

  inline void VBFAnalysisStage2::ArrayMergedCal::
  setEventNumber(unsigned event_num)
  {
    unsigned sup_islice = m_sup->getSliceByEventNum(event_num);
    if(sup_islice!=m_sup_islice)
      {
	m_sup_islice = sup_islice;
	mergeSupSlice();
      }

    if((m_cs==CS_PEDRMS)||(m_pad_ped))
      {
	unsigned ped_islice = 
	  m_ped->getSliceByEventNum(event_num, m_ped_islice);
	if(ped_islice!=m_ped_islice)
	  {
	    m_ped_islice = ped_islice;
	    mergePedSlice();
	  }
      }
  }
  
  // ==========================================================================
  // IED related
  // ==========================================================================

  inline unsigned VBFAnalysisStage2::allocateIED(unsigned nalloc)
  {
    unsigned ialloc=0;
    while((ialloc<nalloc)&&(m_ied_created<m_ied_maximum))
      {
	IED* ied = new IED(m_stage1, m_settings.theta_targets.size()+1,
			   m_cuts_calc);
	m_ied_free.push_front(ied);
	m_ied_created++;
	ialloc++;
      }
    return ialloc;
  }

  // ==========================================================================
  // LOCKING
  // ==========================================================================

  inline void VBFAnalysisStage2::lockIEDList()
  {
#ifndef NOTHREADS
    if(m_threaded)vsassert(pthread_mutex_lock(&m_ied_list_mutex) == 0);
#endif
  }
    
  inline void VBFAnalysisStage2::unlockIEDList()
  {
#ifndef NOTHREADS
    if(m_threaded)vsassert(pthread_mutex_unlock(&m_ied_list_mutex) == 0);
#endif
  }

  inline void VBFAnalysisStage2::lockIEDListAndWait(unsigned seq_number)
  {
#ifndef NOTHREADS
    if(m_threaded)doIEDListWaitForSeqNumber(seq_number);
    else if(m_ied_free.empty())allocateIED();
#else
    if(m_ied_free.empty())allocateIED();
#endif
  }

  inline void VBFAnalysisStage2::unlockIEDListAndBroadcast()
  {
#ifndef NOTHREADS
    if(m_threaded)
      {
	vsassert(pthread_mutex_unlock(&m_ied_list_mutex) == 0);  
	vsassert(pthread_cond_broadcast(&m_ied_list_wait_cond) == 0);
      }
#endif
  }

  inline void VBFAnalysisStage2::lockPointing()
  {
#ifndef NOTHREADS
    if(m_threaded)vsassert(pthread_mutex_lock(&m_pointing_mutex) == 0);
#endif
  }
    
  inline void VBFAnalysisStage2::unlockPointing()
  {
#ifndef NOTHREADS
    if(m_threaded)vsassert(pthread_mutex_unlock(&m_pointing_mutex) == 0);
#endif
  }

  // ==========================================================================
  // Diagnostics
  // ==========================================================================

  inline void VBFAnalysisStage2::PartialScopeDiagnostics::
  recordTriggeredChannel(unsigned ichan, bool l2_trigger_received)
  {
    if(l2_trigger_received)channel_is_triggered[ichan]++;
    else channel_is_triggered_no_l2[ichan]++;
  }
  
  inline void VBFAnalysisStage2::PartialScopeDiagnostics::
  recordHitChannel(bool trigger_l2,
		   unsigned ichan, bool lo_gain, bool trigger, unsigned peak, 
		   double sij, double tij, unsigned nsample, 
		   const uint32_t* samples, bool slow_diagnostics)
  {
    //    channel_is_hit[ichan]++;
    
    if(lo_gain)
      {
	if(trigger_l2)
	  channel_is_lo_gain[ichan]++;
	else
	  channel_is_lo_gain_no_l2[ichan]++;
      }
    else
      {
	if(trigger_l2)
	  channel_is_hi_gain[ichan]++;
	else
	  channel_is_hi_gain_no_l2[ichan]++;
      }
	
    if(slow_diagnostics)
      {
	uint32_t*__restrict__ istat;

	if(lo_gain)
	  {
	    if(trigger_l2)
	      istat = channel_mean_lo_trace + m_nsample*ichan;
	    else 
	      istat = channel_mean_lo_trace_no_l2 + m_nsample*ichan;

	    channel_peak_hist_lowg[peak_hist_iel(ichan,peak)]++;
	  }
	else
	  {
	    if(trigger_l2)
	      istat = channel_mean_hi_trace + m_nsample*ichan;
	    else
	      istat = channel_mean_hi_trace_no_l2 + m_nsample*ichan;

	    if(trigger)
	      channel_peak_hist_trig[peak_hist_iel(ichan,peak)]++;
	    else
	      channel_peak_hist_notr[peak_hist_iel(ichan,peak)]++;
	  }

#ifdef VECTORIZE_FRIENDLY
	for(unsigned i=0;i<nsample;i++)istat[i]+=samples[i];
#else
	const uint32_t* isamp = samples;
	const uint32_t*const nsamp = samples+nsample;
	while(isamp<nsamp)*istat++ += *isamp++;
#endif

	tij*=2.0;
	if((tij>=-15)&&(tij<=30))channel_raw_tij[ichan].accumulate(tij);

	if(sij > 0)
	  channel_raw_log10_sij[ichan].accumulate(log10(sij));
      }
  }

  inline void VBFAnalysisStage2::PartialScopeDiagnostics::
  recordImageChannel(unsigned ichan, bool lo_gain, double sij,
		     double tij, bool slow_diagnostics)
  {
    channel_in_image[ichan]++;

    if(slow_diagnostics)
      {
	if((tij>=-15)&&(tij<=30))channel_image_tij[ichan].accumulate(tij);

	if(sij > 0)
	  {
	    channel_image_log10_sij[ichan].accumulate(log10(sij));
	    if(lo_gain) channel_logain_log10_sij[ichan].accumulate(log10(sij));
	  }
      }
  }

  inline void VBFAnalysisStage2::PartialScopeDiagnostics::
  recordPedestal(unsigned ichan, unsigned nsample, const uint32_t* samples)
  {
    if(nsample>=m_nsample)
      {
	npedestal[ichan]++;
	uint32_t* covariance = channel_pedestal_covariance + ichan*m_nsample;
	for(unsigned isample = 0; isample<m_nsample; isample++)
	  {
	    const uint32_t*__restrict__ si = samples + isample;
	    uint32_t s0 = si[0];
	    unsigned nj = m_nsample-isample;
	    channel_pedestal_mean[ichan*m_nsample+isample] += s0;
	    for(unsigned jsample = 0; jsample<nj; jsample++)
	      covariance[jsample] += s0*si[jsample];
	  }
      }
  }

} // namespace VERITAS

#undef __CPC
#undef __CPA

#endif // not defined VBFANALYSISSTAGE2_HPP
