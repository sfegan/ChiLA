#ifndef VSREVENTVISITOR_HPP
#define VSREVENTVISITOR_HPP

#include <VSSimpleVBF.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSCleaner.hpp>
#include <VSTimingCalc.hpp>
#include <VSChannelMap.hpp>

#include <VSRHistogramCamera.hpp>

namespace VERITAS
{
  class VSRVBFVisitor: public VSSimpleVBFVisitor
  {
  public:

    enum CleaningScale { CS_UNITY, CS_GAIN, CS_PEDRMS };

    struct Options
    {
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
      std::vector<std::string> cleaning_args;
      std::string              cleaning_scale;
      std::string              camera;
      bool                     include_trigger_channels;
      bool                     include_logain_channels;

      std::string              output_file;
    };


    VSRVBFVisitor(const VSAnalysisStage1Data& stage1);
    virtual ~VSRVBFVisitor();

    virtual void visitPacket(bool& veto_packet, 
			     void*& user_data,
			     uint32_t  seq_packet_number,
			     uint32_t  vbf_packet_number,
			     bool      has_array_event,
			     bool      has_sim_header,
			     bool      has_sim_event,
			     uint32_t  num_overflow_datum,
			     const VBFPacket* packet);

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

    virtual void leaveScopeEvent(bool veto_scope_event, 
				 void* user_data);

    virtual void visitChannel(bool& veto_channel, 
			      void* user_data,
			      uint32_t channel_num, 
			      bool hit, 
			      bool trigger);

    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);

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

    void draw();

  private:
    std::ostream& m_stream;

    const int (*m_camera_neighbors)[NUM_NEIGHBORS];

    const VSAnalysisStage1Data*        m_stage1;
    const VSTimeDepSupData*            m_sup;
    const VSTimeDepPedData*            m_ped;
    VSConstantFractionCalc*            m_qcalc;
    VSCleaner*                         m_clean;
    VSChannelMap*                      m_channel_map;
    unsigned                           m_event_num;
    unsigned                           m_iscope;
    unsigned                           m_scope_nchan_total;
    unsigned                           m_scope_nchan_hit;
    std::vector<unsigned>              m_nchan;

    std::vector<bool>                  m_chan_has_signal;
    std::vector<bool>                  m_chan_hit;
    std::vector<bool>                  m_chan_trigger;
    std::vector<bool>                  m_chan_logain;
    std::vector< std::vector<bool> >   m_chan_l2;
    double*                            m_chan_sij;
    double*                            m_chan_sij_total;
    double*                            m_chan_tij;
    std::vector<double>                m_chan_clean_multiplier;

    std::vector< VSRHistogramCamera* > m_scope_hist;

    static Options s_default_options;
    Options m_options;
  };

} // namespace VERITAS

#endif // VSREVENTVISITOR_HPP
