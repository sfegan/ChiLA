//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFDumper.hpp

  Dump VBF data to a stream

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/20/2005

  $Id: VBFDumper.hpp,v 3.8 2009/11/05 23:02:20 matthew Exp $

*/

#ifndef VBFDUMPER_HPP
#define VBFDUMPER_HPP

#include <iostream>

#include <VSSimpleVBF.hpp>

namespace VERITAS
{

  class VBFDumper: public VSSimpleVBFVisitor
  {
  public:
    struct ChanData
    {
      ChanData(): 
	channel_num(), hit(), trigger(), charge(),
	pedestal(), lo_gain(), nsample(), samples(), integrated()
      { }

      uint32_t                  channel_num;
      bool                      hit; 
      bool                      trigger;
      uint32_t                  charge; 
      uint32_t                  pedestal;
      bool                      lo_gain;
      unsigned                  nsample;
      uint32_t*                 samples;
      uint32_t*                 integrated;
    };

    VBFDumper(std::ostream& stream, 
	      bool lo_gain = false,
	      bool one_scope = false,
	      unsigned scope_id = 0,
	      bool no_channels = false, 
	      bool dump_raw_clock = false):
      VSSimpleVBFVisitor(), m_stream(stream), m_lo_gain(lo_gain),
      m_scope_id(scope_id), m_one_scope(one_scope),
      m_dump_raw_clock(dump_raw_clock), m_no_channels(no_channels),
      m_overflow(false) { }
    virtual ~VBFDumper();

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

    virtual void visitChannel(bool& veto_channel, void* user_data,
			      uint32_t                  channel_num, 
			      bool                      hit, 
			      bool                      trigger);
  
    virtual void leaveChannel(bool& veto_channel, void* user_data);

    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);

    virtual void visitOverflowEvent(void* user_data,
				    unsigned            datum_num,
				    const VDatum*       datum);

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

    static void dumpRawClock(std::ostream& stream, const uint16_t* clock);

  private:
    std::ostream& m_stream;
    bool m_lo_gain;
    unsigned m_scope_id;
    bool m_one_scope;
    bool m_dump_raw_clock;
    bool m_no_channels;
    bool m_overflow;

    ChanData m_chan_data;
  };

  std::ostream& operator<< (std::ostream& stream, const VEventType& et);

} // namespace VERITAS

#endif // not defined VBFDUMPER_HPP
