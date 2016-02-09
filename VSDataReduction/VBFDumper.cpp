//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFDumper.hpp

  Dump VBF data to a stream

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/20/2005

  $Id: VBFDumper.cpp,v 3.9 2009/11/05 23:02:20 matthew Exp $

*/

#include <VBFDumper.hpp>
#include "fast_alloc.hpp"

using namespace VERITAS;

VBFDumper::~VBFDumper()
{
  // nothing to see here
}

#ifndef NOTHREADS
void VBFDumper::usingThreads(unsigned nthreads)
{
  // nothing to see here
}
#endif

void VBFDumper::
visitPacket(bool& veto_packet, void*& user_data,
	    uint32_t                   seq_packet_number,
	    uint32_t                   vbf_packet_number,
	    bool                       has_array_event,
	    bool                       has_sim_header,
	    bool                       has_sim_event,
	    uint32_t                   num_overflow_datum,
	    const VBFPacket*           packet)
{
  m_overflow = false;
  m_stream 
    << "PACKET:   " << seq_packet_number << ' ' << vbf_packet_number 
    << ' ' << has_array_event << ' ' << has_sim_header << ' ' << has_sim_event 
    << ' ' << num_overflow_datum
    << '\n';
}

void VBFDumper::
visitArrayEvent(bool& veto_array_event, void* user_data,
		uint32_t               num_scope_events,
		bool                   has_array_trigger,
		bool                   has_event_number,
		uint32_t               event_num,
		bool                   has_good_event_time,
		const VSTime&          best_event_time,
		EventType              event_type,
		uint32_t               l2_trigger_mask,
		const VArrayEvent*     array_event)
{
  m_stream 
    << "EVENT:    " << has_array_trigger << ' ' << has_event_number << ' '
    << num_scope_events << ' ' << event_num << ' '
    << has_good_event_time << ' ' << best_event_time << ' ';

  switch(event_type)
    {
    case ET_UNKNOWN:
      m_stream << "UNKNOWN ";
      break;
    case ET_PED:
      m_stream << "PED ";
      break;
    case ET_L2:
      m_stream << "L2 ";
      break;
    }

  for(unsigned iscope=4;iscope>0;iscope--)
    m_stream << ((l2_trigger_mask&(1<<(iscope-1)))?1:0);
    
  m_stream << '\n';
}

void VBFDumper::
visitArrayTrigger(bool& veto_array_event, void* user_data,
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
		  const VArrayTrigger* trigger)
{
  if(m_overflow)m_stream << "OVERFLOW ";
  m_stream 
    << "TRIGGER:  " << event_num << ' ' << event_type << ' '
    << trigger_mask << ' ' << flags << ' ';
  if(m_dump_raw_clock)
    {
      dumpRawClock(m_stream,trigger->getGPSTime());
      m_stream << ' ';
    }
  m_stream << raw_time << ' ';
  if(m_dump_raw_clock)m_stream << raw_time.getDayNS() << ' ';
  m_stream
    << (raw_time.isFlyWheeling()?1:0) << ' '
    << (raw_time.isLargeTimeOffset()?1:0) << ' '
    << (raw_time.isLargeFreqOffset()?1:0) << ' '
    << (raw_time.isBatteryFailed()?1:0) << ' '
    << at_flags << ' ' << config_mask << ' ' 
    << num_telescopes << ' ' << num_trigger_telescopes << ' ' 
    << run_number << ' '
    << ten_mhz_clocks[0] << ' ' << ten_mhz_clocks[1] << ' '
    << ten_mhz_clocks[2] << ' ' << ten_mhz_clocks[3] << ' '
    << cal_count[0] << ' ' << cal_count[1] << ' ' << cal_count[2] << ' '
    << ped_count[0] << ' ' << ped_count[1] << ' ' << ped_count[2] << ' '
    << '\n';
}

void VBFDumper::
visitArrayTriggerScope(bool& veto_array_event, void* user_data,
		       uint32_t        telescope_num,
		       bool            triggered,
		       uint32_t        event_type,
		       float           altitude,
		       float           azimuth,
		       uint32_t        tdc,
		       uint32_t        shower_delay,
		       uint32_t        comp_delay,
		       const uint32_t* l2_counts,
		       const uint32_t* cal_counts)
{
  if(m_overflow)m_stream << "OVERFLOW ";
  m_stream 
    << "TRSCOPE:  " << telescope_num << ' ' << triggered << ' ' 
    << event_type << ' ' << altitude << ' ' << azimuth << ' ' << tdc << ' '
    << shower_delay << ' ' << comp_delay << ' ' 
    << l2_counts[0] << ' ' << l2_counts[1] << ' ' << l2_counts[2] << ' '
    << cal_counts[0] << ' ' << cal_counts[1] << ' ' << cal_counts[2]
    << '\n';
}

void VBFDumper::dumpRawClock(std::ostream& stream, const uint16_t* clock)
{
  static const char HEX[] = { '0','1','2','3','4','5','6','7',
			      '8','9','a','b','c','d','e','f' };
  stream 
    << HEX[clock[0]>>12&15] 
    << HEX[clock[0]>>8&15]
    << HEX[clock[0]>>4&15] << ':'
    << HEX[clock[0]>>0&15]
    << HEX[clock[1]>>12&15]
    << HEX[clock[1]>>8&15] << ':'
    << HEX[clock[1]>>4&15]
    << HEX[clock[1]>>0&15] << ':'
    << HEX[clock[2]>>12&15]
    << HEX[clock[2]>>8&15] << ':'
    << HEX[clock[2]>>4&15]
    << HEX[clock[2]>>0&15] << '.'
    << HEX[clock[3]>>12&15]
    << HEX[clock[3]>>8&15]
    << HEX[clock[3]>>4&15]
    << HEX[clock[3]>>0&15]
    << HEX[clock[4]>>12&15]
    << HEX[clock[4]>>8&15]
    << HEX[clock[4]>>4&15] << ':'
    << HEX[clock[4]>>0&15];
};

void VBFDumper::
visitScopeEvent(bool& veto_scope_event, void* user_data,
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
		const VEvent*          event)
{
  if(m_one_scope && telescope_num != m_scope_id) 
    {
      veto_scope_event=true;
      return;
    }

  if(m_overflow)m_stream << "OVERFLOW ";
  m_stream 
    << "SCOPE:    " << event_num << ' ' 
    << (uint16_t)telescope_num << ' ' << event_type << ' ' 
    << (uint16_t)trigger_mask << ' ' << flags << ' ';
  if(m_dump_raw_clock)
    {
      dumpRawClock(m_stream,event->getGPSTime());
      m_stream << ' ';
    }
  m_stream << raw_time << ' ';
  if(m_dump_raw_clock)m_stream << raw_time.getDayNS() << ' ';
  m_stream
    << (raw_time.isFlyWheeling()?1:0) << ' '
    << (raw_time.isLargeTimeOffset()?1:0) << ' '
    << (raw_time.isLargeFreqOffset()?1:0) << ' '
    << (raw_time.isBatteryFailed()?1:0) << ' '
    << num_samples << ' ' << num_channels_saved << ' ' 
    << num_channels_total << ' ' << num_clock_trigger << '\n';
  if(m_no_channels)veto_scope_event=true;
}

void VBFDumper::
visitChannel(bool& veto_channel, void* user_data,
	     uint32_t                  channel_num, 
	     bool                      hit, 
	     bool                      trigger)
{
  m_chan_data.channel_num = channel_num;
  m_chan_data.hit = hit;
  m_chan_data.trigger = trigger;
}

void VBFDumper::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		unsigned               nsample,
		const uint32_t*        samples,
		const uint32_t*        integrated)
{
  m_chan_data.charge = charge;
  m_chan_data.pedestal = pedestal;
  m_chan_data.lo_gain = lo_gain;
  
  if(m_chan_data.nsample != nsample)
    {
      if(m_chan_data.nsample)
	{
	  delete[] m_chan_data.samples;
	  delete[] m_chan_data.integrated;
	}

      std::cout << "ALLOCATING" << std::endl;

      m_chan_data.samples = new uint32_t[nsample];
      m_chan_data.integrated = new uint32_t[nsample];
    }

  m_chan_data.nsample = nsample;

  for(unsigned isample=0;isample<nsample;isample++)
    {
      m_chan_data.samples[isample] = samples[isample];
      m_chan_data.integrated[isample] = integrated[isample];
    }

  const ChanData& cd = m_chan_data;
  if(m_lo_gain && (!cd.hit || !cd.lo_gain)) return;

  if(m_overflow)m_stream << "OVERFLOW ";
  m_stream 
    << "CHAN:     " 
    << std::setw(4) << cd.channel_num << ' ' << cd.hit << ' ' << cd.trigger
    << '\n';

  if(cd.hit)
    {
      m_stream 
	<< "HITCHAN:  " 
	<< std::setw(4) << cd.channel_num 
	<< ' ' << cd.charge << ' ' << cd.pedestal << ' ' << cd.lo_gain;
      
      for(unsigned isample=0;isample<cd.nsample;isample++)
	m_stream << ' ' << cd.samples[isample];
      m_stream << ' ' << cd.integrated[cd.nsample-1] << '\n';
    }
}

void VBFDumper::leaveChannel(bool& veto_channel, void* user_data)
{
  std::cout << "VBFDumper::leaveChannel() " << std::endl;

  const ChanData& cd = m_chan_data;

  if(m_lo_gain && (!cd.hit || !cd.lo_gain)) return;

  if(m_overflow)m_stream << "OVERFLOW ";
  m_stream 
    << "CHAN:     " 
    << std::setw(4) << cd.channel_num << ' ' << cd.hit << ' ' << cd.trigger
    << '\n';

  if(cd.hit)
    {
      m_stream 
	<< "HITCHAN:  " 
	<< std::setw(4) << cd.channel_num 
	<< ' ' << cd.charge << ' ' << cd.pedestal << ' ' << cd.lo_gain;
      
      for(unsigned isample=0;isample<cd.nsample;isample++)
	m_stream << ' ' << cd.samples[isample];
      m_stream << ' ' << cd.integrated[cd.nsample-1] << '\n';
    }
}

void VBFDumper::
visitOverflowEvent(void* user_data,
		   unsigned            datum_num,
		   const VDatum*       datum)
{
  if(!datum)return;
  m_overflow = true;
  const VArrayTrigger* at = dynamic_cast<const VArrayTrigger*>(datum);
  if(at)VSSimpleVBFDispatcher::dispatchArrayTrigger(this, user_data, at);

  const VEvent* ev = dynamic_cast<const VEvent*>(datum);
  if(ev)VSSimpleVBFDispatcher::dispatchScopeEvent(this,user_data,datum_num,
						  ev,!m_no_channels);
  m_overflow = false;
}


void VBFDumper::
visitSimulationHeader(void* user_data,
		      uint32_t         run_number,
		      const VSTime&    date_of_sims,
		      uint32_t         simulation_package,
		      const std::string& simulator_name,
		      const VSTime&    date_of_array,
		      uint32_t         corsika_atm_model,
		      float            obs_altitude_m,
		      const std::vector<VArrayConfiguration>& array,
		      const std::string& stringified_config)
{
  m_stream 
    << "SIMHEAD:  " << run_number << ' ' << date_of_sims << ' '
    << simulation_package << ' ' << simulator_name << ' ' 
    << date_of_array << ' ' << corsika_atm_model << ' ' << array.size() << ' ';
  if(stringified_config.size() > 20)
    m_stream << stringified_config.substr(0,20) << "...";
  else
    m_stream << stringified_config;
  m_stream << '\n';
}
  
void VBFDumper::
visitSimulationEvent(bool& veto_packet, void* user_data,
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
		     float             core_elevation_asl_m)
{
  m_stream 
    << "SIMEVENT: " << run_number << ' ' << event_num << ' '
    << corsika_particle_id << ' ' << energy_gev << ' ' 
    << obs_zenith_deg << ' ' << obs_azimuth_deg << ' '
    << primary_zenith_deg << ' ' << primary_azimuth_deg << ' '
    << core_east_m << ' ' << core_south_m << ' ' << core_elevation_asl_m << ' '
    << ref_zenith_deg << ' ' << ref_azimuth_deg << ' ' 
    << ref_position_angle_deg
    << '\n';
}

void VBFDumper::
visitChiLAHeader(void* user_data,
		 const std::string&    database_name,
		 const std::vector<VChiLASimParamTable>&
                                       sim_param_tables,
		 const std::vector<VChiLAOpticsConfiguration>&
		                       optics_configurations,
		 const std::vector<VChiLAElectronicsConfiguration>& 
		                       electronics_configurations)
{
  m_stream 
    << "CHILAHEAD: " << database_name << ' ' << sim_param_tables.size() << ' '
    << optics_configurations.size() << ' '
    << electronics_configurations.size() << '\n';
  for(unsigned itable=0;itable<sim_param_tables.size();itable++)
    {
      const VChiLASimParamTable& spt(sim_param_tables[itable]);
      m_stream
	<< "CHILATABLE: " << spt.fTableID << ' ' 
	<< spt.fPrimaryID << ' ' << spt.fEnergyTeV << ' ' 
	<< spt.fZenithMinRad << ' ' << spt.fZenithMaxRad << ' '
	<< spt.fAzimuthMinRad << ' ' << spt.fAzimuthMaxRad << ' '
	<< spt.fOpticsID << ' ' << spt.fSamplingRadiusM << ' ' 
	<< spt.fEventCount << ' ' << spt.fTableName << '\n';
    }
}

void VBFDumper::
visitChiLAEvent(bool& veto_packet, void* user_data,
		uint32_t               table_index,
		uint32_t               electronics_id,
		uint32_t               table_event_index,
		float                  a_omega,
		float                  a_omega_var)
{
  m_stream 
    << "CHILAEVENT: " << table_index << ' ' << electronics_id << ' '
    << table_event_index << ' ' << a_omega << ' ' << a_omega_var << '\n';
}

void VBFDumper::
visitKascadeHeader(void* user_data,
		   uint32_t            corsika_particle_id,
		   float               energy_gev,
		   uint32_t            shower_id,
		   float               x_area_width_m,
		   float               y_area_width_m,
		   uint32_t            north_south_grid)
{
  m_stream 
    << "KASHEAD: " << corsika_particle_id << ' ' << energy_gev << ' '
    << shower_id << ' ' << x_area_width_m << ' ' << y_area_width_m << ' '
    << north_south_grid << '\n';
}
 
void VBFDumper::
visitKascadeEvent(bool& veto_packet, void* user_data,
		  int32_t              nx,
		  int32_t              ny,
		  int32_t              direction_index,
		  float                emission_altitude_m,
		  float                emission_altitude_sigma,
		  float                muon_ratio,
		  float                a_omega,
		  float                rel_tel_trigger_time_ns,
		  float                differential_rate_per_event_hz,
		  float                integral_rate_per_event_hz)
{
  m_stream 
    << "KASEVENT: " << nx << ' ' << ny << ' ' << direction_index << ' '
    << emission_altitude_m << ' ' << emission_altitude_sigma << ' '
    << muon_ratio << ' ' << a_omega << ' ' << rel_tel_trigger_time_ns << ' '
    << differential_rate_per_event_hz << ' ' << integral_rate_per_event_hz
    << '\n';
}

std::ostream& VERITAS::operator<< (std::ostream& stream, const VEventType& et)
{
  switch(et.trigger)
    {
    case VEventType::L2_TRIGGER:
      stream << "L2";
      break;
    case VEventType::HIGH_MULT_TRIGGER:
      stream << "HIGH_MULT";
      break;
    case VEventType::NEW_PHYS_TRIGGER:
      stream << "NEW_PHYS";
      break;
    case VEventType::CAL_TRIGGER:
      stream << "CAL/";
      switch(et.calibration)
	{
	case VEventType::NOT_CALIBRATION:
	  stream << "NOT";
	  break;
	case VEventType::OPTICAL_CALIBRATION:
	  stream << "OPTICAL";
	  break;
	case VEventType::CHARGE_CALIBRATION:
	  stream << "CHARGE";
	  break;
	default:
	  stream << "UNKNOWN";
	  break;
	}
      break;
    case VEventType::PED_TRIGGER:
      stream << "PED";
      break;
    default:
      stream << "UNKNOWN";
      break;      
    };
  
  if(et.xtra_samples)stream << ",XTRA";
  if(et.force_full_mode)stream << ",FULL";
  //if(et.has_error)stream << ",ERROR";
  
  return stream;
}

