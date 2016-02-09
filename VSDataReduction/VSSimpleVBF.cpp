//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleVBF.cpp
  Simple VBF Dispatcher and Vistor

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/07/2005

  $Id: VSSimpleVBF.cpp,v 3.12 2008/03/05 22:38:10 sfegan Exp $

*/

#include <cstdlib>
#include <time.h>

#include <fast_alloc.hpp>
#include <VSSimpleVBF.hpp>

using namespace VERITAS;

// ============================================================================
// 
// VBFPacket
// 
// ============================================================================

VBFPacket* VBFPacket::deepCopy(bool omit_vbf_packet) const
{
  VPacket* v = 0;
  if(!omit_vbf_packet)
    {
      const VArrayEvent* ae = vbf->getArrayEvent();
#warning DO WE NEED TO SET RUN HERE
      VArrayEvent* nae = new VArrayEvent();
      if(ae->hasTrigger())nae->setTrigger(ae->getTrigger()->copyAT());
      for(unsigned ievent=0;ievent<ae->getNumEvents();ievent++)
	nae->addEvent(ae->getEventAt(ievent)->copyEvent());
      v = new VPacket;
      v->putArrayEvent(nae);
    }

  VBFPacket* p = new VBFPacket;
  p->vbf_num             = vbf_num;
  p->vbf                 = v;
  p->has_best_event_time = has_best_event_time;
  p->best_event_time     = best_event_time;
  return p;
}

// ============================================================================
// 
// VSSimpleVBFVisitor
// 
// ============================================================================

VSSimpleVBFVisitor::~VSSimpleVBFVisitor()
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
registerDispatcher(VSSimpleVBFDispatcherStop* dispatch_stop)
{
  fDispatchStop=dispatch_stop;
}

void VSSimpleVBFVisitor::
visitFile(const char* filename, unsigned npacket)
{
  // nothing to see here
}

#ifndef NOTHREADS
void VSSimpleVBFVisitor::usingThreads(unsigned nthreads)
{
  vsassert(nthreads==1);
}
#endif

void VSSimpleVBFVisitor::
leaveFile()
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitPacket(bool& veto_packet, void*& user_data,
	    uint32_t                   seq_packet_number,
	    uint32_t                   vbf_packet_number,
	    bool                       has_array_event,
	    bool                       has_sim_header,
	    bool                       has_sim_event,
	    uint32_t                   num_overflow_datum,
	    const VBFPacket*           packet)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
leavePacket(bool veto_packet, void* user_data)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
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
  // nothing to see here
}

void VSSimpleVBFVisitor::
leaveArrayEvent(bool veto_array_event, void* user_data)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
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
  // nothing to see here
}


void VSSimpleVBFVisitor::
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
  // nothing to see here
}

void VSSimpleVBFVisitor::
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
  // nothing to see here
}

void VSSimpleVBFVisitor::
leaveScopeEvent(bool veto_scope_event, void* user_data)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitChannel(bool& veto_channel, void* user_data,
	     uint32_t                  channel_num, 
	     bool                      hit, 
	     bool                      trigger)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
leaveChannel(bool veto_channel, void* user_data)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitOverflowEvent(void* user_data,
		   unsigned            datum_num,
		   const VDatum*       datum)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		unsigned               nsample,
		const uint32_t*        samples,
		const uint32_t*        integrated)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitSimulationHeader(void* user_data,
		      uint32_t         run_number,
		      const VSTime&    date_of_sims,
		      uint32_t         simulation_package,
		      const std::string& simulator_name,
		      const VSTime&    date_of_array,
		      uint32_t         corsika_atm_model,
		      float            obs_altitude_m,
		      const std::vector<VArrayConfiguration>&
		                       array,
		      const std::string& stringified_config)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
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
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitChiLAHeader(void* user_data,
		 const std::string&    database_name,
		 const std::vector<VChiLASimParamTable>&
		                       sim_param_tables,
		 const std::vector<VChiLAOpticsConfiguration>&
		                       optics_configurations,
		 const std::vector<VChiLAElectronicsConfiguration>& 
		                       electronics_configurations)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitChiLAEvent(bool& veto_packet, void* user_data,
		uint32_t               table_index,
		uint32_t               electronics_id,
		uint32_t               table_event_index,
		float                  a_omega,
		float                  a_omega_var)
{
  // nothing to see here
}

void VSSimpleVBFVisitor::
visitKascadeHeader(void* user_data,
		   uint32_t            corsika_particle_id,
		   float               energy_gev,
		   uint32_t            shower_id,
		   float               x_area_width_m,
		   float               y_area_width_m,
		   uint32_t            north_south_grid)
{
  // nothing to see here
}
 
void VSSimpleVBFVisitor::
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
  // nothing to see here
}

// ============================================================================
//
// VBFBroadcastVisitor
//
// ============================================================================

VBFBroadcastVisitor::VBFBroadcastVisitor(): 
  VSSimpleVBFVisitor(), 
#ifndef NOTHREADS
  m_visitors_mutex(),
#endif
  m_visitors()
{
#ifndef NOTHREADS
  vsassert(pthread_mutex_init(&m_visitors_mutex, NULL) == 0);
#endif
}

VBFBroadcastVisitor::~VBFBroadcastVisitor()
{
  // nothing to see here
}

void VBFBroadcastVisitor::addVisitor(VSSimpleVBFVisitor* visitor)
{
  lock();
  m_visitors.push_back(Visitor(visitor,this));
  unlock();
}

void VBFBroadcastVisitor::delVisitor(VSSimpleVBFVisitor* visitor)
{
  lock();
  for(std::list<Visitor>::iterator ivisitor=m_visitors.begin(); 
      ivisitor!=m_visitors.end(); ivisitor++)
    if(ivisitor->visitor == visitor)
      {
	m_visitors.erase(ivisitor);
	break;
      }
  unlock();
}

void VBFBroadcastVisitor::stopDispatcherForVisitor(Visitor* visitor)
{
  lock();
  visitor->stop_processing=true;
  unsigned nstopped = 0;
  for(std::list<Visitor>::iterator ivisitor=m_visitors.begin(); 
      ivisitor!=m_visitors.end(); ivisitor++)
    if(ivisitor->stop_processing)nstopped++;
  if(nstopped == m_visitors.size())fDispatchStop->stopProcessingFile();
  unlock();
}

void VBFBroadcastVisitor::
registerDispatcher(VSSimpleVBFDispatcherStop* dispatch_stop)
{
  VSSimpleVBFVisitor::registerDispatcher(dispatch_stop);
  for(std::list<Visitor>::iterator ivisitor=m_visitors.begin(); 
      ivisitor!=m_visitors.end(); ivisitor++)
    ivisitor->visitor->registerDispatcher(&*ivisitor);
}

void VBFBroadcastVisitor::visitFile(const char* filename, unsigned npacket)
{
  lock();
  for(std::list<Visitor>::iterator ivisitor=m_visitors.begin(); 
      ivisitor!=m_visitors.end(); ivisitor++)
    ivisitor->visitor->visitFile(filename, npacket);
  unlock();
}

#ifndef NOTHREADS
void VBFBroadcastVisitor::usingThreads(unsigned nthreads)
{
  for(std::list<Visitor>::iterator ivisitor=m_visitors.begin(); 
      ivisitor!=m_visitors.end(); ivisitor++)
    ivisitor->visitor->usingThreads(nthreads);
}
#endif

void VBFBroadcastVisitor::leaveFile()
{
  for(std::list<Visitor>::iterator ivisitor=m_visitors.begin(); 
      ivisitor!=m_visitors.end(); ivisitor++)
    ivisitor->visitor->leaveFile();
}

void VBFBroadcastVisitor::
visitPacket(bool& veto_packet, void*& user_data,
	    uint32_t                   seq_packet_number,
	    uint32_t                   vbf_packet_number,
	    bool                       has_array_event,
	    bool                       has_sim_header,
	    bool                       has_sim_event,
	    uint32_t                   num_overflow_datum,
	    const VBFPacket*           packet)
{
  EventData* ed = new EventData;
  user_data = static_cast<void*>(ed);

  lock();
  for(std::list<Visitor>::iterator ivisitor=m_visitors.begin(); 
      ivisitor!=m_visitors.end(); ivisitor++)
    ed->push_back(VisitorEventData(*ivisitor));
  unlock();

  unsigned nveto = 0;
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      ivisitor->veto_packet = ivisitor->veto_master;
      if(!ivisitor->veto_master)
	ivisitor->visitor->visitPacket(ivisitor->veto_packet,
				       ivisitor->user_data,
				       seq_packet_number,
				       vbf_packet_number,
				       has_array_event,
				       has_sim_header,
				       has_sim_event,
				       num_overflow_datum,
				       packet);
      if(ivisitor->veto_packet)nveto++;
    }

  if(nveto==m_visitors.size())veto_packet=true;
}

void VBFBroadcastVisitor::leavePacket(bool veto_packet, void* user_data)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    if(!ivisitor->veto_master)
      ivisitor->visitor->leavePacket(ivisitor->veto_packet,
				     ivisitor->user_data);
  delete ed;
}

void VBFBroadcastVisitor::
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
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      ivisitor->veto_array_event_primary = ivisitor->veto_packet;
      if(!ivisitor->veto_packet)
	ivisitor->visitor->visitArrayEvent(ivisitor->veto_array_event_primary,
					   ivisitor->user_data,
					   num_scope_events,
					   has_array_trigger,
					   has_event_number,
					   event_num,
					   has_good_event_time,
					   best_event_time,
					   event_type,
					   l2_trigger_mask,
					   array_event);
      if(ivisitor->veto_array_event_primary)nveto++;
    }

  if(nveto==m_visitors.size())veto_array_event=true;
}

void VBFBroadcastVisitor::
leaveArrayEvent(bool veto_array_event, void* user_data)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    if(!ivisitor->veto_packet)
      ivisitor->visitor->leaveArrayEvent(ivisitor->veto_array_event_primary||
					 ivisitor->veto_array_event_secondary,
					 ivisitor->user_data);
}

void VBFBroadcastVisitor::
visitArrayTrigger(bool& veto_array_event, void* user_data,
		  uint32_t             event_num,
		  const VEventType&    event_type,
		  uint32_t             trigger_mask,
		  uint32_t             flags,
		  const VSTime&        raw_time,
		  uint32_t             at_flags,
		  uint32_t             config_mask,
		  uint32_t             num_array_telescopes,
		  uint32_t             num_trigger_telescopes,
		  uint32_t             run_number,
		  const uint32_t*      ten_mhz_clocks,
		  const uint32_t*      cal_count,
		  const uint32_t*      ped_count,
		  const VArrayTrigger* trigger)
{
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      ivisitor->veto_array_event_secondary = 
	ivisitor->veto_array_event_primary;
      if(!ivisitor->veto_array_event_primary)
	ivisitor->visitor->
	  visitArrayTrigger(ivisitor->veto_array_event_secondary,
			    ivisitor->user_data,
			    event_num,
			    event_type,
			    trigger_mask,
			    flags,
			    raw_time,
			    at_flags,
			    config_mask,
			    num_array_telescopes,
			    num_trigger_telescopes,
			    run_number,
			    ten_mhz_clocks,
			    cal_count,
			    ped_count,
			    trigger);
      if(ivisitor->veto_array_event_secondary)nveto++;
    }

  if(nveto==m_visitors.size())veto_array_event=true;
}
    
void VBFBroadcastVisitor::
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
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_array_event_primary)
	ivisitor->visitor->
	  visitArrayTriggerScope(ivisitor->veto_array_event_secondary,
				 ivisitor->user_data,
				 telescope_num,
				 triggered,
				 event_type,
				 altitude,
				 azimuth,
				 tdc,
				 shower_delay,
				 comp_delay,
				 l2_counts,
				 cal_counts);
      if(ivisitor->veto_array_event_secondary)nveto++;
    }

  if(nveto==m_visitors.size())veto_array_event=true;
}

void VBFBroadcastVisitor::
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
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      ivisitor->veto_scope_event = 
	ivisitor->veto_array_event_primary
	||ivisitor->veto_array_event_secondary;
      if((!ivisitor->veto_array_event_primary)
	 &&(!ivisitor->veto_array_event_secondary))
	ivisitor->visitor->visitScopeEvent(ivisitor->veto_scope_event,
					   ivisitor->user_data,
					   event_num,
					   telescope_num, 
					   event_type,
					   trigger_mask,
					   flags,
					   raw_time,
					   num_samples,
					   num_channels_saved,
					   num_channels_total,
					   num_clock_trigger,
					   event);
      if(ivisitor->veto_scope_event)nveto++;
    }

  if(nveto==m_visitors.size())veto_scope_event=true;
}

void VBFBroadcastVisitor::
leaveScopeEvent(bool veto_scope_event, void* user_data)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    if((!ivisitor->veto_array_event_primary)
       &&(!ivisitor->veto_array_event_secondary))
      ivisitor->visitor->leaveScopeEvent(ivisitor->veto_scope_event,
					 ivisitor->user_data);
}

void VBFBroadcastVisitor::
visitChannel(bool& veto_channel, void* user_data,
	     uint32_t                  channel_num, 
	     bool                      hit, 
	     bool                      trigger)
{
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      ivisitor->veto_channel = ivisitor->veto_scope_event;
      if(!ivisitor->veto_scope_event)
	ivisitor->visitor->visitChannel(ivisitor->veto_channel,
					ivisitor->user_data,
					channel_num, 
					hit, 
					trigger);
      if(ivisitor->veto_channel)nveto++;
    }

  if(nveto==m_visitors.size())veto_channel=true;
}

void VBFBroadcastVisitor::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		unsigned               nsample,
		const uint32_t*        samples,
		const uint32_t*        integrated)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_channel)
	ivisitor->visitor->visitHitChannel(ivisitor->user_data,
					   channel_num,
					   charge, 
					   pedestal,
					   lo_gain,
					   nsample,
					   samples,
					   integrated);
    }
}

void VBFBroadcastVisitor::
leaveChannel(bool veto_channel, void* user_data)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    if(!ivisitor->veto_scope_event)
      ivisitor->visitor->leaveChannel(veto_channel,
				      ivisitor->user_data);
}

void VBFBroadcastVisitor::
visitOverflowEvent(void* user_data,
		   unsigned            datum_num,
		   const VDatum*       datum)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    if(!ivisitor->veto_scope_event)
      ivisitor->visitor->visitOverflowEvent(ivisitor->user_data,
					    datum_num, datum);
}

void VBFBroadcastVisitor::
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
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_packet)
	ivisitor->visitor->visitSimulationHeader(ivisitor->user_data,
						 run_number,
						 date_of_sims,
						 simulation_package,
						 simulator_name,
						 date_of_array,
						 corsika_atm_model,
						 obs_altitude_m,
						 array,
						 stringified_config);
    }
}

void VBFBroadcastVisitor::
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
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_packet)
	ivisitor->visitor->visitSimulationEvent(ivisitor->veto_packet,
						ivisitor->user_data,
						run_number,
						event_num,
						corsika_particle_id,
						energy_gev,
						obs_zenith_deg,
						obs_azimuth_deg,
						primary_zenith_deg,
						primary_azimuth_deg,
						ref_zenith_deg,
						ref_azimuth_deg,
						ref_position_angle_deg,
						core_east_m,
						core_south_m,
						core_elevation_asl_m);
      if(ivisitor->veto_packet)nveto++;
    }

  if(nveto==m_visitors.size())veto_packet=true;
}

void VBFBroadcastVisitor::
visitChiLAHeader(void* user_data,
		 const std::string&    database_name,
		 const std::vector<VChiLASimParamTable>&
		                       sim_param_tables,
		 const std::vector<VChiLAOpticsConfiguration>&
		                       optics_configurations,
		 const std::vector<VChiLAElectronicsConfiguration>& 
		                       electronics_configurations)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_packet)
	ivisitor->visitor->visitChiLAHeader(ivisitor->user_data,
					    database_name,
					    sim_param_tables,
					    optics_configurations,
					    electronics_configurations);
    }
}

void VBFBroadcastVisitor::
visitChiLAEvent(bool& veto_packet, void* user_data,
		uint32_t               table_index,
		uint32_t               electronics_id,
		uint32_t               table_event_index,
		float                  a_omega,
		float                  a_omega_var)
{
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_packet)
	ivisitor->visitor->visitChiLAEvent(ivisitor->veto_packet,
					   ivisitor->user_data,
					   table_index,
					   electronics_id,
					   table_event_index,
					   a_omega,
					   a_omega_var);
      if(ivisitor->veto_packet)nveto++;
    }

  if(nveto==m_visitors.size())veto_packet=true;
}

void VBFBroadcastVisitor::
visitKascadeHeader(void* user_data,
		   uint32_t            corsika_particle_id,
		   float               energy_gev,
		   uint32_t            shower_id,
		   float               x_area_width_m,
		   float               y_area_width_m,
		   uint32_t            north_south_grid)
{
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_packet)
	ivisitor->visitor->visitKascadeHeader(ivisitor->user_data,
					      corsika_particle_id,
					      energy_gev,
					      shower_id,
					      x_area_width_m,
					      y_area_width_m,
					      north_south_grid);
    }
}
 
void VBFBroadcastVisitor::
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
  unsigned nveto = 0;
  EventData* ed = static_cast<EventData*>(user_data);
  for(EventData::iterator ivisitor=ed->begin(); 
      ivisitor!=ed->end(); ivisitor++)
    {
      if(!ivisitor->veto_packet)
	ivisitor->visitor->visitKascadeEvent(ivisitor->veto_packet,
					     ivisitor->user_data,
					     nx,
					     ny,
					     direction_index,
					     emission_altitude_m,
					     emission_altitude_sigma,
					     muon_ratio,
					     a_omega,
					     rel_tel_trigger_time_ns,
					     differential_rate_per_event_hz,
					     integral_rate_per_event_hz);
      if(ivisitor->veto_packet)nveto++;
    }

  if(nveto==m_visitors.size())veto_packet=true;
}

VBFBroadcastVisitor::Visitor::~Visitor()
{
  // nothing to see here
}

void VBFBroadcastVisitor::Visitor::stopProcessingFile()
{
  m_parent->stopDispatcherForVisitor(this);
}

// ============================================================================
//
// VBFPacketTransform - Arbitrary packet manipulation
//
// ============================================================================

VBFPacketTransform::~VBFPacketTransform()
{
  // nothing to see here
}

VBFIntegrateOverflowPacketTransform::
VBFIntegrateOverflowPacketTransform(bool overwrite, bool discard_vaev, 
				    bool leave_in_ovrf):
  m_overwrite(overwrite), m_discard_vaev(discard_vaev), 
  m_leave_in_ovrf(leave_in_ovrf)
{
  // nothing to see here
}

VBFIntegrateOverflowPacketTransform::
~VBFIntegrateOverflowPacketTransform()
{
  // nothing to see here
}

void VBFIntegrateOverflowPacketTransform::transform(VPacket* packet)
{
  integrateOverflowDatums(packet, 
			  m_overwrite, m_discard_vaev, m_leave_in_ovrf);
}

void VBFIntegrateOverflowPacketTransform::
integrateOverflowDatums(VPacket* packet, 
			bool overwrite, bool discard_vaev, bool leave_in_ovrf)
{
  // Purpose: Transfer datums from overflow bank back to array event bank
  // Author:  SJF 2007-04-20
  // Input:   packet:        VBF packet to process
  //          overwrite:     An OVRF datum should overwrite existing VAEV
  //          discard_vaev:  Discard VAEV if OVRF exists (NOT RECOMMENDED!!)
  //          leave_in_ovrf: VEvents should remain in OVRF even after they have
  //                         been added to VAEV, perhaps for counting
  // Output:  packet:        Modified as appropriate

  if(!packet->hasEventOverflow())return;

  // 1 - Make lists of what datums are available in the overflow bank ---------

  VEventOverflow* ovrf_old = packet->getEventOverflow();
  unsigned ovrf_ndatum = ovrf_old->numDatums();

  std::vector<const VArrayTrigger*> ovrf_at;
  std::vector<const VEvent*> ovrf_ev;
  std::vector<const VEvent*> ovrf_dup_ev;
  
  for(unsigned idatum=0;idatum<ovrf_ndatum;idatum++)
    {
      VDatum* datum = ovrf_old->getDatumAt(idatum);

      const VArrayTrigger* at = dynamic_cast<const VArrayTrigger*>(datum);
      if(at)
	{
	  ovrf_at.push_back(at);
	  continue;
	}

      const VEvent* ev = dynamic_cast<const VEvent*>(datum);
      if(ev)
	{
	  unsigned inode = ev->getNodeNumber();
	  if(inode>=ovrf_ev.size())ovrf_ev.resize(inode+1);
	  if(ovrf_ev[inode]==0)ovrf_ev[inode]=ev;
	  else ovrf_dup_ev.push_back(ev);
	  continue;
	}
    }

  // 2 - Set up the array event -----------------------------------------------

  VArrayEvent* ae = 0;
  if(packet->hasArrayEvent() && !discard_vaev)
    ae = packet->getArrayEvent();
  else
    {
      ae = new VArrayEvent;
      packet->putArrayEvent(ae);
    }

  // 3 - Make new overflow bank to copy unusable datums into ------------------

  VEventOverflow* ovrf_new = leave_in_ovrf?0:new VEventOverflow;

  // 4 - Process VEvents - add them to the array event or new overflow --------

  for(unsigned ievent=0;ievent<ae->getNumEvents();ievent++)
    {
      VEvent* ev = ae->getEventAt(ievent);
      if(ev==NULL)continue;
      unsigned inode = ev->getNodeNumber();
      if(inode<ovrf_ev.size() && ovrf_ev[inode]!=0)
	{
	  const VEvent* matched_ev = ovrf_ev[inode];
	  if(overwrite)ae->replaceEventAt(ievent,matched_ev->copyEvent());
	  else if(ovrf_new)ovrf_new->addDatum(matched_ev->copyEvent());
	  ovrf_ev[inode] = 0;
	}
    }

  for(unsigned ievent=0;ievent<ovrf_ev.size();ievent++)
    if(ovrf_ev[ievent])ae->addEvent(ovrf_ev[ievent]->copyEvent());

  if(ovrf_new)
    for(unsigned ievent=0;ievent<ovrf_dup_ev.size();ievent++)
      ovrf_new->addDatum(ovrf_dup_ev[ievent]->copyEvent());

  // 5 - Process the VArrayTrigger --------------------------------------------

  if((!ovrf_at.empty()) && (!ae->hasTrigger() || overwrite))
    {
      ae->setTrigger(ovrf_at.front()->copyAT());
      ovrf_at.front()=0;
    }


  if(ovrf_new)
    for(unsigned iat=0;iat<ovrf_at.size();iat++)
      if(ovrf_at[iat])ovrf_new->addDatum(ovrf_at[iat]->copyAT());

  // 6 - Replace the old overflow with the new one ----------------------------

  if(ovrf_new)
    if(ovrf_new->numDatums())packet->putEventOverflow(ovrf_new);
    else 
      {
	packet->remove(VGetEventOverflowBankName());
	delete ovrf_new;
      }
}

// ============================================================================
//
// VBFRationalizedClock
//
// ============================================================================

VBFRationalizedClock::~VBFRationalizedClock()
{
  // nothing to see here
}

VBFStandardRationalizedClock::~VBFStandardRationalizedClock()
{
  // nothing to see here
}

bool VBFStandardRationalizedClock::
bestEventTime(const VBFPacket* packet, VSTime& best_event_time)
{
  const VPacket* vbf = packet->vbf;
  if((vbf == 0)||(!vbf->hasArrayEvent()))
    {
      m_last_event_num = 0;
      return false;
    }

  VArrayEvent* ae = vbf->getArrayEvent();

  std::vector<Clock> clocks;

  bool ae_has_trigger = ae->hasTrigger();
  unsigned max_clock_num = 0;
  if((!m_dont_use_l3)&&(ae_has_trigger))
    {
      VArrayTrigger* at = 0;
      at = ae->getTrigger();

      uword16 at_gps_year = at->getGPSYear();

      VSTime at_event_time;
      at_event_time.setFromVBF(at_gps_year, 
			       at->getGPSTimeNumElements(),
			       at->getGPSTime());

      if((at_event_time.decodedOK())
	 &&(at_event_time>=m_time0)
	 &&(at_event_time<=m_time1))
	{
	  Clock at_clock;
	  at_clock.clock      = at_event_time;
	  at_clock.clock_good = at_event_time.statusOK();
	  at_clock.clock_num  = 0;
	  at_clock.priority   = 0;
	  clocks.push_back(at_clock);
	  max_clock_num       = 1;
	}
    }

  unsigned num_events = ae->getNumEvents();
  for(unsigned telescope_num=0; telescope_num<num_events; 
      telescope_num++)
    {
      VEvent* ev = ae->getEvent(telescope_num);
      if(ev!=NULL)
	{
	  uword16 ev_gps_year = ev->getGPSYear();

	  unsigned iscope = ev->getNodeNumber();

	  VSTime scope_event_time;
	  scope_event_time.setFromVBF(ev_gps_year,
				      ev->getGPSTimeNumElements(),
				      ev->getGPSTime());

	  if((scope_event_time.decodedOK())
	     &&(scope_event_time>=m_time0)
	     &&(scope_event_time<=m_time1))
	    {
	      Clock scope_clock;
	      scope_clock.clock      = scope_event_time;
	      scope_clock.clock_good = scope_event_time.statusOK();
	      scope_clock.clock_num  = iscope+1;
	      scope_clock.priority   = 0;
	      clocks.push_back(scope_clock);
	      if(max_clock_num<iscope+2)
		max_clock_num        = iscope+2;
	    }
	}
    }

  if(m_clock_count.size()<max_clock_num)
    {
      m_clock_count.resize(max_clock_num);
      m_clock_seen.resize(max_clock_num);
    }

  unsigned nclock = clocks.size();
  if(nclock==0)return false;

  for(unsigned iclock=0;iclock<nclock;iclock++)
    if(m_clock_seen[clocks[iclock].clock_num]==false)
      {
	m_nclock_seen++;
	m_clock_seen[clocks[iclock].clock_num]=true;
      }

  unsigned clock_chosen = 0;

  if(nclock==2)
    {
      int64_t dclock0 = clocks[0].clock-m_last_event_time;
      if(dclock0<0)dclock0=-dclock0;
      int64_t dclock1 = clocks[1].clock-m_last_event_time;
      if(dclock1<0)dclock1=-dclock1;
      if(dclock1<dclock0)clock_chosen=1;
    }
  else if(nclock>2)
    {
      unsigned iclock = 0;
      unsigned jclock = 1;

      unsigned min_iclock = iclock;
      unsigned min_jclock = jclock;
      int64_t min_dclock = clocks[iclock].clock-clocks[jclock].clock;
      if(min_dclock<0)min_dclock=-min_dclock;

      for(unsigned iclock=0;iclock<nclock;iclock++)
	for(unsigned jclock=iclock+1;jclock<nclock;jclock++)
	  {
	    int64_t dclock = clocks[iclock].clock-clocks[jclock].clock;
	    if(dclock<0)dclock=-dclock;
	    if(dclock<min_dclock)
	      min_dclock=dclock,min_iclock=iclock,min_jclock=jclock;
	  }

      iclock = min_iclock;
      jclock = min_jclock;

      unsigned icount = m_clock_count[clocks[iclock].clock_num];
      unsigned jcount = m_clock_count[clocks[jclock].clock_num];

      if((clocks[iclock].clock_good)&&(!clocks[jclock].clock_good))
	clock_chosen = iclock;
      else if((clocks[jclock].clock_good)&&(!clocks[iclock].clock_good))
	clock_chosen = jclock;
      else if(clocks[iclock].priority < clocks[jclock].priority)
	clock_chosen = iclock;
      else if(clocks[jclock].priority < clocks[iclock].priority)
	clock_chosen = jclock;
      else if(icount >= jcount)
	clock_chosen = iclock;
      else
	clock_chosen = jclock;
    }

  m_clock_count[clocks[clock_chosen].clock_num]++;
  best_event_time = clocks[clock_chosen].clock;

#if 0
  std::cout << nclock << ' ' << clock_chosen << ' ' 
	    << clocks[clock_chosen].clock_num << ' '
	    << best_event_time.toString();
  for(unsigned iclock=0;iclock<m_clock_count.size();iclock++)
    std::cout << ' ' << m_clock_count[iclock];
  std::cout << std::endl;
#endif  

  m_last_event_num = packet->vbf_num;
  m_last_event_time = best_event_time;
  return true;
}

// ============================================================================
// 
// VASequentialPacketFetcher
// 
// ============================================================================

VBFSequentialPacketFetcher::~VBFSequentialPacketFetcher()
{
  // nothing to see here
}

VBFSequentialPacketFetcherChained::~VBFSequentialPacketFetcherChained()
{
  // nothing to see here
}

VBFPacket* VBFSequentialPacketFetcherChained::getPacket()
{
  return fFetcher->getPacket();
}

// ----------------------------------------------------------------------------
// VBFSimpleFetcher
// ----------------------------------------------------------------------------

VBFSimpleFetcher::
VBFSimpleFetcher(VBankFileReader* reader, VBFRationalizedClock* clock,
		 VBFPacketTransform* transform):
  VBFSequentialPacketFetcher(),
  fReader(reader), fClock(clock), fTransform(transform),
  fIPacket(0), fNPacket(fReader->numPackets())
{
  // nothing to see here
}

VBFSimpleFetcher::~VBFSimpleFetcher()
{
  // nothing to see here
}

VBFPacket* VBFSimpleFetcher::getPacket()
{
  if(fIPacket==fNPacket)return 0;
  VBFPacket* p = new VBFPacket;
  p->vbf_num               = fIPacket;
  p->vbf                   = fReader->readPacket(fIPacket++);
  if(fTransform)fTransform->transform(p->vbf);
  p->has_best_event_time   = fClock->bestEventTime(p, p->best_event_time);
  return p;
}

// ----------------------------------------------------------------------------
// VBFLinearizedFetcher
// ----------------------------------------------------------------------------

VBFLinearizedFetcher::
VBFLinearizedFetcher(VBFSequentialPacketFetcher* fetcher,
		     unsigned buffer_size,
		     bool discard_if_out_of_sequence,
		     bool verbose):
  VBFSequentialPacketFetcherChained(fetcher),
  fDiscard(discard_if_out_of_sequence), fVerbose(verbose),
  fIPacket(0), fBufferSize(buffer_size), fBuffer(), fNextEvent(0)
{
  // nothing to see here
}

VBFLinearizedFetcher::~VBFLinearizedFetcher()
{
  while(!fBuffer.empty()) 
    {
      delete fBuffer.front().packet;
      fBuffer.pop_front(); 
    }
}

VBFPacket* VBFLinearizedFetcher::getPacket()
{
  unsigned next_event_num_in_buffer = 0;

  unsigned iloop = 0;
  while(1)
    {
      if(!fBuffer.empty())
	{
	  // Always send an invalid node when it reaches the front

	  if(!fBuffer.front().valid)
	    {
	      VBFPacket* packet = fBuffer.front().packet;
	      fBuffer.pop_front();
	      return packet;
	    }

	  // Search list for next packet and simultaneously find lowest
	  // event number in the list for later

	  std::list<BufNode>::iterator inode = fBuffer.begin();
	  while(inode != fBuffer.end())
	    {
	      if(inode->valid)
		{
		  const unsigned inode_event_num = inode->event_num;
		  if(inode_event_num == fNextEvent)
		    {
		      // Found the packet we want, take it off the list and
		      // advance the event number

		      VBFPacket* packet = inode->packet;
		      fBuffer.erase(inode);
		      fNextEvent++;
		      return packet;
		    }

		  if((next_event_num_in_buffer==0)||
		     (inode_event_num<next_event_num_in_buffer))
		    next_event_num_in_buffer = inode_event_num;
		}
	      inode++;
	    }
	}

      vsassert(iloop == 0);

      // We are here when either the buffer is empty OR the buffer has
      // at least one valid packet (but not the one we want) and
      // next_event_num_in_buffer has the lowest event number on the
      // list

      while(fBuffer.size() < fBufferSize)
	{
	  // We only try to read a packet if there is space on the
	  // buffer to put it

	  VBFPacket* packet = VBFSequentialPacketFetcherChained::getPacket();
	  if(packet==0)
	    {
	      // If there are no more packets then either return or
	      // advance the event to the next one on the list

	      if(fBuffer.empty())return 0;
	      else break;
	    }

	  fIPacket++;

	  bool valid = packet->vbf->hasArrayEvent();
	  VArrayEvent* ae = 0;
	  unsigned event_num = 0;
	  if(valid)
	    {
	      ae = packet->vbf->getArrayEvent();
	      valid = ae->hasEventNumber();
	      if(valid)
		{
		  event_num = ae->getEventNumber();
		  valid = (event_num >= fNextEvent);
		}
	    }

	  // If packet is invalid we discard it, if requested. Otherwise 
	  // we return it immediately if the buffer is empty or add it
	  // at the end of the buffer otherwise

	  if(!valid)
	    {
	      if(fDiscard)
		{
		  if((fVerbose)&&(packet->vbf_num))
		    std::cout << __PRETTY_FUNCTION__ << ": discarding packet "
			      << fIPacket-1 << std::endl;
		  delete packet;
		  continue;
		}
	      else if(fBuffer.empty())
		{
		  return packet;
		}
	      else 
		{
		  BufNode node;
		  node.valid          = valid;
		  node.event_num      = event_num;
		  node.packet         = packet;
		  fBuffer.push_back(node);
		  continue;
		}
	    }

	  // If it is the one we are looking for then do a fast return

	  if(event_num == fNextEvent)
	    {
	      fNextEvent++;
	      return packet;
	    }

	  // Put the packet onto the buffer

	  BufNode node;
	  node.valid          = valid;
	  node.event_num      = event_num;
	  node.packet         = packet;
	  fBuffer.push_back(node);

	  // Maintain the invariant that: either the buffer is empty OR
	  // the buffer has at least one valid packet (but not the one
	  // we want) and next_event_num_in_buffer has the lowest event
	  // number on the list

	  if((next_event_num_in_buffer==0)||
	     (event_num<next_event_num_in_buffer))
	    next_event_num_in_buffer = event_num;
	}

      // Advance next to the lowest event on the buffer and go around again.
      // Next time around the loop we are guarenteed to find the packet

      if((fVerbose)&&(fNextEvent))
	std::cout << __PRETTY_FUNCTION__ << ": event " << fNextEvent 
		  << " not found, skipping to "
		  << next_event_num_in_buffer << std::endl;
      fNextEvent = next_event_num_in_buffer;
      iloop++;
    }

  vsassert(0);
  return 0;
}

// ----------------------------------------------------------------------------
// VBFRandomizedFetcher
// ----------------------------------------------------------------------------

VBFRandomizedFetcher::
VBFRandomizedFetcher(VBFSequentialPacketFetcher* fetcher,
		     unsigned buffer_size,
		     float discard_probability,
		     float long_term_keep_probability,
		     bool verbose):
  VBFSequentialPacketFetcherChained(fetcher), 
  fDiscardFate(unsigned(discard_probability*float(RAND_MAX))),
  fLongTermFate(unsigned((discard_probability+long_term_keep_probability)
			 * float(RAND_MAX))),
  fVerbose(verbose), fIPacket(0),
  fBufferSize(buffer_size), fBuffer(buffer_size,0), fLongTermBuffer(0)
{
  time_t t;
  time(&t);
  srand(t);
}

VBFRandomizedFetcher::~VBFRandomizedFetcher()
{
  for(std::vector<VBFPacket*>::iterator ipacket=fBuffer.begin();
      ipacket!=fBuffer.end();ipacket++)delete *ipacket;
  delete fLongTermBuffer;
}

VBFPacket* VBFRandomizedFetcher::getPacket()
{
  VBFPacket* r_packet = VBFSequentialPacketFetcherChained::getPacket();
  while(r_packet)
    {
      unsigned fate = rand();
      if(fate<fDiscardFate)
	{
	  if(fVerbose)
	    std::cout << __PRETTY_FUNCTION__ << ": discarding packet"
		      << std::endl;
	  delete r_packet;
	}
      else if(fate<fLongTermFate)
	{
	  VBFPacket* s_packet = fLongTermBuffer;
	  fLongTermBuffer = r_packet;
	  if(s_packet)
	    {
	      fIPacket++;
	      return s_packet;
	    }
	}
      else
	{
	  unsigned index = rand()%fBufferSize;
	  VBFPacket* s_packet = fBuffer[index];
	  fBuffer[index] = r_packet;
	  if(s_packet)
	    {
	      fIPacket++;
	      return s_packet;
	    }
	}
      r_packet = VBFSequentialPacketFetcherChained::getPacket();
    }

  for(std::vector<VBFPacket*>::iterator ipacket=fBuffer.begin();
      ipacket!=fBuffer.end();ipacket++)
    {
      if(*ipacket)
	{
	  VBFPacket* s_packet = *ipacket;
	  *ipacket = 0;
	  fIPacket++;
	  return s_packet;
	}
    }
  
  if(fLongTermBuffer)
    {
      VBFPacket* s_packet = fLongTermBuffer;
      fLongTermBuffer = 0;
      fIPacket++;
      return s_packet;
    }

  return 0;
}

// ----------------------------------------------------------------------------
// VBFThreadedFetcher
// ----------------------------------------------------------------------------

#ifndef NOTHREADS
VBFThreadedFetcher::
VBFThreadedFetcher(VBFSequentialPacketFetcher* fetcher,
		   unsigned buffer_size):
  VBFSequentialPacketFetcherChained(fetcher), 
  fFinished(false), fBufferSize(buffer_size), fBuffer(),
  fThread(), fMutex(), fConditionNotifyConsumer(), fConditionNotifyProducer()
{
  vsassert(pthread_mutex_init(&fMutex, NULL) == 0);
  vsassert(pthread_cond_init(&fConditionNotifyConsumer, NULL) == 0);
  vsassert(pthread_cond_init(&fConditionNotifyProducer, NULL) == 0);
  vsassert(pthread_create(&fThread, 0, &threadStart, this) == 0);
}
 
VBFThreadedFetcher::~VBFThreadedFetcher()
{
  vsassert(pthread_mutex_lock(&fMutex) == 0);
  fFinished=true;
  vsassert(pthread_mutex_unlock(&fMutex) == 0);
  vsassert(pthread_cond_broadcast(&fConditionNotifyProducer) == 0);
  vsassert(pthread_cond_broadcast(&fConditionNotifyConsumer) == 0);
  vsassert(pthread_join(fThread,NULL) == 0);
  vsassert(pthread_cond_destroy(&fConditionNotifyConsumer) == 0);
  vsassert(pthread_cond_destroy(&fConditionNotifyProducer) == 0);
  vsassert(pthread_mutex_destroy(&fMutex) == 0);
  while(!fBuffer.empty()) 
    {
      delete fBuffer.front().packet;
      fBuffer.pop_front(); 
    }
}

VBFPacket* VBFThreadedFetcher::getPacket()
{
  vsassert(pthread_mutex_lock(&fMutex) == 0);

  while(fBuffer.empty())
    {
      if(fFinished)
	{
	  vsassert(pthread_mutex_unlock(&fMutex) == 0);
	  return 0;
	}

      vsassert(pthread_cond_wait(&fConditionNotifyConsumer,&fMutex) == 0);
    }

  VBFPacket* packet = fBuffer.front().packet;
  fBuffer.pop_front();
  vsassert(pthread_cond_signal(&fConditionNotifyProducer) == 0);
  vsassert(pthread_mutex_unlock(&fMutex) == 0);
  return packet;
}

void* VBFThreadedFetcher::threadStart(void *object)
{
  static_cast<VBFThreadedFetcher*>(object)->run();
  return 0;
}

void VBFThreadedFetcher::run()
{
  vsassert(pthread_mutex_lock(&fMutex) == 0);

  while(!fFinished)
    {
      if(fBuffer.size() == fBufferSize)
	{
	  vsassert(pthread_cond_wait(&fConditionNotifyProducer,&fMutex) == 0);
	}
      else
	{
	  vsassert(pthread_mutex_unlock(&fMutex) == 0);  
	  BufNode node;
	  node.packet = VBFSequentialPacketFetcherChained::getPacket();
	  vsassert(pthread_mutex_lock(&fMutex) == 0);
	  if(node.packet == 0)
	    fFinished = true;
	  else
	    fBuffer.push_back(node);
	  vsassert(pthread_cond_broadcast(&fConditionNotifyConsumer) == 0);
	}
    }

  vsassert(pthread_mutex_unlock(&fMutex) == 0);
}
#endif

// ============================================================================
// 
// VSSimpleVBFDispatcherStop
// 
// ============================================================================

VSSimpleVBFDispatcherStop::~VSSimpleVBFDispatcherStop()
{
  // nothing to see here
}

// ============================================================================
// 
// VSSimpleVBFDispatcher
// 
// ============================================================================

VSSimpleVBFDispatcher::VSSimpleVBFDispatcher(VSSimpleVBFVisitor* visitor,
					     VBFRationalizedClock* clock):
  VSSimpleVBFDispatcherStop(), 
#ifndef NOTHREADS
  fProcessingMutex(),
#endif
  fMyClock(), fClock(clock), fVisitor(visitor), fFilename(), fReader(), 
  fStopProcessingFlag(false), fVisitChannels(true), fDontUseL3Time(false)
{ 
  if(fClock==0)
    {
      fMyClock = new VBFStandardRationalizedClock(fDontUseL3Time);
      fClock = fMyClock;
    }
  fVisitor->registerDispatcher(this);
}

void VSSimpleVBFDispatcher::resetClock(VBFRationalizedClock* clock)
{
  if(clock != 0)
    {
      delete fMyClock;
      fMyClock = 0;
      fClock = clock;
    }
  else
    {
      delete fMyClock;
      fMyClock = new VBFStandardRationalizedClock(fDontUseL3Time);
      fClock = fMyClock;
    }
}

void VSSimpleVBFDispatcher::resetVisitor(VSSimpleVBFVisitor* visitor)
{
  if((fReader)&&(fVisitor))fVisitor->leaveFile();
  fVisitor = visitor;
  if(fVisitor)
    {
      fVisitor->registerDispatcher(this);
      if(fReader)fVisitor->visitFile(fFilename.c_str(),fReader->numPackets());
    }
}


VSSimpleVBFDispatcher::~VSSimpleVBFDispatcher()
{
  closeFile();
  delete fMyClock;
  fMyClock = 0;
}

bool VSSimpleVBFDispatcher::openFile(const char* filename)
{
  if(fReader)closeFile();
  fReader = new VBankFileReader(filename);
  fFilename = filename;
  if(fVisitor)fVisitor->visitFile(filename, fReader->numPackets());  
  return true;
}

void VSSimpleVBFDispatcher::closeFile()
{
  if(fReader)
    {
      if(fVisitor)fVisitor->leaveFile();
      delete fReader;
      fReader = 0;
      fFilename = std::string();
    }
}

unsigned VSSimpleVBFDispatcher::numPackets() const
{
  if(fReader)return fReader->numPackets();
  else return 0;
}

void VSSimpleVBFDispatcher::
dispatchPacket(unsigned vbf_packet_num, unsigned seq_packet_num,
	       uint32_t flags)
{
  VBFPacket* p = new VBFPacket;
  p->vbf_num               = vbf_packet_num;
  p->vbf                   = fReader->readPacket(vbf_packet_num);
  if(flags&DISPATCH_INTEGRATE_OVERFLOW)
    VBFIntegrateOverflowPacketTransform::integrateOverflowDatums(p->vbf,false,false,true);
  p->has_best_event_time   = fClock->bestEventTime(p, p->best_event_time);
  dispatchPacket(p, seq_packet_num);
  delete p;
}

void VSSimpleVBFDispatcher::
dispatchPacket(const VBFPacket* packet, unsigned seq_packet_num)
{
  bool veto_packet = false;

  VPacket* vbf = packet->vbf;

  bool has_array_event = vbf->hasArrayEvent();
  bool has_sim_header = vbf->hasSimulationHeader();
  bool has_sim_event = vbf->hasSimulationData();

  unsigned num_overflow_datum = 
    vbf->hasEventOverflow()?vbf->getEventOverflow()->numDatums():0;

  void* user_data = 0;

  fVisitor->visitPacket(veto_packet, user_data, 
			seq_packet_num, packet->vbf_num,
			has_array_event, has_sim_header, has_sim_event,
			num_overflow_datum,
			packet);
  
  if(num_overflow_datum)
    {
      VEventOverflow* overflow = vbf->getEventOverflow();
      for(unsigned idatum=0;idatum<num_overflow_datum;idatum++)
	fVisitor->visitOverflowEvent(user_data,
				     idatum, overflow->getDatumAt(idatum));
    }

  if((!veto_packet)&&(has_sim_header))
    {
      VSimulationHeader *sh = vbf->getSimulationHeader();
      
      VSTime date_of_sims;
      VSTime date_of_array;
      date_of_sims.setFromMSTimeStamp(uint64_t(sh->fDateOfSimsUTC)
				      * UINT64_C(1000000000));
      date_of_array.setFromMSTimeStamp(uint64_t(sh->fDateOfArrayForSims)
				       * UINT64_C(1000000000));
      
      fVisitor->visitSimulationHeader(user_data, 
				      sh->fRunNumber,
				      date_of_sims,
				      sh->fSimulationPackage,
				      sh->fSimulator,
				      date_of_array,
				      sh->fAtmosphericModel,
				      sh->fObsAltitudeM,
				      sh->fArray,
				      sh->fSimConfigFile);
    }

  if((!veto_packet)&&(vbf->has(VGetChiLASimulationHeaderBankName())))
    {
      VChiLASimulationHeader *ch = 
	vbf->get<VChiLASimulationHeader>(VGetChiLASimulationHeaderBankName());

      fVisitor->visitChiLAHeader(user_data, 
				 ch->fDatabaseName,
				 ch->fSimParamTables,
				 ch->fOpticsConfigurations,
				 ch->fElectronicsConfigurations);
    }

  if((!veto_packet)&&(vbf->has(VGetKascadeSimulationHeaderBankName())))
    {
      VKascadeSimulationHeader *kh = 
    vbf->get<VKascadeSimulationHeader>(VGetKascadeSimulationHeaderBankName());

      fVisitor->visitKascadeHeader(user_data, 
				   kh->fCORSIKAParticleID,
				   kh->fEnergyGeV,
				   kh->fShowerID,
				   kh->fXAreaWidthM,
				   kh->fYAreaWidthM,
				   kh->fNorthSouthGrid);
    }

  if((!veto_packet)&&(has_sim_event))
    {
      VSimulationData* se = vbf->getSimulationData();
     
      fVisitor->visitSimulationEvent(veto_packet, user_data,
				     se->fRunNumber,
				     se->fEventNumber,
				     se->fCORSIKAParticleID,  
				     se->fEnergyGeV,
				     se->fObservationZenithDeg,
				     se->fObservationAzimuthDeg,
				     se->fPrimaryZenithDeg,
				     se->fPrimaryAzimuthDeg,
#if((defined(V_SIMULATION_DATA_VERSION))&&(V_SIMULATION_DATA_VERSION>=2))
				     se->fRefZenithDeg,
				     se->fRefAzimuthDeg,
				     se->fRefPositionAngleDeg,
#else
				     0.0,
				     0.0,
				     0.0,
#endif
				     se->fCoreEastM,
				     se->fCoreSouthM,
				     se->fCoreElevationMASL);
    }

  if((!veto_packet)&&(vbf->has(VGetChiLASimulationDataBankName())))
    {
      VChiLASimulationData *cd =
	vbf->get<VChiLASimulationData>(VGetChiLASimulationDataBankName());

      fVisitor->visitChiLAEvent(veto_packet, user_data,
				cd->fSimParamTableID,
				cd->fElectronicsID,
#if((defined(V_CHILA_SIMULATION_DATA_VERSION))&&(V_CHILA_SIMULATION_DATA_VERSION>=1))
				cd->fTableIndex,
				cd->fAomega,
				cd->fAomegaVar);
#else
                                0,
				0.0,
				0.0);
#endif
    }

  if((!veto_packet)&&(vbf->has(VGetKascadeSimulationDataBankName())))
    {
      VKascadeSimulationData *kd =
	vbf->get<VKascadeSimulationData>(VGetKascadeSimulationDataBankName());

      fVisitor->visitKascadeEvent(veto_packet, user_data,
				  kd->fNx,
				  kd->fNy,
				  kd->fDirectionIndex,
				  kd->fEmissionAltitudeM,
				  kd->fEmissionAltitudeSigma,
				  kd->fMuonRatio,
				  kd->fAomega,
				  kd->fRelTelTriggerTimeNS,
				  kd->fDifferentialRatePerEventHz,
				  kd->fIntegralRatePerEventHz);
    }

  if((!veto_packet)&&(has_array_event))
    {
      bool veto_array_event = veto_packet;

      VArrayEvent* ae = vbf->getArrayEvent();
      bool ae_has_event_number = ae->hasEventNumber();
      bool ae_has_trigger = ae->hasTrigger();

      VSSimpleVBFVisitor::EventType event_type;
      uint32_t l2_trigger_mask;
      getArrayEventData(event_type, l2_trigger_mask, ae);

      fVisitor->visitArrayEvent(veto_array_event, user_data,
				ae->getNumEvents(),
				ae_has_trigger,
				ae_has_event_number, 
				ae_has_event_number?ae->getEventNumber():0,
				packet->has_best_event_time, 
				packet->best_event_time,
				event_type,
				l2_trigger_mask,
				ae);

      if((!veto_array_event)&&(ae_has_trigger))
	{
	  VArrayTrigger* at = ae->getTrigger();
	  veto_array_event = dispatchArrayTrigger(fVisitor, user_data, at);
	}

      if(!veto_array_event)      
	{
	  unsigned nevent = ae->getNumEvents();
	  for(unsigned ievent=0; ievent<nevent; ievent++)
	    {
	      VEvent* ev = ae->getEvent(ievent);
	      if(ev!=NULL)
		dispatchScopeEvent(fVisitor, user_data,
				   ievent, ev, fVisitChannels);
	    }
	}
      fVisitor->leaveArrayEvent(veto_array_event, user_data);
    }

  fVisitor->leavePacket(veto_packet, user_data);
}

void VSSimpleVBFDispatcher::
getArrayEventData(VSSimpleVBFVisitor::EventType& event_type,
		  uint32_t& l2_trigger_mask,
		  const VArrayEvent* ae)
{
  event_type = VSSimpleVBFVisitor::ET_UNKNOWN;
  l2_trigger_mask = 0;

  if(ae->hasTrigger())
    {
      // Use array trigger to determine event type and mask

      VArrayTrigger* at = ae->getTrigger();
      switch(at->getEventType().trigger)
	{
	case VEventType::L2_TRIGGER:
	  event_type = VSSimpleVBFVisitor::ET_L2;
	  break;
	case VEventType::PED_TRIGGER:
	  event_type = VSSimpleVBFVisitor::ET_PED;
	  break;
	case VEventType::HIGH_MULT_TRIGGER:
	case VEventType::NEW_PHYS_TRIGGER:
	case VEventType::CAL_TRIGGER:
	default:
	  event_type = VSSimpleVBFVisitor::ET_UNKNOWN;
	  break;
	}

      if(event_type == VSSimpleVBFVisitor::ET_L2)
	{
	  for(unsigned i=0; i<at->getNumTriggerTelescopes(); i++)
	    {
	      uint32_t scope_id = at->getTriggerTelescopeId(i);
	      l2_trigger_mask |= 0x00000001 << scope_id;
	    }
	}
    }
  else
    {
      // Array trigger missing, try to determine type from telescope 
      // datums

      unsigned nevent = ae->getNumEvents();
      for(unsigned ievent=0; ievent<nevent; ievent++)
	{
	  VEvent* ev = ae->getEvent(ievent);
	  if(ev!=NULL)
	    {
	      switch(ev->getEventType().trigger)
		{
		case VEventType::L2_TRIGGER:
		  if(event_type == VSSimpleVBFVisitor::ET_UNKNOWN)
		    event_type = VSSimpleVBFVisitor::ET_L2;
		  break;
		case VEventType::PED_TRIGGER:
		  event_type = VSSimpleVBFVisitor::ET_PED;
		  break;
		case VEventType::HIGH_MULT_TRIGGER:
		case VEventType::NEW_PHYS_TRIGGER:
		case VEventType::CAL_TRIGGER:
		default:
		  break;
		}
	    }
	}
    }

  // Use per-telescope event type to determine trigger mask
  unsigned nevent = ae->getNumEvents();
  for(unsigned ievent=0; ievent<nevent; ievent++)
    {
      VEvent* ev = ae->getEvent(ievent);
      if(ev!=NULL)
	{
	  if((ev->getEventType().trigger == VEventType::L2_TRIGGER)
	     &&(ev->getEventType().xtra_samples == false))
	    l2_trigger_mask |= 0x00000001 << ev->getNodeNumber();
	}
    }
}

bool VSSimpleVBFDispatcher::
dispatchArrayTrigger(VSSimpleVBFVisitor* visitor, void* user_data, 
		     const VArrayTrigger* at)
{
  bool veto_array_event = false;

  uword16 at_gps_year = at->getGPSYear();
  VSTime gps_time(at_gps_year, 
		  at->getGPSTimeNumElements(),
		  at->getGPSTime());

  uint32_t op_trigger_mask=at->getTriggerMask();

  visitor->visitArrayTrigger(veto_array_event, user_data,
			     at->getEventNumber(),
			     at->getEventType(),
			     op_trigger_mask,
			     at->getFlags(),
			     gps_time,
			     at->getATFlags(),
			     at->getConfigMask(),
			     at->getNumSubarrayTelescopes(),
			     at->getNumTriggerTelescopes(),
			     at->getRunNumber(),
			     at->getTenMHzClockArray(),
			     at->getOptCalCountArray(),
			     at->getPedCountArray(),
			     at);

  uint32_t ip_trigger_mask=0;
  for(uword32 i=0; i<at->getNumTriggerTelescopes(); i++)
    {
      uint32_t scope_id = at->getTriggerTelescopeId(i);
      ip_trigger_mask |= 0x00000001 << scope_id;
    }
	  
  for(uword32 i=0; i<at->getNumSubarrayTelescopes(); i++)
    {
      unsigned scope_id = at->getSubarrayTelescopeId(i);
      visitor->
	visitArrayTriggerScope(veto_array_event, user_data,
			       scope_id,
			       ip_trigger_mask&(1<<scope_id),
			       at->getSpecificRawEventTypeCode(i),
			       at->getAltitude(i),
			       at->getAzimuth(i),
			       at->getTDCTime(i),
			       at->getShowerDelay(i),
			       at->getCompDelay(i),
			       at->getL2CountsArray(i),
			       at->getCalCountsArray(i));
    }

  return veto_array_event;
}

bool VSSimpleVBFDispatcher::
dispatchScopeEvent(VSSimpleVBFVisitor* visitor, void* user_data, 
		   unsigned ievent, const VEvent* ev,
		   bool visit_channels)
{
  bool veto_scope_event = false;

#warning Why do we need to get GPS year from Array Trigger 
  uword16 ev_gps_year = ev->getGPSYear();	      
		  
  VSTime gps_time(ev_gps_year,
		  ev->getGPSTimeNumElements(),
		  ev->getGPSTime());

  unsigned num_samples = ev->getNumSamples();
  unsigned num_channels_stored = ev->getNumChannels();
  unsigned num_channels_total = ev->getMaxNumChannels();

  visitor->visitScopeEvent(veto_scope_event, user_data,
			   ev->getEventNumber(),
			   ev->getNodeNumber(),
			   ev->getEventType(),
			   ev->getTriggerMask(),
			   ev->getFlags(),
			   gps_time,
			   num_samples,
			   num_channels_stored,
			   num_channels_total,
			   ev->getNumClockTrigBoards(),
			   ev);
		
  if((!veto_scope_event)&&(visit_channels))
    {
      uint32_t*__restrict__ q = FASTCALLOC(uint32_t,num_samples);
      uint32_t*__restrict__ Q = FASTCALLOC(uint32_t,num_samples);
    
#ifndef VECTORIZE_FRIENDLY
      uint32_t*const nq = q+num_samples;
#endif
		    
      for(unsigned channel_num=0, stored_channel_num=0;
	  channel_num < num_channels_total; channel_num++)
	{
	  bool veto_channel = veto_scope_event;

	  bool hit = ev->getHitBit(channel_num);
	  bool trigger = ev->getTriggerBit(channel_num);

	  visitor->visitChannel(veto_channel, user_data,
				channel_num, hit, trigger);
	  if(hit)
	    {
	      if(!veto_channel)
		{
		  uint32_t Qdaq = 
		    ev->getChargeWithoutVerify(stored_channel_num);
		  uint32_t P = 
		    ev->getPedestalAndHiLoWithoutVerify(stored_channel_num);

#ifdef VECTORIZE_FRIENDLY
		  const ubyte*__restrict__ samples = 
		    ev->getSamplePtr(stored_channel_num, 0);
		  //ev->getSamplePtrWithoutVerify(stored_channel_num, 0);

		  unsigned Qsum = Q[0] = q[0] = *samples;
		  for(unsigned isample=1;isample<num_samples;isample++)
		    Qsum += q[isample] = samples[isample], Q[isample] = Qsum;
#else
		  uint32_t* iq = q;
		  uint32_t* iQ = Q;

		  const ubyte* isample = 
		    ev->getSamplePtr(stored_channel_num, 0);

		  uint32_t Q1 = 0;
		  uint32_t q1;
		  while(iq<nq)
		    {
		      q1 = *(isample++);
		      Q1 += q1;
		      *(iq++) = q1;
		      *(iQ++) = Q1;
		    }
#endif

		  visitor->
		    visitHitChannel(user_data, channel_num,
				    Qdaq, P&32767, P&32768,
				    num_samples, q, Q);
		}
	      stored_channel_num++;
	    }
			  
	  visitor->leaveChannel(veto_channel, user_data);
	}
      FASTFREE(Q);
      FASTFREE(q);
    }
		  
  visitor->leaveScopeEvent(veto_scope_event, user_data);

  return veto_scope_event;
}

unsigned VSSimpleVBFDispatcher::
dispatchAllPackets(unsigned num_packets_max, uint32_t flags)
{
  sGlobalStopProcessingFlag = fStopProcessingFlag = false;
  unsigned seq_packet_num;

  std::list<VBFSequentialPacketFetcher*> fetcher_stack;

  VBFPacketTransform* transform = 0;
  if(flags&DISPATCH_INTEGRATE_OVERFLOW)
    transform = new VBFIntegrateOverflowPacketTransform;  

  fetcher_stack.push_front(new VBFSimpleFetcher(fReader, fClock, transform));

  if(flags&DISPATCH_SHUFFLE)
    {
      float discard_prob = 0.001;
      if(flags&DISPATCH_SHUFFLE_LOSSLESS)discard_prob = 0;
      VBFSequentialPacketFetcher* fetcher = 
	new VBFRandomizedFetcher(fetcher_stack.front(),
				 5, discard_prob, 0.01,
				 (flags&DISPATCH_VERBOSE)?true:false);
      fetcher_stack.push_front(fetcher);
    }

  if(flags&DISPATCH_IN_ORDER)
    {
      VBFSequentialPacketFetcher* fetcher = 
	new VBFLinearizedFetcher(fetcher_stack.front(),
				 50,
				 (flags&DISPATCH_DISCARD_INVALID)?true:false,
				 (flags&DISPATCH_VERBOSE)?true:false);
      fetcher_stack.push_front(fetcher);
    }

  #ifndef NOTHREADS
  if(flags&DISPATCH_THREADED)
    {
      VBFSequentialPacketFetcher* fetcher = 
	new VBFThreadedFetcher(fetcher_stack.front(),50);
      fetcher_stack.push_front(fetcher);
    }
  #endif

  VBFSequentialPacketFetcher* fetcher = fetcher_stack.front();

#ifndef NOTHREADS
  unsigned nthread = DISPATCH_NTHREAD_GET(flags);
  if(nthread)
    {
      fVisitor->usingThreads(nthread);
      seq_packet_num = 0;
      pthread_mutex_t mutex;
      vsassert(pthread_mutex_init(&mutex, NULL) == 0);
      std::list<DispatcherThreadAssist*> threads;
      for(unsigned ithread=0;ithread<nthread;ithread++)
	{
	  DispatcherThreadAssist* thread =
	    new DispatcherThreadAssist(this, fetcher, mutex, seq_packet_num,
				       num_packets_max, fStopProcessingFlag,
				       sGlobalStopProcessingFlag);
	  threads.push_back(thread);
	}
      for(unsigned ithread=0;ithread<nthread;ithread++)
	{
	  threads.front()->join();
	  delete(threads.front());
	  threads.pop_front();
	}
      vsassert(pthread_mutex_destroy(&mutex) == 0);
    }
  else
#else
  if(1)
#endif
    {
      for(seq_packet_num=0; 
	  ((num_packets_max==0)||(seq_packet_num<num_packets_max))
	    &&(!fStopProcessingFlag)
	    &&(!sGlobalStopProcessingFlag); seq_packet_num++)
	{
	  VBFPacket* packet = fetcher->getPacket();
	  if(packet==0)break;
	  dispatchPacket(packet, seq_packet_num);
	  delete packet;
	}
    }
  
  while(!fetcher_stack.empty())
    {
      delete fetcher_stack.front();
      fetcher_stack.pop_front();
    }

  delete transform;

  // Set this to false before we leave in case we are processing nested files
  sGlobalStopProcessingFlag = fStopProcessingFlag = false;
  return seq_packet_num;
}

unsigned VSSimpleVBFDispatcher::
processFile(const char* filename, unsigned num_packets_max, uint32_t flags)
{
  if(!openFile(filename))return false;
  unsigned packet_num = dispatchAllPackets(num_packets_max, flags);
  closeFile();
  return packet_num;
}

void VSSimpleVBFDispatcher::stopProcessingFile()
{
#ifndef NOTHREADS
  if(fProcessingMutex)vsassert(pthread_mutex_lock(fProcessingMutex) == 0);
#endif

  fStopProcessingFlag = true;

#ifndef NOTHREADS
  if(fProcessingMutex)vsassert(pthread_mutex_unlock(fProcessingMutex) == 0);
#endif
}

void VSSimpleVBFDispatcher::catchSignalAndStopProcessingFile(int signum)
{
  signal(signum,&catchHandler);
}

void VSSimpleVBFDispatcher::catchHandler(int signum)
{
  psignal(signum,"VSSimpleVBFDispatcher::catchHandler caught signal");
  sGlobalStopProcessingFlag=true;
  signal(signum,SIG_DFL); // so that second signal kills (or whatever)
}

bool VSSimpleVBFDispatcher::sGlobalStopProcessingFlag = false;

#ifndef NOTHREADS

VSSimpleVBFDispatcher::DispatcherThreadAssist::
DispatcherThreadAssist(VSSimpleVBFDispatcher* dispatcher,
		       VBFSequentialPacketFetcher* fetcher,
		       pthread_mutex_t& fetcher_mutex,
		       unsigned& seq_packet_num,
		       unsigned num_packets_max,
		       bool& stop_processing_flag,
		       bool& global_stop_processing_flag):
  fDispatcher(dispatcher), fFetcher(fetcher), fFetcherMutex(fetcher_mutex),
  fSeqPacketNum(seq_packet_num), fNumPacketsMax(num_packets_max), 
  fStopProcessingFlag(stop_processing_flag), 
  fGlobalStopProcessingFlag(global_stop_processing_flag), 
  fThread(), fFinishedMutex(), fFinishedCond(), fFinished()
{
  vsassert(pthread_cond_init(&fFinishedCond, NULL) == 0);
  vsassert(pthread_mutex_init(&fFinishedMutex, NULL) == 0);
  vsassert(pthread_create(&fThread, 0, &threadStart, this) == 0);  
}

VSSimpleVBFDispatcher::DispatcherThreadAssist::~DispatcherThreadAssist()
{
  vsassert(fFinished);
  vsassert(pthread_cond_destroy(&fFinishedCond) == 0); 
  vsassert(pthread_mutex_destroy(&fFinishedMutex) == 0);
  vsassert(pthread_join(fThread,NULL) == 0);
}

void VSSimpleVBFDispatcher::DispatcherThreadAssist::join()
{
  vsassert(pthread_mutex_lock(&fFinishedMutex) == 0);
  if(!fFinished)
    vsassert(pthread_cond_wait(&fFinishedCond,&fFinishedMutex) == 0);
  vsassert(pthread_mutex_unlock(&fFinishedMutex) == 0);
}

void* VSSimpleVBFDispatcher::DispatcherThreadAssist::threadStart(void *object)
{
  static_cast<DispatcherThreadAssist*>(object)->run();
  return 0;
}

void VSSimpleVBFDispatcher::DispatcherThreadAssist::run()
{
  vsassert(pthread_mutex_lock(&fFetcherMutex) == 0);

  while((!fStopProcessingFlag)&&(!fGlobalStopProcessingFlag)
	&&((fNumPacketsMax==0)||(fSeqPacketNum<fNumPacketsMax)))
    {
      VBFPacket* packet = fFetcher->getPacket();
      if(!packet)break;
      unsigned old_seq_packet_num = fSeqPacketNum;
      fSeqPacketNum++;
      vsassert(pthread_mutex_unlock(&fFetcherMutex) == 0);
#if 0
      std::cout << "F: " << old_seq_packet_num << " T"
		<< pthread_self() << std::endl;
#endif
      fDispatcher->dispatchPacket(packet, old_seq_packet_num);
      delete packet;
      vsassert(pthread_mutex_lock(&fFetcherMutex) == 0);
    }

  vsassert(pthread_mutex_unlock(&fFetcherMutex) == 0);
  vsassert(pthread_mutex_lock(&fFinishedMutex) == 0);
  fFinished = true;
  vsassert(pthread_mutex_unlock(&fFinishedMutex) == 0);
  vsassert(pthread_cond_broadcast(&fFinishedCond) == 0);
}

#endif
