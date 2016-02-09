//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFRunInfo.cpp

  Extract some information from run

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VBFRunInfo.cpp,v 3.40 2008/09/24 20:03:42 matthew Exp $

*/

#include <Astro.h>
#include <VSPointing.hpp>
#include <VBFRunInfo.hpp>
#include <VSSimCoordTransform.hpp>
#include <VSAAlgebra.hpp>

using namespace VERITAS;
using namespace SEphem;

VBFRunInfo::VBFRunInfo(bool copy_pedestal_packets,
		       double sim_wobble_phi_rad)
  : VSSimpleVBFVisitor(),
    m_copy_pedestal_packets(copy_pedestal_packets),
    m_sim_wobble_phi_rad(sim_wobble_phi_rad),
    m_run_number(), m_nevents(),
    m_first_event_time(), m_first_event_vbf_packet_num(),
    m_lo_event_num(), m_hi_event_num(), 
    m_lo_event_time(), m_hi_event_time(),
    m_config_mask(), m_nsample(), m_nchan(),
    m_event_time(), m_pedestals(), 
    m_sim_energy_limits(), m_sim_energy(), m_sim_data(),
    m_got_run_number(0), m_got_first_event(), m_got_first_event_time(),
    m_vbf_packet_num(), m_has_good_event_time(), 
    m_sim_src_radec(), m_sim_obs_radec(), 
    m_sim_nevent(),
    m_best_event_time(), 
    m_event_num(), m_is_ped_event(), m_packet(0)
{
  m_sim_transform = new VSSimCoordTransform;
}

VBFRunInfo::~VBFRunInfo()
{
  for(std::list<VBFPacket*>::iterator iped=m_pedestals.begin();
      iped!=m_pedestals.end(); iped++)delete *iped;
  delete m_sim_data;
  delete m_sim_transform;
}

void VBFRunInfo::visitFile(const char* filename, unsigned npacket)
{
  m_event_time.reserve(npacket);
}

void VBFRunInfo::leaveFile()
{
  if(m_sim_data)
    {
      if(m_sim_nevent)
	{
	  m_sim_src_radec /= (double)m_sim_nevent;
	  m_sim_obs_radec /= (double)m_sim_nevent;
	}	  
      else
	{
	  m_sim_src_radec = VSAAlgebra::Vec3D::makePolar(M_PI/2.,0);
	  m_sim_obs_radec = VSAAlgebra::Vec3D::makePolar(M_PI/2.,0);
	}

      SphericalCoords src_radec(m_sim_src_radec.theta(),
				m_sim_src_radec.phi());
      SphericalCoords obs_radec(m_sim_obs_radec.theta(),
				m_sim_obs_radec.phi());
      
      m_sim_data->src_ra_deg  = src_radec.phiDeg();
      m_sim_data->src_dec_deg = src_radec.latitudeDeg();
      m_sim_data->obs_ra_deg  = obs_radec.phiDeg();
      m_sim_data->obs_dec_deg = obs_radec.latitudeDeg();
      m_sim_data->wobble_theta_deg = src_radec.separation(obs_radec).deg();
      m_sim_data->wobble_phi_deg = 
	src_radec.compassDirectionTo(obs_radec).deg();
    }

  if(m_sim_data)return;

  std::vector<VSTime> times;
  unsigned nevent = m_event_time.size();
  times.reserve(nevent);
  for(unsigned ievent=0;ievent<nevent;ievent++)
    {
      if((m_event_time[ievent].event_found)
	 &&(m_event_time[ievent].event_time.isGood()))
	times.push_back(m_event_time[ievent].event_time);
    }

  unsigned ntime = times.size();
  if(ntime==0)return;
  
  std::sort(times.begin(), times.end());
 
  VSTime mid_time;
  if(ntime%2 == 1)
    mid_time = times[ntime/2];
  else
    mid_time = times[ntime/2-1]+(times[ntime/2]-times[ntime/2-1])/2;

  const VSTime min_time(mid_time-INT64_C(3600000000000));
  const VSTime max_time(mid_time+INT64_C(3600000000000));

  m_lo_event_time = mid_time;
  m_hi_event_time = mid_time;

  unsigned itime = 0;
  while(times[itime]<min_time)itime++;
  m_lo_event_time = times[itime];

  itime = ntime-1;
  while(times[itime]>max_time)itime--;
  m_hi_event_time = times[itime];
}

void VBFRunInfo::
visitPacket(bool& veto_packet, void*& user_data,
	    uint32_t                   seq_packet_number,
	    uint32_t                   vbf_packet_number,
	    bool                       has_array_event,
	    bool                       has_sim_header,
	    bool                       has_sim_event,
	    uint32_t                   num_overflow_datum,
	    const VBFPacket*           packet)
{
  m_vbf_packet_num=vbf_packet_number;
  m_packet=packet;
}

void VBFRunInfo::
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
  m_has_good_event_time = has_good_event_time;
  m_best_event_time     = best_event_time;
  m_event_num           = event_num;
  m_is_ped_event        = false;
  m_nevents++;

  if(has_good_event_time)
    {
      if(!m_got_first_event_time)
	{
	  m_first_event_time          = m_best_event_time;
	  m_lo_event_time             = m_best_event_time;
	  m_hi_event_time             = m_best_event_time;
	  m_got_first_event_time      = true;
	}
      else
	{
	  if(m_lo_event_time > m_best_event_time)
	    m_lo_event_time = m_best_event_time;
	  if(m_hi_event_time < m_best_event_time)
	    m_hi_event_time = m_best_event_time;
	}

      if(event_num >= m_event_time.size())m_event_time.resize(event_num+1);
      m_event_time[event_num] = EventTime(m_best_event_time);
    }


  if(!m_got_first_event)
    {
      m_first_event_vbf_packet_num    = m_vbf_packet_num;
      m_lo_event_num                  = event_num;
      m_hi_event_num                  = event_num;
      m_got_first_event               = true;
    }
  else
    {
      if(event_num < m_lo_event_num)m_lo_event_num=event_num;
      if(event_num > m_hi_event_num)m_hi_event_num=event_num;
    }
}

void VBFRunInfo::
leaveArrayEvent(bool veto_array_event, void* user_data)
{
  if(m_is_ped_event)
    {
      VBFPacket* packet = m_packet->deepCopy(!m_copy_pedestal_packets);
      m_pedestals.push_back(packet);
    }
}

void VBFRunInfo::
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
  if((m_got_run_number<10)
     &&((m_got_run_number==0)||(m_run_number!=run_number)))
    {
      m_run_number=run_number;
      m_got_run_number=1;
    }
  else m_got_run_number++;
  if(event_type.trigger == VEventType::PED_TRIGGER)m_is_ped_event=true;
  for(unsigned iscope=0; iscope<8; iscope++)
    if(config_mask & (0x0001 << iscope))
      {
	if(m_config_mask.size()<=iscope)m_config_mask.resize(iscope+1);
	m_config_mask[iscope] = true;
      }
}

void VBFRunInfo::
visitArrayTriggerScope(bool& veto_array_event,
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
		       const uint32_t* cal_counts)
{
  if(m_config_mask.size()<=telescope_num)
    m_config_mask.resize(telescope_num+1);
  m_config_mask[telescope_num] = true;
}

void VBFRunInfo::
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
  veto_scope_event = true;

  if(m_nsample.size() <= telescope_num)
    m_nsample.resize(telescope_num+1);
  if(m_nsample[telescope_num]<num_samples)
    m_nsample[telescope_num]=num_samples;

  if(m_nchan.size() <= telescope_num)m_nchan.resize(telescope_num+1);
  if(m_nchan[telescope_num]<num_channels_total)
    m_nchan[telescope_num]=num_channels_total;

  if(event_type.trigger == VEventType::PED_TRIGGER)m_is_ped_event=true;
}

void VBFRunInfo::
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
  m_sim_data = new SimData;
  unsigned nscope = array.size();
  m_sim_data->scope_positions.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      m_sim_data->scope_positions[iscope].first = 
	array[iscope].fRelTelLocEastM;
      m_sim_data->scope_positions[iscope].second = 
	-array[iscope].fRelTelLocSouthM;
      m_sim_data->scope_positions[iscope].third = 
	array[iscope].fRelTelLocUpM;
    }
}

void VBFRunInfo::
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
  vsassert(m_sim_data);

  if(corsika_particle_id == 1)
    {
      SphericalCoords primary_azzn =
	SEphem::SphericalCoords::makeDeg(primary_zenith_deg, 
					 primary_azimuth_deg);
      SphericalCoords obs_azzn =
	SEphem::SphericalCoords::makeDeg(obs_zenith_deg, 
					 obs_azimuth_deg);


      SphericalCoords primary_radec = primary_azzn;
      SphericalCoords obs_radec = obs_azzn;

      m_sim_transform->transform(primary_radec,corsika_particle_id,
				 primary_azzn,obs_azzn);

      m_sim_transform->transform(obs_radec,corsika_particle_id,
				 primary_azzn,obs_azzn);

      m_sim_src_radec += 
	VSAAlgebra::Vec3D::makePolar(primary_radec.thetaRad(),
				     primary_radec.phiRad());
      m_sim_obs_radec += 
	VSAAlgebra::Vec3D::makePolar(obs_radec.thetaRad(),obs_radec.phiRad());

      m_sim_nevent++;
    }

  if(energy_gev>0)
    {
      float etev = energy_gev*0.001;

      m_sim_energy[corsika_particle_id].insert(std::log10(etev));

      if((m_sim_energy_limits[corsika_particle_id].first == 0)
	 ||(etev < m_sim_energy_limits[corsika_particle_id].first))
	m_sim_energy_limits[corsika_particle_id].first = etev;

      if((m_sim_energy_limits[corsika_particle_id].second == 0)
	 ||(etev > m_sim_energy_limits[corsika_particle_id].second))
	m_sim_energy_limits[corsika_particle_id].second = etev;
    }
}

void VBFRunInfo::
visitChiLAHeader(void* user_data,
		 const std::string&    database_name,
		 const std::vector<VChiLASimParamTable>&
		                       sim_param_tables,
		 const std::vector<VChiLAOpticsConfiguration>&
		                       optics_configurations,
		 const std::vector<VChiLAElectronicsConfiguration>& 
		                       electronics_configurations)
{
  vsassert(m_sim_data);
  m_sim_data->package = SimData::SP_CHILA;
}

void VBFRunInfo::
visitChiLAEvent(bool& veto_packet, void* user_data,
		uint32_t               table_index,
		uint32_t               electronics_id,
		uint32_t               table_event_index,
		float                  a_omega,
		float                  a_omega_var)
{
  // nothing t see here
}

void VBFRunInfo::
visitKascadeHeader(void* user_data,
		   uint32_t            corsika_particle_id,
		   float               energy_gev,
		   uint32_t            shower_id,
		   float               x_area_width_m,
		   float               y_area_width_m,
		   uint32_t            north_south_grid)
{
  vsassert(m_sim_data);
  m_sim_data->package = SimData::SP_KASCADE;
}
 
void VBFRunInfo::
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
  // nothing t see here
}

void VBFRunInfo::getData(Data& data) const
{
  data.got_run_number              = (m_got_run_number>=10);
  data.run_number                  = m_run_number;
  data.nevents                     = m_nevents;
  data.first_event_time            = m_first_event_time;
  data.first_event_vbf_packet_num  = m_first_event_vbf_packet_num;
  data.lo_event_num                = m_lo_event_num;
  data.hi_event_num                = m_hi_event_num;
  data.lo_event_time               = m_lo_event_time;
  data.hi_event_time               = m_hi_event_time;
  data.config_mask                 = m_config_mask;
  data.nsample                     = m_nsample;
  data.nchan                       = m_nchan;
  if(data.nsample.size() < data.config_mask.size())
    data.nsample.resize(data.config_mask.size());
  if(data.nchan.size() < data.config_mask.size())
    data.nchan.resize(data.config_mask.size());
}

VBFRunInfo::SimData* VBFRunInfo::getSimData() const
{
  if(!m_sim_data)return 0;

  SimData* sd(new SimData(*m_sim_data));

  // Try to estimate the binning of energy. If the spectrum is
  // continuous or not binned in log energy then this returns a zero.
  for(std::map<unsigned,std::pair<float,float> >::const_iterator ipart
	= m_sim_energy_limits.begin(); ipart != m_sim_energy_limits.end();
      ipart++)
    {
      unsigned cpid = ipart->first;
      std::map<unsigned,std::set<float> >::const_iterator 
	eset(m_sim_energy.find(cpid));
      assert(eset != m_sim_energy.end());

      float dloge = 0.0;
      float bin0 = 0.0;
      if(eset->second.size()>1)
	{
	  std::set<float>::const_iterator ienergy = eset->second.begin();
	  float loge = *ienergy++;

	  std::vector<float> all_dloge;
	  while(ienergy != eset->second.end())
	    {
	      all_dloge.push_back(*ienergy - loge);
	      loge = *ienergy++;
	    }
	  dloge = median(all_dloge);

	  for(std::vector<float>::const_iterator idloge = all_dloge.begin();
	      idloge != all_dloge.end(); idloge++)
	    {
	      float x = *idloge/dloge;
	      if(std::abs(x-roundf(x))>0.025)
		{
		  dloge = 0;
		  break;
		}
	    }

	  float logeref = std::log10(ipart->second.first);
	  all_dloge.clear();
	  while(ienergy != eset->second.end())
	    {
	      float x = (*ienergy - logeref)/dloge;
	      all_dloge.push_back(std::abs(x-roundf(x)));
	    }
	  bin0 = median(all_dloge)*dloge + logeref;
	}

      if(cpid == 1 && m_sim_data->package == SimData::SP_CHILA && dloge == 0)
	dloge = 0.0625;

      SimData::ParticleSpectrum ps;
      ps.corsika_particle_id    = cpid;
      ps.min_energy_tev         = ipart->second.first;
      ps.max_energy_tev         = ipart->second.second;
      ps.dlog_energy            = dloge;
      ps.discrete_binning_start = bin0;
      sd->particle_spectrum.push_back(ps);
    }

  return sd;
}

double VBFRunInfo::getSliceWidthForNSlice(unsigned nslice) const
{
  int64_t dt = m_hi_event_time-m_lo_event_time+INT64_C(1000000000);
  return double(dt)/1e9/double(nslice);
}

double VBFRunInfo::getSliceWidthForApproxWidth(double approx_slice_width) const
{
  int64_t dt = (m_hi_event_time-m_lo_event_time)+INT64_C(1000000000);
  unsigned nslice = unsigned(round(double(dt)/1e9/approx_slice_width));
  if(nslice==0)nslice=1;
  return double(dt)/1e9/double(nslice);
}

std::vector<VBFRunInfo::Slice> 
VBFRunInfo::getTimeSlices(double slice_width) const
{
  double dt = double(m_hi_event_time-m_lo_event_time)/1e9;
  unsigned nslice = slice_width>0 ? unsigned(ceil(dt/slice_width)) : 0;

  std::vector<VBFRunInfo::Slice> slice(nslice);
  if(nslice==0)
    {
      if(m_hi_event_num)
	{
	  slice.resize(1);
	  slice[0].event_time_lo = VSTime::perpetual_past();
	  slice[0].event_time_hi = VSTime::perpetual_future();
	  slice[0].event_num_lo = 0;
	  slice[0].event_num_hi = m_hi_event_num+1;
	}
      return slice;
    }

  for(unsigned islice=0;islice<nslice;islice++)
    {
      slice[islice].event_time_lo = 
	m_lo_event_time+int64_t(slice_width*1e9*islice);
      slice[islice].event_time_hi = 
	m_lo_event_time+int64_t(slice_width*1e9*(islice+1));
      slice[islice].event_num_lo = 0;
      slice[islice].event_num_hi = 0;
    }

  unsigned islice = 0;
  unsigned nevent=m_hi_event_num+1;//m_event_time.size();
  for(unsigned ievent=0;ievent<nevent;ievent++)
    if((ievent<m_event_time.size())
       &&(m_event_time[ievent].event_found)
       &&(m_event_time[ievent].event_time.isGood())
       &&(m_event_time[ievent].event_time>=m_lo_event_time)
       &&(m_event_time[ievent].event_time<=m_hi_event_time))
      {
	while((m_event_time[ievent].event_time>=slice[islice].event_time_hi)
	      &&(++islice<nslice))
	  slice[islice].event_num_lo = ievent;
	if(islice==nslice)break;
      }

  for(unsigned islice=0;islice<nslice-1;islice++)
    slice[islice].event_num_hi = slice[islice+1].event_num_lo;
  slice[nslice-1].event_num_hi = m_hi_event_num+1;

  return slice;
}

std::vector<VBFRunInfo::Slice> 
VBFRunInfo::getEventSlices(unsigned num_events) const
{
  unsigned nevent = m_hi_event_num+1;
  unsigned nslice = (nevent+num_events-1)/num_events;

  std::vector<VBFRunInfo::Slice> slice(nslice);
  if(nslice==0)return slice;

  for(unsigned islice=0;islice<nslice;islice++)
    {
      slice[islice].event_num_lo = islice*num_events;
      slice[islice].event_num_hi = (islice+1)*num_events;
      slice[islice].event_time_lo = m_lo_event_time;
      slice[islice].event_time_hi = VSTime();
    }
  
  VSTime last_start_time(VSTime::perpetual_past());
  unsigned islice = 0;
  for(unsigned ievent=0;ievent<nevent;ievent++)
    if((ievent<m_event_time.size())
       &&(m_event_time[ievent].event_found)
       &&(m_event_time[ievent].event_time.isGood())
       &&(m_event_time[ievent].event_time>=m_lo_event_time)
       &&(m_event_time[ievent].event_time<=m_hi_event_time))       
      {
	while((ievent>=slice[islice].event_num_hi)
	      &&(++islice<nslice))
	  {
	    last_start_time = slice[islice-1].event_time_lo;
	    slice[islice].event_time_lo = m_event_time[ievent].event_time;
	  }
	if(islice==nslice)break;
	if((m_event_time[ievent].event_time<slice[islice].event_time_lo)
	   &&(m_event_time[ievent].event_time>last_start_time))
	  slice[islice].event_time_lo = m_event_time[ievent].event_time;
      }
  
  for(unsigned islice=0;islice<nslice-1;islice++)
    slice[islice].event_time_hi = slice[islice+1].event_time_lo;
  slice[nslice-1].event_time_hi = m_hi_event_time;
  slice[nslice-1].event_time_hi += INT64_C(1000000000);

  return slice;
}

const std::vector<VBFRunInfo::EventTime>& VBFRunInfo::getEventTimes() const
{
  return m_event_time;
}

unsigned VBFRunInfo::dispatchPedestals(VSSimpleVBFDispatcher& dispatcher) const
{
  unsigned ipacket = 0;
  for(std::list<VBFPacket*>::const_iterator iped=m_pedestals.begin();
      iped!=m_pedestals.end(); iped++)
    {
      if((*iped)->vbf)
	dispatcher.dispatchPacket(*iped, ipacket);
      else
	dispatcher.dispatchPacket((*iped)->vbf_num, ipacket);
      ipacket++;
    }
  return ipacket;
}

void VBFRunInfo::Slice::
load(VSOctaveH5ReaderStruct* reader, std::vector<Slice>& slice)
{
  std::vector<unsigned> _event_num_lo;
  std::vector<unsigned> _event_num_hi;
  std::vector<uint32_t> _event_mjd_lo;
  std::vector<uint64_t> _event_dns_lo;
  std::vector<uint32_t> _event_mjd_hi;
  std::vector<uint64_t> _event_dns_hi;
  reader->readVector("event_num_lo",_event_num_lo);
  reader->readVector("event_num_hi",_event_num_hi);
  reader->readVector("event_mjd_lo",_event_mjd_lo);
  reader->readVector("event_dns_lo",_event_dns_lo);
  reader->readVector("event_mjd_hi",_event_mjd_hi);
  reader->readVector("event_dns_hi",_event_dns_hi);
  unsigned nslice=_event_num_lo.size();
  slice.resize(nslice);
  for(unsigned islice=0;islice<nslice;islice++)
    {
      slice[islice].event_num_lo = _event_num_lo[islice];
      slice[islice].event_num_hi = _event_num_hi[islice];
      slice[islice].event_time_lo.setFromMJDIntAndNS(_event_mjd_lo[islice],
						     _event_dns_lo[islice]);
      slice[islice].event_time_hi.setFromMJDIntAndNS(_event_mjd_hi[islice],
						     _event_dns_hi[islice]);
    }
}

void VBFRunInfo::Slice::
save(VSOctaveH5WriterStruct* writer, const std::vector<Slice>& slice)
{
  unsigned nslice=slice.size();
  std::vector<unsigned> _event_num_lo(nslice);
  std::vector<unsigned> _event_num_hi(nslice);
  std::vector<uint32_t> _event_mjd_lo(nslice);
  std::vector<uint64_t> _event_dns_lo(nslice);
  std::vector<uint32_t> _event_mjd_hi(nslice);
  std::vector<uint64_t> _event_dns_hi(nslice);
  for(unsigned islice=0;islice<nslice;islice++)
    {
      _event_num_lo[islice] = slice[islice].event_num_lo;
      _event_num_hi[islice] = slice[islice].event_num_hi;
      _event_mjd_lo[islice] = slice[islice].event_time_lo.getMJDInt();
      _event_dns_lo[islice] = slice[islice].event_time_lo.getDayNS();
      _event_mjd_hi[islice] = slice[islice].event_time_hi.getMJDInt();
      _event_dns_hi[islice] = slice[islice].event_time_hi.getDayNS();
    }
  writer->writeVector("event_num_lo",_event_num_lo);
  writer->writeVector("event_num_hi",_event_num_hi);
  writer->writeVector("event_mjd_lo",_event_mjd_lo);
  writer->writeVector("event_dns_lo",_event_dns_lo);
  writer->writeVector("event_mjd_hi",_event_mjd_hi);
  writer->writeVector("event_dns_hi",_event_dns_hi);
}

void VBFRunInfo::Data::clear()
{
  got_run_number             = false;
  run_number                 = 0;
  nevents                    = 0;
  first_event_time           = VSTime();
  first_event_vbf_packet_num = 0;
  lo_event_num               = 0;
  hi_event_num               = 0;
  lo_event_time              = VSTime();
  hi_event_time              = VSTime();
  ra_mean_deg                = 0;
  dec_mean_deg               = 0;
  ra_rms_deg                 = 0;
  dec_rms_deg                = 0;
  zn_mean_deg                = 0;
  az_mean_deg                = 0;
  zn_rms_deg                 = 0;
  az_rms_deg                 = 0;
  config_mask.clear();
  nsample.clear();
  nchan.clear();
}

void VBFRunInfo::Data::load(VSOctaveH5ReaderStruct* reader)
{
  clear();
  reader->readCompositeHere(*this);
  reader->readVector("nsample",nsample);
  reader->readVector("nchan",nchan);
  if(reader->isValid("config_mask"))
    reader->readVector("config_mask",config_mask);
  else
    config_mask.resize(nchan.size());
}

void VBFRunInfo::Data::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);
  writer->writeVector("config_mask",config_mask);
  writer->writeVector("nsample",nsample);
  writer->writeVector("nchan",nchan);
}

void VBFRunInfo::Data::
calculateMeanPointing(VSPointing& pointing,
		      const SEphem::SphericalCoords& earth_position)
{
  SphericalCoords mean_radec;
  double rms_ra_rad = 0;
  double rms_dec_rad = 0;
  pointing.getMeanTarget(mean_radec, rms_ra_rad, rms_dec_rad,
			 lo_event_time, hi_event_time, earth_position);
      
  SphericalCoords mean_azzn;
  double rms_az_rad = 0;
  double rms_zn_rad = 0;
  double mean_zn = 0;
  pointing.getMeanAzZn(mean_azzn, rms_az_rad, rms_zn_rad, mean_zn,
		       lo_event_time, hi_event_time);

  ra_mean_deg      = mean_radec.longitudeDeg();
  dec_mean_deg     = mean_radec.latitudeDeg();
  ra_rms_deg       = Angle::toDeg(rms_ra_rad);
  dec_rms_deg      = Angle::toDeg(rms_dec_rad);
  zn_mean_deg      = mean_azzn.thetaDeg();
  az_mean_deg      = mean_azzn.phiDeg();
  zn_rms_deg       = Angle::toDeg(rms_zn_rad);
  az_rms_deg       = Angle::toDeg(rms_az_rad);
  mean_zn_1dim_deg = Angle::toDeg(mean_zn);
}

void VBFRunInfo::SimData::clear()
{
  src_ra_deg       = 0;
  src_dec_deg      = 0;
  obs_ra_deg       = 0;
  obs_dec_deg      = 0;
  wobble_theta_deg = 0;
  wobble_phi_deg   = 0;
}

void VBFRunInfo::SimData::load(VSOctaveH5ReaderStruct* reader)
{
  reader->readCompositeHere(*this);
  if(reader->isStruct("particle_spectrum"))
    reader->readCompositeVector("particle_spectrum",particle_spectrum);
  if(reader->isStruct("scope_positions"))
    reader->readCompositeVector("scope_positions",scope_positions);
}

void VBFRunInfo::SimData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);
  writer->writeCompositeVector("particle_spectrum",particle_spectrum);
  writer->writeCompositeVector("scope_positions",scope_positions);
}
