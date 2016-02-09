//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFRunInfo.hpp

  Extract some information from run

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VBFRunInfo.hpp,v 3.19 2008/09/24 20:03:42 matthew Exp $

*/

#ifndef VBFRUNINFO_HPP
#define VBFRUNINFO_HPP

#include <vector>

#include <SphericalCoords.h>
#include <VSSimpleVBF.hpp>
#include <VSOctaveIO.hpp>
#include <VSDataReductionTypes.hpp>
#include <VSSimCoordTransform.hpp>
#include <VSSimpleStat.hpp>
#include <VSAAlgebra.hpp>

namespace VERITAS
{

  class VSPointing;

  class VBFRunInfo: public VSSimpleVBFVisitor
  {
  public:
    VBFRunInfo(bool copy_pedestal_packets,
	       double sim_wobble_phi_rad);
    virtual ~VBFRunInfo();

    virtual void visitFile(const char* filename, unsigned npacket);

    virtual void leaveFile();

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

    class Data
    {
    public:
      Data():
        got_run_number(), run_number(), nevents(), first_event_time(),
        first_event_vbf_packet_num(), lo_event_num(), hi_event_num(),
        lo_event_time(), hi_event_time(),
	ra_mean_deg(), dec_mean_deg(), ra_rms_deg(), dec_rms_deg(),
	zn_mean_deg(), az_mean_deg(), zn_rms_deg(), az_rms_deg(),
	mean_zn_1dim_deg(),
	config_mask(), nsample(), nchan()
      { /* nothing to see here */ }

      bool                    got_run_number;
      unsigned                run_number;
      unsigned                nevents;
      VSTime                  first_event_time;
      unsigned                first_event_vbf_packet_num;
      unsigned                lo_event_num;
      unsigned                hi_event_num;
      VSTime                  lo_event_time;
      VSTime                  hi_event_time;
      double                  ra_mean_deg;
      double                  dec_mean_deg;
      double                  ra_rms_deg;
      double                  dec_rms_deg;
      double                  zn_mean_deg;
      double                  az_mean_deg;
      double                  zn_rms_deg;
      double                  az_rms_deg;
      double                  mean_zn_1dim_deg;
      std::vector<bool>       config_mask;
      std::vector<unsigned>   nsample;
      std::vector<unsigned>   nchan;

      static void _compose(VSOctaveH5CompositeDefinition& c) 
      {
	H5_ADDMEMBER(c,Data,got_run_number);
	H5_ADDMEMBER(c,Data,run_number);
	H5_ADDMEMBER(c,Data,nevents);
	H5_ADDSIMPLECOMPOSITE(c,Data,first_event_time);
	H5_ADDMEMBER(c,Data,first_event_vbf_packet_num);
	H5_ADDMEMBER(c,Data,lo_event_num);
	H5_ADDMEMBER(c,Data,hi_event_num);
	H5_ADDSIMPLECOMPOSITE(c,Data,lo_event_time);
	H5_ADDSIMPLECOMPOSITE(c,Data,hi_event_time);
	H5_ADDMEMBER(c,Data,ra_mean_deg);
	H5_ADDMEMBER(c,Data,dec_mean_deg);
	H5_ADDMEMBER(c,Data,ra_rms_deg);
	H5_ADDMEMBER(c,Data,dec_rms_deg);
	H5_ADDMEMBER(c,Data,zn_mean_deg);
	H5_ADDMEMBER(c,Data,az_mean_deg);
	H5_ADDMEMBER(c,Data,zn_rms_deg);
	H5_ADDMEMBER(c,Data,az_rms_deg);
	H5_ADDMEMBER(c,Data,mean_zn_1dim_deg);
      }

      void clear();
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;

      void calculateMeanPointing(VSPointing& pointing,
				const SEphem::SphericalCoords& earth_position);
    
    };

    class SimData
    {
    public:
      SimData():
        package(SP_UNKNOWN), 
	src_ra_deg(), src_dec_deg(), obs_ra_deg(), obs_dec_deg(),
        wobble_theta_deg(), wobble_phi_deg(), 
	particle_spectrum(), scope_positions()
      { /* nothing to see here */ }

      enum SimPackage { SP_UNKNOWN, SP_CHILA, SP_KASCADE, SP_GRISU,
			SP_MCGILL };

      class ParticleSpectrum
      {
      public:
	ParticleSpectrum(): 
	  corsika_particle_id(), min_energy_tev(), max_energy_tev(),
	  dlog_energy(), discrete_binning_start() { /* nothing to see here */ }

	uint32_t              corsika_particle_id;
	double                min_energy_tev;
	double                max_energy_tev;      
	double                dlog_energy;
	double                discrete_binning_start;

	static void _compose(VSOctaveH5CompositeDefinition& c)
	{
	  H5_ADDMEMBER(c,SimData::ParticleSpectrum,corsika_particle_id);
	  H5_ADDMEMBER(c,SimData::ParticleSpectrum,min_energy_tev);
	  H5_ADDMEMBER(c,SimData::ParticleSpectrum,max_energy_tev);
	  H5_ADDMEMBER(c,SimData::ParticleSpectrum,dlog_energy);
	  H5_ADDMEMBER(c,SimData::ParticleSpectrum,discrete_binning_start);
	}
      };

      SimPackage              package;
      double                  src_ra_deg;
      double                  src_dec_deg;
      double                  obs_ra_deg;
      double                  obs_dec_deg;
      double                  wobble_theta_deg;
      double                  wobble_phi_deg;
      
      std::vector<ParticleSpectrum> particle_spectrum;
      std::vector<pos_type>         scope_positions;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,SimData,package);
	H5_ADDMEMBER(c,SimData,src_ra_deg);
	H5_ADDMEMBER(c,SimData,src_dec_deg);
	H5_ADDMEMBER(c,SimData,obs_ra_deg);
	H5_ADDMEMBER(c,SimData,obs_dec_deg);
	H5_ADDMEMBER(c,SimData,wobble_theta_deg);
	H5_ADDMEMBER(c,SimData,wobble_phi_deg);
      }

      void clear();
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;
    };
  
    struct EventTime
    {
      EventTime(): event_time(), event_found(false) { }
      EventTime(const VSTime& _time): event_time(_time), event_found(true) { }
      VSTime    event_time;
      bool      event_found;
    };

    struct Slice
    {
      Slice(): event_num_lo(), event_num_hi(), 
	       event_time_lo(), event_time_hi() { }
      unsigned event_num_lo;
      unsigned event_num_hi;
      VSTime event_time_lo;
      VSTime event_time_hi;

      static void load(VSOctaveH5ReaderStruct* reader, 
		       std::vector<Slice>& slice);
      static void save(VSOctaveH5WriterStruct* writer, 
		       const std::vector<Slice>& slice);
    };

    void getData(Data& data) const;
    SimData* getSimData() const;
    double getSliceWidthForNSlice(unsigned nslice) const;
    double getSliceWidthForApproxWidth(double approx_slice_width) const;
    std::vector<Slice> getTimeSlices(double slice_width) const;
    std::vector<Slice> getEventSlices(unsigned num_events) const;
    const std::vector<EventTime>& getEventTimes() const;

    unsigned dispatchPedestals(VSSimpleVBFDispatcher& dispatcher) const;

  private:
    VBFRunInfo(const VBFRunInfo&);
    VBFRunInfo& operator= (const VBFRunInfo&);

    // SETTINGS ---------------------------------------------------------------
    bool                    m_copy_pedestal_packets;
    double                  m_sim_wobble_phi_rad;

    // DATA -------------------------------------------------------------------
    unsigned                m_run_number;
    unsigned                m_nevents;
    VSTime                  m_first_event_time;
    unsigned                m_first_event_vbf_packet_num;
    unsigned                m_lo_event_num;
    unsigned                m_hi_event_num;
    VSTime                  m_lo_event_time;
    VSTime                  m_hi_event_time;
    std::vector<bool>       m_config_mask;
    std::vector<unsigned>   m_nsample;
    std::vector<unsigned>   m_nchan;
    std::vector<EventTime>  m_event_time;
    std::list<VBFPacket*>   m_pedestals;

    std::map<unsigned,std::pair<float,float> > m_sim_energy_limits;
    std::map<unsigned,std::set<float> > m_sim_energy;
    SimData*                m_sim_data;

    // INTERMEDIATE -----------------------------------------------------------
    unsigned                m_got_run_number;
    bool                    m_got_first_event;
    bool                    m_got_first_event_time;
    unsigned                m_vbf_packet_num;
    bool                    m_has_good_event_time;
    VSAAlgebra::Vec3D       m_sim_src_radec;
    VSAAlgebra::Vec3D       m_sim_obs_radec;
    unsigned                m_sim_nevent;
    VSSimCoordTransform*    m_sim_transform;
    VSTime                  m_best_event_time;
    unsigned                m_event_num;
    bool                    m_is_ped_event;
    const VBFPacket*        m_packet;
  };
 
  DEFINE_H5ENUM(VBFRunInfo::SimData::SimPackage);

} // namespace VERITAS

#endif // not defined VBFRUNINFO_HPP
