//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimulationData.hpp

  Data structures for array and scope simulation data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSSimulationData.hpp,v 3.13 2008/09/24 20:20:21 matthew Exp $

*/

#ifndef VSSIMULATIONDATA_HPP
#define VSSIMULATIONDATA_HPP

#include<vector>

#include<VSTime.hpp>
#include<VSOctaveIO.hpp>

namespace VERITAS
{
  class VSTableSimulationDatum
  {
  public:
    VSTableSimulationDatum():
      table_index(), table_db_id(), primary_id(), energy_tev(), 
      zenith_min_deg(), zenith_max_deg(), azimuth_min_deg(), azimuth_max_deg(),
      optics_id(), sampling_radius_m(), event_count(), table_name(),
      num_events_read(), num_events_processed(), num_events_written()
    { /* nothing to see here */ }

    ~VSTableSimulationDatum() { /* nothing to see here */ }

    uint32_t      table_index;
    uint32_t      table_db_id;
    uint32_t      primary_id;
    double        energy_tev;
    double        zenith_min_deg;
    double        zenith_max_deg;
    double        azimuth_min_deg;
    double        azimuth_max_deg;
    uint32_t      optics_id;
    double        sampling_radius_m;
    uint32_t      event_count;
    std::string   table_name;
    uint32_t      num_events_read;
    uint32_t      num_events_processed;
    uint32_t      num_events_written;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSTableSimulationDatum,table_index);
      H5_ADDMEMBER(c,VSTableSimulationDatum,table_db_id);
      H5_ADDMEMBER(c,VSTableSimulationDatum,primary_id);
      H5_ADDMEMBER(c,VSTableSimulationDatum,energy_tev);
      H5_ADDMEMBER(c,VSTableSimulationDatum,zenith_min_deg);
      H5_ADDMEMBER(c,VSTableSimulationDatum,zenith_max_deg);
      H5_ADDMEMBER(c,VSTableSimulationDatum,azimuth_min_deg);
      H5_ADDMEMBER(c,VSTableSimulationDatum,azimuth_max_deg);
      H5_ADDMEMBER(c,VSTableSimulationDatum,optics_id);
      H5_ADDMEMBER(c,VSTableSimulationDatum,sampling_radius_m);
      H5_ADDMEMBER(c,VSTableSimulationDatum,event_count);
      H5_ADDMEMBER(c,VSTableSimulationDatum,table_name);
      H5_ADDMEMBER(c,VSTableSimulationDatum,num_events_read);
      H5_ADDMEMBER(c,VSTableSimulationDatum,num_events_processed);
      H5_ADDMEMBER(c,VSTableSimulationDatum,num_events_written);
    }
  };

  class VSHeaderSimulationDatum
  {
  public:
    VSHeaderSimulationDatum():
      run_number(), date_of_sims(), simulation_package(),
      simulator_name(), date_of_array(), corsika_atm_model(), obs_altitude_m(),
      sim_config(), database_name(), tables() 
    { /* nothing to see here */ }
    ~VSHeaderSimulationDatum() {}

    // Generic Simulation Parameters ------------------------------------------
    uint32_t    run_number;
    VSTime      date_of_sims;
    uint32_t    simulation_package;
    std::string simulator_name;
    VSTime      date_of_array;
    uint32_t    corsika_atm_model;
    double      obs_altitude_m;
    std::string sim_config;

    // ChiLA Parameters -------------------------------------------------------
    std::string database_name;
    std::vector< VSTableSimulationDatum > tables;

    void clear();
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSHeaderSimulationDatum,run_number);
      H5_ADDSIMPLECOMPOSITE(c,VSHeaderSimulationDatum,date_of_sims);
      H5_ADDMEMBER(c,VSHeaderSimulationDatum,simulation_package);
      H5_ADDMEMBER(c,VSHeaderSimulationDatum,simulator_name);
      H5_ADDSIMPLECOMPOSITE(c,VSHeaderSimulationDatum,date_of_array);
      H5_ADDMEMBER(c,VSHeaderSimulationDatum,corsika_atm_model);
      H5_ADDMEMBER(c,VSHeaderSimulationDatum,obs_altitude_m);
      H5_ADDMEMBER(c,VSHeaderSimulationDatum,sim_config);
      H5_ADDMEMBER(c,VSHeaderSimulationDatum,database_name);
    }
  };

#if 0
  class VSScopeSimulationDatum
  {
  public:
    VSScopeSimulationDatum() {}
    ~VSScopeSimulationDatum() {}

    static void _compose(VSOctaveH5CompositeDefinition& c)
    { }
  };
#endif

  class VSArraySimulationDatum
  {
  public:
    VSArraySimulationDatum():
      event_num(),
      corsika_particle_id(), energy_tev(), obs_zenith_deg(), obs_azimuth_deg(),
      primary_zenith_deg(), primary_azimuth_deg(), 
      ref_zenith_deg(), ref_azimuth_deg(), ref_position_angle_deg(),
      primary_dec_deg(), primary_ra_deg(), 
      core_east_m(), core_north_m(), core_elevation_asl_m(),
      core_x_m(), core_y_m(), core_R_m(),
      table_index(), table_event_index(), electronics_id()
#if 0
      ,scope()
#endif
    {}

    ~VSArraySimulationDatum()
    {
#if 0
      for(unsigned iscope=0;iscope<scope.size();iscope++)delete scope[iscope];
#endif
    }

    // Generic Simulation Parameters ------------------------------------------
    uint32_t event_num;
    uint32_t corsika_particle_id;
    double   energy_tev;
    double   obs_zenith_deg;
    double   obs_azimuth_deg;
    double   primary_zenith_deg;
    double   primary_azimuth_deg;
    double   ref_zenith_deg;
    double   ref_azimuth_deg;
    double   ref_position_angle_deg;
    double   primary_dec_deg;
    double   primary_ra_deg;
    double   core_east_m;
    double   core_north_m;
    double   core_elevation_asl_m;
    double   core_x_m;
    double   core_y_m;
    double   core_R_m;
    bool     has_array_event;

    // ChiLA Parameters -------------------------------------------------------
    uint32_t table_index;
    uint32_t table_event_index;
    uint32_t electronics_id;

#if 0
    std::vector<VSScopeSimulationDatum*> scope;
#endif

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSArraySimulationDatum,event_num);
      H5_ADDMEMBER(c,VSArraySimulationDatum,corsika_particle_id);
      H5_ADDMEMBER(c,VSArraySimulationDatum,energy_tev);
      H5_ADDMEMBER(c,VSArraySimulationDatum,obs_zenith_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,obs_azimuth_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,primary_zenith_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,primary_azimuth_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,ref_zenith_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,ref_azimuth_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,ref_position_angle_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,primary_dec_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,primary_ra_deg);
      H5_ADDMEMBER(c,VSArraySimulationDatum,core_east_m);
      H5_ADDMEMBER(c,VSArraySimulationDatum,core_north_m);
      H5_ADDMEMBER(c,VSArraySimulationDatum,core_elevation_asl_m);
      H5_ADDMEMBER(c,VSArraySimulationDatum,core_x_m);
      H5_ADDMEMBER(c,VSArraySimulationDatum,core_y_m);
      H5_ADDMEMBER(c,VSArraySimulationDatum,core_R_m);
      H5_ADDMEMBER(c,VSArraySimulationDatum,has_array_event);
      H5_ADDMEMBER(c,VSArraySimulationDatum,table_index);
      H5_ADDMEMBER(c,VSArraySimulationDatum,table_event_index);
      H5_ADDMEMBER(c,VSArraySimulationDatum,electronics_id);
    }
  };

  class VSArraySimulationWriter
  {
  public:
    VSArraySimulationWriter(VSOctaveH5WriterStruct* s, unsigned ntel);
    ~VSArraySimulationWriter();

    bool append(const VSArraySimulationDatum& x);

  private:
    VSArraySimulationWriter(const VSArraySimulationWriter&);
    VSArraySimulationWriter& operator=(const VSArraySimulationWriter&);

    typedef 
    VSOctaveH5WriterCompositeVector<VSArraySimulationDatum> ArrayDataWriter;

#if 0
    typedef 
    VSOctaveH5WriterCompositeVector<VSScopeSimulationDatum> ScopeDataWriter;
#endif

    ArrayDataWriter*                              m_array_writer;
#if 0
    std::vector<ScopeDataWriter*>                 m_scope_writer;
#endif
  };

  class VSArraySimulationReader
  {
  public:
    typedef VSOctaveH5ReaderBase::MemberSubset MemberSubset;

    VSArraySimulationReader(VSOctaveH5ReaderStruct* s,
			    const MemberSubset& array_subset = MemberSubset());
    ~VSArraySimulationReader();

    bool element(VSArraySimulationDatum& x, unsigned index);

    inline unsigned rows() const { return m_array_reader->rows(); }

    VSArraySimulationDatum at(unsigned index) 
    { 
      VSArraySimulationDatum x; 
      if(!element(x,index))throw std::out_of_range(__PRETTY_FUNCTION__); 
      return x; 
    }

    VSArraySimulationDatum operator[] (unsigned index) 
    {
      VSArraySimulationDatum x; 
      element(x,index); 
      return x;
    }

    static bool loadAllEvents(VSOctaveH5ReaderStruct* s, 
			      std::vector<VSArraySimulationDatum>& x);

  private:
    VSArraySimulationReader(const VSArraySimulationReader&);
    VSArraySimulationReader& operator=(const VSArraySimulationReader&);

    typedef VSOctaveH5ReaderCompositeVector<VSArraySimulationDatum> 
    ArrayDataReader;
#if 0
    typedef VSOctaveH5ReaderCompositeVector<VSScopeSimulationDatum> 
    ScopeDataReader;
#endif

    ArrayDataReader*                              m_array_reader;
#if 0
    std::vector<ScopeDataReader*>                 m_scope_reader;
#endif
  };
}

#endif // VSSIMULATIONDATA_HPP
