//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleVBF.hpp
  Simple VBF Dispatcher and Vistor

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/07/2005

  $Id: VSSimpleVBF.hpp,v 3.9 2008/03/05 22:38:10 sfegan Exp $

*/

#ifndef VSSIMPLEVBF_HPP
#define VSSIMPLEVBF_HPP

#include <vsassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <list>

#include <signal.h>

#ifndef NOTHREADS
#include <pthread.h>
#endif

// Filip says to include...
#include <VBF/VBankFileReader.h>
#include <VBF/VPacket.h>
#include <VBF/VArrayEvent.h>
#include <VBF/VDatum.h>

#include "VSTime.hpp"

namespace VERITAS
{

  class VSSimpleVBFDispatcherStop
  {
  public:
    VSSimpleVBFDispatcherStop() { /* nothing to see here */ }
    virtual ~VSSimpleVBFDispatcherStop();
    virtual void stopProcessingFile() = 0;
  };

  // ==========================================================================
  // VBFPacket - Used internally to carry extra information about packet
  // ==========================================================================

  class VBFPacket
  {
  public:
    VBFPacket(): 
      vbf_num(), vbf(), has_best_event_time(), best_event_time() 
    { /* nothing to see here */ }
    ~VBFPacket() { delete vbf; }

    VBFPacket* deepCopy(bool omit_vbf_packet = false) const;

    unsigned   vbf_num;
    VPacket*   vbf;
    bool       has_best_event_time;
    VSTime     best_event_time;
      
  private:
    VBFPacket(const VBFPacket& o);
    VBFPacket& operator= (const VBFPacket& o);
  };

  // ==========================================================================
  // VSSimpleVBFVisitor
  // ==========================================================================

  // Class which receives VBF data. Inherited classes should implement
  // some useful functionality. The dispatcher/visitor pattern makes
  // for an event-loop-like driven program.

  class VSSimpleVBFVisitor
  {
  public:
    VSSimpleVBFVisitor(): fDispatchStop(0) { }
    virtual ~VSSimpleVBFVisitor();

    // Registration
    virtual void registerDispatcher(VSSimpleVBFDispatcherStop* dispatch_stop);

    // File
    virtual void visitFile(const char* filename, unsigned npacket);
#ifndef NOTHREADS
    virtual void usingThreads(unsigned nthreads);
#endif
    virtual void leaveFile();

    // Packet
    virtual void visitPacket(bool& veto_packet, void*& user_data,
			     uint32_t                   seq_packet_number,
			     uint32_t                   vbf_packet_number,
			     bool                       has_array_event,
			     bool                       has_sim_header,
			     bool                       has_sim_event,
			     uint32_t                   num_overflow_datum,
			     const VBFPacket*           packet);

    virtual void leavePacket(bool veto_packet, void* user_data);

    // Array Event

    enum EventType { ET_UNKNOWN, ET_PED, ET_L2 };

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
				   uint32_t             num_array_telescopes,
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

    // Telescope Event
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
#if 0
#warning What do I do with clock trigger data ?
#endif
				 const VEvent*          event);

    virtual void leaveScopeEvent(bool veto_scope_event, void* user_data);

    // Channel
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

    virtual void leaveChannel(bool veto_channel, void* user_data);

    virtual void visitOverflowEvent(void* user_data,
				    unsigned            datum_num,
				    const VDatum*       datum);

    // Simulations
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

  protected:
    VSSimpleVBFDispatcherStop* fDispatchStop;

  private:
    VSSimpleVBFVisitor(const VSSimpleVBFVisitor&);
    VSSimpleVBFVisitor& operator= (const VSSimpleVBFVisitor&);
  };

  // ==========================================================================
  // VBFBroadcastVisitor
  // ==========================================================================

  // VBF visitor class which distributes events to other visitors

  class VBFBroadcastVisitor: public VSSimpleVBFVisitor
  {
  private:
    class Visitor;

  public:
    VBFBroadcastVisitor();
    virtual ~VBFBroadcastVisitor();

    void addVisitor(VSSimpleVBFVisitor* visitor);
    void delVisitor(VSSimpleVBFVisitor* visitor);
    void stopDispatcherForVisitor(Visitor* visitor);

    // Registration
    virtual void registerDispatcher(VSSimpleVBFDispatcherStop* dispatch_stop);

    // File
    virtual void visitFile(const char* filename, unsigned npacket);
#ifndef NOTHREADS
    virtual void usingThreads(unsigned nthreads);
#endif
    virtual void leaveFile();

    // Packet
    virtual void visitPacket(bool& veto_packet, void*& user_data,
			     uint32_t                   seq_packet_number,
			     uint32_t                   vbf_packet_number,
			     bool                       has_array_event,
			     bool                       has_sim_header,
			     bool                       has_sim_event,
			     uint32_t                   num_overflow_datum,
			     const VBFPacket*           packet);

    virtual void leavePacket(bool veto_packet, void* user_data);

    // Array Event
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
				   uint32_t             num_array_telescopes,
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

    // Telescope Event
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
#if 0
#warning What do I do with clock trigger data ?
#endif
				 const VEvent*          event);

    virtual void leaveScopeEvent(bool veto_scope_event, void* user_data);

    // Channel
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

    virtual void leaveChannel(bool veto_channel, void* user_data);

    virtual void visitOverflowEvent(void* user_data,
				    unsigned            datum_num,
				    const VDatum*       datum);

    // Simulations
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
    VBFBroadcastVisitor(const VBFBroadcastVisitor&);
    VBFBroadcastVisitor& operator= (const VBFBroadcastVisitor&);

    class Visitor: public VSSimpleVBFDispatcherStop
    {
    public: 
      Visitor(VSSimpleVBFVisitor* _visitor, VBFBroadcastVisitor* _parent)
	: VSSimpleVBFDispatcherStop(),
	  visitor(_visitor), user_data(0),
	  stop_processing(false), m_parent(_parent) 
      { /* nothing to see here */ }

#define __CPC(x) x(o.x)
#define __CPA(x) x = o.x;

      Visitor(const Visitor& o)
	: VSSimpleVBFDispatcherStop(o),
	  __CPC(visitor), 
	  __CPC(user_data), 
	  __CPC(stop_processing),
	  __CPC(m_parent)
      { /* nothing to see here */ }
	
      Visitor& operator= (const Visitor& o)
      {
	*static_cast<VSSimpleVBFDispatcherStop*>(this) = o;
	__CPA(visitor);
	__CPA(user_data);
	__CPA(stop_processing);
	__CPA(m_parent);
	return *this;
      };

#undef __CPC
#undef __CPA

      virtual ~Visitor();
      virtual void stopProcessingFile();
      VSSimpleVBFVisitor*       visitor;
      void*                     user_data;
      bool                      stop_processing;
    private:
      VBFBroadcastVisitor*      m_parent;
    };

    class VisitorEventData
    {
    public:
      VisitorEventData(const Visitor& _visitor)
	: visitor(_visitor.visitor), user_data(_visitor.user_data),
	  veto_master(_visitor.stop_processing), veto_packet(false),
	  veto_array_event_primary(false), veto_array_event_secondary(false),
	  veto_scope_event(false), veto_channel(false) 
      { /* nothing to see here */ }
      VSSimpleVBFVisitor*       visitor;
      void*                     user_data;
      bool                      veto_master;
      bool                      veto_packet;
      bool                      veto_array_event_primary;
      bool                      veto_array_event_secondary;
      bool                      veto_scope_event;
      bool                      veto_channel;      
    };
    
    typedef std::list<VisitorEventData> EventData;

#ifndef NOTHREADS
    pthread_mutex_t             m_visitors_mutex;
    void lock() { vsassert(pthread_mutex_lock(&m_visitors_mutex) == 0); }
    void unlock() { vsassert(pthread_mutex_unlock(&m_visitors_mutex) == 0); }
#else
    void lock() { }
    void unlock() { }
#endif

    std::list<Visitor>          m_visitors;
  };

  // ==========================================================================
  // VBFPacketTransform - Arbitrary packet manipulation
  // ==========================================================================

  class VBFPacketTransform
  {
  public:
    virtual ~VBFPacketTransform();
    virtual void transform(VPacket* packet) = 0;
  };

  class VBFIntegrateOverflowPacketTransform: public VBFPacketTransform
  {
  public:
    VBFIntegrateOverflowPacketTransform(bool overwrite = false,
					bool discard_vaev = false,
					bool leave_in_ovrf = false);
    virtual ~VBFIntegrateOverflowPacketTransform();

    virtual void transform(VPacket* packet);
    
    static void integrateOverflowDatums(VPacket* packet,
					bool overwrite = false,
					bool discard_vaev = false,
					bool leave_in_ovrf = false);

  private:
    bool m_overwrite;
    bool m_discard_vaev;
    bool m_leave_in_ovrf;
  };

  // ==========================================================================
  // VBFRationalizedClock - Determine event time from VBF packet
  // ==========================================================================

  class VBFRationalizedClock
  {
  public:
    virtual ~VBFRationalizedClock();
    virtual bool bestEventTime(const VBFPacket* packet,
			       VSTime& best_event_time) = 0;
  };

  class VBFStandardRationalizedClock: public VBFRationalizedClock
  {
  public:
    VBFStandardRationalizedClock(bool no_l3 = false,
				 const VSTime& t0 = VSTime::perpetual_past(),
				 const VSTime& t1 = VSTime::perpetual_future()):
      m_dont_use_l3(no_l3),
      m_last_event_num(), m_last_event_time(),
      m_clock_count(), m_clock_seen(), m_nclock_seen(),
      m_time0(t0), m_time1(t1) { /* nothing to see here */ }
    virtual ~VBFStandardRationalizedClock();
    virtual bool bestEventTime(const VBFPacket* packet, 
			       VSTime& best_event_time);
  private:

    struct Clock
    {
      Clock(): clock(), clock_good(), clock_num(), priority() 
      { /* nothing to see here */ }
      VSTime     clock;
      bool       clock_good;
      unsigned   clock_num;
      unsigned   priority;
    };

    bool                        m_dont_use_l3;
    unsigned                    m_last_event_num;
    VSTime                      m_last_event_time;
    std::vector<unsigned>       m_clock_count;
    std::vector<bool>           m_clock_seen;
    unsigned                    m_nclock_seen;
    VSTime                      m_time0;
    VSTime                      m_time1;
  };

  // ==========================================================================
  // VBFSequentialPacketFetcher - Class to supply packets to dispatcher
  // ==========================================================================

  class VBFSequentialPacketFetcher
  {
  public:
    VBFSequentialPacketFetcher() { /* nothing to see here */ }
    virtual ~VBFSequentialPacketFetcher();
    virtual VBFPacket* getPacket() = 0;

  private:  
    VBFSequentialPacketFetcher(const VBFSequentialPacketFetcher&);
    VBFSequentialPacketFetcher& operator= (const VBFSequentialPacketFetcher&);
  };

  class VBFSequentialPacketFetcherChained: public VBFSequentialPacketFetcher
  {
  public:
    VBFSequentialPacketFetcherChained(VBFSequentialPacketFetcher* fetcher)
      : fFetcher(fetcher) { /* nothing to see here */ }
    virtual ~VBFSequentialPacketFetcherChained();
    virtual VBFPacket* getPacket();

  protected:
    VBFSequentialPacketFetcher* fFetcher;

  private:
    VBFSequentialPacketFetcherChained
    (const VBFSequentialPacketFetcherChained&);
    VBFSequentialPacketFetcherChained& operator=
    (const VBFSequentialPacketFetcherChained&);
  };

  // --------------------------------------------------------------------------
  // VBFSimpleFetcher - Retrieve packets sequentially from VBF file
  // --------------------------------------------------------------------------

  class VBFSimpleFetcher: public VBFSequentialPacketFetcher 
  {
  public:
    VBFSimpleFetcher(VBankFileReader* reader, VBFRationalizedClock* clock,
		     VBFPacketTransform* transform);
    virtual ~VBFSimpleFetcher();
    virtual VBFPacket* getPacket();

  private:
    VBFSimpleFetcher(const VBFSimpleFetcher&);
    VBFSimpleFetcher& operator= (const VBFSimpleFetcher&);

    VBankFileReader*            fReader;
    VBFRationalizedClock*       fClock;
    VBFPacketTransform*         fTransform;
    unsigned                    fIPacket;
    unsigned                    fNPacket;
  };

  // --------------------------------------------------------------------------
  // VBFLinearizedFetcher - Reorder packets by event number
  // --------------------------------------------------------------------------

  class VBFLinearizedFetcher: public VBFSequentialPacketFetcherChained 
  {
  public:
    VBFLinearizedFetcher(VBFSequentialPacketFetcher* fetcher,
			 unsigned buffer_size,
			 bool discard_if_out_of_sequence,
			 bool verbose);
    virtual ~VBFLinearizedFetcher();
    virtual VBFPacket* getPacket();

  private:
    VBFLinearizedFetcher(const VBFLinearizedFetcher&);
    VBFLinearizedFetcher& operator= (const VBFLinearizedFetcher&);

    struct BufNode
    {
      VBFPacket*   packet;
      bool         valid;
      unsigned     event_num;
    };

    bool                        fDiscard;
    bool                        fVerbose;
    unsigned                    fIPacket;
    
    // Buffer
    unsigned                    fBufferSize;
    std::list<BufNode>          fBuffer;
    unsigned                    fNextEvent;
  };

  // --------------------------------------------------------------------------
  // VBFRandomizedFetcher - Randomize packet order discarding some fraction
  // --------------------------------------------------------------------------

  class VBFRandomizedFetcher: public VBFSequentialPacketFetcherChained
  {
  public:
    VBFRandomizedFetcher(VBFSequentialPacketFetcher* fetcher,
			 unsigned buffer_size,
			 float discard_probability,
			 float long_term_keep_probability,
			 bool verbose);
    virtual ~VBFRandomizedFetcher();
    virtual VBFPacket* getPacket();

  private:
    VBFRandomizedFetcher(const VBFRandomizedFetcher&);
    VBFRandomizedFetcher& operator= (const VBFRandomizedFetcher&);

    unsigned                    fDiscardFate;
    unsigned                    fLongTermFate;
    bool                        fVerbose;
    unsigned                    fIPacket;

    // Buffer
    unsigned                    fBufferSize;
    std::vector<VBFPacket*>     fBuffer;
    VBFPacket*                  fLongTermBuffer;
  };

  // --------------------------------------------------------------------------
  // VBFThreadedFetcher - Perform packet decoding in seperate thread
  // --------------------------------------------------------------------------

#ifndef NOTHREADS
  class VBFThreadedFetcher: public VBFSequentialPacketFetcherChained
  {
  public:
    VBFThreadedFetcher(VBFSequentialPacketFetcher* fetcher,
		       unsigned buffer_size);
    virtual ~VBFThreadedFetcher();
    virtual VBFPacket* getPacket();

  private:
    VBFThreadedFetcher(const VBFThreadedFetcher&);
    VBFThreadedFetcher& operator= (const VBFThreadedFetcher&);

    struct BufNode
    {
      BufNode(): packet() {}
      VBFPacket*   packet;
    };

    static void* threadStart(void *object);
    void run();

    bool                        fFinished;
    
    // Buffer
    unsigned                    fBufferSize;
    std::list<BufNode>          fBuffer;
    
    // Thread
    pthread_t                   fThread;
    pthread_mutex_t             fMutex;
    pthread_cond_t              fConditionNotifyConsumer;
    pthread_cond_t              fConditionNotifyProducer;
  };
#endif

  // ==========================================================================
  // VSSimpleVBFDispatcher - Parse VBF file and dispatch information to visitor
  // ==========================================================================

  class VSSimpleVBFDispatcher: private VSSimpleVBFDispatcherStop
  {
  public:
    static const uint32_t DISPATCH_NONE                = 0x00000000;
    static const uint32_t DISPATCH_THREADED            = 0x00000001;
    static const uint32_t DISPATCH_IN_ORDER            = 0x00000002;
    static const uint32_t DISPATCH_DISCARD_INVALID     = 0x00000004;
    static const uint32_t DISPATCH_VERBOSE             = 0x00000008;
    static const uint32_t DISPATCH_SHUFFLE             = 0x00000010;
    static const uint32_t DISPATCH_SHUFFLE_LOSSLESS    = 0x00000020;
    static const uint32_t DISPATCH_INTEGRATE_OVERFLOW  = 0x00000040;
    static const uint32_t DISPATCH_NTHREAD_MASK        = 0x00000f00;
    static const uint32_t DISPATCH_NTHREAD_SHIFT       = 8;
    
    static uint32_t DISPATCH_NTHREAD_GET(uint32_t flags)
    { return (flags&DISPATCH_NTHREAD_MASK)>>DISPATCH_NTHREAD_SHIFT; }
    static uint32_t DISPATCH_NTHREAD_SET(uint32_t nthread)
    { if(nthread>(DISPATCH_NTHREAD_MASK>>DISPATCH_NTHREAD_SHIFT))
	return DISPATCH_NTHREAD_MASK;
      else return nthread<<DISPATCH_NTHREAD_SHIFT; }
    
    VSSimpleVBFDispatcher(VSSimpleVBFVisitor* visitor,
			  VBFRationalizedClock* clock = 0);
    virtual ~VSSimpleVBFDispatcher();

    void resetClock(VBFRationalizedClock* clock);
    void resetVisitor(VSSimpleVBFVisitor* visitor);

    void setVisitChannels(bool visit) { fVisitChannels=visit; }
    void setDontUseL3Time(bool dont) { fDontUseL3Time=dont; }

    // Open and close file

    bool openFile(const char* filename);
    void closeFile();
    unsigned numPackets() const;

    // Dispatch a single packet

    void dispatchPacket(const VBFPacket* packet, unsigned seq_packet_num);
    void dispatchPacket(unsigned vbf_packet_num, unsigned seq_packet_num = 0,
			uint32_t flags = DISPATCH_NONE);

    // Dispatch full file - provides "value added" functionality over
    // dispatching individual packets, for example (1) using threaded
    // VBF reader to perform packet decompression on different CPU,
    // (2) reorder packets by event number to simplify analysis code
    // and (3) shuffle the packet order to test the robustness of code

    unsigned processFile(const char* filename, 
			 unsigned num_packets_max = 0,
			 uint32_t flags = DISPATCH_THREADED);
    virtual void stopProcessingFile();
    unsigned dispatchAllPackets(unsigned num_packets_max = 0, 
				uint32_t flags = DISPATCH_THREADED);

    static void getArrayEventData(VSSimpleVBFVisitor::EventType& event_type,
				  uint32_t& l2_trigger_mask,
				  const VArrayEvent* ae);

    static bool dispatchArrayTrigger(VSSimpleVBFVisitor* visitor,
				     void* user_data, 
				     const VArrayTrigger* at);
				     
    static bool dispatchScopeEvent(VSSimpleVBFVisitor* visitor,
				   void* user_data, 
				   unsigned telescope_num, const VEvent* ev,
				   bool visit_channels = true);

    // Configure dispatcher class to catch signal and stop processing file

    static void catchSignalAndStopProcessingFile(int signum);

  private:
    VSSimpleVBFDispatcher(const VSSimpleVBFDispatcher&);
    VSSimpleVBFDispatcher& operator= (const VSSimpleVBFDispatcher&);

#ifndef NOTHREADS
    class DispatcherThreadAssist
    {
    public:
      DispatcherThreadAssist(VSSimpleVBFDispatcher* dispatcher,
			     VBFSequentialPacketFetcher* fetcher,
			     pthread_mutex_t& fetcher_mutex,
			     unsigned& seq_packet_num,
			     unsigned num_packets_max,
			     bool& stop_processing_flag,
			     bool& global_stop_processing_flag);
      ~DispatcherThreadAssist();
      void join();
    private:
      DispatcherThreadAssist(const DispatcherThreadAssist&);
      DispatcherThreadAssist& operator= (const DispatcherThreadAssist&);

      static void* threadStart(void *object);
      void run();

      VSSimpleVBFDispatcher*      fDispatcher;
      VBFSequentialPacketFetcher* fFetcher;
      pthread_mutex_t&            fFetcherMutex;
      unsigned&                   fSeqPacketNum;
      unsigned                    fNumPacketsMax;
      bool&                       fStopProcessingFlag;
      bool&                       fGlobalStopProcessingFlag;
      pthread_t                   fThread;
      pthread_mutex_t             fFinishedMutex;
      pthread_cond_t              fFinishedCond;
      bool                        fFinished;
    };
#endif

    static void catchHandler(int signum);

#ifndef NOTHREADS
    pthread_mutex_t*       fProcessingMutex;
#endif

    VBFRationalizedClock*  fMyClock;
    VBFRationalizedClock*  fClock;
    VSSimpleVBFVisitor*    fVisitor;
    std::string            fFilename;
    VBankFileReader*       fReader;
    bool                   fStopProcessingFlag;
    bool                   fVisitChannels;
    bool                   fDontUseL3Time;
    static bool            sGlobalStopProcessingFlag;
  };

} // namespace VERITAS

#endif // VSSIMPLEVBF_HPP
