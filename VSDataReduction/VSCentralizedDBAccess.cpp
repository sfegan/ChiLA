//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCentralizedDBAccess.cpp

  Centralize all access to the database and provide asynchnonous data
  retrieval through thread

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       12/12/2006

  $Id: VSCentralizedDBAccess.cpp,v 3.15 2008/12/08 00:58:43 matthew Exp $

*/

#define vstream std::cout

#include <VSCentralizedDBAccess.hpp>
#include <VSChannelMap.hpp>

using namespace VERITAS;

std::auto_ptr<VSCentralizedDBAccess> VSCentralizedDBAccess::s_instance;

VSCentralizedDBAccess::~VSCentralizedDBAccess()
{
#ifndef NOTHREADS
  if(m_threaded)
    {
      wait();
      vsassert(pthread_cond_destroy(&m_cond) == 0); 
      vsassert(pthread_mutex_destroy(&m_mutex) == 0);
    }
#endif
}

void VSCentralizedDBAccess::
prefetchFromDB(unsigned run_no, const VSTime& start_time,
	       const std::vector<unsigned>& telescopes, bool no_threads)
{
#ifndef NOTHREADS
  if(no_threads)doPrefetch(run_no,start_time,telescopes,true);
  else
    {
      m_run_no           = run_no;
      m_start_time       = start_time;
      m_telescopes       = telescopes;
      m_threaded         = true;
      m_threaded_active  = true;
      vsassert(pthread_cond_init(&m_cond, NULL) == 0);
      vsassert(pthread_mutex_init(&m_mutex, NULL) == 0);
      vsassert(pthread_create(&m_thread, 0, &threadStart, this) == 0);
      vsassert(pthread_detach(m_thread) == 0);
    }
#else
  doPrefetch(run_no,start_time,telescopes,true);
#endif  
}

bool VSCentralizedDBAccess::getRunInfo(unsigned run_no, RunInfoDatum& datum)
{
  wait();
  if(m_run_info_cache.run_id == run_no)
    {
      datum = m_run_info_cache;
      return true;
    }
  return doFetchRunInfo(run_no,datum);
}

bool VSCentralizedDBAccess::
getTrackingTarget(unsigned iscope, 
		  const VSTime& start, const VSTime& stop,
		  std::vector<TrackingTargetDatum>& data)
{
  wait();
  for(unsigned icache=0;icache<m_tracking_target_cache.size();icache++)
    {
      if((m_tracking_target_cache[icache].scope_num == iscope)
	 &&(m_tracking_target_cache[icache].start_time <= start)
	 &&(m_tracking_target_cache[icache].stop_time >= stop))
	{
	  data.clear();

	  for(unsigned idatum=0;
	      idatum<m_tracking_target_cache[icache].data.size();
	      idatum++)
	    {
	      VSTime& 
		ts(m_tracking_target_cache[icache].data[idatum].db_start_time);
	      VSTime& 
		te(m_tracking_target_cache[icache].data[idatum].db_end_time);
	      bool te_null = 
		m_tracking_target_cache[icache]
		.data[idatum].db_end_time_is_null;

	      if((te>=start || te_null) && (ts < stop))
		data.push_back(m_tracking_target_cache[icache].data[idatum]);
	    }
	  if(m_stream)
	    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
			<< " tracking target records\n";
	  return true;
	}
    }

  bool _ret = doFetchTrackingTarget(iscope,start,stop,data);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
		<< " tracking target records\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getCorrectionParameters(unsigned iscope, 
			std::vector<CorrectionParametersDatum>& data)
{
  wait();
  bool _ret = true;
  if((iscope<m_correction_parameters.size())
     &&(m_correction_parameters[iscope].first))
    data = m_correction_parameters[iscope].second;
  else _ret = doFetchCorrectionParameters(iscope, data);
  
  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
		<< " correction parameters data sets\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getPointing(unsigned iscope, const VSTime& start, const VSTime& stop,
	    std::vector<PointingDatum>& data)
{
  wait();
  for(unsigned icache=0;icache<m_pointing_cache.size();icache++)
    {
      if((m_pointing_cache[icache].scope_num == iscope)
	 &&(m_pointing_cache[icache].start_time <= start)
	 &&(m_pointing_cache[icache].stop_time >= stop))
	{
	  data.clear();
	  for(unsigned idatum=0;idatum<m_pointing_cache[icache].data.size();
	      idatum++)
	    {
	      VSTime& 
		timestamp(m_pointing_cache[icache].data[idatum].timestamp);
	      if(timestamp>=start && timestamp<stop)
		data.push_back(m_pointing_cache[icache].data[idatum]);
	    }
	  if(m_stream)
	    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
			<< " pointing records\n";
	  return true;
	}
    }

  bool _ret = doFetchPointing(iscope,start,stop,data);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
		<< " pointing records\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getPower(unsigned iscope, const VSTime& start, const VSTime& stop,
	 std::vector<PowerDatum>& data)
{
  wait();
  for(unsigned icache=0;icache<m_power_cache.size();icache++)
    {
      if((m_power_cache[icache].scope_num == iscope)
	 &&(m_power_cache[icache].start_time <= start)
	 &&(m_power_cache[icache].stop_time >= stop))
	{
	  data.clear();
	  for(unsigned idatum=0;idatum<m_power_cache[icache].data.size();
	      idatum++)
	    {
	      VSTime& 
		db_start(m_power_cache[icache].data[idatum].db_start_time);
	      VSTime& 
		db_end(m_power_cache[icache].data[idatum].db_end_time);
	      if((db_start<=start && db_end>start)
		 ||(db_start>start && db_start<stop))
		data.push_back(m_power_cache[icache].data[idatum]);
	    }
	  if(m_stream)
	    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
			<< " HV power records\n";
	  return true;
	}
    }

  bool _ret = doFetchPower(iscope,start,stop,data);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
		<< " HV power records\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getHVStatus(unsigned iscope, const VSTime& start, const VSTime& stop,
	    std::vector<HVStatusDatum>& data)
{
  wait();
  for(unsigned icache=0;icache<m_hv_status_cache.size();icache++)
    {
      if((m_hv_status_cache[icache].scope_num == iscope)
	 &&(m_hv_status_cache[icache].start_time <= start)
	 &&(m_hv_status_cache[icache].stop_time >= stop))
	{
	  data.clear();
	  for(unsigned idatum=0;
	      idatum<m_hv_status_cache[icache].data.size();idatum++)
	    {
	      VSTime& dbt(m_hv_status_cache[icache].data[idatum].timestamp);

	      if(dbt>=start && dbt<stop)
		data.push_back(m_hv_status_cache[icache].data[idatum]);
	    }
	  if(m_stream)
	    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
			<< " HV current records\n";
	  return true;
	}
    }

  bool _ret = doFetchHVStatus(iscope,start,stop,data);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
		<< " HV current records\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getL1Rate(unsigned iscope, const VSTime& start, const VSTime& stop,
	  std::vector<L1RateDatum>& data)
{
  wait();
  for(unsigned icache=0;icache<m_l1_rate_cache.size();icache++)
    {
      if((m_l1_rate_cache[icache].scope_num == iscope)
	 &&(m_l1_rate_cache[icache].start_time <= start)
	 &&(m_l1_rate_cache[icache].stop_time >= stop))
	{
	  data.clear();
	  for(unsigned idatum=0;
	      idatum<m_l1_rate_cache[icache].data.size();idatum++)
	    {
	      VSTime& db_time(m_l1_rate_cache[icache].data[idatum].timestamp);

	      if(db_time>=start && db_time<stop)
		data.push_back(m_l1_rate_cache[icache].data[idatum]);
	    }
	  if(m_stream)
	    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
			<< " L1 rate records\n";
	  return true;
	}
    }

  bool _ret = doFetchL1Rate(iscope,start,stop,data);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
		<< " L1 rate records\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getFIR(unsigned iscope, const VSTime& start, const VSTime& stop,
       std::vector<FIRDatum>& data)
{
  wait();

  // --------------------------------------------------------------------------
  // Use scope index 4 for FIR0 (index 0 in the DB)
  // --------------------------------------------------------------------------
  unsigned fir_scope = 0;
  if(iscope != 4) fir_scope = iscope + 1;

  for(unsigned icache=0;icache<m_fir_cache.size();icache++)
    {
      if((m_fir_cache[icache].scope_num == iscope)
	 &&(m_fir_cache[icache].start_time <= start)
	 &&(m_fir_cache[icache].stop_time >= stop))
	{
	  data.clear();
	  for(unsigned idatum=0;
	      idatum<m_fir_cache[icache].data.size();idatum++)
	    {
	      VSTime& db_time(m_fir_cache[icache].data[idatum].timestamp);

	      if(db_time>=start && db_time<stop)
		data.push_back(m_fir_cache[icache].data[idatum]);
	    }

	  if(m_stream)
	    (*m_stream) << "DB: loaded " << data.size() << " FIR " 
			<< fir_scope << " records\n";
	  return true;
	}
    }

  bool _ret = doFetchFIR(iscope,start,stop,data);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << fir_scope
		<< " FIR records\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getL3Scope(unsigned iscope, const VSTime& start, const VSTime& stop,
	   std::vector<L3ScopeDatum>& data)
{
  wait();
  for(unsigned icache=0;icache<m_l3_scope_cache.size();icache++)
    {
      if((m_l3_scope_cache[icache].scope_num == iscope)
	 &&(m_l3_scope_cache[icache].start_time <= start)
	 &&(m_l3_scope_cache[icache].stop_time >= stop))
	{
	  data.clear();
	  for(unsigned idatum=0;
	      idatum<m_l3_scope_cache[icache].data.size();idatum++)
	    {
	      VSTime& db_time(m_l3_scope_cache[icache].data[idatum].timestamp);

	      if(db_time>=start && db_time<stop)
		data.push_back(m_l3_scope_cache[icache].data[idatum]);
	    }
	  if(m_stream)
	    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
			<< " L3 telescope records\n";
	  return true;
	}
    }

  bool _ret = doFetchL3Scope(iscope,start,stop,data);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << data.size() << " T" << iscope+1
		<< " L3 telescope records\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getFADC(const VSTime& time,
	std::vector<FADCChannelSettingsDatum>& chan_settings,
	std::vector<FADCBoardSettingsDatum>& board_settings,
	std::vector<FADCChannelRelation>& chan_relation,
	std::vector<FADCSlotRelation>& slot_relation)
{
  wait();
  bool _ret = false;

  if((m_fadc_chan_settings.config_time == time)
     &&(m_fadc_board_settings.config_time == time)
     &&(m_fadc_chan_relation.config_time == time)
     &&(m_fadc_slot_relation.config_time == time))
    {
      chan_settings    = m_fadc_chan_settings.data;
      board_settings   = m_fadc_board_settings.data;
      chan_relation    = m_fadc_chan_relation.data;
      slot_relation    = m_fadc_slot_relation.data;
      _ret = true;      
    }
  else
    {
      _ret = 
	doFetchFADC(time,
		    chan_settings,board_settings,chan_relation,slot_relation);
    }

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << chan_settings.size() 
		<< " FADC channel settings\n"
		<< "DB: loaded " << board_settings.size() 
		<< " FADC board settings\n"
		<< "DB: loaded " << chan_relation.size() 
		<< " FADC channel relations\n"
		<< "DB: loaded " << slot_relation.size() 
		<< " FADC slot relations\n";

  return _ret;
}

bool VSCentralizedDBAccess::
getTargetTable(std::vector<TargetTableCoord>& target_table)
{
  wait();
  bool _ret = true;
  if(m_target_table.first)target_table = m_target_table.second;
  else _ret = doFetchTargetTable(target_table);

  if((_ret)&&(m_stream))
    (*m_stream) << "DB: loaded " << target_table.size() 
		<< " target table records\n";

  return _ret;
}

VSCentralizedDBAccess* VSCentralizedDBAccess::getInstance()
{
  if(s_instance.get() == 0)s_instance.reset(new VSCentralizedDBAccess);
  return s_instance.get();
}

// ============================================================================
//
// Internal functions
//
// ============================================================================

VSCentralizedDBAccess::VSCentralizedDBAccess():
  m_stream(&vstream),
  m_run_info_cache(), m_tracking_target_cache(), m_correction_parameters(),
  m_pointing_cache(), m_power_cache(), m_hv_status_cache(), 
  m_l1_rate_cache(), m_fir_cache(), 
  m_l3_scope_cache(), m_vpm_stars_scope_cache(), m_vpm_leds_scope_cache(),
  m_fadc_chan_settings(), m_fadc_board_settings(), 
  m_fadc_chan_relation(), m_fadc_slot_relation(), m_target_table()
#ifndef NOTHREADS
  , m_run_no(), m_start_time(),  m_telescopes(),
  m_threaded(false), m_threaded_active(false), m_mutex(), m_cond(), m_thread()
#endif
{
  // nothing to see here
}

VSDatabase* VSCentralizedDBAccess::
constructDBIfRequired(VSDatabase*& supplied_db)
{
  VSDatabase* my_db = 0;
  if(supplied_db == 0)
    {
      supplied_db = my_db = VSDBFactory::getInstance()->createVSDB();

      if(supplied_db->useDatabase("VERITAS") < 0)
	{
	  std::cerr << __PRETTY_FUNCTION__ 
		    << ": could not connect to database: VERITAS " 
		    << std::endl;
	  exit(EXIT_FAILURE);
	}
    }
  return my_db;
}

void VSCentralizedDBAccess::
doPrefetch(unsigned run_no, const VSTime& start_time,
	   const std::vector<unsigned>& telescopes, bool no_real_fetch)
{
  VSChannelMap channel_map(start_time);
  VSDatabase* db = 0;
  constructDBIfRequired(db);

  bool got_run_info = doFetchRunInfo(run_no, m_run_info_cache, db);

  VSTime fetch_start_time = start_time;
  VSTime fetch_stop_time = start_time;
  fetch_stop_time += INT64_C(1800000000000);
  
  if(got_run_info)
    {
      if(m_run_info_cache.db_start_time < fetch_start_time)
	fetch_start_time = m_run_info_cache.db_start_time;
      if(m_run_info_cache.db_end_time.isGood())
	fetch_stop_time = m_run_info_cache.db_end_time;
      if(fetch_stop_time - fetch_start_time > INT64_C(3600000000000))
	{
	  fetch_stop_time = start_time;
	  fetch_stop_time += INT64_C(3600000000000);
	}
    }

  fetch_start_time -= INT64_C(60000000000);
  fetch_stop_time += INT64_C(60000000000);;

  m_correction_parameters.clear();
  m_correction_parameters.resize(telescopes.size());
  for(unsigned iscope = 0; iscope<telescopes.size(); iscope++)
    m_correction_parameters[iscope].first =
      doFetchCorrectionParameters(iscope,
				  m_correction_parameters[iscope].second,
				  db);

  const unsigned nfir = channel_map.nfir();
  m_fir_cache.clear();
  m_fir_cache.resize(nfir);
  for(unsigned ifir = 0; ifir<nfir; ifir++)
    {
      m_fir_cache[ifir].start_time = fetch_start_time;
      m_fir_cache[ifir].stop_time  = fetch_stop_time;
      m_fir_cache[ifir].scope_num  = ifir;

      doFetchFIR(m_fir_cache[ifir].scope_num,
		 m_fir_cache[ifir].start_time,
		 m_fir_cache[ifir].stop_time, 
		 m_fir_cache[ifir].data, db);
    }

  m_fadc_chan_settings.config_time  = start_time;
  m_fadc_board_settings.config_time = start_time;
  m_fadc_chan_relation.config_time  = start_time;
  m_fadc_slot_relation.config_time  = start_time;
  doFetchFADC(start_time,
	      m_fadc_chan_settings.data, m_fadc_board_settings.data,
	      m_fadc_chan_relation.data, m_fadc_slot_relation.data,
	      db);

  m_target_table.first = doFetchTargetTable(m_target_table.second, db);

#define FILLSCOPECACHE(cache,fillfn)					\
  cache.clear();							\
  cache.resize(telescopes.size());					\
  for(unsigned iscope = 0; iscope<telescopes.size(); iscope++)		\
    {									\
      cache[iscope].start_time = fetch_start_time;			\
      cache[iscope].stop_time  = fetch_stop_time;			\
      cache[iscope].scope_num  = telescopes[iscope];			\
      fillfn(cache[iscope].scope_num,					\
	     cache[iscope].start_time,					\
	     cache[iscope].stop_time,					\
	     cache[iscope].data, db);					\
    }
  
  FILLSCOPECACHE(m_tracking_target_cache, doFetchTrackingTarget);
  FILLSCOPECACHE(m_pointing_cache,        doFetchPointing);
  FILLSCOPECACHE(m_power_cache,           doFetchPower);
  FILLSCOPECACHE(m_hv_status_cache,       doFetchHVStatus);
  FILLSCOPECACHE(m_l1_rate_cache,         doFetchL1Rate);
  FILLSCOPECACHE(m_l3_scope_cache,        doFetchL3Scope);
  FILLSCOPECACHE(m_vpm_stars_scope_cache, doFetchVPMStarsScope);
  FILLSCOPECACHE(m_vpm_leds_scope_cache,  doFetchVPMLEDsScope);
  
  delete db;
}

static inline void setVSTime(VSTime& vs_time, const VSDBDateTime& db_time)
{
  vs_time.
    setFromCalendarDateAndTime(db_time.year,    db_time.month,
			       db_time.day,     db_time.hour,
			       db_time.minute,  db_time.second,
			       db_time.microsecond * 1000000);
}
  
bool VSCentralizedDBAccess::
doFetchRunInfo(unsigned run_no, RunInfoDatum& datum, VSDatabase* db)
{
  VSDatabase* my_db = constructDBIfRequired(db);

  VSDBStatement* stmt = 
    db->createSelectQuery("tblRun_Info", "run_id=?", 
			  "run_id,run_type,observing_mode,run_status,"
			  "db_start_time,db_end_time,data_start_time,"
			  "data_end_time,duration,weather,config_mask,"
			  "pointing_mode,trigger_config,trigger_multiplicity,"
			  "trigger_coincidence,offset_distance,offset_angle,"
			  "source_id");

  vsassert(stmt);
  
  stmt->bindToParam(run_no);

  int nresults = stmt->execute();
  vsassert(nresults>=0);

  bool null_db_end_time;
  bool null_data_end_time;
  
  VSDBDateTime db_start_time;
  VSDBDateTime db_end_time;
  VSDBDateTime data_start_time;
  VSDBDateTime data_end_time;
  VSDBTime duration;

  stmt->bindToResult(datum.run_id);
  stmt->bindToResult(datum.run_type);
  stmt->bindToResult(datum.observing_mode);
  stmt->bindToResult(datum.run_status);
  stmt->bindToResult(db_start_time);
  stmt->bindToResult(db_end_time,&null_db_end_time);
  stmt->bindToResult(data_start_time);
  stmt->bindToResult(data_end_time,&null_data_end_time);
  stmt->bindToResult(duration);
  stmt->bindToResult(datum.weather);
  stmt->bindToResult(datum.config_mask);
  stmt->bindToResult(datum.pointing_mode);
  stmt->bindToResult(datum.trigger_config);
  stmt->bindToResult(datum.trigger_multiplicity);
  stmt->bindToResult(datum.trigger_coincidence);
  stmt->bindToResult(datum.offset_distance);
  stmt->bindToResult(datum.offset_angle);
  stmt->bindToResult(datum.source_id);

  bool retval = false;
  if(stmt->retrieveNextRow() == 1)
    {
      setVSTime(datum.db_start_time, db_start_time);
      if(!null_db_end_time)setVSTime(datum.db_end_time, db_end_time);
      else datum.db_end_time = VSTime();

      setVSTime(datum.data_start_time, data_start_time);
      if(!null_data_end_time)setVSTime(datum.data_end_time, data_end_time);
      else datum.data_end_time = VSTime();

      datum.duration_msec = 
	duration.hour*3660000 + duration.minute*60000 + duration.second*1000 
	+ duration.microsecond;

      vsassert(stmt->retrieveNextRow() != 1);

      retval = true;
    }

  delete stmt;
  delete my_db;

  return retval;
}

bool VSCentralizedDBAccess::
doFetchTrackingTarget(unsigned iscope, 
		      const VSTime& start, const VSTime& stop,
		      std::vector<TrackingTargetDatum>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);
	  
  std::ostringstream query_table;
  query_table << "tblPositioner_Telescope" << iscope << "_Targets";

  VSDBStatement* stmt = 
    db->createSelectQuery(query_table.str(),
			  "( db_end_time IS NULL OR db_end_time>=? ) AND "
			  "db_start_time<=? ORDER BY db_start_time",
			  "db_start_time,db_end_time,mode,angle1,angle2,epoch,"
			  "source_id,target_mode,target_param1,target_param2,"
			  "config_mask,pointing_mode,pointing_param");
  vsassert(stmt);

  uint64_t starttime = start.getDBTimeStamp();
  uint64_t stoptime = stop.getDBTimeStamp();
  
  stmt->bindToParam(starttime);
  stmt->bindToParam(stoptime);

  int nresults = stmt->execute();
  vsassert(nresults>=0);

  data.reserve(nresults);

  TrackingTargetDatum target;
  VSDBDateTime db_start_time;
  VSDBDateTime db_end_time;
  bool source_id_is_null;

  stmt->bindToResult(db_start_time);
  stmt->bindToResult(db_end_time,&target.db_end_time_is_null);
  stmt->bindToResult(target.target_type);
  stmt->bindToResult(target.target_angle1_rad);
  stmt->bindToResult(target.target_angle2_rad);
  stmt->bindToResult(target.target_epoch);
  stmt->bindToResult(target.target_name,&source_id_is_null);
  stmt->bindToResult(target.tracking_mode);
  stmt->bindToResult(target.tracking_param1);
  stmt->bindToResult(target.tracking_param2);
  stmt->bindToResult(target.pointing_array_mask);
  stmt->bindToResult(target.pointing_mode);
  stmt->bindToResult(target.pointing_param);

  while(stmt->retrieveNextRow() == 1)
    {
      setVSTime(target.db_start_time, db_start_time);
      if(!target.db_end_time_is_null)setVSTime(target.db_end_time,db_end_time);
      else target.db_end_time = VSTime();

      target.has_extended_target_info = !source_id_is_null;
      data.push_back(target);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchCorrectionParameters(unsigned iscope, 
			    std::vector<CorrectionParametersDatum>& data,
			    VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);
	  
  std::ostringstream query_table;
  query_table << "tblPositioner_Telescope" << iscope << "_Corrections";

  VSDBStatement* stmt = 
    db->createSelectQuery(query_table.str(), "",
			  "db_start_time,db_end_time,"
			  "enable_offsets,enable_corrections,"
			  "az_gear_ratio,el_gear_ratio,"
			  "az_offset_rad,el_offset_rad,"
			  "azimuth_axisNS_rad,azimuth_axisEW_rad,"
			  "perpendicularity_rad,collimation_rad,"
			  "flexure_cosEL_rad,flexure_sin2EL_rad,"
			  "enable_velocity_feed_forward,"
			  "vff_pos_el_slope_sec,"
			  "vff_pos_el_threshold_rad_per_sec,"
			  "vff_neg_el_slope_sec,"
			  "vff_neg_el_threshold_rad_per_sec,"
			  "vff_pos_az_slope_sec,"
			  "vff_pos_az_threshold_rad_per_sec,"
			  "vff_neg_az_slope_sec,"
			  "vff_neg_az_threshold_rad_per_sec," 
			  "comment");
  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults>=0);

  data.reserve(nresults);

  CorrectionParametersDatum cpd;
  VSDBDateTime db_start_time;
  VSDBDateTime db_end_time;

  std::string str_enable_offsets;
  std::string str_enable_corrections;
  std::string str_enable_vff;

  stmt->bindToResult(db_start_time);
  stmt->bindToResult(db_end_time,&cpd.db_end_time_is_null);
  stmt->bindToResult(str_enable_offsets);
  stmt->bindToResult(str_enable_corrections);
  stmt->bindToResult(cpd.az_ratio);
  stmt->bindToResult(cpd.el_ratio);
  stmt->bindToResult(cpd.az_offset);
  stmt->bindToResult(cpd.el_offset);
  stmt->bindToResult(cpd.az_ns);
  stmt->bindToResult(cpd.az_ew);
  stmt->bindToResult(cpd.el_udew);
  stmt->bindToResult(cpd.fp_az);
  stmt->bindToResult(cpd.flex_el_A);
  stmt->bindToResult(cpd.flex_el_B);
  stmt->bindToResult(str_enable_vff);
  stmt->bindToResult(cpd.el_pos_vff_s);
  stmt->bindToResult(cpd.el_pos_vff_t);
  stmt->bindToResult(cpd.el_neg_vff_s);
  stmt->bindToResult(cpd.el_neg_vff_t);
  stmt->bindToResult(cpd.az_pos_vff_s);
  stmt->bindToResult(cpd.az_pos_vff_t);
  stmt->bindToResult(cpd.az_neg_vff_s);
  stmt->bindToResult(cpd.az_neg_vff_t);
  stmt->bindToResult(cpd.comment);

  while(stmt->retrieveNextRow() == 1)
    {
      setVSTime(cpd.db_start_time, db_start_time);
      if(!cpd.db_end_time_is_null)setVSTime(cpd.db_end_time,db_end_time);
      else cpd.db_end_time = VSTime();
      cpd.enable_offsets     = (str_enable_offsets == "enabled");
      cpd.enable_corrections = (str_enable_corrections == "enabled");
      cpd.enable_vff         = (str_enable_vff == "enabled");
      data.push_back(cpd);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchPointing(unsigned iscope, const VSTime& start, const VSTime& stop,
		std::vector<PointingDatum>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);
	  
  std::ostringstream query_table;
  query_table << "tblPositioner_Telescope" << iscope << "_Status";

  VSDBStatement* stmt = 
    db->createSelectQuery(query_table.str(),
			  "timestamp>=? AND timestamp<? ORDER BY timestamp",
			  "timestamp,elevation_raw,azimuth_raw,"
			  "elevation_meas,azimuth_meas,"
			  "elevation_target,azimuth_target");
  vsassert(stmt);

  uint64_t starttime = start.getMSTimeStamp();
  uint64_t stoptime = stop.getMSTimeStamp();
  
  stmt->bindToParam(starttime);
  stmt->bindToParam(stoptime);

  int nresults = stmt->execute();
  vsassert(nresults>=0);

  data.reserve(nresults);

  PointingDatum pointing;
  uint64_t ts;
  stmt->bindToResult(ts);
  stmt->bindToResult(pointing.elevation_raw);
  stmt->bindToResult(pointing.azimuth_raw);
  stmt->bindToResult(pointing.elevation_meas);
  stmt->bindToResult(pointing.azimuth_meas);
  stmt->bindToResult(pointing.elevation_target);
  stmt->bindToResult(pointing.azimuth_target);

  while(stmt->retrieveNextRow() == 1)
    {
      pointing.timestamp.setFromMSTimeStamp(ts);
      data.push_back(pointing);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchPower(unsigned iscope, 
	     const VSTime& start, const VSTime& stop,
	     std::vector<PowerDatum>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);

  std::ostringstream query_table;
  query_table << "tblHV_Telescope" << iscope << "_Power";

  const uint64_t lo_dbtime = start.getDBTimeStamp();
  const uint64_t hi_dbtime = stop.getDBTimeStamp();
      
  std::ostringstream query_cond;
  query_cond << "( ( ( db_start_time<=" 
	     << lo_dbtime << " ) AND ( db_end_time>=" << hi_dbtime
	     << " OR db_end_time IS NULL ) ) OR ( db_start_time>="
	     << lo_dbtime << " AND db_start_time<=" << hi_dbtime << " ) )"
	     << " ORDER BY db_start_time";

  VSDBStatement* stmt = 
    db->createSelectQuery(query_table.str(),
			  query_cond.str(),
			  "db_start_time, db_end_time, channel, status");
  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults >= 0);

  PowerDatum datum;
  bool null_db_end_time;
  VSDBDateTime db_start_time;
  VSDBDateTime db_end_time;
  std::string status;

  stmt->bindToResult(db_start_time);
  stmt->bindToResult(db_end_time,&null_db_end_time);
  stmt->bindToResult(datum.ichan);
  stmt->bindToResult(status);

  while(stmt->retrieveNextRow() == 1)
    {
      if(datum.ichan == 0)continue;
      datum.ichan--;

      setVSTime(datum.db_start_time,db_start_time);
      if(!null_db_end_time)setVSTime(datum.db_end_time,db_end_time);
      else datum.db_end_time = VSTime::perpetual_future();

      if(status == "ON")datum.power_on=true;
      else datum.power_on=false;

      data.push_back(datum);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchHVStatus(unsigned iscope, const VSTime& start, const VSTime& stop,
		std::vector<HVStatusDatum>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);

  std::ostringstream query_table;
  query_table << "tblHV_Telescope" << iscope << "_Status";

  const uint64_t lo_dbtime = start.getDBTimeStamp();
  const uint64_t hi_dbtime = stop.getDBTimeStamp();
      
  std::ostringstream query_cond;
  query_cond << "( db_start_time>=" << lo_dbtime 
	     << " AND db_start_time<=" << hi_dbtime << " )"
	     << " ORDER BY db_start_time";

  VSDBStatement* stmt = 
    db->createSelectQuery(query_table.str(),
			  query_cond.str(),
			  "db_start_time,channel,voltage_meas,current_meas");
  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults >= 0);

  VSDBDateTime db_start_time;
  HVStatusDatum datum;

  stmt->bindToResult(db_start_time);
  stmt->bindToResult(datum.ichan);
  stmt->bindToResult(datum.voltage);
  stmt->bindToResult(datum.current);

  while(stmt->retrieveNextRow() == 1)
    {
      if(datum.ichan == 0)continue;
      datum.ichan--;
      setVSTime(datum.timestamp,db_start_time);
      data.push_back(datum);
    }

  delete stmt;
  delete my_db;

  return true;
}


bool VSCentralizedDBAccess::
doFetchL1Rate(unsigned iscope, const VSTime& start, const VSTime& stop,
	      std::vector<L1RateDatum>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);

  const uint64_t lo_dbtime = start.getDBTimeStamp();
  const uint64_t hi_dbtime = stop.getDBTimeStamp();
      
  std::ostringstream query_cond;
  query_cond << "( telescope_id=" << iscope 
	     << " AND timestamp>=" << lo_dbtime 
	     << " AND timestamp<=" << hi_dbtime << " )"
	     << " ORDER BY timestamp";

  VSDBStatement* stmt = 
    db->createSelectQuery("tblL1_TriggerInfo", query_cond.str(),
			  "timestamp,pixel_id,rate");
  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults >= 0);

  L1RateDatum datum;
  VSDBDateTime timestamp;

  stmt->bindToResult(timestamp);
  stmt->bindToResult(datum.ichan);
  stmt->bindToResult(datum.rate);

  while(stmt->retrieveNextRow() == 1)
    {
      setVSTime(datum.timestamp,timestamp);
      data.push_back(datum);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchFIR(unsigned iscope, const VSTime& start, const VSTime& stop,
	   std::vector<FIRDatum>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);

  const uint64_t lo_dbtime = start.getDBTimeStamp();
  const uint64_t hi_dbtime = stop.getDBTimeStamp();
      
  // --------------------------------------------------------------------------
  // Use scope index 4 for FIR0 (index 0 in the DB)
  // --------------------------------------------------------------------------
  unsigned fir_scope = 0;
  if(iscope != 4) fir_scope = iscope+1;

  std::ostringstream query_cond;
  query_cond << "( telescope_id=" << fir_scope
	     << " AND timestamp>=" << lo_dbtime 
	     << " AND timestamp<=" << hi_dbtime << " )"
	     << " ORDER BY timestamp";

  VSDBStatement* stmt = 
    db->createSelectQuery("tblFIR_Pyrometer_Info", query_cond.str(),
			  "timestamp,ambient_temp,radiant_sky_temp");
  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults >= 0);

  FIRDatum datum;
  VSDBDateTime timestamp;

  stmt->bindToResult(timestamp);
  stmt->bindToResult(datum.ambient);
  stmt->bindToResult(datum.sky);

  while(stmt->retrieveNextRow() == 1)
    {
      setVSTime(datum.timestamp,timestamp);
      data.push_back(datum);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchL3Scope(unsigned iscope, const VSTime& start, const VSTime& stop,
	       std::vector<L3ScopeDatum>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);

  const uint64_t lo_dbtime = start.getMSTimeStamp();
  const uint64_t hi_dbtime = stop.getMSTimeStamp();
      
  std::ostringstream query_cond;
  query_cond << "( telescope_id=" << iscope 
	     << " AND timestamp>=" << lo_dbtime 
	     << " AND timestamp<=" << hi_dbtime << " )"
	     << " ORDER BY timestamp";

  VSDBStatement* stmt = 
    db->createSelectQuery("tblL3_Telescope_TriggerInfo", query_cond.str(),
			  "timestamp,L2,QI,HM,NP,L2LL3,L3,"
			  "VDAQBusyScaler,TenMHzScaler");

  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults >= 0);

  L3ScopeDatum datum;
  uint64_t timestamp;

  stmt->bindToResult(timestamp);
  stmt->bindToResult(datum.scaler_l2);
  stmt->bindToResult(datum.scaler_qi);
  stmt->bindToResult(datum.scaler_hm);
  stmt->bindToResult(datum.scaler_np);
  stmt->bindToResult(datum.scaler_l2_less_l3);
  stmt->bindToResult(datum.scaler_l3);
  stmt->bindToResult(datum.scaler_vdacq_busy);
  stmt->bindToResult(datum.scaler_ten_mhz);

  while(stmt->retrieveNextRow() == 1)
    {
      datum.timestamp.setFromMSTimeStamp(timestamp);
      data.push_back(datum);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchVPMStarsScope(unsigned iscope, 
		     const VSTime& start, const VSTime& stop,
		     std::vector<VPMCentroids>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);

  const uint64_t lo_dbtime = llround(start.getMJDDbl()*1e8);
  const uint64_t hi_dbtime = llround(stop.getMJDDbl()*1e8);
      
  std::ostringstream query_table;
  query_table << "tblPointing_Monitor_Telescope" << iscope << "_Centroids";

  std::ostringstream query_cond;
  query_cond << "( mjd>=" << lo_dbtime << " AND mjd<=" << hi_dbtime << " )"
	     << " ORDER BY mjd";

  VSDBStatement* stmt = 
    db->createSelectQuery(query_table.str(), query_cond.str(),
			  "mjd,"
			  "x01,y01,b01,x02,y02,b02,x03,y03,b03,"
			  "x04,y04,b04,x05,y05,b05,x06,y06,b06,"
			  "x07,y07,b07,x08,y08,b08,x09,y09,b09,"
			  "x10,y10,b10,x11,y11,b11,x12,y12,b12,"
			  "x13,y13,b13,x14,y14,b14,x15,y15,b15,"
			  "x16,y16,b16,x17,y17,b17,x18,y18,b18,"
			  "x19,y19,b19,x20,y20,b20,x21,y21,b21,"
			  "x22,y22,b22,x23,y23,b23,x24,y24,b24,"
			  "x25,y25,b25,x26,y26,b26,x27,y27,b27,"
			  "x28,y28,b28,x29,y29,b29,x30,y30,b30");

  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults >= 0);

  VPMCentroids datum;
  uint64_t mjd;

  stmt->bindToResult(mjd);
#define BINDXYB(x,y,b) stmt->bindToResult(datum.x);   \
  stmt->bindToResult(datum.y);			      \
  stmt->bindToResult(datum.b)			      \

  BINDXYB(x01, y01, b01); BINDXYB(x02, y02, b02); BINDXYB(x03, y03, b03);
  BINDXYB(x04, y04, b04); BINDXYB(x05, y05, b05); BINDXYB(x06, y06, b06);
  BINDXYB(x07, y07, b07); BINDXYB(x08, y08, b08); BINDXYB(x09, y09, b09);
  BINDXYB(x10, y10, b10); BINDXYB(x11, y11, b11); BINDXYB(x12, y12, b12);
  BINDXYB(x13, y13, b13); BINDXYB(x14, y14, b14); BINDXYB(x15, y15, b15);
  BINDXYB(x16, y16, b16); BINDXYB(x17, y17, b17); BINDXYB(x18, y18, b18);
  BINDXYB(x19, y19, b19); BINDXYB(x20, y20, b20); BINDXYB(x21, y21, b21);
  BINDXYB(x22, y22, b22); BINDXYB(x23, y23, b23); BINDXYB(x24, y24, b24);
  BINDXYB(x25, y25, b25); BINDXYB(x26, y26, b26); BINDXYB(x27, y27, b27);
  BINDXYB(x28, y28, b28); BINDXYB(x29, y29, b29); BINDXYB(x30, y30, b30);

#undef BINDXYB

  while(stmt->retrieveNextRow() == 1)
    {
      datum.timestamp.setFromMJDDbl(double(mjd)*1e-8);
      data.push_back(datum);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchVPMLEDsScope(unsigned iscope, 
		    const VSTime& start, const VSTime& stop,
		    std::vector<VPMLEDs>& data, VSDatabase* db)
{
  data.clear();
  
  VSDatabase* my_db = constructDBIfRequired(db);

  const uint64_t lo_dbtime = llround(start.getMJDDbl()*1e8);
  const uint64_t hi_dbtime = llround(stop.getMJDDbl()*1e8);
      
  std::ostringstream query_table;
  query_table << "tblPointing_Monitor_Telescope" << iscope << "_LEDs";

  std::ostringstream query_cond;
  query_cond << "( mjd>=" << lo_dbtime << " AND mjd<=" << hi_dbtime << " )"
	     << " ORDER BY mjd";

  VSDBStatement* stmt = 
    db->createSelectQuery(query_table.str(), query_cond.str(),
			  "mjd,n1,x1,y1,n2,x2,y2,n3,x3,y3,n4,x4,y4");

  vsassert(stmt);

  int nresults = stmt->execute();
  vsassert(nresults >= 0);

  VPMLEDs datum;
  uint64_t mjd;

  stmt->bindToResult(mjd);

  #define BINDNXY(n,x,y) stmt->bindToResult(datum.n);	      \
  stmt->bindToResult(datum.x);			      \
  stmt->bindToResult(datum.y)			      \

  BINDNXY(n1, x1, y1);
  BINDNXY(n2, x2, y2);
  BINDNXY(n3, x3, y3);
  BINDNXY(n4, x4, y4);

#undef BINDNYB

  while(stmt->retrieveNextRow() == 1)
    {
      datum.timestamp.setFromMJDDbl(double(mjd)*1e-8);
      data.push_back(datum);
    }

  delete stmt;
  delete my_db;

  return true;
}

bool VSCentralizedDBAccess::
doFetchFADC(const VSTime& time,
	    std::vector<FADCChannelSettingsDatum>& chan_settings,
	    std::vector<FADCBoardSettingsDatum>& board_settings,
	    std::vector<FADCChannelRelation>& chan_relation,
	    std::vector<FADCSlotRelation>& slot_relation,
	    VSDatabase* db)
{
  chan_settings.clear();
  board_settings.clear();
  chan_relation.clear();
  slot_relation.clear();

  VSDatabase* my_db = constructDBIfRequired(db);

  const uint64_t dbtime = time.getDBTimeStamp();
      
  std::ostringstream query_cond;
  query_cond << "( db_start_time <= " << dbtime
	     << " AND ( db_end_time IS NULL OR db_end_time > " << dbtime 
	     << " ) )";
  
  VSDBStatement* stmt;
  int nresults;

  // FADC Channel Settings ---------------------------------------------------

  stmt=
    db->createSelectQuery("tblFADC_Channel_Settings", query_cond.str(),
			  "fadc_id,fadc_channel,LookbackTime,AreaWidth,"
			  "AreaOffset,HiLoOffset,RereadOffset,PedOffset,"
			  "IsEnabled,ZeroSupThresh,delay");
  vsassert(stmt);
  nresults = stmt->execute();
  vsassert(nresults >= 0);

  FADCChannelSettingsDatum cs_datum;
  std::string is_enabled;

  stmt->bindToResult(cs_datum.fadc_id);
  stmt->bindToResult(cs_datum.fadc_channel);
  stmt->bindToResult(cs_datum.lookback_time);
  stmt->bindToResult(cs_datum.area_width);
  stmt->bindToResult(cs_datum.area_offset);
  stmt->bindToResult(cs_datum.hi_lo_offset);
  stmt->bindToResult(cs_datum.reread_offset);
  stmt->bindToResult(cs_datum.ped_offset);
  stmt->bindToResult(is_enabled);
  stmt->bindToResult(cs_datum.zero_sup_thresh);
  stmt->bindToResult(cs_datum.delay);

  while(stmt->retrieveNextRow() == 1)
    {
      if(is_enabled=="on")cs_datum.is_enabled=true;
      else cs_datum.is_enabled=false;
      chan_settings.push_back(cs_datum);
    }

  delete stmt;

  // FADC Board Settings ------------------------------------------------------

  stmt=
    db->createSelectQuery("tblFADC_Board_Settings", query_cond.str(),
			  "fadc_id,Mode,EventCounter,DWordsPerChan,IsEnabled,"
			  "HiLoWindowSize,RereadWindowSize,PedWindowSize");
  vsassert(stmt);
  nresults = stmt->execute();
  vsassert(nresults >= 0);

  FADCBoardSettingsDatum bs_datum;
  std::string mode;

  stmt->bindToResult(bs_datum.fadc_id);
  stmt->bindToResult(mode);
  stmt->bindToResult(bs_datum.event_counter);
  stmt->bindToResult(bs_datum.dwords_per_chan);
  stmt->bindToResult(is_enabled);
  stmt->bindToResult(bs_datum.hi_lo_window_size);
  stmt->bindToResult(bs_datum.reread_window_size);
  stmt->bindToResult(bs_datum.ped_window_size);

  while(stmt->retrieveNextRow() == 1)
    {
      if(is_enabled=="on")bs_datum.is_enabled=true;
      else bs_datum.is_enabled=false;
      if(mode=="FADC")
	bs_datum.mode=FADCBoardSettingsDatum::M_FADC;
      else if(mode=="FULL")
	bs_datum.mode=FADCBoardSettingsDatum::M_FULL;
      else if(mode=="QADC")
	bs_datum.mode=FADCBoardSettingsDatum::M_QADC;
      else if(mode=="WORD")
	bs_datum.mode=FADCBoardSettingsDatum::M_WORD;
      else vsassert(0);
      board_settings.push_back(bs_datum);
    }

  delete stmt;

  // FADC Channel Relation ----------------------------------------------------

  stmt=
    db->createSelectQuery("tblFADC_Channel_Relation", query_cond.str(),
			  "telescope_id,fadc_crate,fadc_slot,fadc_channel,"
			  "pixel_id,chan_type,channel_id");
  vsassert(stmt);
  nresults = stmt->execute();
  vsassert(nresults >= 0);

  FADCChannelRelation cr_datum;
  std::string chan_type;

  stmt->bindToResult(cr_datum.telescope_id);
  stmt->bindToResult(cr_datum.fadc_crate);
  stmt->bindToResult(cr_datum.fadc_slot);
  stmt->bindToResult(cr_datum.fadc_channel);
  stmt->bindToResult(cr_datum.pixel_id);
  stmt->bindToResult(chan_type);
  stmt->bindToResult(cr_datum.channel_id);

  while(stmt->retrieveNextRow() == 1)
    {
      if(chan_type == "PMT")
	cr_datum.chan_type=FADCChannelRelation::T_PMT;
      else if(chan_type == "L2")
	cr_datum.chan_type=FADCChannelRelation::T_L2;
      else if(chan_type == "LASER")
	cr_datum.chan_type=FADCChannelRelation::T_LASER;
      else if(chan_type == "EMPTY")
	cr_datum.chan_type=FADCChannelRelation::T_EMPTY;
      else vsassert(0);
      chan_relation.push_back(cr_datum);
    }

  delete stmt;

  // FADC Slot Relation -------------------------------------------------------

  stmt=
    db->createSelectQuery("tblFADC_Slot_Relation", query_cond.str(),
			  "telescope_id,fadc_crate,fadc_slot,fadc_id");
  vsassert(stmt);
  nresults = stmt->execute();
  vsassert(nresults >= 0);

  FADCSlotRelation sr_datum;

  stmt->bindToResult(sr_datum.telescope_id);
  stmt->bindToResult(sr_datum.fadc_crate);
  stmt->bindToResult(sr_datum.fadc_slot);
  stmt->bindToResult(sr_datum.fadc_id);

  while(stmt->retrieveNextRow() == 1)
    slot_relation.push_back(sr_datum);

  delete stmt;
  delete my_db;

  return true;  
}

bool VSCentralizedDBAccess::
doFetchTargetTable(std::vector<TargetTableCoord>& target_table,
		   VSDatabase* db)
{
  target_table.clear();

  VSDatabase* my_db = constructDBIfRequired(db);
  VSDBStatement* stmt = 
    db->createSelectQuery("tblObserving_Sources",
			  "source_id NOT IN ( SELECT source_id FROM "
			  "tblObserving_Collection WHERE "
			  "collection_id='yale_bright_star' )",
			  "source_id,ra,decl,epoch");

  vsassert(stmt);
  vsassert(stmt->execute() >= 0);

  TargetTableCoord ttc_datum;

  stmt->bindToResult(ttc_datum.name);
  stmt->bindToResult(ttc_datum.ra_rad);
  stmt->bindToResult(ttc_datum.dec_rad);
  stmt->bindToResult(ttc_datum.epoch);

  while(stmt->retrieveNextRow() == 1)
    target_table.push_back(ttc_datum);

  delete stmt;
  delete my_db;

  return true;  
}

void VSCentralizedDBAccess::wait()
{
#ifndef NOTHREADS
  if(!m_threaded)return;
  vsassert(pthread_mutex_lock(&m_mutex) == 0);
  while(m_threaded_active)vsassert(pthread_cond_wait(&m_cond, &m_mutex) == 0);
  vsassert(pthread_mutex_unlock(&m_mutex) == 0);
#endif
}

#ifndef NOTHREADS
void* VSCentralizedDBAccess::threadStart(void *object)
{
  static_cast<VSCentralizedDBAccess*>(object)->run();
  return 0;
}

void VSCentralizedDBAccess::run()
{
  doPrefetch(m_run_no, m_start_time, m_telescopes,false);
  vsassert(pthread_mutex_lock(&m_mutex) == 0);
  m_threaded_active=false;
  vsassert(pthread_mutex_unlock(&m_mutex) == 0);
  vsassert(pthread_cond_broadcast(&m_cond) == 0);
}
#endif
