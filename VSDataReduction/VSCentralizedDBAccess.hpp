//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCentralizedDBAccess.hpp

  Centralize all access to the database and provide asynchnonous data
  retrieval through thread

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       12/12/2006

  $Id: VSCentralizedDBAccess.hpp,v 3.5 2008/04/14 22:29:59 sfegan Exp $

*/

#ifndef VSCENTRALIZEDDBACCESS_HPP
#define VSCENTRALIZEDDBACCESS_HPP

#include<memory>

#ifndef NOTHREADS
#include<pthread.h>
#endif

#include<VSDBFactory.hpp>
#include<VSPointing.hpp>
#include<VSCentralizedDBAccessData.hpp>

namespace VERITAS
{

  class VSCentralizedDBAccess
  {
  public:
    ~VSCentralizedDBAccess();

    typedef VSDBPointingData::Datum PointingDatum;

    typedef VSCentralizedDBAccessData::RunInfoDatum RunInfoDatum;
    typedef VSCentralizedDBAccessData::TrackingTargetDatum TrackingTargetDatum;
    typedef VSCentralizedDBAccessData::CorrectionParametersDatum CorrectionParametersDatum;
    typedef VSCentralizedDBAccessData::VPMCentroids VPMCentroids;
    typedef VSCentralizedDBAccessData::VPMLEDs VPMLEDs;
    typedef VSCentralizedDBAccessData::PowerDatum PowerDatum;
    typedef VSCentralizedDBAccessData::HVStatusDatum HVStatusDatum;
    typedef VSCentralizedDBAccessData::L1RateDatum L1RateDatum;
    typedef VSCentralizedDBAccessData::FIRDatum FIRDatum;
    typedef VSCentralizedDBAccessData::L3ScopeDatum L3ScopeDatum;
    typedef VSCentralizedDBAccessData::FADCChannelSettingsDatum FADCChannelSettingsDatum;
    typedef VSCentralizedDBAccessData::FADCBoardSettingsDatum FADCBoardSettingsDatum;
    typedef VSCentralizedDBAccessData::FADCChannelRelation FADCChannelRelation;
    typedef VSCentralizedDBAccessData::FADCSlotRelation FADCSlotRelation;
    typedef VSCentralizedDBAccessData::TargetTableCoord TargetTableCoord;

    void prefetchFromDB(unsigned run_no, const VSTime& start_time,
			const std::vector<unsigned>& telescopes, 
			bool no_threads = false);

    bool getRunInfo(unsigned run_no, RunInfoDatum& datum);

    bool getTrackingTarget(unsigned iscope, 
			   const VSTime& start, const VSTime& stop,
			   std::vector<TrackingTargetDatum>& data);

    bool getCorrectionParameters(unsigned iscope, 
				 std::vector<CorrectionParametersDatum>& data);

    bool getPointing(unsigned iscope, const VSTime& start, const VSTime& stop,
		     std::vector<PointingDatum>& data);

    bool getPower(unsigned iscope, 
		  const VSTime& start, const VSTime& stop,
		  std::vector<PowerDatum>& data);

    bool getHVStatus(unsigned iscope, 
		     const VSTime& start, const VSTime& stop,
		     std::vector<HVStatusDatum>& data);

    bool getL1Rate(unsigned iscope, 
		   const VSTime& start, const VSTime& stop,
		   std::vector<L1RateDatum>& data);

    bool getFIR(unsigned iscope,
		const VSTime& start, const VSTime& stop,
		std::vector<FIRDatum>& data);

    bool getL3Scope(unsigned iscope, 
		    const VSTime& start, const VSTime& stop,
		    std::vector<L3ScopeDatum>& data);

    bool getFADC(const VSTime& time,
		 std::vector<FADCChannelSettingsDatum>& chan_settings,
		 std::vector<FADCBoardSettingsDatum>& board_settings,
		 std::vector<FADCChannelRelation>& chan_relation,
		 std::vector<FADCSlotRelation>& slot_relation);

    bool getTargetTable(std::vector<TargetTableCoord>& target_table);

    static VSCentralizedDBAccess* getInstance();

  protected:
    VSCentralizedDBAccess();

  private:
    VSCentralizedDBAccess(const VSCentralizedDBAccess&);
    VSCentralizedDBAccess& operator=(const VSCentralizedDBAccess&);

    VSDatabase* constructDBIfRequired(VSDatabase*& supplied_db);

    void doPrefetch(unsigned run_no, const VSTime& start_time,
		    const std::vector<unsigned>& telescopes, 
		    bool no_real_fetch = false);

    bool doFetchRunInfo(unsigned run_no, RunInfoDatum& datum,
			VSDatabase* db = 0);

    bool doFetchTrackingTarget(unsigned iscope, 
			       const VSTime& start, const VSTime& stop,
			       std::vector<TrackingTargetDatum>& data,
			       VSDatabase* db = 0);

    bool doFetchCorrectionParameters(unsigned iscope, 
				 std::vector<CorrectionParametersDatum>& data,
				     VSDatabase* db = 0);

    bool doFetchPointing(unsigned iscope, 
			 const VSTime& start, const VSTime& stop,
			 std::vector<PointingDatum>& data, VSDatabase* db = 0);

    bool doFetchPower(unsigned iscope, 
		      const VSTime& start, const VSTime& stop,
		      std::vector<PowerDatum>& data, VSDatabase* db = 0);

    bool doFetchHVStatus(unsigned iscope, 
			 const VSTime& start, const VSTime& stop,
			 std::vector<HVStatusDatum>& data, VSDatabase* db = 0);

    bool doFetchL1Rate(unsigned iscope, 
		       const VSTime& start, const VSTime& stop,
		       std::vector<L1RateDatum>& data, VSDatabase* db = 0);

    bool doFetchFIR(unsigned iscope, 
		    const VSTime& start, const VSTime& stop,
		    std::vector<FIRDatum>& data, VSDatabase* db = 0);

    bool doFetchL3Scope(unsigned iscope, 
			const VSTime& start, const VSTime& stop,
			std::vector<L3ScopeDatum>& data, VSDatabase* db = 0);

    bool doFetchVPMStarsScope(unsigned iscope, 
			      const VSTime& start, const VSTime& stop,
			      std::vector<VPMCentroids>& data, 
			      VSDatabase* db = 0);

    bool doFetchVPMLEDsScope(unsigned iscope, 
			     const VSTime& start, const VSTime& stop,
			     std::vector<VPMLEDs>& data, 
			     VSDatabase* db = 0);

    bool doFetchFADC(const VSTime& time,
		     std::vector<FADCChannelSettingsDatum>& chan_settings,
		     std::vector<FADCBoardSettingsDatum>& board_settings,
		     std::vector<FADCChannelRelation>& chan_relation,
		     std::vector<FADCSlotRelation>& slot_relation,
		     VSDatabase* db = 0);

    bool doFetchTargetTable(std::vector<TargetTableCoord>& target_table,
			    VSDatabase* db = 0);

    void wait();

    template<typename T> struct TSCache
    {
      TSCache(): start_time(), stop_time(), scope_num(), data()
      { /* nothing to see here */ }
      VSTime start_time;
      VSTime stop_time;
      unsigned scope_num;
      std::vector<T> data;
    };

    template<typename T> struct ACCache
    {
      ACCache(): config_time(), data()
      { /* nothing to see here */ }
      VSTime config_time;
      std::vector<T> data;
    };

    static std::auto_ptr<VSCentralizedDBAccess>      s_instance;

    std::ostream*                                    m_stream;

    RunInfoDatum                                     m_run_info_cache;
    std::vector<TSCache<TrackingTargetDatum> >       m_tracking_target_cache;
    std::vector<std::pair<bool, std::vector<CorrectionParametersDatum> > >
                                                     m_correction_parameters;
    std::vector<TSCache<PointingDatum> >             m_pointing_cache;
    std::vector<TSCache<PowerDatum> >                m_power_cache;
    std::vector<TSCache<HVStatusDatum> >             m_hv_status_cache;
    std::vector<TSCache<L1RateDatum> >               m_l1_rate_cache;
    std::vector<TSCache<FIRDatum> >                  m_fir_cache;
    std::vector<TSCache<L3ScopeDatum> >              m_l3_scope_cache;
    std::vector<TSCache<VPMCentroids> >              m_vpm_stars_scope_cache;
    std::vector<TSCache<VPMLEDs> >                   m_vpm_leds_scope_cache;
    ACCache<FADCChannelSettingsDatum>                m_fadc_chan_settings;
    ACCache<FADCBoardSettingsDatum>                  m_fadc_board_settings;
    ACCache<FADCChannelRelation>                     m_fadc_chan_relation;
    ACCache<FADCSlotRelation>                        m_fadc_slot_relation;
    std::pair<bool, std::vector<TargetTableCoord> >  m_target_table;

#ifndef NOTHREADS
    static void* threadStart(void *object);
    void run();

    unsigned                                         m_run_no;
    VSTime                                           m_start_time;
    std::vector<unsigned>                            m_telescopes;

    bool                                             m_threaded;
    bool                                             m_threaded_active;
    pthread_mutex_t                                  m_mutex;
    pthread_cond_t                                   m_cond;
    pthread_t                                        m_thread;
#endif
  };

}

#endif // defined VSCENTRALIZEDDBACCESS_HPP
