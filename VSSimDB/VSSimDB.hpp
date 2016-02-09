//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimDB.hpp

  Classes for input and output to simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       27/02/2005
  \note
*/

#ifndef VSSIMDB_HPP
#define VSSIMDB_HPP

#include <string>
#include <map>
#include <memory>
#include <stdint.h>

#include <VSDatabase.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSSimDBPEData
  {
  public:
    VSSimDBPEData():
      fEventID(), fScopeID(), fPixelID(), fPixelTimeNS() 
    { /* nothing to see here */ }

    uint32_t    fEventID;
    uint16_t    fScopeID;
    uint32_t    fPixelID;
    float       fPixelTimeNS;
  };

  class VSSimDBScopeData
  {
  public:
    VSSimDBScopeData():
      fEventID(), fScopeID(), fScopeZenithRad(), fScopeAzimuthRad(), 
      fNumHitPixels() 
    { /* nothing to see here */ }

    uint32_t    fEventID;
    uint16_t    fScopeID;
    float       fScopeZenithRad;
    float       fScopeAzimuthRad;
    uint32_t    fNumHitPixels;
  };

  class VSSimDBEventData
  {
  public:
    VSSimDBEventData(): 
      fEventID(), fWorkunitRunID(), fEnergyGeV(),
      fTargetZenithRad(), fTargetAzimuthRad(),
      fPrimaryZenithRad(), fPrimaryAzimuthRad(), 
      fPrimaryCoreEastM(), fPrimaryCoreNorthM(), fPrimaryCoreUpASLM(),
      fNumHitScopes(), fEventComplete()
    { /* nothing to see here */ }

    uint32_t    fEventID;
    uint32_t    fWorkunitRunID;
    float       fEnergyGeV;
    float       fTargetZenithRad;
    float       fTargetAzimuthRad;
    float       fPrimaryZenithRad;
    float       fPrimaryAzimuthRad;
    float       fPrimaryCoreEastM;
    float       fPrimaryCoreNorthM;
    float       fPrimaryCoreUpASLM;
    uint16_t    fNumHitScopes;
    bool        fEventComplete;
  };

  class VSSimDBTableParam
  {
  public:
    VSSimDBTableParam():
      fTableID(), fPrimaryID(), fEnergyGeV(),
      fZenithMinRad(), fZenithMaxRad(), fAzimuthMinRad(), fAzimuthMaxRad(),
      fOpticsID(), fSamplingRadiusM(), fTargetEventCount(), 
      fWorkunitEventCount(), fWorkunitPriority(), fTableName()
    { /* nothing to see here */ }

    uint32_t    fTableID;
    uint32_t    fPrimaryID;
    float       fEnergyGeV;
    float       fZenithMinRad;
    float       fZenithMaxRad;
    float       fAzimuthMinRad;
    float       fAzimuthMaxRad;
    uint32_t    fOpticsID;
    float       fSamplingRadiusM;
    uint32_t    fTargetEventCount;
    uint32_t    fWorkunitEventCount;
    int32_t     fWorkunitPriority;
    std::string fTableName;
  };

  class VSSimDBWorkunitRun
  {
  public:
    VSSimDBWorkunitRun():
      fWorkunitRunID(), fTableID(), fHostName(), fJobID(), fStartTime(),
      fFinishTime(), fWorkunitRunComplete()
    { /* nothing to see here */ }

    uint64_t     fWorkunitRunID;
    uint32_t     fTableID;
    std::string  fHostName;
    uint32_t     fJobID;
    VSDBDateTime fStartTime;
    VSDBDateTime fFinishTime;
    bool         fWorkunitRunComplete;
  };

  class VSSimDBEventStore
  {
  public:
    virtual ~VSSimDBEventStore();
    virtual int insertEvent(VSSimDBEventData& event) = 0;
    virtual int insertScope(const VSSimDBScopeData& scope) = 0;
    virtual int insertPE(const VSSimDBPEData& pe) = 0;
    virtual int updateHitScopes(const VSSimDBEventData& event) = 0;
    virtual int updateHitPixels(const VSSimDBScopeData& scope) = 0;
  };

  class VSSimDB: public VSSimDBEventStore
  {
  public:
    VSSimDB(VSDatabase* db): 
      VSSimDBEventStore(), fDB(db), 
      fStmtInsertEvent(), fStmtInsertScope(), fStmtInsertPE(),
      fStmtUpdateHitScopes(), fStmtUpdateHitPixels(), fStmtRetrieval(),
      fInsertBoundEvent(), fInsertBoundScope(), fInsertBoundPE(),
      fUpdateBoundEvent(), fUpdateBoundScope() { }
    virtual ~VSSimDB();

    std::string tableName(const VSSimDBTableParam& param);
    VSDatabase* db() const { return fDB; }

    // ------------------------------------------------------------------------
    // DB structure creation / query
    // ------------------------------------------------------------------------

    void createInfrastructureTables();
    void createDataTable(VSSimDBTableParam& table, bool event_only = false);

    VSSimDBTableParam* getDataTableByName(const std::string& prefix);
    std::vector<VSSimDBTableParam> getAllDataTables();    

    uint32_t getEventCountByName(const std::string& prefix,
				 bool include_incomplete=false);
    uint32_t getMaxEventNumByName(const std::string& prefix);
    uint32_t getCompleteEventCountByName(const std::string& prefix,
					 uint32_t max_event_num = 0);
    uint32_t getTableIDByName(const std::string& tablename);

    std::vector<uint32_t>
    getEventCountOfAllTables(const std::vector<VSSimDBTableParam>& param,
			     bool include_incomplete=false);

    uint64_t registerWorkunitRunStart(uint32_t table_id,
				      const std::string& host_name,
				      uint32_t job_id);
    void registerWorkunitRunFinish(uint32_t workunit_run_id);

    // ------------------------------------------------------------------------
    // Data retrieval
    // ------------------------------------------------------------------------

    int getAllCompleteEvents(const std::string& tablename, 
			     std::vector<VSSimDBEventData>& events);
    int getLimitedCompleteEventNums(const std::string& tablename,
				    uint32_t min_event_num,
				    uint32_t max_event_num,
				    std::vector<uint32_t>& event_nums);
    int getEventByNum(const std::string& tablename,  uint32_t event_id,
		      VSSimDBEventData& event);
    int getAllScopeByNum(const std::string& tablename,  uint32_t event_id,
			 std::vector<VSSimDBScopeData>& scopes);
    int getAllPEByNum(const std::string& tablename, uint32_t event_id,
		      std::vector<VSSimDBPEData>& pes);
    int getWorkunitIDs(const std::string& tablename, 
		       std::vector<unsigned>& workunit_ids);

    // ------------------------------------------------------------------------
    // Data insertion
    // ------------------------------------------------------------------------

    void setInsertTableName(const std::string& prefix, bool event_only=false);
    virtual int insertEvent(VSSimDBEventData& event);
    virtual int insertScope(const VSSimDBScopeData& scope);
    virtual int insertPE(const VSSimDBPEData& pe);
    virtual int updateHitScopes(const VSSimDBEventData& event);
    virtual int updateHitPixels(const VSSimDBScopeData& scope);

  private:
    VSSimDB(const VSSimDB&);
    VSSimDB& operator= (const VSSimDB&);

    struct RetrievalStatements
    {
      RetrievalStatements():
	all_complete_events(), limited_complete_eventnums(), event_by_num(), 
	all_scope_by_num(), all_pe_by_num() { }

      RetrievalStatements(const RetrievalStatements& o):
	all_complete_events(o.all_complete_events), 
	limited_complete_eventnums(o.limited_complete_eventnums), 
	event_by_num(o.event_by_num), 
	all_scope_by_num(o.all_scope_by_num), 
	all_pe_by_num(o.all_pe_by_num) { }

      RetrievalStatements& operator=(const RetrievalStatements& o)
      {
	all_complete_events        = o.all_complete_events;
	limited_complete_eventnums = o.limited_complete_eventnums;
	event_by_num               = o.event_by_num;
	all_scope_by_num           = o.all_scope_by_num;
	all_pe_by_num              = o.all_pe_by_num;
	return *this;
      }

      VSDBStatement* all_complete_events;
      VSDBStatement* limited_complete_eventnums;
      VSDBStatement* event_by_num;
      VSDBStatement* all_scope_by_num;
      VSDBStatement* all_pe_by_num;
    };

    VSDatabase*                    fDB;

    std::auto_ptr<VSDBStatement>   fStmtInsertEvent; 
    std::auto_ptr<VSDBStatement>   fStmtInsertScope; 
    std::auto_ptr<VSDBStatement>   fStmtInsertPE; 
    std::auto_ptr<VSDBStatement>   fStmtUpdateHitScopes; 
    std::auto_ptr<VSDBStatement>   fStmtUpdateHitPixels; 

    std::map<std::string,RetrievalStatements>  fStmtRetrieval;

    VSSimDBEventData*              fInsertBoundEvent;
    const VSSimDBScopeData*        fInsertBoundScope;
    const VSSimDBPEData*           fInsertBoundPE;
    const VSSimDBEventData*        fUpdateBoundEvent;
    const VSSimDBScopeData*        fUpdateBoundScope;
  };

}

#endif // VSSIMDB_HPP
