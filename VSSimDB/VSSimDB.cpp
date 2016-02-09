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

#include <cmath>

#include <VSDBMySQLBase.hpp>
#include <VSDBParameterTable.hpp>

#include "VSSimDB.hpp"
#include "VSSimDBTables.hpp"

using namespace VERITAS;

VSSimDBEventStore::~VSSimDBEventStore()
{
  // nothing to see here
}

VSSimDB::~VSSimDB()
{
  for(std::map<std::string,RetrievalStatements>::iterator istmt =
	fStmtRetrieval.begin(); istmt != fStmtRetrieval.end(); istmt++)
    {
      delete istmt->second.all_complete_events;
      delete istmt->second.limited_complete_eventnums;
      delete istmt->second.event_by_num;
      delete istmt->second.all_scope_by_num;
      delete istmt->second.all_pe_by_num;
    }
}

std::string VSSimDB::tableName(const VSSimDBTableParam& param)
{
  std::ostringstream stream;
      
  stream << "P";
  switch(param.fPrimaryID)
    {
    case 1:    stream << "g";  break;
    case 3:    stream << "e";  break;
    case 6:    stream << "m";  break;
    case 13:   stream << "n";  break;
    case 14:   stream << "p";  break;
    case 402:  stream << "He"; break;
    case 5626: stream << "Fe"; break;
    default:   stream << param.fPrimaryID; break;
    }
  
  unsigned e = unsigned(round(param.fEnergyGeV*10));

  stream << "_En" 
	 << std::setfill('0') << std::setw(6) << std::setprecision(5)
	 << e/10 << 'd' << e%10

	 << "_Zn"
	 << std::setfill('0') << std::setw(2) << std::setprecision(2)
	 << unsigned(round(param.fZenithMinRad/M_PI*180.0))
	 << std::setfill('0') << std::setw(2) << std::setprecision(2)
	 << unsigned(round(param.fZenithMaxRad/M_PI*180.0))

	 << "_Az"
	 << std::setfill('0') << std::setw(3) << std::setprecision(3)
	 << unsigned(round(fmod(fmod(param.fAzimuthMinRad/M_PI*180.0,360)+360,360)))
	 << std::setfill('0') << std::setw(3) << std::setprecision(3)
	 << unsigned(round(fmod(fmod(param.fAzimuthMaxRad/M_PI*180.0,360)+360,360)))

	 << "_Opt" << param.fOpticsID;
  
  return stream.str();
}

// ============================================================================
// DB structure creation / query
// ============================================================================

void VSSimDB::createInfrastructureTables()
{
  VSDBParameterTable parameter_table(fDB);
  parameter_table.createParameterTable();

  // stupid to create this --- C++ should have a "typeof" operator
  // (g++ does but it is not portable)
  VSSimDBTableParam temp; 

  fDB->createTable(VSIMDB_TABLE_NAME_DIRECTORY,

		   fDB->sqlSpecOf("TableID",
				  temp.fTableID,true,
				  "NOT NULL AUTO_INCREMENT")+

		   fDB->sqlSpecOf("PrimaryID",
				  temp.fPrimaryID,false,"NOT NULL")+
		   fDB->sqlSpecOf("EnergyGeV",
				  temp.fEnergyGeV,false,"NOT NULL")+
		   fDB->sqlSpecOf("ZenithMinRad",
				  temp.fZenithMinRad,false,"NOT NULL")+
		   fDB->sqlSpecOf("ZenithMaxRad",
				  temp.fZenithMaxRad,false,"NOT NULL")+
		   fDB->sqlSpecOf("AzimuthMinRad",
				  temp.fAzimuthMinRad,false,"NOT NULL")+

		   fDB->sqlSpecOf("AzimuthMaxRad",
				  temp.fAzimuthMaxRad,false,"NOT NULL")+
		   fDB->sqlSpecOf("OpticsID",
				  temp.fOpticsID,false,"NOT NULL")+
		   fDB->sqlSpecOf("SamplingRadiusM",
				  temp.fSamplingRadiusM,false,"NOT NULL")+
		   fDB->sqlSpecOf("TargetEventCount",
				  temp.fTargetEventCount,false,"NOT NULL")+
		   fDB->sqlSpecOf("WorkunitEventCount",
				  temp.fWorkunitEventCount,false,"NOT NULL")+

		   fDB->sqlSpecOf("WorkunitPriority",
				  temp.fWorkunitPriority,false,"NOT NULL")+
 		   fDB->sqlSpecOf("TableName",
				  temp.fTableName,false,"NOT NULL")+

		   ", PRIMARY KEY ( TableID ), KEY ( TableName )",

		   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);  


  VSSimDBWorkunitRun temp2;

  fDB->createTable(VSIMDB_TABLE_NAME_WORKUNIT,

		   fDB->sqlSpecOf("WorkunitRunID",
				  temp2.fWorkunitRunID,true,
				  "NOT NULL AUTO_INCREMENT")+
		   fDB->sqlSpecOf("TableID",
				  temp2.fTableID,false,"NOT NULL")+
		   fDB->sqlSpecOf("HostName",
				  temp2.fHostName,false,"NOT NULL")+
		   fDB->sqlSpecOf("JobID",
				  temp2.fJobID,false,"NOT NULL")+
		   fDB->sqlSpecOf("StartTime",
				  temp2.fStartTime,false,"NOT NULL")+

		   fDB->sqlSpecOf("FinishTime",
				  temp2.fFinishTime,false)+
		   fDB->sqlSpecOf("WorkunitRunComplete",
				  temp2.fWorkunitRunComplete,false,
				  "NOT NULL DEFAULT 0")+
		   
		   ", Primary Key ( WorkunitRunID )");
}

void VSSimDB::createDataTable(VSSimDBTableParam& table, bool event_only)
{
  VSDBStatement* stmt = 
    fDB->createInsertQuery(VSIMDB_TABLE_NAME_DIRECTORY,12,
			   "PrimaryID,EnergyGeV,ZenithMinRad,ZenithMaxRad,"
			   "AzimuthMinRad,AzimuthMaxRad,OpticsID,"
			   "SamplingRadiusM,TargetEventCount,"
			   "WorkunitEventCount,WorkunitPriority,TableName",
			   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);

  stmt->bindToParam(table.fPrimaryID);
  stmt->bindToParam(table.fEnergyGeV);
  stmt->bindToParam(table.fZenithMinRad);
  stmt->bindToParam(table.fZenithMaxRad);
  stmt->bindToParam(table.fAzimuthMinRad);

  stmt->bindToParam(table.fAzimuthMaxRad);
  stmt->bindToParam(table.fOpticsID);
  stmt->bindToParam(table.fSamplingRadiusM);
  stmt->bindToParam(table.fTargetEventCount);
  stmt->bindToParam(table.fWorkunitEventCount);

  stmt->bindToParam(table.fWorkunitPriority);
  stmt->bindToParam(table.fTableName);

  unsigned result = stmt->execute(); 
  if(result)table.fTableID = stmt->getInsertID();
  delete stmt;

  if(1 /* result */)
    {
      std::string table_name = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ table.fTableName 
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

      // stupid to create this --- C++ should have a "typeof" operator
      // (g++ does but it is not portable)
      VSSimDBEventData temp;

      fDB->createTable(table_name,

		       fDB->sqlSpecOf("EventID",
			       temp.fEventID,true,"NOT NULL AUTO_INCREMENT")+
		       fDB->sqlSpecOf("WorkunitRunID",
			       temp.fWorkunitRunID,false,"NOT NULL")+
		       fDB->sqlSpecOf("TargetZenithRad",
			       temp.fTargetZenithRad,false,"NOT NULL")+
 		       fDB->sqlSpecOf("TargetAzimuthRad",
			       temp.fTargetAzimuthRad,false,"NOT NULL")+
		       fDB->sqlSpecOf("PrimaryZenithRad",
			       temp.fPrimaryZenithRad,false,"NOT NULL")+
 		       fDB->sqlSpecOf("PrimaryAzimuthRad",
			       temp.fPrimaryAzimuthRad,false,"NOT NULL")+

  		       fDB->sqlSpecOf("PrimaryCoreEastM",
			       temp.fPrimaryCoreEastM,false,"NOT NULL")+
  		       fDB->sqlSpecOf("PrimaryCoreNorthM",
			       temp.fPrimaryCoreNorthM,false,"NOT NULL")+
 		       fDB->sqlSpecOf("PrimaryCoreUpASLM",
			       temp.fPrimaryCoreUpASLM,false,"NOT NULL")+
  		       fDB->sqlSpecOf("NumHitScopes",
			       temp.fNumHitScopes,false,"NOT NULL DEFAULT 0")+
  		       fDB->sqlSpecOf("EventComplete",
			       temp.fEventComplete,false,"NOT NULL DEFAULT 0")+

		       ", PRIMARY KEY ( EventID ), INDEX ( EventComplete )",

		       VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);  
    }     

  if(!event_only)
    {
      std::string table_name = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ table.fTableName 
	+ std::string(VSIMDB_TABLE_POSTFIX_SCOPES);

      // stupid to create this --- C++ should have a "typeof" operator
      // (g++ does but it is not portable)
      VSSimDBScopeData temp;
      
      fDB->createTable(table_name,

		       fDB->sqlSpecOf("EventID",
				temp.fEventID,true,"NOT NULL")+
		       fDB->sqlSpecOf("ScopeID",
				temp.fScopeID,false,"NOT NULL")+
		       fDB->sqlSpecOf("ScopeZenithRad",
				temp.fScopeZenithRad,false,"NOT NULL")+
		       fDB->sqlSpecOf("ScopeAzimuthRad",
				temp.fScopeAzimuthRad,false,"NOT NULL")+
		       fDB->sqlSpecOf("NumHitPixels",
				temp.fNumHitPixels,false,"NOT NULL DEFAULT 0")+

		       ", PRIMARY KEY ( EventID, ScopeID )",

		       VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
    }

  if(!event_only)
    {
      std::string table_name = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ table.fTableName 
	+ std::string(VSIMDB_TABLE_POSTFIX_PES);

      // stupid to create this --- C++ should have a "typeof" operator
      // (g++ does but it is not portable)
      VSSimDBPEData temp;

      fDB->createTable(table_name,

		       fDB->sqlSpecOf("EventID",
				      temp.fEventID,true,"NOT NULL")+
		       fDB->sqlSpecOf("ScopeID",
				      temp.fScopeID,false,"NOT NULL")+
		       fDB->sqlSpecOf("PixelID",
				      temp.fPixelID,false,"NOT NULL")+
		       fDB->sqlSpecOf("PixelTimeNS",
				      temp.fPixelTimeNS,false,"NOT NULL")+

		       ", INDEX ( EventID, ScopeID, PixelID )",

		       VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
    }
}

VSSimDBTableParam* VSSimDB::getDataTableByName(const std::string& prefix)
{
  std::string cond = std::string("TableName=") + MySQLValueToString(prefix);

  VSDBStatement* stmt = 
    fDB->createSelectQuery(VSIMDB_TABLE_NAME_DIRECTORY, cond, "*",
			   VSDatabase::FLAG_NO_SERVER_PS);

  assert(stmt);

  int result;
  result = stmt->execute();
  assert(result>=0);
  assert(result<2);
  if(result == 0){ return 0; }
  
  VSSimDBTableParam* param = new VSSimDBTableParam;

  stmt->bindToResult(param->fTableID);

  stmt->bindToResult(param->fPrimaryID);
  stmt->bindToResult(param->fEnergyGeV);
  stmt->bindToResult(param->fZenithMinRad);
  stmt->bindToResult(param->fZenithMaxRad);
  stmt->bindToResult(param->fAzimuthMinRad);

  stmt->bindToResult(param->fAzimuthMaxRad);
  stmt->bindToResult(param->fOpticsID);
  stmt->bindToResult(param->fSamplingRadiusM);
  stmt->bindToResult(param->fTargetEventCount);
  stmt->bindToResult(param->fWorkunitEventCount);

  stmt->bindToResult(param->fWorkunitPriority);
  stmt->bindToResult(param->fTableName);
  
  assert(stmt->retrieveNextRow() == 1);
  assert(stmt->retrieveNextRow() == 0);  

  delete(stmt);

  return param;
}

std::vector<VSSimDBTableParam> VSSimDB::getAllDataTables()
{
  std::vector<VSSimDBTableParam> params;

  VSDBStatement* stmt = 
    fDB->createSelectQuery(VSIMDB_TABLE_NAME_DIRECTORY, "", "*",
			   VSDatabase::FLAG_NO_SERVER_PS);

  assert(stmt);
  assert(stmt->execute());
  
  VSSimDBTableParam param;

  stmt->bindToResult(param.fTableID);

  stmt->bindToResult(param.fPrimaryID);
  stmt->bindToResult(param.fEnergyGeV);
  stmt->bindToResult(param.fZenithMinRad);
  stmt->bindToResult(param.fZenithMaxRad);
  stmt->bindToResult(param.fAzimuthMinRad);

  stmt->bindToResult(param.fAzimuthMaxRad);
  stmt->bindToResult(param.fOpticsID);
  stmt->bindToResult(param.fSamplingRadiusM);
  stmt->bindToResult(param.fTargetEventCount);
  stmt->bindToResult(param.fWorkunitEventCount);

  stmt->bindToResult(param.fWorkunitPriority);
  stmt->bindToResult(param.fTableName);
  
  while(stmt->retrieveNextRow() == 1)params.push_back(param);

  delete(stmt);

  return params;
}

uint32_t VSSimDB::getEventCountByName(const std::string& prefix,
				       bool include_incomplete)
{
  std::string table_name = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    + prefix
    + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

  std::string cond;
  if(!include_incomplete)cond="EventComplete=1";

  VSDBStatement* stmt = 
    fDB->createSelectQuery(table_name, cond, "COUNT(*)",
			   VSDatabase::FLAG_NO_SERVER_PS);
  assert(stmt);
  assert(stmt->execute());
  
  uint32_t count = 0;
  stmt->bindToResult(count);

  assert(stmt->retrieveNextRow() == 1);
  assert(stmt->retrieveNextRow() == 0);
  
  delete(stmt);

  return count;
}

uint32_t VSSimDB::getMaxEventNumByName(const std::string& prefix)
{
  std::string table_name = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    + prefix
    + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);


  VSDBStatement* stmt = 
    fDB->createSelectQuery(table_name, "", "MAX(EventID)",
			   VSDatabase::FLAG_NO_SERVER_PS);
  assert(stmt);
  assert(stmt->execute());
  
  uint32_t max_event_num = 0;
  stmt->bindToResult(max_event_num);

  assert(stmt->retrieveNextRow() == 1);
  assert(stmt->retrieveNextRow() == 0);
  
  delete(stmt);

  return max_event_num;
}

uint32_t VSSimDB::getTableIDByName(const std::string& prefix)
{
  std::string table_name = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    + prefix
    + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

  std::string cond = "TableName='" + prefix + "'";

  VSDBStatement* stmt = 
    fDB->createSelectQuery(VSIMDB_TABLE_NAME_DIRECTORY, cond, "TableID",
			   VSDatabase::FLAG_NO_SERVER_PS);

  vsassert(stmt);
  vsassert(stmt->execute());

  uint32_t table_id = 0;
  stmt->bindToResult(table_id);

  vsassert(stmt->retrieveNextRow() == 1);
  vsassert(stmt->retrieveNextRow() == 0);

  delete stmt;

  return table_id;
}

uint32_t VSSimDB::getCompleteEventCountByName(const std::string& prefix,
					      uint32_t max_event_num)
{
  std::string table_name = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    + prefix
    + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

  std::string cond; // = "EventComplete=1";
  if(max_event_num)cond += " AND EventID<=?";

  VSDBStatement* stmt = 
    fDB->createSelectQuery(table_name, cond, "COUNT(*)",
			   VSDatabase::FLAG_NO_SERVER_PS);

  if(max_event_num)stmt->bindToParam(max_event_num);
  assert(stmt);
  assert(stmt->execute());
  
  uint32_t count = 0;
  stmt->bindToResult(count);

  assert(stmt->retrieveNextRow() == 1);
  assert(stmt->retrieveNextRow() == 0);
  
  delete(stmt);

  return count;
}

std::vector<uint32_t> VSSimDB::
getEventCountOfAllTables(const std::vector<VSSimDBTableParam>& param,
			 bool include_incomplete)
{
  std::vector<uint32_t> count;
  for(std::vector<VSSimDBTableParam>::const_iterator iparam=param.begin();
      iparam!=param.end(); iparam++)
    {
      count.push_back(getEventCountByName(iparam->fTableName,
					  include_incomplete));
    }
  return count;
}

uint64_t VSSimDB::registerWorkunitRunStart(uint32_t table_id,
					   const std::string& host_name,
					   uint32_t job_id)
{
  VSDBStatement* stmt = 
    fDB->createQuery("INSERT INTO " VSIMDB_TABLE_NAME_WORKUNIT " "
		     "(TableID,HostName,JobID,StartTime) "
		     "VALUES (?,?,?,now())");
  assert(stmt);

  stmt->bindToParam(table_id);
  stmt->bindToParam(host_name);
  stmt->bindToParam(job_id);
 
  assert(stmt->execute());
 
  return stmt->getInsertID();
}

void VSSimDB::registerWorkunitRunFinish(uint32_t workunit_run_id)
{
  VSDBStatement* stmt = 
    fDB->createQuery("UPDATE " VSIMDB_TABLE_NAME_WORKUNIT 
		     " SET FinishTime=now(),WorkunitRunComplete=1"
		     " WHERE WorkunitRunID=?");
  assert(stmt);

  stmt->bindToParam(workunit_run_id);
 
  assert(stmt->execute());
}

// ----------------------------------------------------------------------------
// Data retrieval
// ----------------------------------------------------------------------------

int VSSimDB::getAllCompleteEvents(const std::string& tablename, 
				  std::vector<VSSimDBEventData>& events)
{
  VSDBStatement* stmt = fStmtRetrieval[tablename].all_complete_events;

  if(stmt == 0)
    {
      std::string tablename_ev = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);  

      stmt = fDB->createSelectQuery(tablename_ev, "EventComplete=1",
				    "EventID, "
				    "WorkunitRunID, "
				    "TargetZenithRad, "
				    "TargetAzimuthRad, "
				    "PrimaryZenithRad, "
				    "PrimaryAzimuthRad, "
				    "PrimaryCoreEastM, "
				    "PrimaryCoreNorthM, "
				    "PrimaryCoreUpASLM, "
				    "NumHitScopes, "
				    "EventComplete",
				    VSDatabase::FLAG_NO_BUFFER);
      
      if(stmt == 0)
	{
	  std::cerr << __PRETTY_FUNCTION__ << ": could not query table: " 
		    << tablename_ev << std::endl;
	  return EXIT_FAILURE;
	}

      fStmtRetrieval[tablename].all_complete_events = stmt;
    }

  VSSimDBEventData event;

  stmt->bindToResult(event.fEventID);
  stmt->bindToResult(event.fWorkunitRunID);
  stmt->bindToResult(event.fTargetZenithRad);
  stmt->bindToResult(event.fTargetAzimuthRad);
  stmt->bindToResult(event.fPrimaryZenithRad);
  stmt->bindToResult(event.fPrimaryAzimuthRad);
  stmt->bindToResult(event.fPrimaryCoreEastM);
  stmt->bindToResult(event.fPrimaryCoreNorthM);
  stmt->bindToResult(event.fPrimaryCoreUpASLM);
  stmt->bindToResult(event.fNumHitScopes);
  stmt->bindToResult(event.fEventComplete);

  int ret = stmt->execute();
  if(ret < 0)
    {
      std::string tablename_ev = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);  

      std::cerr << __PRETTY_FUNCTION__ << ": error executing query on: " 
		<< tablename_ev << std::endl;
      return ret;
    }

  events.clear();
  while(stmt->retrieveNextRow())events.push_back(event);

  return ret;
}

int VSSimDB::getLimitedCompleteEventNums(const std::string& tablename,
					 uint32_t min_event_num,
					 uint32_t max_event_num,
					 std::vector<uint32_t>& event_nums)
{
  VSDBStatement* stmt = fStmtRetrieval[tablename].limited_complete_eventnums;

  if(stmt == 0)
    {
      std::string tablename_ev = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);  

      stmt = fDB->createSelectQuery(tablename_ev, 
				    "EventComplete=1 AND "
				    "EventID>? AND EventID<=?",
				    "EventID",
				    VSDatabase::FLAG_NO_BUFFER);
      
      if(stmt == 0)
	{
	  std::cerr << __PRETTY_FUNCTION__ << ": could not query table: " 
		    << tablename_ev << std::endl;
	  return EXIT_FAILURE;
	}

      fStmtRetrieval[tablename].limited_complete_eventnums = stmt;
    }

  unsigned event_num;

  stmt->bindToParam(min_event_num);
  stmt->bindToParam(max_event_num);

  stmt->bindToResult(event_num);

  int ret = stmt->execute();
  if(ret < 0)
    {
      std::string tablename_ev = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);  

      std::cerr << __PRETTY_FUNCTION__ << ": error executing query on: " 
		<< tablename_ev << std::endl;
      return ret;
    }

  event_nums.clear();
  while(stmt->retrieveNextRow())event_nums.push_back(event_num);

  return ret;
}

int VSSimDB::getEventByNum(const std::string& tablename,  uint32_t event_id,
			   VSSimDBEventData& event)
{
  VSDBStatement* stmt = fStmtRetrieval[tablename].event_by_num;

  if(stmt == 0)
    {
      std::string tablename_ev = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

      stmt = fDB->createSelectQuery(tablename_ev, "EventID=?",
				    "WorkunitRunID, "
				    "TargetZenithRad, "
				    "TargetAzimuthRad, "
				    "PrimaryZenithRad, "
				    "PrimaryAzimuthRad, "
				    "PrimaryCoreEastM, "
				    "PrimaryCoreNorthM, "
				    "PrimaryCoreUpASLM, "
				    "NumHitScopes, "
				    "EventComplete",
				    VSDatabase::FLAG_NO_BUFFER);
      
      if(stmt == 0)
	{
	  std::cerr << __PRETTY_FUNCTION__ << ": could not query table: " 
		    << tablename_ev << std::endl;
	  return EXIT_FAILURE;
	}

      fStmtRetrieval[tablename].event_by_num = stmt;
    }

  stmt->bindToParam(event_id);

  stmt->bindToResult(event.fWorkunitRunID);
  stmt->bindToResult(event.fTargetZenithRad);
  stmt->bindToResult(event.fTargetAzimuthRad);
  stmt->bindToResult(event.fPrimaryZenithRad);
  stmt->bindToResult(event.fPrimaryAzimuthRad);
  stmt->bindToResult(event.fPrimaryCoreEastM);
  stmt->bindToResult(event.fPrimaryCoreNorthM);
  stmt->bindToResult(event.fPrimaryCoreUpASLM);
  stmt->bindToResult(event.fNumHitScopes);
  stmt->bindToResult(event.fEventComplete);

  int ret = stmt->execute();
  if(ret < 0)
    {
      std::string tablename_ev = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

      std::cerr << __PRETTY_FUNCTION__ << ": error executing query on: "
		<< tablename_ev << std::endl;
      return ret;
    }

  if(stmt->retrieveNextRow())
    {
      assert(stmt->retrieveNextRow() == 0);
      event.fEventID = event_id;
      return 1;
    }

  return 0;
}

int VSSimDB::getAllScopeByNum(const std::string& tablename, 
			      uint32_t event_id,
			      std::vector<VSSimDBScopeData>& scopes)
{
  VSDBStatement* stmt = fStmtRetrieval[tablename].all_scope_by_num;

  if(stmt == 0)
    {
      std::string tablename_tel = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_SCOPES);

      stmt = fDB->createSelectQuery(tablename_tel, "EventID=?",
				    "ScopeID, "
				    "ScopeZenithRad, "
				    "ScopeAzimuthRad, "
				    "NumHitPixels",
				    VSDatabase::FLAG_NO_BUFFER);
      
      if(stmt == 0)
	{
	  std::cerr << __PRETTY_FUNCTION__ << ": could not query table: " 
		    << tablename_tel << std::endl;
	  return EXIT_FAILURE;
	}

      fStmtRetrieval[tablename].all_scope_by_num = stmt;
    }

  stmt->bindToParam(event_id);

  VSSimDBScopeData scope;
  scope.fEventID = event_id;

  stmt->bindToResult(scope.fScopeID);
  stmt->bindToResult(scope.fScopeZenithRad);
  stmt->bindToResult(scope.fScopeAzimuthRad);
  stmt->bindToResult(scope.fNumHitPixels);

  int ret = stmt->execute();
  if(ret < 0)
    {
      std::string tablename_tel = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_SCOPES);

      std::cerr << __PRETTY_FUNCTION__ << ": error executing query on: "
		<< tablename_tel << std::endl;
      return ret;
    }

  scopes.clear();
  while(stmt->retrieveNextRow())scopes.push_back(scope);

  return ret;
}

int VSSimDB::getAllPEByNum(const std::string& tablename, 
			   uint32_t event_id,
			   std::vector<VSSimDBPEData>& pes)
{
  VSDBStatement* stmt = fStmtRetrieval[tablename].all_pe_by_num;

  if(stmt == 0)
    {
      std::string tablename_pe = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_PES);

      stmt = fDB->createSelectQuery(tablename_pe, "EventID=?",
				    "ScopeID, "
				    "PixelID, "
				    "PixelTimeNS",
				    VSDatabase::FLAG_NO_BUFFER);
      if(stmt == 0)
	{
	  std::cerr << __PRETTY_FUNCTION__ << ": could not query table: " 
		    << tablename_pe << std::endl;
	  return EXIT_FAILURE;
	}

      fStmtRetrieval[tablename].all_pe_by_num = stmt;
    }

  stmt->bindToParam(event_id);

  VSSimDBPEData pe;
  pe.fEventID = event_id;

  stmt->bindToResult(pe.fScopeID);
  stmt->bindToResult(pe.fPixelID);
  stmt->bindToResult(pe.fPixelTimeNS);

  int ret = stmt->execute();
  if(ret < 0)
    {
      std::string tablename_pe = 
	std::string(VSIMDB_TABLE_PREFIX_DATA)
	+ tablename
	+ std::string(VSIMDB_TABLE_POSTFIX_PES);

      std::cerr << __PRETTY_FUNCTION__ << ": error executing query on: "
		<< tablename_pe << std::endl;
      return ret;
    }

  pes.clear();
  while(stmt->retrieveNextRow())pes.push_back(pe);

  return ret;
}

int VSSimDB::getWorkunitIDs(const std::string& tablename, 
			    std::vector<unsigned>& workunit_ids)
{
  std::string tablename_ev = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    + tablename
    + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);
  
  uint32_t table_id = getTableIDByName(tablename);

  VSDBStatement* stmt = 
    fDB->createSelectQuery(VSIMDB_TABLE_NAME_WORKUNIT, 
			   "WorkunitRunComplete=1 AND TableID=?",
  			   "DISTINCT WorkunitRunID",
  			   VSDatabase::FLAG_NO_BUFFER);
  
  stmt->bindToParam(table_id);

//   VSDBStatement* stmt = 
//     fDB->createSelectQuery(tablename_ev, "EventComplete=1",
// 			   "DISTINCT WorkunitRunID",
// 			   VSDatabase::FLAG_NO_BUFFER);
  if(stmt == 0)
    {
      std::cerr << __PRETTY_FUNCTION__ << ": could not query table: " 
		<< tablename_ev << std::endl;
      return EXIT_FAILURE;
    }

  unsigned workunit_id;
  stmt->bindToResult(workunit_id);

  int ret = stmt->execute();
  if(ret < 0)
    {
      std::cerr << __PRETTY_FUNCTION__ << ": error executing query on: "
		<< tablename_ev << std::endl;
      return ret;
    }

  workunit_ids.clear();
  while(stmt->retrieveNextRow())workunit_ids.push_back(workunit_id);

  delete stmt;

  return ret;
}

// ============================================================================
// Data insertion
// ============================================================================

void VSSimDB::setInsertTableName(const std::string& prefix, bool event_only)
{
  VSDBStatement* stmt;
  std::string table_name;

  table_name = std::string(VSIMDB_TABLE_PREFIX_DATA)
    + prefix + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);
  stmt = fDB->createInsertQuery(table_name,9 /* 10 */,
				"WorkunitRunID, "
				"TargetZenithRad, "
				"TargetAzimuthRad, "
				"PrimaryZenithRad, "
				"PrimaryAzimuthRad, "
				"PrimaryCoreEastM, "
				"PrimaryCoreNorthM, "
				"PrimaryCoreUpASLM, "
				"EventComplete");
  fStmtInsertEvent.reset(stmt);
  fInsertBoundEvent = NULL;

  stmt = fDB->createQuery(std::string("UPDATE ")+table_name+
			  std::string(" SET NumHitScopes=?, EventComplete=1"
				      " WHERE EventID=?"));
  fStmtUpdateHitScopes.reset(stmt);
  fUpdateBoundEvent = NULL;

  if(event_only)return;

  table_name = std::string(VSIMDB_TABLE_PREFIX_DATA)
    + prefix + std::string(VSIMDB_TABLE_POSTFIX_SCOPES);
  stmt = fDB->createInsertQuery(table_name,4 /* 5 */,
				"EventID, "
				"ScopeID, "
				"ScopeZenithRad, "
				"ScopeAzimuthRad");
  fStmtInsertScope.reset(stmt);
  fInsertBoundScope = NULL;

  stmt = fDB->createQuery(std::string("UPDATE ")+table_name+
			  std::string(" SET NumHitPixels=?"
				      " WHERE EventID=? AND ScopeID=?"));
  fStmtUpdateHitPixels.reset(stmt);
  fUpdateBoundScope = NULL;

  table_name = std::string(VSIMDB_TABLE_PREFIX_DATA)
    + prefix + std::string(VSIMDB_TABLE_POSTFIX_PES);
  stmt = fDB->createInsertQuery(table_name,4);
  fStmtInsertPE.reset(stmt);
  fInsertBoundPE = NULL;
}

int VSSimDB::insertEvent(VSSimDBEventData& event)
{
  if(fInsertBoundEvent != &event)
    {
      fInsertBoundEvent = &event;
      fStmtInsertEvent->bindToParam(event.fWorkunitRunID);
      fStmtInsertEvent->bindToParam(event.fTargetZenithRad);
      fStmtInsertEvent->bindToParam(event.fTargetAzimuthRad);
      fStmtInsertEvent->bindToParam(event.fPrimaryZenithRad);
      fStmtInsertEvent->bindToParam(event.fPrimaryAzimuthRad);
      fStmtInsertEvent->bindToParam(event.fPrimaryCoreEastM);
      fStmtInsertEvent->bindToParam(event.fPrimaryCoreNorthM);
      fStmtInsertEvent->bindToParam(event.fPrimaryCoreUpASLM);
      fStmtInsertEvent->bindToParam(event.fEventComplete);
    }
  int ret = fStmtInsertEvent->execute();
  if(ret==1)event.fEventID=fStmtInsertEvent->getInsertID();
  return ret;
}

int VSSimDB::insertScope(const VSSimDBScopeData& scope)
{
  if(fInsertBoundScope != &scope)
    {
      fInsertBoundScope = &scope;
      fStmtInsertScope->bindToParam(scope.fEventID);
      fStmtInsertScope->bindToParam(scope.fScopeID);
      fStmtInsertScope->bindToParam(scope.fScopeZenithRad);
      fStmtInsertScope->bindToParam(scope.fScopeAzimuthRad);
    }
  return fStmtInsertScope->execute();
}

int VSSimDB::insertPE(const VSSimDBPEData& pe)
{
  if(fInsertBoundPE != &pe)
    {
      fInsertBoundPE = &pe;
      fStmtInsertPE->bindToParam(pe.fEventID);
      fStmtInsertPE->bindToParam(pe.fScopeID);
      fStmtInsertPE->bindToParam(pe.fPixelID);
      fStmtInsertPE->bindToParam(pe.fPixelTimeNS);
    }
  return fStmtInsertPE->execute();  
}

int VSSimDB::updateHitScopes(const VSSimDBEventData& event)
{
  if(fUpdateBoundEvent != &event)
    {
      fUpdateBoundEvent = &event;
      fStmtUpdateHitScopes->bindToParam(event.fNumHitScopes);
      fStmtUpdateHitScopes->bindToParam(event.fEventID);
    }
  return fStmtUpdateHitScopes->execute();    
}

int VSSimDB::updateHitPixels(const VSSimDBScopeData& scope)
{
  if(fUpdateBoundScope != &scope)
    {
      fUpdateBoundScope = &scope;
      fStmtUpdateHitPixels->bindToParam(scope.fNumHitPixels);
      fStmtUpdateHitPixels->bindToParam(scope.fEventID);
      fStmtUpdateHitPixels->bindToParam(scope.fScopeID);
    }
  return fStmtUpdateHitPixels->execute();    
}





// ============================================================================
// TEST MAIN FUNCTION
// ============================================================================

#ifdef TEST_MAIN

#include <VSOptions.hpp>
#include <VSDBFactory.hpp>

#include "RandomNumbers.hpp"

static inline float d2r(const float& x)
{
  return x/180.0*M_PI;
}

int main(int argc, char** argv)
{
  VSOptions options(argc,argv);
  VSDBFactory::configure(&options);

  char *progname = *argv;
  argv++, argc--;

  // --------------------------------------------------------------------------
  // Get the database name from the command line arguments
  // --------------------------------------------------------------------------
  if(argc<1)
    {
      std::cerr << "Usage: " << progname << " database" 
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  std::string database(*argv);
  argc--; argv++;

  // --------------------------------------------------------------------------
  // Create the database
  // --------------------------------------------------------------------------
  
  VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
  db->createDatabase(database,VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  db->useDatabase(database);
 
  // --------------------------------------------------------------------------
  // Simulation database
  // --------------------------------------------------------------------------
  VSSimDB sim_db(db);
  sim_db.createInfrastructureTables();

  // --------------------------------------------------------------------------
  // Create a table
  // --------------------------------------------------------------------------
  VSSimDBTableParam table;
  table.fPrimaryID          = 1;
  table.fEnergyGeV          = 10.0;
  table.fZenithMinRad       = d2r(0);
  table.fZenithMaxRad       = d2r(20);
  table.fAzimuthMinRad      = d2r(0);
  table.fAzimuthMaxRad      = d2r(360);
  table.fSamplingRadiusM    = 300.0;
  table.fTargetEventCount   = 1000;
  table.fWorkunitEventCount = 100;
  table.fWorkunitPriority   = 0;
  table.fOpticsID           = 0;
  table.fTableName          = sim_db.tableName(table);
  sim_db.createDataTable(table);

  // --------------------------------------------------------------------------
  // INSERT some data
  // --------------------------------------------------------------------------
  RandomNumbers rng("/tmp/seeds.dat");
  sim_db.setInsertTableName(table.fTableName);
  for(unsigned ievent=0; ievent<20; ievent++)
    {
      VSSimDBEventData event;
      event.fTargetZenithRad   = rng.Uniform()*M_PI/2.0;
      event.fTargetAzimuthRad  = rng.Uniform()*M_PI*2.0;
      event.fPrimaryZenithRad  = event.fTargetZenithRad;
      event.fPrimaryAzimuthRad = event.fTargetAzimuthRad;
      event.fPrimaryCoreEastM  = (rng.Uniform()-0.5)*1000;
      event.fPrimaryCoreNorthM = (rng.Uniform()-0.5)*1000;
      event.fPrimaryCoreUpASLM = 0;
      event.fNumHitScopes      = 0;
      sim_db.insertEvent(event);

      for(unsigned iscope=0; iscope<7; iscope++)
	if(rng.Uniform() < 0.5)
	  {
	    VSSimDBScopeData scope;
	    scope.fEventID         = event.fEventID;
	    scope.fScopeID         = iscope;
	    scope.fScopeZenithRad  = event.fTargetZenithRad;
	    scope.fScopeAzimuthRad = event.fTargetAzimuthRad;
	    scope.fNumHitPixels    = 0;
	    sim_db.insertScope(scope);
	    
	    bool scope_hit=false;
	    
	    for(unsigned ipix=0; ipix<1000; ipix++)
	      {
		unsigned npe = unsigned(floor(rng.Uniform()*10));
		for(unsigned ipe=0; ipe<npe; ipe++)
		  {
		    VSSimDBPEData pe;
		    pe.fEventID     = event.fEventID;
		    pe.fScopeID     = scope.fScopeID;
		    pe.fPixelID     = ipix;
		    pe.fPixelTimeNS = ipe;
		    sim_db.insertPE(pe);
		  }
		if(npe)scope.fNumHitPixels++;
	      }
  
	    sim_db.updateHitPixels(scope);
	    if(scope.fNumHitPixels)event.fNumHitScopes++;
	  }

      sim_db.updateHitScopes(event);

      std::cerr << "Event " << event.fEventID << " complete" << std::endl;
    }
}

#endif
