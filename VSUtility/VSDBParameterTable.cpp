//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBParameterTable.cpp
  Database access to parameter table

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       05/18/2005
*/

#include <VSSimDBTables.hpp>

#include "VSDBParameterTable.hpp"

using namespace VERITAS;

VSDBParameterTable::~VSDBParameterTable()
{
  // nothing to see here
}

void VSDBParameterTable::
createParameterTable()
{
  fDB->createTable(VSIMDB_TABLE_NAME_PARAMETERS,
		   "Collection varchar(40) NOT NULL,"
		   "Parameter varchar(40) NOT NULL,"
		   "Value varchar(255), PRIMARY KEY (Collection, Parameter)",
		   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
}

int VSDBParameterTable::
deleteParameterSet(const VSDBCollection& collection)
{
  return fDB->deleteFromTable(VSIMDB_TABLE_NAME_PARAMETERS,
		std::string("Collection='")+collection+std::string("'"));
}

int VSDBParameterTable::
storeParameterSet(const VSDBCollection& collection,
		  const VSDBParameterSet& parameters)
{
  VSDBStatement* stmt = 
    fDB->createInsertQuery(VSIMDB_TABLE_NAME_PARAMETERS, 3, "",
			   VSDatabase::FLAG_NO_SERVER_PS|
			   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);

  int insert_count=0;
  for(VSDBParameterSet::const_iterator ips = parameters.begin(); 
      ips != parameters.end(); ips++)
    {
      stmt->bindToParam(collection);
      stmt->bindToParam(ips->first);
      stmt->bindToParam(ips->second);
      int c = stmt->execute();
      if(c>0)insert_count += c;
    }

  delete stmt;
  return insert_count;
}

int VSDBParameterTable::
retrieveParameterSet(const VSDBCollection& collection,
		     VSDBParameterSet& parameters)
{
  VSDBStatement* stmt = 
    fDB->createSelectQuery(VSIMDB_TABLE_NAME_PARAMETERS, 
		   std::string("Collection='")+collection+std::string("'"),
			   "Parameter, Value", VSDatabase::FLAG_NO_SERVER_PS);

  std::string parameter;
  std::string value;

  int select_count = stmt->execute();
  stmt->bindToResult(parameter);
  stmt->bindToResult(value);  

  while(stmt->retrieveNextRow())parameters[parameter] = value;

  delete stmt;
  return select_count;
}

int VSDBParameterTable::
retrieveAllParameterSets(VSDBCollectionSet& collections)
{
  VSDBStatement* stmt = 
    fDB->createSelectQuery(VSIMDB_TABLE_NAME_PARAMETERS, "", "",
			   VSDatabase::FLAG_NO_SERVER_PS);

  std::string collection;
  std::string parameter;
  std::string value;

  int select_count = stmt->execute();
  stmt->bindToResult(collection);
  stmt->bindToResult(parameter);
  stmt->bindToResult(value);  

  while(stmt->retrieveNextRow())
    collections[collection][parameter] = value;

  delete stmt;
  return select_count;
}

int VSDBParameterTable::
setOrUpdateParameterValue(const VSDBCollection& collection,
			  const VSDBParameter& parameter,
			  const VSDBValue& value)
{
  VSDBStatement* stmt = 
    fDB->createInsertQuery(VSIMDB_TABLE_NAME_PARAMETERS, 3, "",
			   VSDatabase::FLAG_NO_SERVER_PS|
			   VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(stmt==0)return -1;
  
  stmt->bindToParam(collection);
  stmt->bindToParam(parameter);
  stmt->bindToParam(value);
  int result = stmt->execute();
  delete stmt;

  if(result==1)return result;

  stmt = 
    fDB->createQuery("UPDATE " VSIMDB_TABLE_NAME_PARAMETERS 
		     " SET Value=? WHERE Collection=? AND PARAMETER=?",
		     VSDatabase::FLAG_NO_SERVER_PS|
		     VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(stmt==0)return -1;
  
  stmt->bindToParam(value);
  stmt->bindToParam(collection);
  stmt->bindToParam(parameter);
  result = stmt->execute();
  delete stmt;

  return result;
}

int VSDBParameterTable::deleteParameter(const VSDBCollection& collection,
					const VSDBParameter& parameter)
{
  VSDBStatement* stmt = 
    fDB->createQuery("DELETE FROM " VSIMDB_TABLE_NAME_PARAMETERS 
  		     " WHERE Collection=? AND Parameter=?",
		     VSDatabase::FLAG_NO_SERVER_PS|
		     VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
  if(stmt==0)return -1;

  stmt->bindToParam(collection);
  stmt->bindToParam(parameter);
  int result = stmt->execute();
  
  delete stmt;
  return result;
}
