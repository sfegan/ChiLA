//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBMySQL3X.cpp
  MySQL 3.X database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
*/

#ifndef VSCONFIG_NO_MYSQL

#include "VSDBMySQL3X.hpp"

// ============================================================================
// ============================================================================

// VSDBMySQL3X

// ============================================================================
// ============================================================================

/*! \class VERITAS::VSDBMySQL3X
  Blah
*/

/*!
  See documentation for VSDBMySQLBase::VSDBMySQLBase().
*/
VERITAS::VSDBMySQL3X::
VSDBMySQL3X(bool loud, bool ignore_error,
	    const std::string& mysql_hostname, 
	    const std::string& mysql_username,
	    const std::string& mysql_password,
	    unsigned mysql_port, 
	    const std::string& mysql_unixpath, 
	    unsigned long mysql_client_flag,
	    bool parameters_override_all) throw(std::string)
  : VSDBMySQLBase(loud,ignore_error,
		  mysql_hostname,mysql_username,mysql_password, 
		  mysql_port,mysql_unixpath,mysql_client_flag,
		  parameters_override_all),
    fStatements()
{
  // nothing to see here
}

/*!
  Delete all buffered data which have not been retrieved. Delete all
  prepared statements.
*/
VERITAS::VSDBMySQL3X::
~VSDBMySQL3X()
{
  for(std::set<VSDBMySQL3XPreparedStatement*>::iterator i=fStatements.begin();
      i!=fStatements.end();i++)
    delete *i;
}
  
/*!
  The client side prepared statement class, VSDBMySQL3XPreparedStatement, 
  emulates the functionality of a server side prepared statement by creating
  individual queries to the server.
  \param query The SQL for the query with optional parameter markers.
  \param flags Flags which modify the behavior of the prepared statement.
  VSDBMySQL3X recognises and honors only FLAG_NO_BUFFER.
  \return An instance of VSDBMySQL3XPreparedStatement which embodies the
  query. The instance is created on the heap and should be deleted when
  it is no longer needed.
*/
  
VERITAS::VSDBStatement* VERITAS::VSDBMySQL3X::
createQuery(const std::string& query, uint32_t flags) throw()
{
  VSDBMySQL3XPreparedStatement* ps = 
    new VSDBMySQL3XPreparedStatement(this, query, flags);
  fStatements.insert(ps);
  return ps;
}
    
int VERITAS::VSDBMySQL3X::
psExecute(VSDBMySQL3XPreparedStatement* ps, const std::string& query)
{
  mysql_free_result(ps->fResultSet);
  ps->fResultSet=0;
  ps->fResultSetColumns=0;

  fQuery=query;

  int rows = executeQuery(ps->fQueryFlags,&ps->fResultSet);

#if(MYSQL_VERSION_ID<32224)
  ps->fResultSetColumns = mysql_num_fields(fDB);
#else
  ps->fResultSetColumns = mysql_field_count(fDB);
#endif

  return rows;
}

MYSQL_ROW VERITAS::VSDBMySQL3X::
psFetch(VSDBMySQL3XPreparedStatement* ps)
{
  MYSQL_ROW row = mysql_fetch_row(ps->fResultSet);

  if(row==0)
    {
      if(ps->fQueryFlags&FLAG_NO_BUFFER)
	{
	  if(mysql_errno(fDB)!=0)
	    {
	      std::ostringstream stream;
	      stream << "mysql_fetch_row: " << mysql_error(fDB) << std::endl;
	      fError = stream.str();
	      doError();
	    }
	  return 0;
	}
      mysql_free_result(ps->fResultSet);
      ps->fResultSet=0;
      ps->fResultSetColumns=0;
    }
  return row;
}

uint64_t VERITAS::VSDBMySQL3X::
psGetInsertID(VSDBMySQL3XPreparedStatement* ps)
{
  return mysql_insert_id(fDB);  
}

int VERITAS::VSDBMySQL3X::
psDelete(VSDBMySQL3XPreparedStatement* ps)
{
  if(ps->fResultSet)mysql_free_result(ps->fResultSet);
  fStatements.erase(ps);
  return 0;
}

// ----------------------------------------------------------------------------
// VSDBMySQL3XPreparedStatement
// ----------------------------------------------------------------------------

VERITAS::VSDBMySQL3XPreparedStatement::
~VSDBMySQL3XPreparedStatement()
{
  fMySQL->psDelete(this);
}

std::string VERITAS::VSDBMySQL3XPreparedStatement::
getErrorMessage() throw()
{
  return fMySQL->getErrorMessage();
}

template<typename T> inline std::string str_set(const T* val, const bool* null)
{
  if((null)&&(*null))return "NULL";
  else return VERITAS::MySQLValueToString(*val);
}

int VERITAS::VSDBMySQL3XPreparedStatement::
execute() throw()
{
  if(fParamReset==false)
    {
      vsassert(fParam.size() == fQueryBits.size()-1);
      fParamReset=true;
    }

  std::string query;
  std::vector<std::string>::iterator query_bit=fQueryBits.begin();
  query+=*query_bit;
  query_bit++;
  std::vector<BoundParam>::iterator param=fParam.begin();
  while(query_bit!=fQueryBits.end())
    {
      std::string param_str;
      switch(param->fDataType)
	{
	case VSDatabase::DT_BOOL:
	  param_str=str_set(param->fDataPtr.bool_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_INT8:
	  param_str=str_set(param->fDataPtr.i8_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_UINT8:
	  param_str=str_set(param->fDataPtr.u8_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_INT16:
	  param_str=str_set(param->fDataPtr.i16_ptr, param->fNullPtr);
	  break;
 	case VSDatabase::DT_UINT16:
	  param_str=str_set(param->fDataPtr.u16_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_INT32: 
	  param_str=str_set(param->fDataPtr.i32_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_UINT32:
	  param_str=str_set(param->fDataPtr.u32_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_INT64:
	  param_str=str_set(param->fDataPtr.i64_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_UINT64:
	  param_str=str_set(param->fDataPtr.u64_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_FLOAT: 
	  param_str=str_set(param->fDataPtr.flt_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_DOUBLE:
	  param_str=str_set(param->fDataPtr.dbl_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_STLSTRING:
	  param_str=str_set(param->fDataPtr.stl_str_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_CSTRING:
	  if((param->fNullPtr)&&(*param->fNullPtr))param_str="NULL";
	  else param_str=MySQLValueToString(param->fDataPtr.c_str_ptr,
					  *param->fDataCount);
	  break;
	case VSDatabase::DT_DATE:
	  param_str=str_set(param->fDataPtr.date_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_TIME:
	  param_str=str_set(param->fDataPtr.time_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_DATETIME:
	  param_str=str_set(param->fDataPtr.datetime_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_TIMESTAMP:
	  param_str=str_set(param->fDataPtr.timestamp_ptr, param->fNullPtr);
	  break;
	case VSDatabase::DT_UNKNOWN:
	  vsassert(0);
	  break;
	}
      query+=param_str;
      query+=*query_bit;

      param++;
      query_bit++;
    }
  return(fMySQL->psExecute(this,query));
}

template<typename T> void inline val_set(const char* col, T* val, bool* null)
{
  if(col)
    {
      if(null)*null = false;
      VERITAS::MySQLValueFromString(col,*val);
    }
  else
    {
      if(null)*null = true;
      *val = T();
    }
}

int VERITAS::VSDBMySQL3XPreparedStatement::
retrieveNextRow() throw()
{
  vsassert(fResultSet);

  MYSQL_ROW row = fMySQL->psFetch(this);
  if(!row)
    {
      fValReset=true;
      return 0;
    }

  if(fValReset==false)
    {
      vsassert(fVal.size()==fResultSetColumns);
      fValReset=true;
    }

  for(unsigned i=0;i<fResultSetColumns;i++)
    {
      switch(fVal[i].fDataType)
	{
	case VSDatabase::DT_BOOL:
	  val_set(row[i], fVal[i].fDataPtr.bool_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_INT8:
	  val_set(row[i], fVal[i].fDataPtr.i8_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_UINT8:
	  val_set(row[i], fVal[i].fDataPtr.u8_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_INT16:
	  val_set(row[i], fVal[i].fDataPtr.i16_ptr, fVal[i].fNullPtr);
	  break;
 	case VSDatabase::DT_UINT16:
	  val_set(row[i], fVal[i].fDataPtr.u16_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_INT32: 
	  val_set(row[i], fVal[i].fDataPtr.i32_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_UINT32:
	  val_set(row[i], fVal[i].fDataPtr.u32_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_INT64:
	  val_set(row[i], fVal[i].fDataPtr.i64_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_UINT64:
	  val_set(row[i], fVal[i].fDataPtr.u64_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_FLOAT: 
	  val_set(row[i], fVal[i].fDataPtr.flt_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_DOUBLE:
	  val_set(row[i], fVal[i].fDataPtr.dbl_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_STLSTRING:
	  val_set(row[i], fVal[i].fDataPtr.stl_str_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_CSTRING:
	  if(row[i])
	    {
	      MySQLValueFromString(row[i],fVal[i].fDataPtr.c_str_ptr,
				   fVal[i].fDataBufferCount, 
				   fVal[i].fDataCount);
	      if(fVal[i].fNullPtr)*fVal[i].fNullPtr = true;
	    }
	  else
	    {
	      if(fVal[i].fDataBufferCount>0)*fVal[i].fDataPtr.c_str_ptr = '\0';
	      if(fVal[i].fNullPtr)*fVal[i].fNullPtr = true;
	    }
	  break;
	case VSDatabase::DT_DATE:
	  val_set(row[i], fVal[i].fDataPtr.date_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_TIME:
	  val_set(row[i], fVal[i].fDataPtr.time_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_DATETIME:
	  val_set(row[i], fVal[i].fDataPtr.datetime_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_TIMESTAMP:
	  val_set(row[i], fVal[i].fDataPtr.timestamp_ptr, fVal[i].fNullPtr);
	  break;
	case VSDatabase::DT_UNKNOWN:
	  vsassert(0);
	  break;
	}
    }

  return 1;
}

uint64_t VERITAS::VSDBMySQL3XPreparedStatement::
getInsertID() throw()
{
  return fMySQL->psGetInsertID(this);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
clearBoundParams() throw()
{
  fParam.clear();
  fParamReset=false;
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const bool& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_BOOL;
  bind.fDataPtr.bool_ptr        = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const int8_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_INT8;
  bind.fDataPtr.i8_ptr          = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const uint8_t& val, const bool* null) throw ()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_UINT8;
  bind.fDataPtr.u8_ptr          = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const int16_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_INT16;
  bind.fDataPtr.i16_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const uint16_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_UINT16;
  bind.fDataPtr.u16_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const int32_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_INT32;
  bind.fDataPtr.i32_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const uint32_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_UINT32;
  bind.fDataPtr.u32_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const int64_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_INT64;
  bind.fDataPtr.i64_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
} 

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const uint64_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_UINT64;
  bind.fDataPtr.u64_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const float& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_FLOAT;
  bind.fDataPtr.flt_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const double& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_DOUBLE;
  bind.fDataPtr.dbl_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const std::string& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_STLSTRING;
  bind.fDataPtr.stl_str_ptr     = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const char* val, const unsigned long* count,
	    const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_CSTRING;
  bind.fDataPtr.c_str_ptr       = val;
  bind.fDataCount               = count;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}
    
void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const VSDBDate& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_DATE;
  bind.fDataPtr.date_ptr        = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const VSDBTime& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_TIME;
  bind.fDataPtr.time_ptr        = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const VSDBDateTime& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_DATETIME;
  bind.fDataPtr.datetime_ptr    = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToParam(const VSDBTimestamp& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_TIMESTAMP;
  bind.fDataPtr.timestamp_ptr   = &val;
  bind.fNullPtr                 = null;
  if(fParamReset)VSDBMySQL3XPreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::clearBoundResults() throw()
{
  fVal.clear();
  fValReset=false;
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(bool& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_BOOL;
  bind.fDataPtr.bool_ptr        = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(int8_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_INT8;
  bind.fDataPtr.i8_ptr          = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(uint8_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_UINT8;
  bind.fDataPtr.u8_ptr          = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(int16_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_INT16;
  bind.fDataPtr.i16_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(uint16_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_UINT16;
  bind.fDataPtr.u16_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(int32_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_INT32;
  bind.fDataPtr.i32_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(uint32_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_UINT32;
  bind.fDataPtr.u32_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(int64_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_INT64;
  bind.fDataPtr.i64_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(uint64_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_UINT64;
  bind.fDataPtr.u64_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(float& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_FLOAT;
  bind.fDataPtr.flt_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(double& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_DOUBLE;
  bind.fDataPtr.dbl_ptr         = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(std::string& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_STLSTRING;
  bind.fDataPtr.stl_str_ptr     = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(char* val, unsigned long buffer_count, 
	     unsigned long* count, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_CSTRING;
  bind.fDataPtr.c_str_ptr       = val;
  bind.fDataCount               = count;
  bind.fDataBufferCount         = buffer_count;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(VSDBDate& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_DATE;
  bind.fDataPtr.date_ptr        = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(VSDBTime& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_TIME;
  bind.fDataPtr.time_ptr        = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(VSDBDateTime& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_DATETIME;
  bind.fDataPtr.datetime_ptr    = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL3XPreparedStatement::
bindToResult(VSDBTimestamp& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_TIMESTAMP;
  bind.fDataPtr.timestamp_ptr   = &val;
  bind.fNullPtr                 = null;
  if(fValReset)VSDBMySQL3XPreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

VERITAS::VSDBMySQL3XPreparedStatement::
VSDBMySQL3XPreparedStatement(VSDBMySQL3X* mysql, const std::string& query,
			     uint32_t flags)
  : VSDBStatement(),
    fMySQL(mysql), fQueryFlags(flags), fResultSet(), fResultSetColumns(),
    fParam(), fParamReset(), fVal(), fValReset(), fQueryBits(), 
    fQuery(query)
{
  parseSQLForParameters(fQuery,fQueryBits);
}

#endif // VSCONFIG_NO_MYSQL

// ****************************************************************************
// ****************************************************************************
// **                                                                        **
// ** TEST CODE                                                              **
// **                                                                        **
// ****************************************************************************
// ****************************************************************************

#ifdef TEST_MAIN

using namespace VERITAS;

main(int argc, char** argv)
{
  VSOptions options(argc,argv);
  VSDBMySQLBase::configureFromCommandLine(options);
  VSDBMySQL3X test(true);
  std::cout << test.getClientVersion() << std::endl
	    << test.getServerVersion() << std::endl;
  test.useDatabase("array");
  VSDBStatement* statement = 
    test.createQuery("SELECT * FROM Parameters");
  
  std::string c;
  std::string p;
  std::string v;

  statement->bindToResult(c);
  statement->bindToResult(p);
  statement->bindToResult(v);
  statement->execute();

  while(statement->retrieveNextRow() == 1)
    {
      std::cout << std::left << std::setw(15) << c << ' '
		<< std::left << std::setw(35) << p << ' '
		<< std::left << std::setw(35) << v << std::endl;
    }
}
#endif
