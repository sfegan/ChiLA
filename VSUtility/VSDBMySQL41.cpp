//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBMySQL41.cpp
  MySQL 4.1 database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/02/2005
*/

#include "VSDBMySQL41.hpp"

#ifndef VSCONFIG_NO_MYSQL
#if(MYSQL_VERSION_ID>=40100)

/*!
  See documentation for VSDBMySQLBase::VSDBMySQLBase().
*/
VERITAS::VSDBMySQL41::
VSDBMySQL41(bool loud, bool ignore_error,
	    const std::string& mysql_hostname, 
	    const std::string& mysql_username,
	    const std::string& mysql_password,
	    unsigned mysql_port, 
	    const std::string& mysql_unixpath, 
	    unsigned long mysql_client_flag,
	    bool parameters_override_all) throw(std::string)
  : VSDBMySQL3X(loud,ignore_error,
		mysql_hostname,mysql_username,mysql_password, 
		mysql_port,mysql_unixpath,mysql_client_flag,
		parameters_override_all),
    fServerVersion(getServerVersion()), fStatements()
{
  // nothing to see here
}

/*!
  All buffered data which have not been retrieved are destroyed.
*/
VERITAS::VSDBMySQL41::
~VSDBMySQL41()
{
  std::set<VSDBMySQL41PreparedStatement*> delete_us = fStatements;
  for(std::set<VSDBMySQL41PreparedStatement*>::iterator i=delete_us.begin();
      i!=delete_us.end();i++)
    delete *i;
}
  
VERITAS::VSDBStatement* VERITAS::VSDBMySQL41::
createQuery(const std::string& query, uint32_t flags) throw()
{
  // Drop down to the MySQL 3.X class if the server does not support
  // prepared statements or if the user asks us not to use a
  // server-side prepared statement
  if((fServerVersion<40100)||(flags&FLAG_NO_SERVER_PS))
    return VSDBMySQL3X::createQuery(query,flags);
  
  MYSQL_STMT* stmt = mysql_stmt_init(fDB);
  if(!stmt)
    {
      fError = std::string("mysql_stmt_init: ")+mysql_error(fDB);
      doError();
      return 0;
    }
  
  if(fLoud)std::cerr << "mysql_stmt_prepare: " << query << std::endl;
  if(mysql_stmt_prepare(stmt, query.c_str(), query.size()))
    {
      unsigned int error_number = mysql_stmt_errno(stmt);

      if((flags&FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST == 0)
	 ||(!isErrorExistNotExist(error_number)))
	{
	  fQuery = query;
	  fError = std::string("mysql_stmt_prepare: ")+mysql_stmt_error(stmt);
	  doError();
	}
      mysql_stmt_close(stmt);
      return 0;
    }

  MYSQL_RES* meta = mysql_stmt_result_metadata(stmt);  
  
  VSDBMySQL41PreparedStatement* ps = 
    new VSDBMySQL41PreparedStatement(this, query, stmt, meta, flags);

  fStatements.insert(ps);
  return ps;
}
    
int VERITAS::VSDBMySQL41::
psExecute(VSDBMySQL41PreparedStatement* ps, MYSQL_BIND* bind)
{
  if((bind)&&(mysql_stmt_bind_param(ps->fStmt,bind)))
    {
      std::ostringstream stream;
      stream << "mysql_stmt_bind_param: " 
	     << mysql_stmt_error(ps->fStmt) << std::endl
	     << "mysql_stmt_bind_param: SQL: " << ps->fQuery;
      fError = stream.str();
      doError();
      return -1;
    }

  if(fLoud)std::cerr << ps->fQuery << std::endl;
  if(mysql_stmt_execute(ps->fStmt))
    {
      std::ostringstream stream;
      stream << "mysql_stmt_execute: " 
	     << mysql_stmt_error(ps->fStmt) << std::endl
	     << "mysql_stmt_execute: SQL: " << ps->fQuery;
      fError = stream.str();
      doError();
      return -1;
    }

  if(ps->fValCount)
    {
      if(ps->fQueryFlags&FLAG_NO_BUFFER)return 1; // What to return here ?

      if(mysql_stmt_store_result(ps->fStmt))
	{
	  std::ostringstream stream;
	  stream << "mysql_stmt_store_result: " 
		 << mysql_stmt_error(ps->fStmt) << std::endl
		 << "mysql_stmt_store_result: SQL: " << ps->fQuery;
	  fError = stream.str();
	  doError();
	  return -1;
	}

      int rows = mysql_stmt_num_rows(ps->fStmt);
      if((fLoud)&&(rows!=1))std::cerr << rows << " rows in set" << std::endl;
      else if(fLoud)std::cerr << "1 row in set" << std::endl;
      return rows;
    }
  else
    {
      int rows = mysql_stmt_affected_rows(ps->fStmt);
      if((fLoud)&&(rows!=1))std::cerr << rows << " rows affected" << std::endl;
      else if(fLoud)std::cerr << "1 row affected" << std::endl;
      return rows;
    }
}

int VERITAS::VSDBMySQL41::
psFetch(VSDBMySQL41PreparedStatement* ps, MYSQL_BIND* bind)
{
  if((bind)&&(mysql_stmt_bind_result(ps->fStmt, bind)))
    {
      std::ostringstream stream;
      stream << "mysql_stmt_bind_result: " 
	     << mysql_stmt_error(ps->fStmt) << std::endl
	     << "mysql_stmt_bind_result: SQL: " << ps->fQuery;
      fError = stream.str();
      doError();
      return -1;
    }

  int fetch = mysql_stmt_fetch(ps->fStmt);
  if(fetch == 0)return 1;
#ifdef MYSQL_DATA_TRUNCATED
  else if(fetch == MYSQL_DATA_TRUNCATED)return 1;
#endif
  else if(fetch == MYSQL_NO_DATA)
    {
      if((ps->fQueryFlags&FLAG_NO_BUFFER)&&(fLoud))
	{
	  int rows = mysql_stmt_num_rows(ps->fStmt);
	  if(rows!=1)std::cerr << rows << " rows in set" << std::endl;
	  else std::cerr << "1 row in set" << std::endl;
	}
      mysql_stmt_free_result(ps->fStmt);
      return 0;
    }
  else
    {
      std::ostringstream stream;
      stream << "mysql_stmt_fetch: " 
	     << mysql_stmt_error(ps->fStmt) << std::endl
	     << "mysql_stmt_fetch: SQL: " << ps->fQuery;
      fError = stream.str();
      doError();
      return -1;
    }
}

uint64_t VERITAS::VSDBMySQL41::
psGetInsertID(VSDBMySQL41PreparedStatement* ps)
{
  return mysql_stmt_insert_id(ps->fStmt);
}

int VERITAS::VSDBMySQL41::
psDelete(VSDBMySQL41PreparedStatement* ps)
{
  std::set<VSDBMySQL41PreparedStatement*>::iterator psi = fStatements.find(ps);
  if(psi != fStatements.end())fStatements.erase(psi);

  if(ps->fMeta)
    {
      mysql_free_result(ps->fMeta);
      ps->fMeta=0;
    }

  if(mysql_stmt_close(ps->fStmt))
    {
      std::ostringstream stream;
      stream << "mysql_stmt_close: "
	     << mysql_stmt_error(ps->fStmt) << std::endl
	     << "mysql_stmt_close: SQL: " << ps->fQuery;
      fError = stream.str();
      doError();
      return -1;
    }

  return 1;
}

// ----------------------------------------------------------------------------
// VSDBMySQL41PreparedStatement
// ----------------------------------------------------------------------------

VERITAS::VSDBMySQL41PreparedStatement::
~VSDBMySQL41PreparedStatement()
{
  fMySQL->psDelete(this);
  VSDBMySQL41PreparedStatement::clearBoundParams();
  VSDBMySQL41PreparedStatement::clearBoundResults();
}

std::string VERITAS::VSDBMySQL41PreparedStatement::
getErrorMessage() throw()
{
  return fMySQL->getErrorMessage();
}

int VERITAS::VSDBMySQL41PreparedStatement::
execute() throw()
{
  MYSQL_BIND* bind(0);
  if(fParamReset==false)
    {
      vsassert(fParam.size() == fParamCount);
      bind = new MYSQL_BIND[fParam.size()];
      MYSQL_BIND* b = bind;
  
      for(std::vector<BoundParam>::iterator i=fParam.begin(); 
	  i!=fParam.end(); i++, b++)
	{
	  memset(b,0,sizeof(*b));
	  switch(i->fDataType)
	    {
	    case VSDatabase::DT_BOOL:
	      i->fMirrorPtr.u8_ptr        = new uint8_t;
	      b->buffer_type              = MYSQL_TYPE_TINY;
	      b->buffer                   = i->fMirrorPtr.u8_ptr;
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_INT8:
	      b->buffer_type              = MYSQL_TYPE_TINY;
	      b->buffer                   = const_cast<int8_t*>(i->fUserPtr.i8_ptr);
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT8:
	      b->buffer_type              = MYSQL_TYPE_TINY;
	      b->buffer                   = const_cast<uint8_t*>(i->fUserPtr.u8_ptr);
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_INT16:
	      b->buffer_type              = MYSQL_TYPE_SHORT;
	      b->buffer                   = const_cast<int16_t*>(i->fUserPtr.i16_ptr);
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT16:
	      b->buffer_type              = MYSQL_TYPE_SHORT;
	      b->buffer                   = const_cast<uint16_t*>(i->fUserPtr.u16_ptr);
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_INT32:
	      b->buffer_type              = MYSQL_TYPE_LONG;
	      b->buffer                   = const_cast<int32_t*>(i->fUserPtr.i32_ptr);
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT32:
	      b->buffer_type              = MYSQL_TYPE_LONG;
	      b->buffer                   = const_cast<uint32_t*>(i->fUserPtr.u32_ptr);
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_INT64:
	      b->buffer_type              = MYSQL_TYPE_LONGLONG;
	      b->buffer                   = const_cast<int64_t*>(i->fUserPtr.i64_ptr);
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT64:
	      b->buffer_type              = MYSQL_TYPE_LONGLONG;
	      b->buffer                   = const_cast<uint64_t*>(i->fUserPtr.u64_ptr);
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_FLOAT:
	      b->buffer_type              = MYSQL_TYPE_FLOAT;
	      b->buffer                   = const_cast<float*>(i->fUserPtr.flt_ptr);
	      break;
	    case VSDatabase::DT_DOUBLE:
	      b->buffer_type              = MYSQL_TYPE_DOUBLE;
	      b->buffer                   = const_cast<double*>(i->fUserPtr.dbl_ptr);
	      break;
	    case VSDatabase::DT_CSTRING:
	      b->buffer_type              = MYSQL_TYPE_STRING;
	      b->buffer                   = const_cast<char*>(i->fUserPtr.c_str_ptr);
	      b->length                   = const_cast<unsigned long*>(i->fUserCount);
	      break;
	    case VSDatabase::DT_STLSTRING:
#warning STL string limit of 255
	      i->fMirrorBufferCount       = 255;
	      i->fMirrorPtr.c_str_ptr     = new char[i->fMirrorBufferCount];
	      b->buffer_type              = MYSQL_TYPE_STRING;
	      b->buffer                   = i->fMirrorPtr.c_str_ptr;
	      b->length                   = &i->fMirrorCount;
	      break;
	    case VSDatabase::DT_DATE:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_DATE;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_TIME:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_TIME;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_DATETIME:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_DATETIME;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_TIMESTAMP:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_TIMESTAMP;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_UNKNOWN:
	      vsassert(0);
	      break;
	    }
	  if(i->fUserNullPtr)b->is_null   = &i->fMirrorNull;
	}

      fParamReset=true;
    }


  for(std::vector<BoundParam>::iterator i=fParam.begin(); 
      i!=fParam.end(); i++)
    {
      switch(i->fDataType)
	{
	case VSDatabase::DT_BOOL:
	  if((i->fUserNullPtr)&&(*i->fUserNullPtr))
	    {
	      i->fMirrorNull=true;
	      *i->fMirrorPtr.u8_ptr = false;
	    }
	  else
	    {
	      i->fMirrorNull=false;
	      *i->fMirrorPtr.u8_ptr = *i->fUserPtr.bool_ptr;
	    }
	  break;
	case VSDatabase::DT_INT8:
	case VSDatabase::DT_UINT8:
	case VSDatabase::DT_INT16:
	case VSDatabase::DT_UINT16:
	case VSDatabase::DT_INT32:
	case VSDatabase::DT_UINT32:
	case VSDatabase::DT_INT64:
	case VSDatabase::DT_UINT64:
	case VSDatabase::DT_FLOAT:
	case VSDatabase::DT_DOUBLE:
	case VSDatabase::DT_CSTRING:
	  break;
	case VSDatabase::DT_STLSTRING:
	  if((i->fUserNullPtr)&&(*i->fUserNullPtr))
	    {
	      i->fMirrorCount=0;
	      if(i->fMirrorBufferCount)i->fMirrorPtr.c_str_ptr[0]='\0';
	      i->fMirrorNull=true;
	    }
	  else
	    {
	      i->fMirrorCount = i->fUserPtr.stl_str_ptr->size();
	      if(i->fMirrorCount >= i->fMirrorBufferCount)
		i->fMirrorCount = i->fMirrorBufferCount-1;
	      strncpy(i->fMirrorPtr.c_str_ptr, 
		      i->fUserPtr.stl_str_ptr->c_str(), i->fMirrorCount);
	      i->fMirrorPtr.c_str_ptr[i->fMirrorCount]='\0';
	      i->fMirrorNull=false;
	    }
	  break;
	case VSDatabase::DT_DATE:
	  if((i->fUserNullPtr)&&(*i->fUserNullPtr))i->fMirrorNull=true;
	  else
	    {
	      i->fMirrorNull=false;
	      i->fMirrorPtr.my_time_ptr->year   = i->fUserPtr.date_ptr->year;
	      i->fMirrorPtr.my_time_ptr->month  = i->fUserPtr.date_ptr->month;
	      i->fMirrorPtr.my_time_ptr->day    = i->fUserPtr.date_ptr->day;
	      i->fMirrorPtr.my_time_ptr->hour   = 0;
	      i->fMirrorPtr.my_time_ptr->minute = 0;
	      i->fMirrorPtr.my_time_ptr->second = 0;
	      i->fMirrorPtr.my_time_ptr->second_part = 0;
	      i->fMirrorPtr.my_time_ptr->neg    = false;
	    }
	  break;
	case VSDatabase::DT_TIME:
	  if((i->fUserNullPtr)&&(*i->fUserNullPtr))i->fMirrorNull=true;
	  else
	    {
	      i->fMirrorNull=false;
	      i->fMirrorPtr.my_time_ptr->year   = 0;
	      i->fMirrorPtr.my_time_ptr->month  = 0;
	      i->fMirrorPtr.my_time_ptr->day    = 0;
	      i->fMirrorPtr.my_time_ptr->hour   = i->fUserPtr.time_ptr->hour;
	      i->fMirrorPtr.my_time_ptr->minute = i->fUserPtr.time_ptr->minute;
	      i->fMirrorPtr.my_time_ptr->second = i->fUserPtr.time_ptr->second;
	      i->fMirrorPtr.my_time_ptr->second_part = 0;
	      i->fMirrorPtr.my_time_ptr->neg    = false;
	    }
	  break;
	case VSDatabase::DT_DATETIME:
	  if((i->fUserNullPtr)&&(*i->fUserNullPtr))i->fMirrorNull=true;
	  else
	    {
	      i->fMirrorNull=false;
	      i->fMirrorPtr.my_time_ptr->year   = 
		i->fUserPtr.datetime_ptr->year;
	      i->fMirrorPtr.my_time_ptr->month  = 
		i->fUserPtr.datetime_ptr->month;
	      i->fMirrorPtr.my_time_ptr->day    = 
		i->fUserPtr.datetime_ptr->day;
	      i->fMirrorPtr.my_time_ptr->hour   = 
		i->fUserPtr.datetime_ptr->hour;
	      i->fMirrorPtr.my_time_ptr->minute = 
		i->fUserPtr.datetime_ptr->minute;
	      i->fMirrorPtr.my_time_ptr->second = 
		i->fUserPtr.datetime_ptr->second;
	      i->fMirrorPtr.my_time_ptr->second_part = 0;
	      i->fMirrorPtr.my_time_ptr->neg    = false;
	    }
	  break;
	case VSDatabase::DT_TIMESTAMP:
	  if(((i->fUserNullPtr)&&(*i->fUserNullPtr))||
	     (i->fUserPtr.timestamp_ptr->current))i->fMirrorNull=true;
	  else
	    {
	      i->fMirrorNull=false;
	      i->fMirrorPtr.my_time_ptr->year   = 
		i->fUserPtr.timestamp_ptr->year;
	      i->fMirrorPtr.my_time_ptr->month  = 
		i->fUserPtr.timestamp_ptr->month;
	      i->fMirrorPtr.my_time_ptr->day    = 
		i->fUserPtr.timestamp_ptr->day;
	      i->fMirrorPtr.my_time_ptr->hour   = 
		i->fUserPtr.timestamp_ptr->hour;
	      i->fMirrorPtr.my_time_ptr->minute = 
		i->fUserPtr.timestamp_ptr->minute;
	      i->fMirrorPtr.my_time_ptr->second = 
		i->fUserPtr.timestamp_ptr->second;
	      i->fMirrorPtr.my_time_ptr->second_part = 0;
	      i->fMirrorPtr.my_time_ptr->neg    = false;
	    }
	  break;
	case VSDatabase::DT_UNKNOWN:
	  vsassert(0);
	  break;
	}
      if(i->fUserNullPtr)i->fMirrorNull=*i->fUserNullPtr;
    }
  
  int ret = fMySQL->psExecute(this,bind);
  if(bind)delete[] bind;
  return ret;
}

int VERITAS::VSDBMySQL41PreparedStatement::
retrieveNextRow() throw()
{
  MYSQL_BIND* bind(0);
  
  if(fValReset==false)
    {
      vsassert(fVal.size() == fValCount);
      bind = new MYSQL_BIND[fVal.size()];
      MYSQL_BIND* b = bind;
      
      for(std::vector<BoundVal>::iterator i=fVal.begin(); 
	  i!=fVal.end(); i++, b++)
	{
	  memset(b,0,sizeof(*b));
	  switch(i->fDataType)
	    {
	    case VSDatabase::DT_BOOL:
	      i->fMirrorPtr.u8_ptr        = new uint8_t;
	      b->buffer_type              = MYSQL_TYPE_TINY;
	      b->buffer                   = i->fMirrorPtr.u8_ptr;
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_INT8:
	      b->buffer_type              = MYSQL_TYPE_TINY;
	      b->buffer                   = i->fUserPtr.i8_ptr;
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT8:
	      b->buffer_type              = MYSQL_TYPE_TINY;
	      b->buffer                   = i->fUserPtr.u8_ptr;
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_INT16:
	      b->buffer_type              = MYSQL_TYPE_SHORT;
	      b->buffer                   = i->fUserPtr.i16_ptr;
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT16:
	      b->buffer_type              = MYSQL_TYPE_SHORT;
	      b->buffer                   = i->fUserPtr.u16_ptr;
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_INT32:
	      b->buffer_type              = MYSQL_TYPE_LONG;
	      b->buffer                   = i->fUserPtr.i32_ptr;
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT32:
	      b->buffer_type              = MYSQL_TYPE_LONG;
	      b->buffer                   = i->fUserPtr.u32_ptr;
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_INT64:
	      b->buffer_type              = MYSQL_TYPE_LONGLONG;
	      b->buffer                   = i->fUserPtr.i64_ptr;
	      b->is_unsigned              = false;
	      break;
	    case VSDatabase::DT_UINT64:
	      b->buffer_type              = MYSQL_TYPE_LONGLONG;
	      b->buffer                   = i->fUserPtr.u64_ptr;
	      b->is_unsigned              = true;
	      break;
	    case VSDatabase::DT_FLOAT:
	      b->buffer_type              = MYSQL_TYPE_FLOAT;
	      b->buffer                   = i->fUserPtr.flt_ptr;
	      break;
	    case VSDatabase::DT_DOUBLE:
	      b->buffer_type              = MYSQL_TYPE_DOUBLE;
	      b->buffer                   = i->fUserPtr.dbl_ptr;
	      break;
	    case VSDatabase::DT_CSTRING:
	      b->buffer_type              = MYSQL_TYPE_STRING;
	      b->buffer                   = i->fUserPtr.c_str_ptr;
	      b->length                   = i->fUserCount;
	      break;
	    case VSDatabase::DT_STLSTRING:	
#warning STL string limit of 255
	      i->fMirrorBufferCount       = 255;
	      i->fMirrorPtr.c_str_ptr     = new char[i->fMirrorBufferCount];
	      b->buffer_type              = MYSQL_TYPE_STRING;
	      b->buffer                   = i->fMirrorPtr.c_str_ptr;
	      b->buffer_length            = i->fMirrorBufferCount;
	      b->length                   = &i->fMirrorCount;
	      break;
	    case VSDatabase::DT_DATE:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_DATE;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_TIME:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_TIME;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_DATETIME:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_DATETIME;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_TIMESTAMP:
	      i->fMirrorPtr.my_time_ptr   = new MYSQL_TIME;
 	      b->buffer_type              = MYSQL_TYPE_TIMESTAMP;
	      b->buffer                   = i->fMirrorPtr.my_time_ptr;
	      break;
	    case VSDatabase::DT_UNKNOWN:
	      vsassert(0);
	      break;
	    }
	  if(i->fUserNullPtr)b->is_null   = &i->fMirrorNull;
	}

      fValReset=true;
    }

  int row = fMySQL->psFetch(this,bind);
  
  if(row==1)
    {
      for(std::vector<BoundVal>::iterator i=fVal.begin();  i!=fVal.end(); i++)
	{
	  switch(i->fDataType)
	    {
	    case VSDatabase::DT_BOOL:
	      if(i->fMirrorNull)*i->fUserPtr.bool_ptr=false;
	      else *i->fUserPtr.bool_ptr=*i->fMirrorPtr.u8_ptr;
	      break;
	    case VSDatabase::DT_INT8:
	    case VSDatabase::DT_UINT8:
	    case VSDatabase::DT_INT16:
	    case VSDatabase::DT_UINT16:
	    case VSDatabase::DT_INT32:
	    case VSDatabase::DT_UINT32:
	    case VSDatabase::DT_INT64:
	    case VSDatabase::DT_UINT64:
	    case VSDatabase::DT_FLOAT:
	    case VSDatabase::DT_DOUBLE:
	    case VSDatabase::DT_CSTRING:
	      break;
	    case VSDatabase::DT_STLSTRING:
	      if(i->fMirrorNull)i->fUserPtr.stl_str_ptr->erase();
	      else 
		{
		  unsigned count = i->fMirrorCount;
		  if(count > i->fMirrorBufferCount)
		    count = i->fMirrorBufferCount;
		  i->fUserPtr.stl_str_ptr->assign(i->fMirrorPtr.c_str_ptr,
						  count);
		}
	      break;
	    case VSDatabase::DT_DATE:
	      if(!i->fMirrorNull)
		{
		  i->fUserPtr.date_ptr->year   = 
		    i->fMirrorPtr.my_time_ptr->year;
		  i->fUserPtr.date_ptr->month  = 
		    i->fMirrorPtr.my_time_ptr->month;
		  i->fUserPtr.date_ptr->day    = 
		    i->fMirrorPtr.my_time_ptr->day;
		}
	      break;
	    case VSDatabase::DT_TIME:
	      if(!i->fMirrorNull)
		{
		  i->fUserPtr.time_ptr->hour   = 
		    i->fMirrorPtr.my_time_ptr->hour;
		  i->fUserPtr.time_ptr->minute = 
		    i->fMirrorPtr.my_time_ptr->minute;
		  i->fUserPtr.time_ptr->second = 
		    i->fMirrorPtr.my_time_ptr->second;
		  i->fUserPtr.time_ptr->microsecond = 0;
		}
	      break;
	    case VSDatabase::DT_DATETIME:
	      if(!i->fMirrorNull)
		{
		  i->fUserPtr.datetime_ptr->year   = 
		    i->fMirrorPtr.my_time_ptr->year;
		  i->fUserPtr.datetime_ptr->month  = 
		    i->fMirrorPtr.my_time_ptr->month;
		  i->fUserPtr.datetime_ptr->day    = 
		    i->fMirrorPtr.my_time_ptr->day;
		  i->fUserPtr.datetime_ptr->hour   = 
		    i->fMirrorPtr.my_time_ptr->hour;
		  i->fUserPtr.datetime_ptr->minute = 
		    i->fMirrorPtr.my_time_ptr->minute;
		  i->fUserPtr.datetime_ptr->second = 
		    i->fMirrorPtr.my_time_ptr->second;
		  i->fUserPtr.datetime_ptr->microsecond = 0;
		}
	      break;
	    case VSDatabase::DT_TIMESTAMP:
	      if(!i->fMirrorNull)
		{
		  i->fUserPtr.timestamp_ptr->year   = 
		    i->fMirrorPtr.my_time_ptr->year;
		  i->fUserPtr.timestamp_ptr->month  = 
		    i->fMirrorPtr.my_time_ptr->month;
		  i->fUserPtr.timestamp_ptr->day    = 
		    i->fMirrorPtr.my_time_ptr->day;
		  i->fUserPtr.timestamp_ptr->hour   = 
		    i->fMirrorPtr.my_time_ptr->hour;
		  i->fUserPtr.timestamp_ptr->minute = 
		    i->fMirrorPtr.my_time_ptr->minute;
		  i->fUserPtr.timestamp_ptr->second = 
		    i->fMirrorPtr.my_time_ptr->second;
		  i->fUserPtr.timestamp_ptr->microsecond = 0;
		}
	      break;
	    case VSDatabase::DT_UNKNOWN:
	      vsassert(0);
	      break;
	    }
	  if(i->fUserNullPtr)*i->fUserNullPtr=i->fMirrorNull;
	}
    }

  delete[] bind;
  return row;
}

uint64_t VERITAS::VSDBMySQL41PreparedStatement::
getInsertID() throw()
{
  return fMySQL->psGetInsertID(this);
}

void VERITAS::VSDBMySQL41PreparedStatement::
clearBoundParams() throw()
{
  for(std::vector<BoundParam>::iterator i=fParam.begin(); i!=fParam.end(); i++)
    {
      switch(i->fDataType)
	{
	case VSDatabase::DT_BOOL:
	  delete i->fMirrorPtr.u8_ptr;
	  break;
	case VSDatabase::DT_INT8:
	case VSDatabase::DT_UINT8:
	case VSDatabase::DT_INT16:
	case VSDatabase::DT_UINT16:
	case VSDatabase::DT_INT32:
	case VSDatabase::DT_UINT32:
	case VSDatabase::DT_INT64:
	case VSDatabase::DT_UINT64:
	case VSDatabase::DT_FLOAT:
	case VSDatabase::DT_DOUBLE:
	case VSDatabase::DT_CSTRING:
	  break;
	case VSDatabase::DT_STLSTRING:
	  delete[] i->fMirrorPtr.c_str_ptr;
	  break;
	case VSDatabase::DT_DATE:
	case VSDatabase::DT_TIME:
	case VSDatabase::DT_DATETIME:
	case VSDatabase::DT_TIMESTAMP:
	  delete i->fMirrorPtr.my_time_ptr;
	  break;
	case VSDatabase::DT_UNKNOWN:
	  vsassert(0);
	  break;
	}
    }
  fParam.clear();
  fParamReset=false;
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const bool& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_BOOL;
  bind.fUserPtr.bool_ptr    = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const int8_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_INT8;
  bind.fUserPtr.i8_ptr      = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const uint8_t& val, const bool* null) throw ()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_UINT8;
  bind.fUserPtr.u8_ptr      = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const int16_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_INT16;
  bind.fUserPtr.i16_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const uint16_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_UINT16;
  bind.fUserPtr.u16_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const int32_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_INT32;
  bind.fUserPtr.i32_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const uint32_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_UINT32;
  bind.fUserPtr.u32_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const int64_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_INT64;
  bind.fUserPtr.i64_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
} 

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const uint64_t& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_UINT64;
  bind.fUserPtr.u64_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const float& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_FLOAT;
  bind.fUserPtr.flt_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const double& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_DOUBLE;
  bind.fUserPtr.dbl_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const std::string& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_STLSTRING;
  bind.fUserPtr.stl_str_ptr = &val;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const char* val, const unsigned long* count,
	    const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType            = VSDatabase::DT_CSTRING;
  bind.fUserPtr.c_str_ptr   = val;
  bind.fUserCount           = count;
  bind.fUserNullPtr         = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const VSDBDate& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_DATE;
  bind.fUserPtr.date_ptr        = &val;
  bind.fUserNullPtr             = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const VSDBTime& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_TIME;
  bind.fUserPtr.time_ptr        = &val;
  bind.fUserNullPtr             = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const VSDBDateTime& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_DATETIME;
  bind.fUserPtr.datetime_ptr    = &val;
  bind.fUserNullPtr             = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToParam(const VSDBTimestamp& val, const bool* null) throw()
{
  BoundParam bind;
  bind.fDataType                = VSDatabase::DT_TIMESTAMP;
  bind.fUserPtr.timestamp_ptr   = &val;
  bind.fUserNullPtr             = null;
  if(fParamReset)VSDBMySQL41PreparedStatement::clearBoundParams();
  fParam.push_back(bind);
}
    
void VERITAS::VSDBMySQL41PreparedStatement::clearBoundResults() throw()
{
  for(std::vector<BoundVal>::iterator i=fVal.begin();  i!=fVal.end(); i++)
    {
      switch(i->fDataType)
	{
	case VSDatabase::DT_BOOL:
	  delete i->fMirrorPtr.u8_ptr;
	  break;
	case VSDatabase::DT_INT8:
	case VSDatabase::DT_UINT8:
	case VSDatabase::DT_INT16:
	case VSDatabase::DT_UINT16:
	case VSDatabase::DT_INT32:
	case VSDatabase::DT_UINT32:
	case VSDatabase::DT_INT64:
	case VSDatabase::DT_UINT64:
	case VSDatabase::DT_FLOAT:
	case VSDatabase::DT_DOUBLE:
	case VSDatabase::DT_CSTRING:
	  break;
	case VSDatabase::DT_STLSTRING:
	  delete[] i->fMirrorPtr.c_str_ptr;
	  break;
	case VSDatabase::DT_DATE:
	case VSDatabase::DT_TIME:
	case VSDatabase::DT_DATETIME:
	case VSDatabase::DT_TIMESTAMP:
	  delete i->fMirrorPtr.my_time_ptr;
	  break;
	case VSDatabase::DT_UNKNOWN:
	  vsassert(0);
	  break;
	}
    }
  fVal.clear();
  fValReset=false;
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(bool& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_BOOL;
  bind.fUserPtr.bool_ptr    = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(int8_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_INT8;
  bind.fUserPtr.i8_ptr      = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(uint8_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_UINT8;
  bind.fUserPtr.u8_ptr      = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(int16_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_INT16;
  bind.fUserPtr.i16_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(uint16_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_UINT16;
  bind.fUserPtr.u16_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(int32_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_INT32;
  bind.fUserPtr.i32_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(uint32_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_UINT32;
  bind.fUserPtr.u32_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(int64_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_INT64;
  bind.fUserPtr.i64_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(uint64_t& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_UINT64;
  bind.fUserPtr.u64_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(float& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_FLOAT;
  bind.fUserPtr.flt_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(double& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_DOUBLE;
  bind.fUserPtr.dbl_ptr     = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(std::string& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_STLSTRING;
  bind.fUserPtr.stl_str_ptr = &val;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(char* val, unsigned long buffer_count, 
	     unsigned long* count, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType            = VSDatabase::DT_CSTRING;
  bind.fUserPtr.c_str_ptr   = val;
  bind.fUserCount           = count;
  bind.fUserBufferCount     = buffer_count;
  bind.fUserNullPtr         = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(VSDBDate& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_DATE;
  bind.fUserPtr.date_ptr        = &val;
  bind.fUserNullPtr             = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(VSDBTime& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_TIME;
  bind.fUserPtr.time_ptr        = &val;
  bind.fUserNullPtr             = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(VSDBDateTime& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_DATETIME;
  bind.fUserPtr.datetime_ptr    = &val;
  bind.fUserNullPtr             = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

void VERITAS::VSDBMySQL41PreparedStatement::
bindToResult(VSDBTimestamp& val, bool* null) throw()
{
  BoundVal bind;
  bind.fDataType                = VSDatabase::DT_TIMESTAMP;
  bind.fUserPtr.timestamp_ptr   = &val;
  bind.fUserNullPtr             = null;
  if(fValReset)VSDBMySQL41PreparedStatement::clearBoundResults();
  fVal.push_back(bind);
}

VERITAS::VSDBMySQL41PreparedStatement::
VSDBMySQL41PreparedStatement(VSDBMySQL41* mysql, const std::string& query,
			     MYSQL_STMT* stmt, MYSQL_RES* meta, 
			     uint32_t flags)
  : VSDBStatement(),
    fMySQL(mysql), fStmt(stmt), fMeta(meta), fQueryFlags(flags),
    fParam(), fParamCount(), fParamReset(),
    fVal(), fValCount(), fValReset(),
    fQuery(query)
{
  fParamCount = mysql_stmt_param_count(fStmt);
  if(fMeta)fValCount = mysql_num_fields(fMeta);
}

#endif // (MYSQL_VERSION_ID>=40100)
#endif // VSCONFIG_NO_MYSQL

// ****************************************************************************
// ****************************************************************************
// **                                                                        **
// ** TEST CODE                                                              **
// **                                                                        **
// ****************************************************************************
// ****************************************************************************

#ifdef TEST_MAIN

#if(TEST_MAIN == 1)
using namespace VERITAS;

#include <sys/time.h>
#include <sstream>

int main(int argc, char** argv)
{
  unsigned C=100000;

  VSOptions options(argc,argv);
  VSDBMySQLBase::configureFromCommandLine(options);
  VSDBMySQL41 db(false);

  db.dropDatabase("pstest");
  db.createDatabase("pstest");
  db.useDatabase("pstest");
  db.createTable("test","s SMALLINT, i INTEGER, u INTEGER UNSIGNED, "
		 "f FLOAT, d DOUBLE");

  VSDBStatement* ps1 = db.createQuery("INSERT INTO test VALUES (?,?,?,?,?)");
  VSDBStatement* ps2 = db.createQuery("INSERT INTO test VALUES (?,?,?,?,?)",
				      VSDBMySQL41::FLAG_NO_SERVER_PS);
  
  short s;
  int i;
  unsigned int u;
  float f;
  double d;
  ps1->bindToParam(s);
  ps1->bindToParam(i);
  ps1->bindToParam(u);
  ps1->bindToParam(f);
  ps1->bindToParam(d);

  ps2->bindToParam(s);
  ps2->bindToParam(i);
  ps2->bindToParam(u);
  ps2->bindToParam(f);
  ps2->bindToParam(d);
  
  struct timeval tv0;
  struct timeval tv1;
  struct timeval tv2;

  // --------------------------------------------------------------------------
  // WRITING TEST
  // --------------------------------------------------------------------------

  std::cout << "DB writing test:" << std::endl;

  gettimeofday(&tv0,0);
  
  for(unsigned c=0;c<C;c++)
    {
      s = (c%10000)-5000;
      i = c-50000;
      u = c;
      f = float(c)/double(C);
      d = double(c)/double(C) - 0.5;
      ps1->execute();
    }

  gettimeofday(&tv1,0);

  for(unsigned c=0;c<C;c++)
    {
      s = (c%10000)-5000;
      i = c-50000;
      u = c;
      f = float(c)/double(C);
      d = double(c)/double(C) - 0.5;
      ps2->execute();
    }
  
  gettimeofday(&tv2,0);
  
  double t1 = 
    (double(tv1.tv_sec-tv0.tv_sec)*1000.0+
     double(tv1.tv_usec-tv0.tv_usec)/1000.0)/C;
  double t2 = 
    (double(tv2.tv_sec-tv1.tv_sec)*1000.0+
     double(tv2.tv_usec-tv1.tv_usec)/1000.0)/C;
  
  std::cout << std::setprecision(6) << std::fixed
	    << "  Prepared statement: " << C << " entries at " 
	    << t1 << " ms/entry (" << int(1000.0/t1) << " Hz)" << std::endl
	    << "  Multiple INSERT:    " << C << " entries at " 
	    << t2 << " ms/entry (" << int(1000.0/t2) << " Hz)" << std::endl;

  delete ps1;
  delete ps2;

  // --------------------------------------------------------------------------
  // READING TEST
  // --------------------------------------------------------------------------

  std::cout << std::endl
	    << "DB reading test:" << std::endl;
  
  ps1 = db.createQuery("SELECT * FROM test");
  ps2 = db.createQuery("SELECT * FROM test", VSDBMySQL41::FLAG_NO_SERVER_PS);

  ps1->bindToResult(s);
  ps1->bindToResult(i);
  ps1->bindToResult(u);
  ps1->bindToResult(f);
  ps1->bindToResult(d);

  ps2->bindToResult(s);
  ps2->bindToResult(i);
  ps2->bindToResult(u);
  ps2->bindToResult(f);
  ps2->bindToResult(d);

  unsigned n1=0;
  unsigned n2=0;

  gettimeofday(&tv0,0);

  ps1->execute();
  while(ps1->retrieveNextRow() == 1)n1++;
  
  gettimeofday(&tv1,0);

  ps2->execute();
  while(ps2->retrieveNextRow() == 1)n2++; 

  gettimeofday(&tv2,0);
  
  t1 = 
    (double(tv1.tv_sec-tv0.tv_sec)*1000.0+
     double(tv1.tv_usec-tv0.tv_usec)/1000.0)/C;
  t2 = 
    (double(tv2.tv_sec-tv1.tv_sec)*1000.0+
     double(tv2.tv_usec-tv1.tv_usec)/1000.0)/C;

  std::cout << "  Prepared statement: " << n1 << " entries at " 
	    << t1 << " ms/entry (" << int(1000.0/t1) << " Hz)" << std::endl
	    << "  Multiple SELECT:    " << n2 << " entries at " 
	    << t2 << " ms/entry (" << int(1000.0/t2) << " Hz)" << std::endl;
}

#elif(TEST_MAIN == 2)

using namespace VERITAS;

#include <sstream>

struct all_data
{
  bool b;
  int8_t i8;
  uint8_t u8;
  int16_t i16;
  uint16_t u16;
  int32_t i32;
  uint32_t u32;
  int64_t i64;
  uint64_t u64;
  float f;
  double d;
  char* c;
  unsigned long int c_len;
  std::string s;
  VSDBDate date;
  VSDBTime time;
  VSDBDateTime datetime;
  VSDBTimestamp timestamp;
};

int main(int argc, char** argv)
{
  VSOptions options(argc,argv);
  VSDBMySQLBase::configureFromCommandLine(options);
  VSDBMySQL41 db(false);

  db.dropDatabase("pstest");
  db.createDatabase("pstest");
  db.useDatabase("pstest");

  all_data data_out;
  db.createTable("test",
		 db.sqlSpecOf("b",data_out.b,true)+
		 db.sqlSpecOf("i8",data_out.i8,false)+
		 db.sqlSpecOf("u8",data_out.u8,false)+
		 db.sqlSpecOf("i16",data_out.i16,false)+
		 db.sqlSpecOf("u16",data_out.u16,false)+
		 db.sqlSpecOf("i32",data_out.i32,false)+
		 db.sqlSpecOf("u32",data_out.u32,false)+
		 db.sqlSpecOf("i64",data_out.i64,false)+
		 db.sqlSpecOf("u64",data_out.u64,false)+
		 db.sqlSpecOf("f",data_out.f,false)+
		 db.sqlSpecOf("d",data_out.d,false)+
		 db.sqlSpecOf("c",data_out.c,false)+
		 db.sqlSpec("s",db.sqlTypeOf(data_out.s,30000),false)+
		 db.sqlSpecOf("da",data_out.date,false)+
		 db.sqlSpecOf("ti",data_out.time,false)+
		 db.sqlSpecOf("dt",data_out.datetime,false)+
		 db.sqlSpecOf("ts",data_out.timestamp,false));

  bool isnull;
  data_out.b=true;
  data_out.i8=-123;
  data_out.u8=253;
  data_out.i16=-23456;
  data_out.u16=63321;
  data_out.i32=-666252;
  data_out.u32=123456789;
  data_out.i64=-5123456789LL;
  data_out.u64=0xEFFFFFFFFFFFFLL;
  data_out.f=-3.4e-12;
  data_out.d=1.2345678901234567E256;
  data_out.c_len=45;
  data_out.c=new char[data_out.c_len+1];
  memset(data_out.c,'s',data_out.c_len);
  data_out.c[data_out.c_len]='\0';
  data_out.s=std::string(500,'c');
  data_out.c[data_out.c_len-1]='\0';
  data_out.date.year=1974;
  data_out.date.month=2;
  data_out.date.day=16;
  data_out.time.hour=12;
  data_out.time.minute=34;
  data_out.time.second=56;
  data_out.datetime.year=1976;
  data_out.datetime.month=6;
  data_out.datetime.day=29;
  data_out.datetime.hour=01;
  data_out.datetime.minute=23;
  data_out.datetime.second=45;
  data_out.timestamp.current=true;  
  
  VSDBStatement* stmt[2];
  stmt[0] = db.createInsertQuery("test",17,"",VSDatabase::FLAG_NO_SERVER_PS);
  stmt[1] = db.createInsertQuery("test",17,"");

  for(unsigned i=0;i<2;i++)
    {
      stmt[i]->bindToParam(data_out.b,&isnull);
      stmt[i]->bindToParam(data_out.i8,&isnull);
      stmt[i]->bindToParam(data_out.u8,&isnull);
      stmt[i]->bindToParam(data_out.i16,&isnull);
      stmt[i]->bindToParam(data_out.u16,&isnull);
      stmt[i]->bindToParam(data_out.i32,&isnull);
      stmt[i]->bindToParam(data_out.u32,&isnull);
      stmt[i]->bindToParam(data_out.i64,&isnull);
      stmt[i]->bindToParam(data_out.u64,&isnull);
      stmt[i]->bindToParam(data_out.f,&isnull);
      stmt[i]->bindToParam(data_out.d,&isnull);
      stmt[i]->bindToParam(data_out.c,&data_out.c_len,&isnull);
      stmt[i]->bindToParam(data_out.s,&isnull);
      stmt[i]->bindToParam(data_out.date,&isnull);
      stmt[i]->bindToParam(data_out.time,&isnull);
      stmt[i]->bindToParam(data_out.datetime,&isnull);
      stmt[i]->bindToParam(data_out.timestamp,&isnull);
    }

  for(unsigned i=0;i<2;i++)stmt[i]->execute();
  isnull=true;
  for(unsigned i=0;i<2;i++)stmt[i]->execute();

  delete stmt[0];
  delete stmt[1];

  stmt[0] = db.createQuery("SELECT * FROM test LIMIT 1",
			   VSDatabase::FLAG_NO_SERVER_PS);
  stmt[1] = db.createQuery("SELECT * FROM test LIMIT 1");

  all_data data_in;
  data_in.c_len=200;
  data_in.c = new char[data_in.c_len];

  std::string s;
  for(unsigned i=0;i<2;i++)
    {
      stmt[i]->bindToResult(data_in.b,&isnull);
      stmt[i]->bindToResult(data_in.i8,&isnull);
      stmt[i]->bindToResult(data_in.u8,&isnull);
      stmt[i]->bindToResult(data_in.i16,&isnull);
      stmt[i]->bindToResult(data_in.u16,&isnull);
      stmt[i]->bindToResult(data_in.i32,&isnull);
      stmt[i]->bindToResult(data_in.u32,&isnull);
      stmt[i]->bindToResult(data_in.i64,&isnull);
      stmt[i]->bindToResult(data_in.u64,&isnull);
      stmt[i]->bindToResult(data_in.f,&isnull);
      stmt[i]->bindToResult(data_in.d,&isnull);
      stmt[i]->bindToResult(data_in.c,data_in.c_len,&data_in.c_len,&isnull);
      stmt[i]->bindToResult(data_in.s,&isnull);
      stmt[i]->bindToResult(data_in.date,&isnull);
      stmt[i]->bindToResult(data_in.time,&isnull);
      stmt[i]->bindToResult(data_in.datetime,&isnull);
      stmt[i]->bindToResult(data_in.timestamp,&isnull);
      //stmt[i]->bindToResult(s,&isnull);
    }
  
  for(unsigned i=0;i<2;i++)
    {
      stmt[i]->execute();
      stmt[i]->retrieveNextRow();

      std::cout << data_in.b << std::endl;
      std::cout << int16_t(data_in.i8) << std::endl;
      std::cout << uint16_t(data_in.u8) << std::endl;
      std::cout << data_in.i16 << std::endl;
      std::cout << data_in.u16 << std::endl;
      std::cout << data_in.i32 << std::endl;
      std::cout << data_in.u32 << std::endl;
      std::cout << data_in.i64 << std::endl;
      std::cout << data_in.u64 << std::endl;
      std::cout << data_in.f << std::endl;
      std::cout << data_in.d << std::endl;
      std::cout << data_in.c_len << ' ' << data_in.c << std::endl;
      std::cout << data_in.s << std::endl;
      std::cout << data_in.date << std::endl;
      std::cout << data_in.time << std::endl;
      std::cout << data_in.datetime << std::endl;
      std::cout << data_in.timestamp << std::endl << std::endl;
      //std::cout << s << std::endl << std::endl;
    }
}

#else

using namespace VERITAS;

#include <sstream>

struct ConfigCollection
{
  std::string         fName;
  std::string         fTask;
  unsigned int        fID;
};

struct ConfigData
{
  uint32_t            fRunNo;
  uint32_t            fCollectionID;
  std::string         fParameter;
  std::string         fValue;
};

struct PEPerDC
{
  uint32_t            fRunNo;
  uint32_t            fCollectionID;
  uint8_t             fTelescopeID;
  uint16_t            fChannelNum;
  float               fGain;
  float               fDev;
}

struct OffChannels
{
  uint32_t            fRunNo;
  uint8_t             fTelescopeID;
  uint16_t            fChannelNum;
};

struct Pedestal
{
  uint32_t            fRunNo;
  uint32_t            fCollectionID;
  uint8_t             fTelescopeID;
  uint16_t            fChannelNum;
  float               fPed;
  float               fDev;
};

struct RelativeGain
{
  uint32_t            fRunNo;
  uint32_t            fCollectionID;
  uint8_t             fTelescopeID;
  uint16_t            fChannelNum;
  float               fGain;
  float               fDev;
};

struct Pointing
{
  uint32_t            fRunNo;
  uint32_t            fCollectionID;
  uint8_t             fTelescopeID;
  uint16_t            fChannelNum;
  float               fRealAz;
  float               fRealEl;
};

int main(int argc, char** argv)
{
  VSOptions options(argc,argv);
  VSDBMySQLBase::configureFromCommandLine(options);
  VSDBMySQL41 db(false);

  db.dropDatabase("pstest");
  db.createDatabase("pstest");
  db.useDatabase("pstest");

  all_data data_out;
  db.createTable("test",
		 db.sqlSpecOf("b",data_out.b,true)+
		 db.sqlSpecOf("i8",data_out.i8,false)+
		 db.sqlSpecOf("u8",data_out.u8,false)+
		 db.sqlSpecOf("i16",data_out.i16,false)+
		 db.sqlSpecOf("u16",data_out.u16,false)+
		 db.sqlSpecOf("i32",data_out.i32,false)+
		 db.sqlSpecOf("u32",data_out.u32,false)+
		 db.sqlSpecOf("i64",data_out.i64,false)+
		 db.sqlSpecOf("u64",data_out.u64,false)+
		 db.sqlSpecOf("f",data_out.f,false)+
		 db.sqlSpecOf("d",data_out.d,false)+
		 db.sqlSpecOf("c",data_out.c,false)+
		 db.sqlSpec("s",db.sqlTypeOf(data_out.s,30000),false)+
		 db.sqlSpecOf("da",data_out.date,false)+
		 db.sqlSpecOf("ti",data_out.time,false)+
		 db.sqlSpecOf("dt",data_out.datetime,false)+
		 db.sqlSpecOf("ts",data_out.timestamp,false));

#endif

#endif
