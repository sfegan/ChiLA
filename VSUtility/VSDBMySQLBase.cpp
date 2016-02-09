//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBMySQLBase.cpp
  Base class for MySQL database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
*/

#ifndef VSCONFIG_NO_MYSQL

#include <iostream>
#include <iomanip>
#include <sstream>

#include "VSDBMySQLBase.hpp"

#include <mysqld_error.h>

using namespace VERITAS;

std::string VSDBMySQLBase::sDemandHostname("");
std::string VSDBMySQLBase::sDemandUsername("");
std::string VSDBMySQLBase::sDemandPassword("");
std::string VSDBMySQLBase::sDemandPortNStr("");
std::string VSDBMySQLBase::sDemandUnixPath("");

/*! \class VSDBMySQLBase

  The VSDBMySQLBase class provided functionality that is shared
  between the various MySQL database classes. It is not intended to be
  used directly, indeed since the createQuery() function is not
  defined here it cannot be instantiated.
*/

VSDBMySQLBase::~VSDBMySQLBase()
{
  mysql_close(fDB);
  delete fDB;
}

/*!
  \return Message which describes what went wrong with last query
 */
std::string VSDBMySQLBase::getErrorMessage() throw()
{
  return fError;
}

std::string VSDBMySQLBase::getRDBM() throw()
{
  return "MySQL";
}

std::string VSDBMySQLBase::getHostname() throw()
{
  return fHostname;
}
    
std::string VSDBMySQLBase::getUsername() throw()
{
  return fUsername;
}

std::string VSDBMySQLBase::getPassword() throw()
{
  return fPassword;
}

std::string VSDBMySQLBase::getDatabase() throw()
{
  return fDatabase;
}

/*!
  \param dbname The name of the database to create, if it does not
  exist already.

  \return Negative value on failure. Zero or positive on success.

  \note The following SQL query is generated: 
  <b>CREATE DATABASE IF NOT EXISTS <i>dbname</i></b>
*/
int VSDBMySQLBase::createDatabase(const std::string& dbname,
				  uint32_t flags) throw()
{
  if(flags&FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST)
    fQuery = std::string("CREATE DATABASE IF NOT EXISTS ")+dbname;
  else
    fQuery = std::string("CREATE DATABASE ")+dbname;
  int rows = executeQuery(flags);
  if(rows<0)doError();
  return rows;
}

/*!
  \param dbname The name of the database to drop, if it exists.

  \return Negative value on failure. Zero or positive on success.

  \note The following SQL query is generated: 
  <b>DROP DATABASE IF EXISTS <i>dbname</i></b>
*/
int VSDBMySQLBase::dropDatabase(const std::string& dbname,
				uint32_t flags) throw()
{
  if(flags&FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST)
    fQuery = std::string("DROP DATABASE IF EXISTS ")+dbname;
  else
    fQuery = std::string("DROP DATABASE ")+dbname;
  int rows = executeQuery(flags);
  if(rows<0)doError();
  return rows;
}

/*!
  \param dbname The name of the database to select.

  \return Negative value on failure. Zero or positive on success.
*/
int VSDBMySQLBase::useDatabase(const std::string& dbname,
			       uint32_t flags) throw()
{
  fDatabase=dbname;
  if(fLoud)std::cerr << "USE DATABASE " << dbname << std::endl;
  if(mysql_select_db(fDB, dbname.c_str()))
    {
      fError = 
	std::string("mysql_select_db(\"") + dbname +
	std::string("\"): ") + std::string(mysql_error(fDB));
      doError();
      return -1;
    }
  return 1;
}

/*!
  \param table The name of the table to create, if it does not
  exist already.

  \param spec SQL specification of the columns in the table.

  \return Negative value on failure. Zero or positive on success.

  \note The following SQL query is generated: 
  <b>CREATE TABLE IF NOT EXISTS <i>table</i> (<i>spec</i>)</b>
*/
int VSDBMySQLBase::createTable(const std::string& table, 
			       const std::string& spec,
					uint32_t flags) throw()
{
  if(flags&FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST)
    fQuery = std::string("CREATE TABLE IF NOT EXISTS ");
  else
    fQuery = std::string("CREATE TABLE ");
  fQuery += table+std::string(" ( ")+spec+std::string(" )");
  int rows = executeQuery(flags);
  if(rows<0)doError();
  return rows;
}

/*!
  \param table The name of the table to drop, if it exists.

  \return Negative value on failure. Zero or positive on success.

  \note The following SQL query is generated: 
  <b>DROP TABLE IF EXISTS <i>table</i></b>
*/
int VSDBMySQLBase::dropTable(const std::string& table,
			     uint32_t flags) throw()
{
  if(flags&FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST)
    fQuery = std::string("DROP TABLE IF EXISTS ")+table;
  else 
    fQuery = std::string("DROP TABLE ")+table;
  int rows = executeQuery(flags);
  if(rows<0)doError();
  return rows;
}

/*!
  \param table The name of the database to select.

  \param cond SQL condition defining what to delete. If it is
  empty then all rows are deleted.

  \return Number of rows deleted, negative value on failure. 

  \note If <i>cond</i> is not empty, the following SQL query is
  generated: <b>DELETE FROM <i>table</i> WHERE <i>cond</i></b>
*/
int VSDBMySQLBase::deleteFromTable(const std::string& table,
				   const std::string& cond,
				   uint32_t flags) throw ()
{
  if(cond.empty())fQuery = std::string("DELETE FROM ")+table;
  else fQuery = std::string("DELETE FROM ")+table+std::string(" WHERE ")+cond;
  int rows = executeQuery(flags);
  if(rows<0)doError();
  return rows;
}

/*!
  \param table Name of the table to within which to execute the INSERT query.

  \param num_columns Number of columns to INSERT data into, which 
  corresponds to the number of parameters which must be bound to the 
  query.

  \param columns List of column names to insert data into or an empty 
  string, meaning insert into all columns.

  \param Flags which modify the behavior of the query. Multiple
  flags can be or'd together. See, for example, FLAG_NO_BUFFER.

  \return A statement which can be used to execute the INSERT query
*/

VSDBStatement* 
VSDBMySQLBase::createInsertQuery(const std::string& table, 
				 unsigned num_columns,
				 const std::string& columns, 
				 uint32_t flags) throw()
{
  std::ostringstream qstream;
  qstream << "INSERT ";
  if(flags&FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST)qstream << "IGNORE ";
  qstream << "INTO " << table << ' ';
  if(!columns.empty())qstream << "( " << columns << " ) ";
  qstream << "VALUES ( ";
  if(num_columns)qstream << '?';
  for(unsigned i=1;i<num_columns;i++)qstream << ",?";
  qstream << " )";
  return createQuery(qstream.str(), flags);
}

/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"BOOL"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const bool& x) const throw() 
{
  return "BOOL"; 
}

/*!
   \param x Value is unused. The parameter is only for the purpose of
   overloading.
   \return The string <code>"TINYINT"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const int8_t& x) const throw() 
{ 
  return "TINYINT"; 
}

/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"TINYINT UNSIGNED"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const uint8_t& x) const throw() 
{ 
  return "TINYINT UNSIGNED"; 
}
  
/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"SMALLINT"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const int16_t& x) const throw() 
{ 
  return "SMALLINT"; 
}
  
/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"SMALLINT UNSIGNED"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const uint16_t& x) const throw() 
{
  return "SMALLINT UNSIGNED"; 
}
  
/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"INT"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const int32_t& x) const throw() 
{
  return "INT"; 
}
  
/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"INT UNSIGNED"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const uint32_t& x) const throw() 
{
  return "INT UNSIGNED"; 
}
  
/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"BINGINT"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const int64_t& x) const throw() 
{
  return "BIGINT"; 
}
  
/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"BIGINT UNSIGNED"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const uint64_t& x) const throw() 
{
  return "BIGINT UNSIGNED"; 
}  

/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"FLOAT"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const float& x) const throw() 
{
  return "FLOAT"; 
}
  
/*!
  \param x Value is unused. The parameter is only for the purpose of 
  overloading.
   \return The string <code>"DOUBLE"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const double& x) const throw() 
{
  return "DOUBLE"; 
}
    
/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \param length Maximum length of the string. If <i>length</i> is zero then
   the maximum is taken as 255 charcters.
   \return The string <code>"VARCHAR(length)"</code> 
   if <i>length</i> is less than 256. <code>"TEXT"</code>, 
   <code>"MEDIUMTEXT"</code> or <code>"LONGTEXT"</code> otherwise.
*/
std::string VSDBMySQLBase::sqlTypeOf(const char* x, unsigned length)
  const throw() 
{
  if(length==0)length=255;
  if(length<=0x000000FF)
    {
      std::ostringstream stream;
      stream << "VARCHAR(" << length << ")";
      return stream.str();
    }
  else if(length<=0x0000FFFF)return "TEXT";
  else if(length<=0x00FFFFFF)return "MEDIUMTEXT";
  else return "LONGTEXT";
}

/*!
   See documentation for sqlTypeOf(const char* x, unsigned length)
*/
std::string VSDBMySQLBase::sqlTypeOf(const std::string& x,
					      unsigned length) const throw()
{
  return sqlTypeOf((const char*)0,length);
}

/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"DATE"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const VSDBDate& x) const throw()
{
  return "DATE";
}

/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"TIME"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const VSDBTime& x) const throw()
{
  return "TIME";
}

/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"DATETIME"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const VSDBDateTime& x) 
  const throw()
{
  return "DATETIME";
}

/*!
   \param x Value is unused. The parameter is only for the purpose of 
   overloading.
   \return The string <code>"TIMESTAMP"</code>.
*/
std::string VSDBMySQLBase::sqlTypeOf(const VSDBTimestamp& x) 
  const throw()
{
  return "TIMESTAMP";
}

/*! \fn VSDBMySQLBase::mysql() throw()

  \return The MySQL connection handle. 

  This function is available for convenience so that unimplemented
  functionality can be implemented externally to this class. You
  should take care when using the handle as you may interfere with the
  operation of the class.

*/

/*! \fn VSDBMySQLBase::sql() const throw()

  \return The SQL of the last query executed.

*/

/*! \fn VSDBMySQLBase::error() const throw()

  \return The MySQL portion of the error message reported by
  getErrorMessage().

*/

/*!

  Set the database connection options from the command line. The
  following options are recognised:

  - <code>-VSDBHost</code> and <code>--VSDBHost</code>
  - <code>-VSDBUser</code> and <code>--VSDBUser</code>
  - <code>-VSDBPass</code> and <code>--VSDBPass</code>
  - <code>-VSDBPort</code> and <code>--VSDBPort</code>
  - <code>-VSDBUnix</code> and <code>--VSDBUnix</code>

  \param options Options class which has been loaded with the command
  line options.

*/
void VSDBMySQLBase::
configureFromCommandLine(VSOptions& options) throw()
{
  options.findWithValue(std::string("VSDBHost"), sDemandHostname);
  options.findWithValue(std::string("VSDBUser"), sDemandUsername);
  options.findWithValue(std::string("VSDBPass"), sDemandPassword);
  options.findWithValue(std::string("VSDBPort"), sDemandPortNStr);
  options.findWithValue(std::string("VSDBUnix"), sDemandUnixPath);
}

/*!    
  \return Version of client library encoded as 10000*major_version +
  100*minor_version + sub_version
*/
int VSDBMySQLBase::getClientVersion() throw()
{
#if(MYSQL_VERSION_ID<40016)
  return MYSQL_VERSION_ID;
#else
  return mysql_get_client_version();
#endif
}

/*!    
  \return Version of client library encoded as 10000*major_version +
  100*minor_version + sub_version
*/
int VSDBMySQLBase::getServerVersion() throw()
{
#if(MYSQL_VERSION_ID<40100)
  fQuery = std::string("SELECT VERSION()");
  MYSQL_RES* result;
  int rows = executeQuery(0,&result);
  if(rows<0)
    {
      doError();
      return rows;
    }
  else if(rows!=1)
    {
      std::ostringstream stream;
      stream << "VSDBMySQLBase::getServerVersion(): Query returned " 
	     << rows << " rows" << std::endl
	     << "VSDBMySQLBase::getServerVersion(): SQL: " << fQuery;
      fError = stream.str();
      doError();
      return -1;
    }

  mysql_free_result(result);

  MYSQL_ROW row = mysql_fetch_row(result);
  std::string version_string(row[0]);
  std::istringstream version_stream(version_string);

  char c;
  int major_version;
  int minor_version;
  int sub_version;
  version_stream >> major_version >> c >> minor_version >> c >> sub_version;
  
  return major_version*10000 + minor_version*100 + sub_version;
#else
  return mysql_get_server_version(fDB);
#endif
}

static inline std::string 
parEnvCmd(const std::string& par, const std::string& env, 
	  const std::string& cmd) throw()
{
  if(!cmd.empty())return cmd;
  else if(!env.empty())
    {
      char* env_val_ptr = getenv(env.c_str());
      if(env_val_ptr)
	{
	  std::string env_val(env_val_ptr);
	  return env_val.substr(env_val.find('=')+1);
	}
    }
  
  return par;
}

static inline const char* str_or_null(const std::string& s) throw()
{
  if(s.length() == 0)return NULL;
  else return s.c_str();
}

/*!

  Construct the class and create the connection to the database. The
  database connection can be specified in the parameters to the
  constructor or configured from the environment variables or through
  the command line, using the configureFromCommandLine() function. 
  Normally the values for the connection options are chosen in the
  following order:

  -# From the command line using configureFromCommandLine() which 
     recognises the options 
     - <code>-VSDBHost</code> and <code>--VSDBHost</code>
     - <code>-VSDBUser</code> and <code>--VSDBUser</code>
     - <code>-VSDBPass</code> and <code>--VSDBPass</code>
     - <code>-VSDBPort</code> and <code>--VSDBPort</code>
     - <code>-VSDBUnix</code> and <code>--VSDBUnix</code>
  -# From the environmental variables 
     - <code>VSDB_HOST</code>
     - <code>VSDB_USER</code>
     - <code>VSDB_PASS</code>
     - <code>VSDB_POST</code>
     - <code>VSDB_UNIX</code>
  -# Using the parameters passed to the class
  -# As the default value of the MySQL C API

  However, if the <i>parameters_override_all</i> is set, then
  the values passed to the class override all others.

  \param loud If <i>true</i> print the SQL for each query to standard
  error

  \param ignore_error If <i>true</i> do not kill the program when an
  error occurs. If <i>false</i>, a message is printed to standard
  output when an error occurs and the program is killed.

  \param mysql_hostname Name or IP address of database host.

  \param mysql_username Name of user to use when connecting to MySQL.

  \param mysql_password Password to use when connecting to MySQL.

  \param mysql_port Network port to which to connect with database.

  \param mysql_unixpath Path of UNIX domain socket on local host to use
  when connecting with database. This is overridden by 
  <i>mysql_hostname</i>.
  
  \param parameters_override_all If <i>true</i> then non-zero
  parameters passwd to the constructor override those selected in the
  environmental variables or on the command line.

  \param mysql_client_flag See documentation of MySQL C API.

  \param parameters_override_all Override the parameter values set
  from the command line options or environment variables with the
  values passed to the constructor.

 */
VSDBMySQLBase::
VSDBMySQLBase(bool loud, bool ignore_error,
	      const std::string& mysql_hostname, 
	      const std::string& mysql_username,
	      const std::string& mysql_password,
	      unsigned mysql_port, const std::string& mysql_unixpath, 
	      unsigned long mysql_client_flag, bool parameters_override_all) 
  throw(std::string)
  : VSDatabase(),
    fHostname(), fUsername(), fPassword(), fPortNumber(), fUnixPath(),
    fDatabase(),
    fDB(), fQuery(), fError(), fLoud(loud), fIgnoreError(ignore_error)
{
  fDB = new MYSQL;

  // Initialize the DB object
  mysql_init(fDB);
  
  // Figure out where the DBMS is and how to connect to it
  std::ostringstream port_ostream;
  if(mysql_port > 0)port_ostream << mysql_port;
  std::string port_string = port_ostream.str();

  fHostname = parEnvCmd(mysql_hostname, "VSDB_HOST", sDemandHostname);
  fUsername = parEnvCmd(mysql_username, "VSDB_USER", sDemandUsername);
  fPassword = parEnvCmd(mysql_password, "VSDB_PASS", sDemandPassword);
  std::string port(parEnvCmd(port_string, "VSDB_PORT", sDemandPortNStr));
  fUnixPath = parEnvCmd(mysql_unixpath, "VSDB_UNIX", sDemandUnixPath);

  if(parameters_override_all)
    {
      if(!mysql_hostname.empty())fHostname=mysql_hostname;
      if(!mysql_username.empty())fUsername=mysql_username;
      if(!mysql_password.empty())fPassword=mysql_password;
      if(!port_string.empty())port=port_string;
      if(!mysql_unixpath.empty())fUnixPath=mysql_unixpath;
    }
  
  VSDataConverter::fromString(fPortNumber,port);
  
  if(fLoud)
    std::cerr << "mysql_real_connect(\"" 
	      << fHostname << "\",\"" << fUsername << "\",\""
	      << fPassword << "\",\"" << "\"," << fPortNumber << ",\""
	      << fUnixPath << "\"," << mysql_client_flag << ")"
	      << std::endl;
  
  // Connect to the DB
  if(mysql_real_connect(fDB,
			str_or_null(fHostname), str_or_null(fUsername),
			str_or_null(fPassword),
			NULL, fPortNumber, str_or_null(fUnixPath),
			mysql_client_flag) == NULL)
    {
      std::ostringstream stream;
      stream << "mysql_init: " << mysql_error(fDB) << std::endl;
      std::string error = stream.str();
      if(!ignore_error)
	{
	  std::cerr << error;
	  exit(EXIT_FAILURE);
	}
      delete fDB;
      throw error;
    }
}

/*!

  \param result Pointer to MYSQL_RES* in which the result set is
  returned if the query is succesful. The value of <i>result</i> can
  be set to NULL entering the function if the query is known not to
  generate a result set (for example an INSERT).

  \param flags Flags, as defined in VSDatabase.hpp. This function
  honors FLAG_NO_BUFFER and FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST.  The
  former stops the result set from being buffered immediately. If you
  do not buffer the data onto the client then MySQL will not allow
  another query to execute until all the rows are fetched. You should
  read and understand the implications of not buffering the data in
  the MySQL documentation of the <i>mysql_use_results</i> C API
  function. The FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST causes the function
  to ignore existance/non-existance errors.

  \return If the query was successful and generated a result set,
  which was buffered then the number of rows in the set is returned,
  which can be zero or greater. If the query was succesful, generated a
  result set which was not buffered then zero is returned. If the
  result was successful but did not generate a result set then the
  number of affected rows is returned. Finally if the query was not
  successful then -1 is returned

*/
int VSDBMySQLBase::executeQuery(uint32_t flags,
				MYSQL_RES** result) throw()
{
  if(fLoud)std::cerr << fQuery << std::endl;

  if(mysql_real_query(fDB, fQuery.c_str(), fQuery.length()))
    {     
      unsigned int error_number = mysql_errno(fDB);
      if((flags&FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST)
	 &&(isErrorExistNotExist(error_number)))return 0;

      std::ostringstream stream;
      stream << "mysql_real_query: " << mysql_error(fDB) << std::endl
	     << "mysql_real_query: SQL: " << fQuery;
      fError = stream.str();
      doError();
      return -1;
    }
  
  bool buffer_result = !(flags&FLAG_NO_BUFFER);

  if(result)
    {
      if(buffer_result)*result = mysql_store_result(fDB);
      else *result = mysql_use_result(fDB);
    }
  else
    {
      vsassert(mysql_use_result(fDB) == 0);
    }

  if((result)&&(*result))
    {
      int rows = mysql_num_rows(*result);
      if((fLoud)&&(rows!=1))std::cerr << rows << " rows in set" << std::endl;
      else if(fLoud)std::cerr << "1 row in set" << std::endl;
      return rows;
    }
  else
    {
#if(MYSQL_VERSION_ID<32224)
      int field_count = mysql_num_fields(fDB);
#else
      int field_count = mysql_field_count(fDB);
#endif
      if(field_count==0)
	{
	  int rows = mysql_affected_rows(fDB);
	  if((fLoud)&&(rows!=1))
	    std::cerr << rows << " rows affected" << std::endl;
	  else if(fLoud)
	    std::cerr << "1 row affected" << std::endl;
	  return rows;
	}
      else
	{
	  std::string fn = "mysql_use_result";
	  if((result)&&(buffer_result))fn="mysql_store_result";
	  std::ostringstream stream;
	  stream << fn << ": " << mysql_error(fDB) << std::endl
		 << fn << ": SQL: " << fQuery;
	  fError=stream.str();
	  doError();
	  return -1;
	}
    }

  vsassert(0);
}

/*!
  Test whether mysql error code corresponds to an error that could be
  due to the existance or non-existance of tables in the database so they
  can be ignored when the FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST flag is given.

  \param error_number MySQL error number to test
*/

bool VSDBMySQLBase::isErrorExistNotExist(unsigned int error_number)
{
  switch(error_number)
    {
    case ER_DB_CREATE_EXISTS:    // 1007
    case ER_DB_DROP_EXISTS:      // 1008
    case ER_TABLE_EXISTS_ERROR:  // 1050
    case ER_BAD_TABLE_ERROR:     // 1051
    case ER_NO_SUCH_TABLE:       // 1146
      return true;
    default:
      return false;
    }
}


/*!

  If the ignore errors option was set when the class was instantiated then
  this function returns without doing anything. If not then the error is
  printed to standard output and the program dies.

 */
void VSDBMySQLBase::doError() throw()
{
  if(!fIgnoreError)
    {
      std::cerr << fError << std::endl;
      exit(EXIT_FAILURE);
    }
}

#endif // VSCONFIG_NO_MYSQL

// ****************************************************************************
// ****************************************************************************

// TEST CODE

// ****************************************************************************
// ****************************************************************************

#ifdef TEST_MAIN

// Complete the virtual functions so we can instantiate VSDBMySQLBase,
// also make the constructor public for the same reason.
class Test: public VSDBMySQLBase
{
public:
  Test(bool loud=false, bool ignore_error=false,
       const std::string& mysql_hostname="", 
       const std::string& mysql_username="",
       const std::string& mysql_password="",
       unsigned mysql_port=0, 
       const std::string& mysql_unixpath="", 
       unsigned long mysql_client_flag=0,
       bool parameters_override_all=false)
    throw(std::string)
    : VSDBMySQLBase(loud,ignore_error,
		    mysql_hostname,mysql_username,mysql_password,
		    mysql_port, mysql_unixpath, mysql_client_flag, 
		    parameters_override_all) {}
  virtual VSDBStatement* createQuery(const std::string& query, 
				     uint32_t flags = 0) throw();
};

VSDBStatement* Test::createQuery(const std::string& query, uint32_t flags) 
  throw()
{
  return 0;
}

main(int argc, char** argv)
{
  VSOptions options(argc,argv);
  VSDBMySQLBase::configureFromCommandLine(options);
  Test test;  
  std::cout << test.getClientVersion() << std::endl
	    << test.getServerVersion() << std::endl;
}
#endif

