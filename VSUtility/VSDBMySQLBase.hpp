//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBMySQLBase.hpp
  Base class for MySQL database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
*/

#ifndef VSDBMYSQLBASE_HPP
#define VSDBMYSQLBASE_HPP

#ifndef VSCONFIG_NO_MYSQL

#include <memory>

#include <mysql.h>

#include "VSDatabase.hpp"
#include "VSOptions.hpp"
#include "VSDataConverter.hpp"

namespace VERITAS
{
  //! Base class for MySQL databases implementing common fuctionality
  class VSDBMySQLBase: public VSDatabase
  {
  public:
    //! Destructor closes connection to database
    virtual ~VSDBMySQLBase();

    //! Return the error message reported by the last call
    virtual std::string getErrorMessage() throw();

    //! Return the database manager
    virtual std::string getRDBM() throw();

    //! Return the hostname we are using
    virtual std::string getHostname() throw();
    
    //! Return the username we are using
    virtual std::string getUsername() throw();

    //! Return the password we are using
    virtual std::string getPassword() throw();

    //! Return the database we are using
    virtual std::string getDatabase() throw();

    //! Create a database if it does not exist
    virtual int createDatabase(const std::string& dbname,
			       uint32_t flags = 0) throw();
    
    //! Delete a database if it exists
    virtual int dropDatabase(const std::string& dbname,
			     uint32_t flags = 0) throw();
    
    //! Select a database
    virtual int useDatabase(const std::string& dbname,
			    uint32_t flags = 0) throw();
    
    //! Create a table in the database
    virtual int createTable(const std::string& table, const std::string& spec,
			    uint32_t flags = 0) throw();
    
    //! Drop a table from the database
    virtual int dropTable(const std::string& table,
			  uint32_t flags = 0) throw();
    
    //! Delete entries from table
    virtual int deleteFromTable(const std::string& table,
				const std::string& cond,
				uint32_t flags = 0) throw();

    //! Create an INSERT query
    virtual VSDBStatement* createInsertQuery(const std::string& table, 
					     unsigned num_columns,
					     const std::string& columns = "",
					     uint32_t flags = 0) throw();

    //! Return the MySQL DB handle
    MYSQL* mysql() throw() { return fDB; }

    //! Return the SQL for the last executed statement
    const std::string& sql() const throw() { return fQuery; }

    //! Return the MySQL error message
    std::string error() const throw() { return mysql_error(fDB); }

    //! Configure the default database connection from the command line options
    static void configureFromCommandLine(VERITAS::VSOptions& options) throw();

    //! Return the MySQL version of the client library
    int getClientVersion() throw();

    //! Return the MySQL version of the server
    int getServerVersion() throw();

    //! Function to return SQL type corresponding to bool
    virtual std::string sqlTypeOf(const bool& x) const throw();
    //! Function to return SQL type corresponding to int8_t
    virtual std::string sqlTypeOf(const int8_t& x) const throw();
    //! Function to return SQL type corresponding to uint8_t
    virtual std::string sqlTypeOf(const uint8_t& x) const throw();
    //! Function to return SQL type corresponding to int16_t
    virtual std::string sqlTypeOf(const int16_t& x) const throw();
    //! Function to return SQL type corresponding to uint16_t
    virtual std::string sqlTypeOf(const uint16_t& x) const throw();
    //! Function to return SQL type corresponding to int32_t
    virtual std::string sqlTypeOf(const int32_t& x) const throw();
    //! Function to return SQL type corresponding to uint32_t
    virtual std::string sqlTypeOf(const uint32_t& x) const throw();
    //! Function to return SQL type corresponding to int64_t
    virtual std::string sqlTypeOf(const int64_t& x) const throw();
    //! Function to return SQL type corresponding to uint64_t
    virtual std::string sqlTypeOf(const uint64_t& x) const throw();
    //! Function to return SQL type corresponding to float
    virtual std::string sqlTypeOf(const float& x) const throw();
    //! Function to return SQL type corresponding to double
    virtual std::string sqlTypeOf(const double& x) const throw();
    //! Function to return SQL type corresponding to C-string
    virtual std::string sqlTypeOf(const char* x, unsigned length=0) 
      const throw();
    //! Function to return SQL type corresponding to STL-string
    virtual std::string sqlTypeOf(const std::string& x, unsigned length=0) 
      const throw();
    //! Function to return SQL type corresponding to VSDBDate
    virtual std::string sqlTypeOf(const VSDBDate& x) const throw();
    //! Function to return SQL type corresponding to VSDBTime
    virtual std::string sqlTypeOf(const VSDBTime& x) const throw();
    //! Function to return SQL type corresponding to VSDBDateTime
    virtual std::string sqlTypeOf(const VSDBDateTime& x) const throw();
    //! Function to return SQL type corresponding to VSDBTimestamp
    virtual std::string sqlTypeOf(const VSDBTimestamp& x) const throw();

  protected:
    //! Protected constructor
    VSDBMySQLBase(bool loud=false, bool ignore_error=false,
		  const std::string& mysql_hostname="", 
		  const std::string& mysql_username="",
		  const std::string& mysql_password="",
		  unsigned mysql_port=0, 
		  const std::string& mysql_unixpath="", 
		  unsigned long mysql_client_flag=0,
		  bool parameters_override_all=false)
      throw(std::string);

    //! Internal function to execute the query and check for errors
    int executeQuery(uint32_t flags=0, MYSQL_RES** result=0) throw();

    //! Internal function to test for EXIST/NON EXIST error codes
    bool isErrorExistNotExist(unsigned int error_number);

    //! Implement the error handling policy
    void doError() throw();

  protected:
    VSDBMySQLBase(const VSDBMySQLBase&);
    VSDBMySQLBase& operator= (const VSDBMySQLBase&);

    static std::string       sDemandHostname;
    static std::string       sDemandUsername;
    static std::string       sDemandPassword;
    static std::string       sDemandPortNStr;
    static std::string       sDemandUnixPath;

    std::string              fHostname;
    std::string              fUsername;
    std::string              fPassword;
    unsigned                 fPortNumber;
    std::string              fUnixPath;
    std::string              fDatabase;

    MYSQL*                   fDB;
    std::string              fQuery;
    std::string              fError;
    bool                     fLoud;
    bool                     fIgnoreError;
  };

  // ==========================================================================
  // ==========================================================================
  //
  // Data conversion functions
  //
  // ==========================================================================
  // ==========================================================================

  // Ideally these would be members of the VSDBMySQLBase class but there is no
  // such thing as specialized member functions in C++ 

  // --------------------------------------------------------------------------
  // MySQLValueToString
  // --------------------------------------------------------------------------
  
  //! Template function to convert a value to its SQL representation
  template<class T> inline std::string MySQLValueToString(const T& x)
  {
    std::string s;
    VSDatumConverter<T>::toString(s,x);
    return s;
  }
    
  //! Overloaded function to convert a char* string to its SQL representation
  inline std::string MySQLValueToString(const char* x, 
					unsigned int count)
  {
    if(count==0)count=strlen(x);
    std::string s;
    s += '\'';
    const char* y = x;
    for(unsigned i=0; (*y)&&(i<count); i++,y++)
      {
	switch(*y)
	  {
	  case '\\': 
	  case '\'':
	    // escape ' and \ in the string
	    s += '\\';
	    // fall through
	  default:
	    s += *y;
	  }
      }
    s += '\'';
    return s;
  }
  
  //! Specialization of template to convert an STL string to its 
  // SQL representation
  template<> inline std::string MySQLValueToString<>(const std::string& x)
  {
    return MySQLValueToString(x.c_str(), x.size());
  }

  template<> inline std::string MySQLValueToString<>(const VSDBDate& x)
  {
    char buffer[13];
    snprintf(buffer,13,"'%04u-%02u-%02u'",x.year,x.month,x.day);
    return std::string(buffer);
  }

  template<> inline std::string MySQLValueToString<>(const VSDBTime& x)
  {
    char buffer[11];
    snprintf(buffer,11,"'%02u:%02u:%02u'",x.hour,x.minute,x.second);
    return std::string(buffer);
  }

  template<> inline std::string MySQLValueToString<>(const VSDBDateTime& x)
  {
    char buffer[22];
    snprintf(buffer,22,"'%04u-%02u-%02u %02u:%02u:%02u'",
	     x.year,x.month,x.day,x.hour,x.minute,x.second);
    return std::string(buffer);
  }

  template<> inline std::string MySQLValueToString<>(const VSDBTimestamp& x)
  {
    if(x.current)return std::string("NULL");
    char buffer[22];
    snprintf(buffer,22,"'%04u-%02u-%02u %02u:%02u:%02u'",
	     x.year,x.month,x.day,x.hour,x.minute,x.second);
    return std::string(buffer);
  }
  
  // --------------------------------------------------------------------------
  // MySQLValueFromString
  // --------------------------------------------------------------------------
  
  //! Template to extract a value from an SQL representation
  template<class T> inline bool MySQLValueFromString(const char* s, T& x)
  {
    return VSDatumConverter<T>::fromString(x,s);
  }
  
  //! Specialization of template to extract an STL string from 
  // SQL representation
  template<> inline bool MySQLValueFromString<>(const char* s, VSDBDate& x)
  {
    x.year        = strtoul(s,NULL,10);
    x.month       = strtoul(s+5,NULL,10);
    x.day         = strtoul(s+8,NULL,10);
    return true;
  }

  //! Specialization of template to extract an STL string from 
  // SQL representation
  template<> inline bool MySQLValueFromString<>(const char* s, VSDBTime& x)
  {
    x.hour        = strtoul(s,NULL,10);
    x.minute      = strtoul(s+3,NULL,10);
    x.second      = strtoul(s+6,NULL,10);
    x.microsecond = 0;
    return true;
  }

  //! Specialization of template to extract an STL string from 
  // SQL representation
  template<> inline bool MySQLValueFromString<>(const char* s, VSDBDateTime& x)
  {
    x.year        = strtoul(s,NULL,10);
    x.month       = strtoul(s+5,NULL,10);
    x.day         = strtoul(s+8,NULL,10);
    x.hour        = strtoul(s+11,NULL,10);
    x.minute      = strtoul(s+14,NULL,10);
    x.second      = strtoul(s+17,NULL,10);
    x.microsecond = 0;
    return true;
  }

  //! Specialization of template to extract an STL string from 
  // SQL representation
  template<> inline bool MySQLValueFromString<>(const char* s,
						VSDBTimestamp& x)
  {
#if(MYSQL_VERSION_ID<40100)
    sscanf(s,"%4hu%2hu%2hu%2hu%2hu%2hu",
	   &x.year, &x.month, &x.day, &x.hour, &x.minute, &x.second);    
#else
    x.year        = strtoul(s,NULL,10);
    x.month       = strtoul(s+5,NULL,10);
    x.day         = strtoul(s+8,NULL,10);
    x.hour        = strtoul(s+11,NULL,10);
    x.minute      = strtoul(s+14,NULL,10);
    x.second      = strtoul(s+17,NULL,10);
#endif
    x.microsecond = 0;
    x.current     = false;
    return true;
  }
  
  //! Overloaded function to extract a char* string from SQL representation
  inline bool MySQLValueFromString(const char* s, char* x, 
				   unsigned long length, unsigned long* count)
  {
    size_t ss=strlen(s);
    strncpy(x, s, length);
    if((ss>length)&&(count))*count=length;
    else if(count)*count=ss;
    return true;
  }
  
}

#endif // VSCONFIG_NO_MYSQL
#endif // VSDBMYSQLBASE_HPP
