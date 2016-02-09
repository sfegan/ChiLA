//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDatabase.hpp
  Base class for database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
*/

#ifndef VSDATABASE_HPP
#define VSDATABASE_HPP

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <stdint.h>

//! VERITAS namespace
namespace VERITAS
{
  class VSDBStatement;

  //! Time type used by VSDBStatement classes for input/output of date
  class VSDBDate
  {
  public:
    VSDBDate(): year(), month(), day() { }
    ~VSDBDate() { }
    uint16_t year;
    uint16_t month;
    uint16_t day;
  };

  //! Time type used by VSDBStatement classes for input/output of time
  class VSDBTime
  {
  public:
    VSDBTime(): hour(), minute(), second(), microsecond() { }
    ~VSDBTime() { }
    uint16_t hour;
    uint16_t minute;
    uint16_t second;
    uint32_t microsecond;
  };

  //! Time type used by VSDBStatement classes for input/output of date and time
  class VSDBDateTime: public VSDBDate, public VSDBTime
  {
  public:
    VSDBDateTime(): VSDBDate(), VSDBTime() { }
    ~VSDBDateTime() { }
  };

  //! Time type used by VSDBStatement classes for input/output of timestamps
  class VSDBTimestamp: public VSDBDate, public VSDBTime
  {
  public:
    VSDBTimestamp(): VSDBDate(), VSDBTime(), current() { }
    ~VSDBTimestamp() { }
    bool     current;
  };
  
  //! Base for classes providing database access
  class VSDatabase
  {
  public:
    //! Request that results of query not be buffered
    static const uint32_t FLAG_NO_BUFFER                          = 0x00000001;
    //! Request that query not generate server-side prepared statement
    static const uint32_t FLAG_NO_SERVER_PS                       = 0x00000002;
    //! Request that INSERT/CREATE not return error on key/table duplication
    static const uint32_t FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST     = 0x00000004;

    //! Virtual destructor
    virtual ~VSDatabase();

    //! Return the error message reported by the last call
    virtual std::string getErrorMessage() throw() = 0;

    //! Return the database manager
    virtual std::string getRDBM() throw() = 0;

    //! Return the hostname we are using
    virtual std::string getHostname() throw() = 0;
    
    //! Return the username we are using
    virtual std::string getUsername() throw() = 0;

    //! Return the password we are using
    virtual std::string getPassword() throw() = 0;

    //! Return the database we are using
    virtual std::string getDatabase() throw() = 0;

    //! Create a database if it does not exist
    virtual int createDatabase(const std::string& dbname,
			       uint32_t flags = 0) throw () = 0;

    //! Delete a database if it exists
    virtual int dropDatabase(const std::string& dbname,
			     uint32_t flags = 0) throw() = 0;

    //! Select a database
    virtual int useDatabase(const std::string& dbname,
			    uint32_t flags = 0) throw() = 0;

    //! Create a table in the database
    virtual int createTable(const std::string& table, const std::string& spec,
			    uint32_t flags = 0) throw() = 0;

    //! Drop a table from the database
    virtual int dropTable(const std::string& table, 
			  uint32_t flags = 0) throw() = 0;

    //! Delete entries from table
    virtual int deleteFromTable(const std::string& table,
				const std::string& cond,
				uint32_t flags = 0) throw() = 0;
    
    //! Create a prepared statement from a query
    virtual VSDBStatement* createQuery(const std::string& query, 
				       uint32_t flags = 0) throw() = 0;

    //! Utility function to create a prepared statement from a SELECT query
    virtual VSDBStatement* createSelectQuery(const std::string& table, 
					     const std::string& condition = "",
					     const std::string& columns = "*",
					     uint32_t flags = 0) throw();
    
    //! Utility function to create a prepared statement from an INSERT query
    virtual VSDBStatement* createInsertQuery(const std::string& table, 
					     unsigned num_columns,
					     const std::string& columns = "",
					     uint32_t flags = 0) throw() = 0;

    //! Function to return SQL type corresponding to bool
    virtual std::string sqlTypeOf(const bool& x) const throw() = 0;
    //! Function to return SQL type corresponding to int8_t
    virtual std::string sqlTypeOf(const int8_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to uint8_t
    virtual std::string sqlTypeOf(const uint8_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to int16_t
    virtual std::string sqlTypeOf(const int16_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to uint16_t
    virtual std::string sqlTypeOf(const uint16_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to int32_t
    virtual std::string sqlTypeOf(const int32_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to uint32_t
    virtual std::string sqlTypeOf(const uint32_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to int64_t
    virtual std::string sqlTypeOf(const int64_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to uint64_t
    virtual std::string sqlTypeOf(const uint64_t& x) const throw() = 0;
    //! Function to return SQL type corresponding to float
    virtual std::string sqlTypeOf(const float& x) const throw() = 0;
    //! Function to return SQL type corresponding to double
    virtual std::string sqlTypeOf(const double& x) const throw() = 0;
    //! Function to return SQL type corresponding to C-string
    virtual std::string sqlTypeOf(const char* x, unsigned length=0) 
      const throw() = 0;
    //! Function to return SQL type corresponding to STL-string
    virtual std::string sqlTypeOf(const std::string& x, unsigned length=0) 
      const throw() = 0;
    //! Function to return SQL type corresponding to VSDBDate
    virtual std::string sqlTypeOf(const VSDBDate& x) const throw() = 0;
    //! Function to return SQL type corresponding to VSDBTime
    virtual std::string sqlTypeOf(const VSDBTime& x) const throw() = 0;
    //! Function to return SQL type corresponding to VSDBDateTime
    virtual std::string sqlTypeOf(const VSDBDateTime& x) const throw() = 0;
    //! Function to return SQL type corresponding to VSDBTimestamp
    virtual std::string sqlTypeOf(const VSDBTimestamp& x) const throw() = 0;

    //! Function for constructing table specifications
    std::string sqlSpec(const std::string& name, const std::string& type,
			bool first, const std::string& post="") const throw();

    //! Templated convenience function for constructing table specifications
    template<typename T> inline
    std::string sqlSpecOf(const std::string& name, const T& x,
			  bool first, const std::string& post="")const throw();

    enum DataType { DT_BOOL, DT_INT8, DT_UINT8, DT_INT16, DT_UINT16,
		    DT_INT32, DT_UINT32, DT_INT64, DT_UINT64, 
		    DT_FLOAT, DT_DOUBLE, DT_STLSTRING, DT_CSTRING,
		    DT_DATE, DT_TIME, DT_DATETIME, DT_TIMESTAMP,
		    DT_UNKNOWN};

  protected:
    //! Constructor
    VSDatabase() { }    

  private:
    // Copy protection for all derived classes
    VSDatabase(const VSDatabase&);
    const VSDatabase& operator =(const VSDatabase&);
  };

  //! Base class for query statement
  class VSDBStatement
  {
  public:
    //! Virtual destructor
    virtual ~VSDBStatement();

    //! Return the error message reported by the last call
    virtual std::string getErrorMessage() throw() = 0;

    //! Execute the query using the values in the bound parameter variables
    virtual int execute() throw() = 0;

    //! Fetch the next row from the result set into the bound result variables
    virtual int retrieveNextRow() throw() = 0;

    //! Get the value of the AUTOINCREMENT variable in the last query
    virtual uint64_t getInsertID() throw() = 0;

    //! Clear all the bound parameters
    virtual void clearBoundParams() throw() = 0;
    //! Bind char to next parameter
    virtual void bindToParam(const bool& val, const bool* null = 0) 
      throw() = 0;
    virtual void bindToParam(const int8_t& val, const bool* null = 0) 
      throw() = 0;
    //! Bind unsigned char to next parameter
    virtual void bindToParam(const uint8_t& val, const bool* null = 0) 
      throw () = 0;
    //! Bind short to next parameter
    virtual void bindToParam(const int16_t& val, const bool* null = 0)
      throw() = 0;
    //! Bind unsigned short to next parameter
    virtual void bindToParam(const uint16_t& val, const bool* null = 0) 
      throw() = 0;
    //! Bind int to next parameter
    virtual void bindToParam(const int32_t& val, const bool* null = 0)
      throw() = 0;
    //! Bind unsigned int to next parameter
    virtual void bindToParam(const uint32_t& val, const bool* null = 0)
      throw() = 0;
    //! Bind unsigned long long to next parameter
    virtual void bindToParam(const int64_t& val, const bool* null = 0)
      throw() = 0;
    //! Bind unsigned long long to next parameter
    virtual void bindToParam(const uint64_t& val, const bool* null = 0)
      throw() = 0;
    //! Bind float to next parameter
    virtual void bindToParam(const float& val, const bool* null = 0)
      throw() = 0;
    //! Bind double to next parameter
    virtual void bindToParam(const double& val, const bool* null = 0)
      throw() = 0;
    //! Bind STL string to next parameter
    virtual void bindToParam(const std::string& val, const bool* null = 0)
      throw() = 0;
    //! Bind char* string to next parameter
    virtual void bindToParam(const char* val, const unsigned long* count = 0,
			     const bool* null = 0)
      throw() = 0;
    //! Bind VSDBDate to next parameter
    virtual void bindToParam(const VSDBDate& val, const bool* null = 0)
      throw() = 0;
    //! Bind VSDBTime to next parameter
    virtual void bindToParam(const VSDBTime& val, const bool* null = 0)
      throw() = 0;
    //! Bind VSDBDateTime to next parameter
    virtual void bindToParam(const VSDBDateTime& val, const bool* null = 0)
      throw() = 0;
    //! Bind VSDBTimestamp to next parameter
    virtual void bindToParam(const VSDBTimestamp& val, const bool* null = 0)
      throw() = 0;
        
    //! Clear all the bound results
    virtual void clearBoundResults()throw() = 0;
    //! Bind char to next result
    virtual void bindToResult(bool& val, bool* null = 0) 
      throw() = 0;
    virtual void bindToResult(int8_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind unsigned char to next result
    virtual void bindToResult(uint8_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind short to next result
    virtual void bindToResult(int16_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind unsigned short to next result
    virtual void bindToResult(uint16_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind int to next result
    virtual void bindToResult(int32_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind unsigned int to next result
    virtual void bindToResult(uint32_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind long long to next result
    virtual void bindToResult(int64_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind unsigned long long to next result
    virtual void bindToResult(uint64_t& val, bool* null = 0) 
      throw() = 0;
    //! Bind float to next result
    virtual void bindToResult(float& val, bool* null = 0) 
      throw() = 0;
    //! Bind double to next result
    virtual void bindToResult(double& val, bool* null = 0) 
      throw() = 0;
    //! Bind STL string to next result
    virtual void bindToResult(std::string& val, bool* null = 0) 
      throw() = 0;
    //! Bind char* (STRING) to next result
    virtual void bindToResult(char* val, unsigned long buffer_count, 
			      unsigned long* count = 0, bool* null = 0)
      throw() = 0;
    //! Bind VSDBDate to next parameter
    virtual void bindToResult(VSDBDate& val, bool* null = 0)
      throw() = 0;
    //! Bind VSDBTime to next parameter
    virtual void bindToResult(VSDBTime& val, bool* null = 0)
      throw() = 0;
    //! Bind VSDBDateTime to next parameter
    virtual void bindToResult(VSDBDateTime& val, bool* null = 0)
      throw() = 0;
    //! Bind VSDBTimestamp to next parameter
    virtual void bindToResult(VSDBTimestamp& val, bool* null = 0)
      throw() = 0;

  protected:
    //! Internal constructor
    VSDBStatement() { }

    //! Internal union which cuts down on casting
    union MultiPointer
    {
      void*                  v_ptr;
      bool*                  bool_ptr;
      int8_t*                i8_ptr;
      uint8_t*               u8_ptr;
      int16_t*               i16_ptr;
      uint16_t*              u16_ptr;
      int32_t*               i32_ptr;
      uint32_t*              u32_ptr;
      int64_t*               i64_ptr;
      uint64_t*              u64_ptr;
      float*                 flt_ptr;
      double*                dbl_ptr;
      std::string*           stl_str_ptr;
      char*                  c_str_ptr;
      VSDBDate*              date_ptr;
      VSDBTime*              time_ptr;
      VSDBDateTime*          datetime_ptr;
      VSDBTimestamp*         timestamp_ptr;
    };

    //! Internal union which cuts down on casting
    union ConstMultiPointer
    {
      const void*            v_ptr;
      const bool*            bool_ptr;
      const int8_t*          i8_ptr;
      const uint8_t*         u8_ptr;
      const int16_t*         i16_ptr;
      const uint16_t*        u16_ptr;
      const int32_t*         i32_ptr;
      const uint32_t*        u32_ptr;
      const int64_t*         i64_ptr;
      const uint64_t*        u64_ptr;
      const float*           flt_ptr;
      const double*          dbl_ptr;
      const std::string*     stl_str_ptr;
      const char*            c_str_ptr;
      const VSDBDate*        date_ptr;
      const VSDBTime*        time_ptr;
      const VSDBDateTime*    datetime_ptr;
      const VSDBTimestamp*   timestamp_ptr;
    };

  private:
    // Copy protection for all derived classes
    VSDBStatement(const VSDBStatement&);
    const VSDBStatement& operator =(const VSDBStatement&);
  };

  template<typename T> inline
  std::string VSDatabase::sqlSpecOf(const std::string& name, const T& x,
				    bool first, const std::string& post)
    const throw()
  {
    return sqlSpec(name,sqlTypeOf(x),first,post);
  }
    
  // --------------------------------------------------------------------------
  // parseSQLForParameters
  // --------------------------------------------------------------------------

  //! Parse SQL query for paramater markers
  bool parseSQLForParameters(const std::string& query, 
			     std::vector<std::string>& query_bits);


}

//! Write date to stream in format of <code>YYYY-MM-DD</code>
std::ostream& operator<<(std::ostream& stream,const VERITAS::VSDBDate& x);
//! Write time to stream in format of <code>HH:MM:SS</code>
std::ostream& operator<<(std::ostream& stream,const VERITAS::VSDBTime& x);
//! Write date and time to stream in format of <code>YYYY-MM-DD HH:MM:SS</code>
std::ostream& operator<<(std::ostream& stream,const VERITAS::VSDBDateTime& x);
//! Write timestamp to stream in format of <code>YYYY-MM-DD HH:MM:SS</code>
std::ostream& operator<<(std::ostream& stream,const VERITAS::VSDBTimestamp& x);

#endif // VSDATABASE_HPP
