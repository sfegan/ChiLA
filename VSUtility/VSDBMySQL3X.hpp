//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBMySQL3X.hpp
  MySQL 3.X database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
*/

#ifndef VSDBMYSQL3X_HPP
#define VSDBMYSQL3X_HPP

#ifndef VSCONFIG_NO_MYSQL

#include <set>

#include "VSDBMySQLBase.hpp"

namespace VERITAS
{
  class VSDBMySQL3XPreparedStatement;

  //! Input and output access to MySQL 3.XX database.
  class VSDBMySQL3X: public VSDBMySQLBase
  {
  public:
    //! Constructor
    VSDBMySQL3X(bool loud=false, bool ignore_error=false,
		const std::string& mysql_hostname="", 
		const std::string& mysql_username="",
		const std::string& mysql_password="",
		unsigned mysql_port=0, 
		const std::string& mysql_unixpath="", 
		unsigned long mysql_client_flag=0,
		bool parameters_override_all=false) throw(std::string);

    //! Destructor, destroys all prepared statements
    virtual ~VSDBMySQL3X();
    
    //! Create a client side prepared statement from a query
    virtual VSDBStatement* createQuery(const std::string& query, 
				       uint32_t flags = 0) throw();
    
  private:
    VSDBMySQL3X(const VSDBMySQL3X&);
    VSDBMySQL3X& operator= (const VSDBMySQL3X&);

    friend class VSDBMySQL3XPreparedStatement; // ugh

    //! Internal function which executes a prepared statement
    int psExecute(VSDBMySQL3XPreparedStatement* ps, const std::string& query);

    //! Internal function which fetches data for next row
    MYSQL_ROW psFetch(VSDBMySQL3XPreparedStatement* ps);

    //! Internal function which gets the AUTO_INCREMENT id from the last query
    uint64_t psGetInsertID(VSDBMySQL3XPreparedStatement* ps);

    //! Internal function to remove links to a prepared statement
    int psDelete(VSDBMySQL3XPreparedStatement* ps);

    std::set<VSDBMySQL3XPreparedStatement*> fStatements;
  };

  //! Client-side prepared statement class for MySQL 3.XX
  class VSDBMySQL3XPreparedStatement: public VSDBStatement
  {
  public:
    //! Destructor
    virtual ~VSDBMySQL3XPreparedStatement();

    //! Return the error message reported by the last call
    virtual std::string getErrorMessage() throw();

    //! Execute the query using the values in the bound parameter variables
    virtual int execute() throw();

    //! Fetch the next row from the result set into the bound result variables
    virtual int retrieveNextRow() throw();

    //! Get the value of the AUTOINCREMENT variable in the last query
    virtual uint64_t getInsertID() throw();

    //! Clear all the bound parameters
    virtual void clearBoundParams() throw();
    //! Bind char to next parameter
    virtual void bindToParam(const bool& val, const bool* null = 0) 
      throw();
    //! Bind char to next parameter
    virtual void bindToParam(const int8_t& val, const bool* null = 0) 
      throw();
    //! Bind unsigned char to next parameter
    virtual void bindToParam(const uint8_t& val, const bool* null = 0) 
      throw ();
    //! Bind short to next parameter
    virtual void bindToParam(const int16_t& val, const bool* null = 0)
      throw();
    //! Bind unsigned short to next parameter
    virtual void bindToParam(const uint16_t& val, const bool* null = 0) 
      throw();
    //! Bind int to next parameter
    virtual void bindToParam(const int32_t& val, const bool* null = 0)
      throw();
    //! Bind unsigned int to next parameter
    virtual void bindToParam(const uint32_t& val, const bool* null = 0)
      throw();
    //! Bind unsigned long long to next parameter
    virtual void bindToParam(const int64_t& val, const bool* null = 0)
      throw();
    //! Bind unsigned long long to next parameter
    virtual void bindToParam(const uint64_t& val, const bool* null = 0)
      throw();
    //! Bind float to next parameter
    virtual void bindToParam(const float& val, const bool* null = 0)
      throw();
    //! Bind double int to next parameter
    virtual void bindToParam(const double& val, const bool* null = 0)
      throw();
    //! Bind STL string int to next parameter
    virtual void bindToParam(const std::string& val, const bool* null = 0)
      throw();
    //! Bind char* string to next parameter
    virtual void bindToParam(const char* val, const unsigned long* count = 0, 
			     const bool* null = 0)
      throw();
    //! Bind VSDBDate to next parameter
    virtual void bindToParam(const VSDBDate& val, const bool* null = 0)
      throw();
    //! Bind VSDBTime to next parameter
    virtual void bindToParam(const VSDBTime& val, const bool* null = 0)
      throw();
    //! Bind VSDBDateTime to next parameter
    virtual void bindToParam(const VSDBDateTime& val, const bool* null = 0)
      throw();
    //! Bind VSDBTimestamp to next parameter
    virtual void bindToParam(const VSDBTimestamp& val, const bool* null = 0)
      throw();
    
    //! Clear all the bound results
    virtual void clearBoundResults()throw();
    //! Bind char to next result
    virtual void bindToResult(bool& val, bool* null = 0) 
      throw();
    //! Bind char to next result
    virtual void bindToResult(int8_t& val, bool* null = 0) 
      throw();
    //! Bind unsigned char to next result
    virtual void bindToResult(uint8_t& val, bool* null = 0) 
      throw();
    //! Bind short to next result
    virtual void bindToResult(int16_t& val, bool* null = 0) 
      throw();
    //! Bind unsigned short to next result
    virtual void bindToResult(uint16_t& val, bool* null = 0) 
      throw();
    //! Bind int to next result
    virtual void bindToResult(int32_t& val, bool* null = 0) 
      throw();
    //! Bind unsigned int to next result
    virtual void bindToResult(uint32_t& val, bool* null = 0) 
      throw();
    //! Bind long long to next result
    virtual void bindToResult(int64_t& val, bool* null = 0) 
      throw();
    //! Bind unsigned long long to next result
    virtual void bindToResult(uint64_t& val, bool* null = 0) 
      throw();
    //! Bind float to next result
    virtual void bindToResult(float& val, bool* null = 0) 
      throw();
    //! Bind double to next result
    virtual void bindToResult(double& val, bool* null = 0) 
      throw();
    //! Bind STL string to next result
    virtual void bindToResult(std::string& val, bool* null = 0) 
      throw();
    //! Bind char* string to next result
    virtual void bindToResult(char* val, unsigned long buffer_count, 
			      unsigned long* count = 0, bool* null = 0)
      throw();
    //! Bind VSDBDate to next parameter
    virtual void bindToResult(VSDBDate& val, bool* null = 0)
      throw();
    //! Bind VSDBTime to next parameter
    virtual void bindToResult(VSDBTime& val, bool* null = 0)
      throw();
    //! Bind VSDBDateTime to next parameter
    virtual void bindToResult(VSDBDateTime& val, bool* null = 0)
      throw();
    //! Bind VSDBTimestamp to next parameter
    virtual void bindToResult(VSDBTimestamp& val, bool* null = 0)
      throw();

  protected:
    // Protected constructor called from VSDBMySQL3X::psCreate()
    VSDBMySQL3XPreparedStatement(VSDBMySQL3X* mysql, const std::string& query,
				 uint32_t flags);

  private:
    VSDBMySQL3XPreparedStatement(const VSDBMySQL3XPreparedStatement&);
    VSDBMySQL3XPreparedStatement& 
    operator= (const VSDBMySQL3XPreparedStatement&);

    friend class VSDBMySQL3X; // ugh
    
    struct BoundParam
    {
      VSDatabase::DataType fDataType;
      ConstMultiPointer    fDataPtr;
      const unsigned long* fDataCount;
      unsigned long        fDataBufferCount;
      const bool*          fNullPtr;
      BoundParam(): 
	fDataType(), fDataPtr(), 
	fDataCount(), fDataBufferCount(), fNullPtr() {}
    };

    struct BoundVal
    {
      VSDatabase::DataType fDataType;
      MultiPointer         fDataPtr;
      unsigned long*       fDataCount;
      unsigned long        fDataBufferCount;
      bool*                fNullPtr;
      BoundVal(): 
	fDataType(), fDataPtr(), 
	fDataCount(), fDataBufferCount(), fNullPtr() {}
    };
      
    VSDBMySQL3X*              fMySQL;
    uint32_t                  fQueryFlags;
    MYSQL_RES*                fResultSet;
    unsigned                  fResultSetColumns;
    std::vector<BoundParam>   fParam;
    bool                      fParamReset;
    std::vector<BoundVal>     fVal;
    bool                      fValReset;
    std::vector<std::string>  fQueryBits;
    std::string               fQuery;
  };

}

#endif // VSCONFIG_NO_MYSQL
#endif // VSDBMYSQL3X_HPP
