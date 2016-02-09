//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBMySQL41.hpp
  MySQL 3.X database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/02/2005
*/

#ifndef VSDBMYSQL41_HPP
#define VSDBMYSQL41_HPP

#ifndef VSCONFIG_NO_MYSQL

#include <set>
#include <map>

#include "VSDBMySQL3X.hpp"

#if(MYSQL_VERSION_ID>=40100)

namespace VERITAS
{
  class VSDBMySQL41PreparedStatement;

  //! Input and output access to MySQL 4.1 database.
  class VSDBMySQL41: public VSDBMySQL3X
  {
  public:
    //! Constructor
    VSDBMySQL41(bool loud=false, bool ignore_error=false,
		const std::string& mysql_hostname="", 
		const std::string& mysql_username="",
		const std::string& mysql_password="",
		unsigned mysql_port=0, 
		const std::string& mysql_unixpath="", 
		unsigned long mysql_client_flag=0,
		bool parameters_override_all=false) throw(std::string);

    //! Destructor, destroys all prepared statements
    virtual ~VSDBMySQL41();
    
    //! Create a client side prepared statement from a query
    virtual VSDBStatement* createQuery(const std::string& query, 
				       uint32_t flags = 0) throw();
    
  private:
    VSDBMySQL41(const VSDBMySQL41&);
    VSDBMySQL41& operator= (const VSDBMySQL41&);
    
    friend class VSDBMySQL41PreparedStatement; // ugh

    //! Internal function which executes a prepared statement
    int psExecute(VSDBMySQL41PreparedStatement* ps, MYSQL_BIND* bind);

    //! Internal function which fetches data for next row
    int psFetch(VSDBMySQL41PreparedStatement* ps, MYSQL_BIND* bind);

    //! Internal function which gets the AUTO_INCREMENT id from the last query
    uint64_t psGetInsertID(VSDBMySQL41PreparedStatement* ps);

    //! Internal function to remove links to a prepared statement
    int psDelete(VSDBMySQL41PreparedStatement* ps);

    int                                     fServerVersion;
    std::set<VSDBMySQL41PreparedStatement*> fStatements;
  };

  //! Client-side prepared statement class for MySQL 4.XX
  class VSDBMySQL41PreparedStatement: public VSDBStatement
  {
  public:
    //! Destructor
    virtual ~VSDBMySQL41PreparedStatement();

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
    // Protected constructor called from VSDBMySQL41::psCreate()
    VSDBMySQL41PreparedStatement(VSDBMySQL41* mysql, const std::string& query,
				 MYSQL_STMT* stmt, MYSQL_RES* meta, 
				 uint32_t flags);

  private:
    VSDBMySQL41PreparedStatement(const VSDBMySQL41PreparedStatement&);
    VSDBMySQL41PreparedStatement& 
    operator= (const VSDBMySQL41PreparedStatement&);

    friend class VSDBMySQL41; // ugh

    union MirrorPointer
    {
      void*                v_ptr;
      uint8_t*             u8_ptr;
      char*                c_str_ptr;
      MYSQL_TIME*          my_time_ptr;
    };

    struct BoundParam
    {
      VSDatabase::DataType fDataType;
      ConstMultiPointer    fUserPtr;
      const unsigned long* fUserCount;
      unsigned long        fUserBufferCount;
      const bool*          fUserNullPtr;
      
      MirrorPointer        fMirrorPtr;
      unsigned long        fMirrorCount;
      unsigned long        fMirrorBufferCount;
      my_bool              fMirrorNull;
      
      BoundParam(): 
	fDataType(), 
	fUserPtr(), fUserCount(), fUserBufferCount(), fUserNullPtr(),
	fMirrorPtr(), fMirrorCount(), fMirrorBufferCount(), fMirrorNull() {}
    };

    struct BoundVal
    {
      VSDatabase::DataType fDataType;
      MultiPointer         fUserPtr;
      unsigned long*       fUserCount;
      unsigned long        fUserBufferCount;
      bool*                fUserNullPtr;
      
      MirrorPointer        fMirrorPtr;
      unsigned long        fMirrorCount;
      unsigned long        fMirrorBufferCount;
      my_bool              fMirrorNull;
      
      BoundVal(): 
	fDataType(), 
	fUserPtr(), fUserCount(), fUserBufferCount(), fUserNullPtr(),
	fMirrorPtr(), fMirrorCount(), fMirrorBufferCount(), fMirrorNull() {}
    };

    VSDBMySQL41*                         fMySQL;
    MYSQL_STMT*                          fStmt;
    MYSQL_RES*                           fMeta;
    uint32_t                             fQueryFlags;

    std::vector<BoundParam>              fParam;
    unsigned int                         fParamCount;
    bool                                 fParamReset;

    std::vector<BoundVal>                fVal;
    unsigned int                         fValCount;
    bool                                 fValReset;

    std::string                          fQuery;
  };

}

#endif // (MYSQL_VERSION_ID>=40100)
#endif // VSCONFIG_NO_MYSQL
#endif // VSDBMYSQL41_HPP
