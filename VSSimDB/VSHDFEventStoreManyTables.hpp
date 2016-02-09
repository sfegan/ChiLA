//-*-mode:c++; mode:font-lock;-*-

/*! \file VSHDFEventStoreManyTables.hpp

  Classes for writing event data to HDF5 files with DB assistance

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       22/10/2006
  \note
*/

#ifndef VSHDFEVENTSTOREMANYTABLES_HPP
#define VSHDFEVENTSTOREMANYTABLES_HPP

#include <map>
#include <set>
#include <string>

#include <VSSimDB.hpp>
#include <VSHDFEventStore.hpp>
#include <VSOctaveIO.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSHDFEventStoreManyTables: public VSSimDBEventStore
  {
  public:
    VSHDFEventStoreManyTables(VSSimDB* simdb);
    virtual ~VSHDFEventStoreManyTables();
    virtual int insertEvent(VSSimDBEventData& event);
    virtual int insertScope(const VSSimDBScopeData& scope);
    virtual int insertPE(const VSSimDBPEData& pe);
    virtual int updateHitScopes(const VSSimDBEventData& event);
    virtual int updateHitPixels(const VSSimDBScopeData& scope);
  private:
    unsigned findTable(const VSSimDBEventData& event);
    void openFile(unsigned table_id);

    VSSimDB*                                fSimDB;
    std::map< uint32_t, VSHDFEventStore* >  fEventStore;       
    std::map< uint32_t, VSSimDBTableParam > fTableParam;
    uint32_t                                fTableID;
    std::string                             fTableName;
  };

}

#endif // VSHDFEVENTSTOREMANYTABLES_HPP
