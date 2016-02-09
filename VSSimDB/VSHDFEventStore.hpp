//-*-mode:c++; mode:font-lock;-*-

/*! \file VSHDFEventStore.hpp

  Classes for writing event data to HDF5 files with DB assistance

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       22/10/2006
  \note
*/

#ifndef VSHDFEVENTSTORE_HPP
#define VSHDFEVENTSTORE_HPP

#include <map>
#include <set>
#include <string>

#include <VSSimDB.hpp>
#include <VSOctaveIO.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSHDFEventStore: public VSSimDBEventStore
  {
  public:
    VSHDFEventStore(VSSimDB* simdb, uint32_t workunit_run_id,
		    const VSSimDBTableParam* table_param,
		    VSOctaveH5WriterStruct* writer = 0);
    virtual ~VSHDFEventStore();
    virtual int insertEvent(VSSimDBEventData& event);
    virtual int insertScope(const VSSimDBScopeData& scope);
    virtual int insertPE(const VSSimDBPEData& pe);
    virtual int updateHitScopes(const VSSimDBEventData& event);
    virtual int updateHitPixels(const VSSimDBScopeData& scope);

    int getWorkunitRunID() { return fWorkunitRunID; }

    static std::string nameForWorkunitID(unsigned workunit_run_id);
  private:

    // SETTINGS
    VSSimDB*                fSimDB;
    uint32_t                fWorkunitRunID;
    VSOctaveH5WriterStruct* fWriter;
    VSOctaveH5WriterStruct* fMyWriter;

    // INTERIM SCOPE DATA STORE

    struct ChanCache
    {
      ChanCache(): PEs() { /* nothing to see here */ }
      std::multiset<float> PEs;
    };
      
    unsigned                            fScopeCount;
    std::map<uint32_t, ChanCache*>      fScopeCache;

    // WRITERS -- INFRASTRUCTURE
    VSOctaveH5WriterStruct*             fEV;
    VSOctaveH5WriterStruct*             fSC;
    VSOctaveH5WriterStruct*             fPX;
    VSOctaveH5WriterStruct*             fPE;

    // WRITERS -- EVENT
    VSOctaveH5WriterVector<uint32_t>*   fWEventID;
    VSOctaveH5WriterVector<float>*      fWTargetZenithRad;
    VSOctaveH5WriterVector<float>*      fWTargetAzimuthRad;
    VSOctaveH5WriterVector<float>*      fWPrimaryZenithRad;
    VSOctaveH5WriterVector<float>*      fWPrimaryAzimuthRad;
    VSOctaveH5WriterVector<float>*      fWPrimaryCoreEastM;
    VSOctaveH5WriterVector<float>*      fWPrimaryCoreNorthM;
    VSOctaveH5WriterVector<float>*      fWPrimaryCoreUpASLM;
    VSOctaveH5WriterVector<uint16_t>*   fWNumHitScopes; 
    VSOctaveH5WriterVector<uint32_t>*   fWScopeStart;
    VSOctaveH5WriterVector<uint16_t>*   fWScopeCount;

    // WRITERS -- SCOPE
    VSOctaveH5WriterVector<uint16_t>*   fWScopeID;
    VSOctaveH5WriterVector<float>*      fWScopeZenithRad;
    VSOctaveH5WriterVector<float>*      fWScopeAzimuthRad;
    //VSOctaveH5WriterVector<uint32_t>*   fWNumHitPixels;
    VSOctaveH5WriterVector<uint32_t>*   fWPixelStart;
    VSOctaveH5WriterVector<uint32_t>*   fWPixelCount;

    // WRITERS -- PIXEL
    VSOctaveH5WriterVector<uint32_t>*   fWPixelID;
    VSOctaveH5WriterVector<uint32_t>*   fWPEStart;
    VSOctaveH5WriterVector<uint32_t>*   fWPECount;

    // WRITERS -- PE
    VSOctaveH5WriterVector<float>*      fWPixelTimeNS;
 };

}

#endif // VSHDFEVENTSTORE_HPP
