//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierOneTableHDF.hpp

  Input of event from single table from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       07/08/2006
  \note
*/

#ifndef VSEVENTSUPPLIERONETABLEHDF_HPP
#define VSEVENTSUPPLIERONETABLEHDF_HPP

#include<vector>
#include<list>

#include<VSSimDB.hpp>
#include<VSEventSupplier.hpp>
#include<VSEventSupplierOneTableCommon.hpp>
#include<VSOctaveIO.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplierOneTableHDF: public VSEventSupplierOneTableCommon
  {
  public:    
    VSEventSupplierOneTableHDF(const std::string& table_name, VSSimDB* sim_db,
			       std::string directory = "",
			       bool check_hdf_event_count = true);
    VSEventSupplierOneTableHDF(const std::string& table_name, VSSimDB* sim_db,
			       uint32_t nevent_to_supply, 
			       double ievent_to_start_fraction,
			       std::string directory = "",
			       bool check_hdf_event_count = true);
    virtual ~VSEventSupplierOneTableHDF();
    virtual bool getNextEvent(Event& e);

  private:
    bool open(unsigned workunit_run_id);
    void close();
    void initialize(bool check_hdf_event_count);

    std::string                         fDirectory;
    std::list<unsigned>                 fWorkunitIds;

    // READER
    VSOctaveH5ReaderStruct*             fReader;
    std::string                         fWorkunitFile;
    unsigned                            fIDOfOpenWorkunit;
    unsigned                            fNEventInWorkunit;
    unsigned                            fICompleteEventAtWorkunitStart;
    
    // READERS -- INFRASTRUCTURE
    VSOctaveH5ReaderStruct*             fEV;
    VSOctaveH5ReaderStruct*             fSC;
    VSOctaveH5ReaderStruct*             fPX;
    VSOctaveH5ReaderStruct*             fPE;

    // READERS -- EVENT
    VSOctaveH5ReaderVector<uint32_t>*   fREventID;
    VSOctaveH5ReaderVector<float>*      fRTargetZenithRad;
    VSOctaveH5ReaderVector<float>*      fRTargetAzimuthRad;
    VSOctaveH5ReaderVector<float>*      fRPrimaryZenithRad;
    VSOctaveH5ReaderVector<float>*      fRPrimaryAzimuthRad;
    VSOctaveH5ReaderVector<float>*      fRPrimaryCoreEastM;
    VSOctaveH5ReaderVector<float>*      fRPrimaryCoreNorthM;
    VSOctaveH5ReaderVector<float>*      fRPrimaryCoreUpASLM;
    VSOctaveH5ReaderVector<uint16_t>*   fRNumHitScopes; 
    VSOctaveH5ReaderVector<uint32_t>*   fRScopeStart;
    VSOctaveH5ReaderVector<uint16_t>*   fRScopeCount;

    // READERS -- SCOPE
    VSOctaveH5ReaderVector<uint16_t>*   fRScopeID;
    VSOctaveH5ReaderVector<float>*      fRScopeZenithRad;
    VSOctaveH5ReaderVector<float>*      fRScopeAzimuthRad;
    VSOctaveH5ReaderVector<uint32_t>*   fRPixelStart;
    VSOctaveH5ReaderVector<uint32_t>*   fRPixelCount;

    // READERS -- PIXEL
    VSOctaveH5ReaderVector<uint32_t>*   fRPixelID;
    VSOctaveH5ReaderVector<uint32_t>*   fRPEStart;
    VSOctaveH5ReaderVector<uint32_t>*   fRPECount;

    // READERS -- PE
    VSOctaveH5ReaderVector<float>*      fRPixelTimeNS;    
  };

}

#endif // defined VSEVENTSUPPLIERONETABLEHDF_HPP
