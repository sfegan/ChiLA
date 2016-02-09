//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplier.hpp

  Base class for input of Event from simulations database 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       19/07/2006
  \note
*/

#ifndef VSEVENTSUPPLIER_HPP
#define VSEVENTSUPPLIER_HPP

#include<vector>
#include<stdint.h>

#include <VSSimDB.hpp>

//! VERITAS namespace
namespace VERITAS 
{

  class VSEventSupplier
  {
  public:
    struct PE
    {
      PE(): fTimeNS() { /* nothing to see here */ }

      float                fTimeNS;
    };
  
    struct Pixel
    {
      Pixel(): fPixelID(), fPEs() { /* nothing to see here */ }

      uint32_t             fPixelID;
      std::vector<PE>      fPEs;
    };
    
    struct Scope
    {
      Scope(): 
	fScopeID(), fScopeZenithRad(), fScopeAzimuthRad(), fPixels()
      { /* nothing to see here */ }

      uint32_t             fScopeID;
      float                fScopeZenithRad;
      float                fScopeAzimuthRad;
      std::vector<Pixel>   fPixels;
    };

    struct Event
    {
      Event():
	fEventID(), fTargetZenithRad(), fTargetAzimuthRad(),
	fPrimarySpeciesCORSIKA(), fPrimaryEnergyTeV(), 
	fPrimaryZenithRad(), fPrimaryAzimuthRad(), 
	fPrimaryCoreEastM(), fPrimaryCoreNorthM(), fPrimaryCoreUpASLM(),
	fSamplingRadiusM(), fEventComplete(), fTableIndex(), fScopes()
      { /* nothing to see here */ }

      uint32_t             fEventID;
      float                fTargetZenithRad;
      float                fTargetAzimuthRad;
      uint32_t             fPrimarySpeciesCORSIKA;
      float                fPrimaryEnergyTeV;
      float                fPrimaryZenithRad;
      float                fPrimaryAzimuthRad;
      float                fPrimaryCoreEastM;
      float                fPrimaryCoreNorthM;
      float                fPrimaryCoreUpASLM;
      float                fSamplingRadiusM;
      bool                 fEventComplete;
      uint32_t             fTableIndex;
      std::vector<Scope>   fScopes;
    };

    struct SimParam
    {
      SimParam(): 
        fTableID(), fPrimaryID(), fEnergyTeV(), 
	fZenithMinRad(), fZenithMaxRad(), fAzimuthMinRad(), fAzimuthMaxRad(),
        fOpticsID(), fSamplingRadiusM(), fEventCount(), fTableName()
      { /* nothing to see here */ }

      SimParam(const VERITAS::VSSimDBTableParam& p, unsigned event_count = 0): 
        fTableID(p.fTableID),
	fPrimaryID(p.fPrimaryID),
        fEnergyTeV(p.fEnergyGeV*0.001),
        fZenithMinRad(p.fZenithMinRad),
        fZenithMaxRad(p.fZenithMaxRad),
        fAzimuthMinRad(p.fAzimuthMinRad),
        fAzimuthMaxRad(p.fAzimuthMaxRad),
        fOpticsID(p.fOpticsID),
        fSamplingRadiusM(p.fSamplingRadiusM),
        fEventCount(event_count),
        fTableName(p.fTableName) 
      { /* nothing to see here */ }
      
      uint32_t             fTableID;
      uint32_t             fPrimaryID;
      float                fEnergyTeV;
      float                fZenithMinRad;
      float                fZenithMaxRad;
      float                fAzimuthMinRad;
      float                fAzimuthMaxRad;
      uint32_t             fOpticsID;
      float                fSamplingRadiusM;
      uint32_t             fEventCount;
      std::string          fTableName;
    };
    
    virtual ~VSEventSupplier();
    virtual std::vector<SimParam> getSimParam() = 0;
    virtual bool getNextEvent(Event& e) = 0;
    virtual uint32_t getNumEvents() = 0;
  };

} // namespace VERITAS

#endif // VSEVENTSUPPLIER_HPP
