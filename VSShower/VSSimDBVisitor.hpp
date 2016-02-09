//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimDBVisitor.hpp
  CORSIKA visitor stores ray-traced PEs in the database

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/19/2005
*/

#ifndef VSSIMDBVISITOR_HPP
#define VSSIMDBVISITOR_HPP

#include <set>

#include <VSRayTracingVisitor.hpp>
#include <VSSimDBCORSIKADatasets.hpp>

#include <VSSimDB.hpp>

namespace VERITAS
{
  class VSSimDBVisitor: public VSRayTracingVisitor
  {
  public:
    VSSimDBVisitor(uint32_t workunit_runid, VSSimDBEventStore* simdb,
		   bool store_empty_telescopes,
		   VSOTelescopeArray* array, RandomNumbers* rng,
		   VSTargeting* targeting,
		   VSArrayTracking* tracking, 
		   VSPrimaryArrivalDistribution* pad,
		   VSSimDBWavelengthData* quaneff_data,
		   VSSimDBWavelengthData* mirrref_data,
		   VSSimDBWavelengthAltitudeData* atmoabs_data);
    virtual ~VSSimDBVisitor();

    // Event setup
    virtual void visitEvent(const VSCORSIKAEvent& event, bool& veto);

    // Event as seen by sampled array
    virtual void visitEventUse(const VSCORSIKAEventUse& use, bool& veto);
    virtual void leaveEventUse(bool veto);

    // Event as seen from each telescope
    virtual void visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope,
				     bool& veto);
    virtual void leaveTelescopeEvent(bool veto);

    // PE generated in pixel
    virtual void visitPE(const VSORayTracer::TraceInfo& info);

  protected:
    uint32_t                       fWorkunitRunID;
    VSSimDBEventStore*             fDB;
    bool                           fStoreEmptyTelescopes;
    VSOTelescopeArray*             fArray;
    VSSimDBEventData               fDBEventData;
    VSSimDBScopeData               fDBScopeData;
    VSSimDBPEData                  fDBPEData;
    VSSimDBWavelengthData*         fQuanEffData;
    VSSimDBWavelengthData*         fMirrRefData;
    VSSimDBWavelengthAltitudeData* fAtmoAbsData;
    bool                           fHitScope;
    std::vector<bool>              fHitPixels;
    std::set<unsigned>             fExtraHitPixels;
  };

}

#endif // VSSIMDBVISITOR_HPP
