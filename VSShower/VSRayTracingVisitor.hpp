//-*-mode:c++; mode:font-lock;-*-

/*! \file VSRayTracingVisitor.hpp
  CORSIKA visitor which traces photons through the optics

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/19/2005
*/

#ifndef VSRAYTRACINGVISITOR_HPP
#define VSRAYTRACINGVISITOR_HPP

#include <VSCORSIKAEvent.hpp>
#include <VSSimDBCORSIKADatasets.hpp>
#include <VSOTelescopeArray.hpp>
#include <VSOTelescope.hpp>
#include <VSORayTracer.hpp>
#include <VSTargeting.hpp>
#include <RandomNumbers.hpp>

namespace VERITAS
{
  class VSRayTracingVisitor: public VSCORSIKAEventVisitor
  {
  public:
    VSRayTracingVisitor(VSOTelescopeArray* array, RandomNumbers* rng,
			VSTargeting* targeting,
			VSArrayTracking* tracking, 
			VSPrimaryArrivalDistribution* pad,
			VSSimDBWavelengthData* quaneff_data,
			VSSimDBWavelengthData* mirrref_data,
			VSSimDBWavelengthAltitudeData* atmoabs_data);
    virtual ~VSRayTracingVisitor();

    // CORSIKA run info
    virtual void visitRun(const VSCORSIKARunParameters& param);
    virtual void visitRunExtra(const VSCORSIKAExtraRunParameters& param);

    // CORSIKA configuarion -- Telescope Positions
    virtual void visitArraySpec(const VSCORSIKAArraySpec& arrayspec);
    virtual void visitTelescopeSpec(const VSCORSIKATelescopeSpec& scopespec);

    // Event setup
    virtual void visitEvent(const VSCORSIKAEvent& event, bool& veto);

    // Event as seen by sampled array
    virtual void visitEventUse(const VSCORSIKAEventUse& use, bool& veto);

    // Event as seen from each telescope
    virtual void visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope,
				     bool& veto);

    // Photon generated in telescope
    virtual void visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch);

    // EXTENSION OF INTERFACE TO RECORD OUTPUT DATA

    virtual void visitPE(const VSORayTracer::TraceInfo& info);

    // ENABLE "FAKE CAMERA" MODE

    void setFakeCamera(float radius) 
    { fFakeCamera=true; fFakeCameraRadius=radius; }

  protected: 
    VSOTelescopeArray*            fArray;
    RandomNumbers*                fRNG;
    VSORayTracer*                 fRayTracer;
    VSTargeting*                  fTargeting;
    VSArrayTracking*              fArrayTracking;
    VSPrimaryArrivalDistribution* fPrimaryArrivalDistribution;
    VSSimDBWavelengthData*        fQuanEffData;
    VSSimDBWavelengthData*        fMirrRefData;
    VSSimDBWavelengthAltitudeData* fAtmoAbsData;
    float                         fObsLevel;
    VSCORSIKARunParameters        fRunParam;
    VSCORSIKAExtraRunParameters   fExtraRunParam;
    VSCORSIKAEvent                fEvent;
    Physics::Vec3D                fCore;
    float                         fTargetZenithRad;
    float                         fTargetAzimuthRad;
    unsigned                      fNPE;
    const VSOTelescope*           fScope;

    bool                          fFakeCamera;
    float                         fFakeCameraRadius;
    std::vector< VSSimDBWavelengthData > fScopeOpticalDepth;
  private:
    float interpolate(const std::map<unsigned, float>& data,float x);
  };

}

#endif // VSRAYTRACINGVISITOR_HPP
