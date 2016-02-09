//-*-mode:c++; mode:font-lock;-*-

/*! \file VSRayTracingVisitor.hpp
  CORSIKA visitor which traces photons through the optics

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/19/2005
*/

#include <cassert>

#include "VSRayTracingVisitor.hpp"

using namespace VERITAS;
using namespace Physics;

static inline float d2r(const float& x)
{
  return x/180.0*M_PI;
}

VSRayTracingVisitor::
VSRayTracingVisitor(VSOTelescopeArray* array, RandomNumbers* rng,
		    VSTargeting* targeting, VSArrayTracking* tracking, 
		    VSPrimaryArrivalDistribution* pad,
		    VSSimDBWavelengthData* quaneff_data,
		    VSSimDBWavelengthData* mirrref_data,
		    VSSimDBWavelengthAltitudeData* atmoabs_data)
  : fArray(array), fRNG(rng), fRayTracer(), fTargeting(targeting),
    fArrayTracking(tracking), fPrimaryArrivalDistribution(pad),
    fQuanEffData(quaneff_data), fMirrRefData(mirrref_data),
    fAtmoAbsData(atmoabs_data),
    fObsLevel(), fRunParam(), fExtraRunParam(), fEvent(), fCore(),
    fTargetZenithRad(), fTargetAzimuthRad(), fNPE(), fScope(),
    fFakeCamera(false), fFakeCameraRadius()
{
  fRayTracer = new VSORayTracer(*array,*rng);
}

VSRayTracingVisitor::~VSRayTracingVisitor()
{
  delete fRayTracer;
}

void VSRayTracingVisitor::
visitRun(const VSCORSIKARunParameters& param)
{
  fRunParam = param;
  fObsLevel = param.obs_height*100.0;
}

void VSRayTracingVisitor::
visitRunExtra(const VSCORSIKAExtraRunParameters& param)
{
  fExtraRunParam = param;
}

void VSRayTracingVisitor::
visitArraySpec(const VSCORSIKAArraySpec& arrayspec)
{
  fScopeOpticalDepth.resize(fArray->numTelescopes());

  assert(arrayspec.array_num_tel == fArray->numTelescopes());
}

void VSRayTracingVisitor::
visitTelescopeSpec(const VSCORSIKATelescopeSpec& scopespec)
{
  // --------------------------------------------------------------------------
  // Calculate optical depth to each telescope
  // --------------------------------------------------------------------------
  for(VSSimDBWavelengthAltitudeData::iterator itr = 
	fAtmoAbsData->begin(); itr != fAtmoAbsData->end(); ++itr)
    {
      double tau = 
	interpolate(itr->second,
		    fArray->telescope(scopespec.scope_num)->position().z/1E5);
      
      fScopeOpticalDepth[scopespec.scope_num][itr->first] = tau;
    }

  assert(fabs(-scopespec.scope_y*100.0 - 
	      fArray->telescope(scopespec.scope_num)->position().x)<1);
  assert(fabs(scopespec.scope_x*100.0 - 
	      fArray->telescope(scopespec.scope_num)->position().y)<1);
  assert(fabs((fObsLevel+scopespec.scope_z*100.0) -
	      fArray->telescope(scopespec.scope_num)->position().z)<1);
  assert(fabs(scopespec.scope_r*200.0 - 
	      fArray->telescope(scopespec.scope_num)->reflectorIP())<1);
}

void VSRayTracingVisitor::
visitEvent(const VSCORSIKAEvent& event, bool& veto)
{
  if(fExtraRunParam.sample_cone_hi < DBL_EPSILON)
    {
      fTargetZenithRad  = d2r(90.0-event.pri_elevation);
      fTargetAzimuthRad = d2r(event.pri_azimuth);
    }
  else
    {
      fTargetZenithRad  = d2r(fExtraRunParam.sample_theta_lo);
      fTargetAzimuthRad = 
	fmod(5.0*M_PI-d2r(fExtraRunParam.sample_phi_lo),2.0*M_PI);
    }
  if(fPrimaryArrivalDistribution!=NULL)
    fPrimaryArrivalDistribution->setTarget(fTargetZenithRad,fTargetAzimuthRad);
  fTargeting->setTarget(fTargetZenithRad,fTargetAzimuthRad);
  fEvent = event;
}

void VSRayTracingVisitor::
visitEventUse(const VSCORSIKAEventUse& use, bool& veto)
{
  fArrayTracking->pointArrayForEvent(fArray, 
				     fTargetZenithRad, fTargetAzimuthRad);
  fCore = Vec3D(-use.pri_ycore * 100.0, use.pri_xcore * 100.0, fObsLevel);
}

void VSRayTracingVisitor::
visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope, bool& veto)
{
  fScope = fArray->telescope(scope.scope_num);
}

void VSRayTracingVisitor::
visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch)
{
  fNPE = 0;
  float bunch_size = bunch.ph_count;
  float lambda_nm = bunch.ph_lambda;

  // --------------------------------------------------------------------------
  // The ph_lambda parameter can have 4 possible meanings:
  // 
  // ph_lambda > 0  : Photon Bunch with wavelength equal to ph_lambda.
  // ph_lambda = 0  : Photon Bunch for which wavelength has not been assigned.
  // ph_lambda = -1 : PE for which QE, Atmospheric Absorption, etc. 
  //                  have already been applied.
  // ph_lambda < -1 : PE for which QE, Atmospheric Absorption, etc. 
  //                  have already been applied.  The absolute value of 
  //                  ph_lambda is equal to the wavelength of the original 
  //                  photon that produced this PE.
  // --------------------------------------------------------------------------

  // --------------------------------------------------------------------------
  // If lambda = 0 then no wavelength has been assigned to this photon
  // so we set it here according to 1/lam^2 distribution.
  // --------------------------------------------------------------------------
  if(fabs(lambda_nm) < DBL_EPSILON)
    {
      lambda_nm = 1./(1./fExtraRunParam.cerenk_lambda_lo-fRNG->Uniform()*
		      (1./fExtraRunParam.cerenk_lambda_lo-
		       1./fExtraRunParam.cerenk_lambda_hi));
    }

  // Convert to Photo Electrons -----------------------------------------------
  if(lambda_nm > 0)
    {
      VSSimDBWavelengthAltitudeData::iterator itr = 
	fAtmoAbsData->lower_bound((unsigned)lambda_nm);

      if(itr == fAtmoAbsData->begin()) ++itr;
     
      float tau2 = interpolate(itr->second,bunch.ph_height/1E3);
      float lam2 = itr->first;
      --itr;
      float tau1 = interpolate(itr->second,bunch.ph_height/1E3);
      float lam1 = itr->first;

      float tau = tau1 + (lambda_nm-lam1)*(tau2-tau1)/(lam2-lam1);
      float tau_scope = interpolate(fScopeOpticalDepth[fScope->id()],
				    lambda_nm);

      float dtau = fabs(tau-tau_scope);
      float atmoabs = exp(-dtau*fabs(1/bunch.ph_cosine_z));
      float qeff = interpolate(*fQuanEffData,lambda_nm);
      float mirrref = interpolate(*fMirrRefData,lambda_nm);

      bunch_size *= qeff*mirrref*atmoabs;
    }

  for(float npe = bunch_size; npe > 0; npe--)
    {
      if((npe < 0.9999)&&(fRNG->Uniform()>npe)) break;
      fNPE++;
    }

  if(fNPE == 0) return;

  // Compute 4-vector of impact event on the ground
  double ix = -bunch.ph_rel_y*100.0 + fScope->position().x;
  double iy = bunch.ph_rel_x*100.0 + fScope->position().y;
  double iz = fScope->position().z;
  Vec3D impact_point(ix,iy,iz);
  Vec4D impact_event(Constants::cgs_c*bunch.ph_time*1e-9, impact_point);

  // Compute 4-vector velocity of idealized (vacuum) photon
  Vec3D d_hat(-bunch.ph_cosine_y, bunch.ph_cosine_x, bunch.ph_cosine_z);
  Vec4D V(1,d_hat); 
  V *= Constants::cgs_c;
  
  // Compute 4-vector of emission event for idealized (vacuum) photon
  Vec4D emission_event(impact_event);
  emission_event -= V*((bunch.ph_height-iz)/V.r.z);

  // Test that we are dealing with PHOTO ELECTRONS from CORSIKA
  lambda_nm = fabs(lambda_nm);
  if(fabs(lambda_nm - 1.0)<DBL_EPSILON)lambda_nm=300;

  // Compute 4-vector momentum of photon
  Vec4D P(V); // 4-momentum vector of photon
  P *= Constants::cgs_h / ( lambda_nm * 1e-7 );

  // Photon particle
  Particle photon(emission_event,P,0);

  // Ray-trace photon and record it
  VSORayTracer::TraceInfo info;
  bool visit = fRayTracer->trace(photon, info, fScope);

  // If the fake camera option is on then we fill the camera with
  // fake pixels allowing for huge cameras with no database Pixel entries
  if(visit)
    visitPE(info);	  
  else if((fFakeCamera)&&(info.status==VSORayTracer::TS_NO_PIXEL))
    {
      double r2 = info.fplane_x*info.fplane_x+info.fplane_z*info.fplane_z;
      double rmax = fFakeCameraRadius*info.scope->focalPlanePosition().y
	/info.scope->pixelSpacing();
      if(r2<rmax*rmax)
	{
	  VSOPixel fake_pixel(info.scope, info.pixel_hexid, 
			      info.pixel_hexid,
			      false, 
			      Vec3D(info.fplane_x*
				    info.scope->pixelSpacing(),
				    0,
				    info.fplane_z*
				    info.scope->pixelSpacing()));
	  info.pixel = &fake_pixel;
	  visitPE(info);	  
	}
    }
}

void VSRayTracingVisitor::
visitPE(const VSORayTracer::TraceInfo& info)
{
  // nothing to see here
}

float VSRayTracingVisitor::interpolate(const std::map<unsigned, float>& data,
				       float x)
{
  std::map<unsigned, float>::const_iterator itr = 
    data.lower_bound((unsigned)x);

  float x2 = itr->first;
  float y2 = itr->second;
  --itr;
  float x1 = itr->first;
  float y1 = itr->second;
  return y1 + (x-x1)*(y2-y1)/(x2-x1);
}
