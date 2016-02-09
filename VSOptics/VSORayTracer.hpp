//-*-mode:c++; mode:font-lock;-*-

/*! \file VSORayTracer.hpp

  Ray tracing class header file

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \author  Maciej Nicewicz             \n
           UCLA                        \n
	   nicewicz@physics.ucla.edu   \n

  \date    12/05/2004
  \version 0.2
  \note
*/

#ifndef VSORAYTRACER_HPP
#define VSORAYTRACER_HPP

#include <Vec3D.hpp>
#include <Vec4D.hpp>
#include <Particle.hpp>
#include <RandomNumbers.hpp>

#include "VSOTelescopeArray.hpp"

namespace VERITAS
{

  class VSORayTracer
  {
  public:
    VSORayTracer(const VSOTelescopeArray& array, RandomNumbers& rng): 
      fArray(array), fRNG(rng) { }
    virtual ~VSORayTracer();
    
    enum Status { TS_NONE,                               // 0
		  TS_DOES_INTERSECT_GROUND,              // 1
		  TS_NO_SCOPE,                           // 2
		  TS_OUTSIDE_REFLECTOR_IP,               // 3
		  TS_MISSED_REFLECTOR_SPHERE,            // 4
		  TS_NO_MIRROR,                          // 5
		  TS_MIRROR_REMOVED,                     // 6
		  TS_MISSED_MIRROR_SPHERE,               // 7
		  TS_MISSED_MIRROR_EDGE,                 // 8
		  TS_ABSORBED_AT_MIRROR,                 // 9
		  TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE,   // 10
		  TS_NO_PIXEL,                           // 11
		  TS_ABSORBED_AT_CONCENTRATOR,           // 12
		  TS_PE_GENERATED };                     // 13
    
    class TraceInfo
    {
    public:					
      const VSOTelescopeArray* array;
      Status              status;
      double              ground_x;
      double              ground_y;
      double              ground_dx;
      double              ground_dy;
      int                 scope_hexid;
      const VSOTelescope* scope;
      double              reflec_x;
      double              reflec_z;
      double              reflec_dx;
      double              reflec_dz;
      int                 mirror_hexid;
      const VSOMirror*    mirror;
      double              mirror_x;
      double              mirror_y;
      double              mirror_z;
      Physics::Vec3D      mirror_normal;
      double              mirror_normal_dispersion;
      Physics::Vec3D      mirror_scattered;
      double              mirror_reflection_angle;
      double              fplane_x;
      double              fplane_z;
      double              fplane_dx;
      double              fplane_dz;
      double              fplane_t;
      int                 pixel_hexid;
      const VSOPixel*     pixel;
      double              pixel_dist;
      bool                concentrator_hit;
      
      TraceInfo();
      void reset();
      std::ostream& write(std::ostream& stream, bool convert_to_physical_units=false) const;
    };
    
    const VSOPixel* trace(Physics::Particle& ray, TraceInfo& info);
    const VSOPixel* trace(Physics::Particle& ray, TraceInfo& info,
			  const VSOTelescope* scope_hint);
    
    bool beam(Physics::Particle& photon,
	      const Physics::Vec3D& origin, const Physics::Vec3D& direction, 
	      double beam_start, double beam_stop, 
	      double beam_radius_in, double beam_radius_out,
	      double beam_angle_lo, double beam_angle_hi,
	      double lambda_nm = 400);
    
    bool laserBeam(Physics::Particle& photon, 
		   const Physics::Vec3D& center,
		   const Physics::Vec3D& direction, 
		   double d0, double sampling_radius,
		   double lambda_nm = 400);
    
    bool fanBeam(Physics::Particle& photon,
		 const Physics::Vec3D& origin, 
		 const Physics::Vec3D& direction, 
		 double half_angle_spread, double lambda_nm = 400);
    
    bool muonBeam(Physics::Particle& photon,
		  const Physics::Vec3D& origin, 
		  const Physics::Vec3D& direction, 
		  double muon_travel_distance, double opening_angle, 
		  double lambda_nm = 400);
    
  private:
    const VSOTelescopeArray& fArray;
    RandomNumbers&           fRNG;

    const VSOPixel* scope_trace(Physics::Particle& ray, TraceInfo& info);
  };
  
  std::ostream& operator <<(std::ostream& stream, 
			    const VSORayTracer::TraceInfo& o);
  
}

#endif // VSORAYTRACER_HPP
