//-*-mode:c++; mode:font-lock;-*-

/*! \file RayTracer.cpp

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

#include <cmath>
#include <iostream>

#include <Vec3D.hpp>
#include <Vec4D.hpp>

#include <VSDataConverter.hpp>

#include "VSORayTracer.hpp"

using namespace Physics;
using namespace VERITAS;

VSORayTracer::TraceInfo::TraceInfo():
  array(), status(TS_NONE), ground_x(), ground_y(),  ground_dx(), ground_dy(),
  scope_hexid(), scope(), reflec_x(), reflec_z(), reflec_dx(), reflec_dz(),
  mirror_hexid(), mirror(), mirror_x(), mirror_y(), mirror_z(),
  mirror_normal(), mirror_normal_dispersion(), mirror_scattered(), 
  mirror_reflection_angle(), fplane_x(), fplane_z(), fplane_dx(), fplane_dz(),
  fplane_t(), pixel_hexid(), pixel(), pixel_dist(), concentrator_hit()
{
  // nothing to see here
}

void VSORayTracer::TraceInfo::reset()
{
  *this = VSORayTracer::TraceInfo();
}

std::ostream& VSORayTracer::TraceInfo::write(std::ostream& stream, bool cpu) const
{
  stream 
    << VSDataConverter::toString(status) << ' '                                                    // $1  --  AWK column
    << VSDataConverter::toString(ground_x * ( cpu ? array->spacing() : 1.0 )) << ' '               // $2
    << VSDataConverter::toString(ground_y * ( cpu ? array->spacing() : 1.0 )) << ' '               // $3
    << VSDataConverter::toString(ground_dx * ( cpu ? array->spacing() : 1.0 )) << ' '              // $4
    << VSDataConverter::toString(ground_dy * ( cpu ? array->spacing() : 1.0 )) << ' '              // $5
    << VSDataConverter::toString(scope_hexid) << ' '                                               // $6
    << VSDataConverter::toString(( scope ? scope->hexID() : 0 )) << ' '                            // $7
    << VSDataConverter::toString(reflec_x * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '   // $8
    << VSDataConverter::toString(reflec_z * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '   // $9
    << VSDataConverter::toString(reflec_dx * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '  // $10
    << VSDataConverter::toString(reflec_dz * ( scope&&cpu ? scope->facetSpacing() : 1.0 )) << ' '  // $11
    << VSDataConverter::toString(mirror_hexid) << ' '                                              // $12
    << VSDataConverter::toString(( mirror ? mirror->hexID() : 0 )) << ' '                          // $13
    << VSDataConverter::toString(mirror_x) << ' '                                                  // $14
    << VSDataConverter::toString(mirror_y) << ' '                                                  // $15
    << VSDataConverter::toString(mirror_z) << ' '                                                  // $16
    << VSDataConverter::toString(mirror_normal) << ' '                                             // $17 $18 $19 $20 $21
    << VSDataConverter::toString(mirror_normal_dispersion) << ' '                                  // $22
    << VSDataConverter::toString(mirror_scattered) << ' '                                          // $23 $24 $25 $26 $27
    << VSDataConverter::toString(mirror_reflection_angle) << ' '                                   // $28
    << VSDataConverter::toString(fplane_x * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '   // $29
    << VSDataConverter::toString(fplane_z * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '   // $30
    << VSDataConverter::toString(fplane_dx * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '  // $31
    << VSDataConverter::toString(fplane_dz * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' '  // $32
    << VSDataConverter::toString(fplane_t) << ' '                                                  // $33
    << VSDataConverter::toString(pixel_hexid) << ' '                                               // $34
    << VSDataConverter::toString(( pixel ? pixel->hexID() : 0 )) << ' '                            // $35
    << VSDataConverter::toString(pixel_dist * ( scope&&cpu ? scope->pixelSpacing() : 1.0 )) << ' ' // $36
    << VSDataConverter::toString(concentrator_hit);                                                // $37
  return stream;
}


std::ostream& operator << (std::ostream& stream, const VSORayTracer::TraceInfo& o)
{
  return o.write(stream,false);
}


VSORayTracer::~VSORayTracer()
{
  // nothing to see here
}

const VSOPixel* VSORayTracer::trace(Physics::Particle& ray, TraceInfo& info,
				    const VSOTelescope* scope_hint)
{
  // Initialize array
  info.reset();
  info.array = &fArray;

  // Require photon to be down going
  if(ray.Momenta().r.z>0)
    {
      info.status = TS_DOES_INTERSECT_GROUND;
      return 0;
    }

  info.scope = scope_hint;

  return(scope_trace(ray,info));
}

const VSOPixel* VSORayTracer::trace(Physics::Particle& ray, TraceInfo& info)
{
  // Initialize array
  info.reset();
  info.array = &fArray;

  // Require photon to be down going
  if(ray.Momenta().r.z>0)
    {
      info.status = TS_DOES_INTERSECT_GROUND;
      return 0;
    }

  // Propagate to ground
  Particle ray_copy(ray);
  bool good;
  good = ray_copy.PropagateFreeToPlane(Vec3D(0,0,1), -fArray.altitude(), true);
  if(!good)
    {
      info.status = TS_DOES_INTERSECT_GROUND;
      return 0;
    }

  // Array is assumed to be hexagonal use VVV look-up function to find site
  info.ground_x = ray_copy.Position().r.x / fArray.spacing();
  info.ground_y = ray_copy.Position().r.y / fArray.spacing();
  info.ground_dx = info.ground_x;
  info.ground_dy = info.ground_y;
  if(fArray.arrayParity())info.ground_dx = -info.ground_dx;
  xy_to_nh(&info.ground_dx,&info.ground_dy,&info.scope_hexid);
  if(fArray.arrayParity())info.ground_dx = -info.ground_dx;

  // Find telescope (if there is a real telescope at that site)
  info.scope = fArray.telescopeByHexID(info.scope_hexid);
  if(info.scope==0)
    {
      info.status = TS_NO_SCOPE;
      return 0;
    }

#if 0 
  // Propagate to reflector impact sphere -- test whether ray can hit scope
  // This probably does not speed things up much as the propagation to
  // the reflector (two steps down) does a similar thing
  good = ray_copy.PropagateFreeToSphere(info.scope->position(),
					info.scope->reflectorIP(),
					IP_EARLIEST,true);
  if(!good)
    {
      info.status = TS_OUTSIDE_REFLECTOR_IP;
      return 0;
    }
#endif

  return(scope_trace(ray,info));
}

const VSOPixel* 
VSORayTracer::scope_trace(Physics::Particle& ray, TraceInfo& info)
{
  bool good;

  // Transform to reflector co-ordinate system
  info.scope->globalToReflector(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "A: " << ray.Momenta().r/ray.Momenta().r0 << std::endl;
#endif
  // **************************************************************************
  // ****************** RAY IS NOW IN REFLECTOR COORDINATES *******************
  // **************************************************************************

  // Propagate to intersection with the reflector sphere
  good = ray.PropagateFreeToSphere(Vec3D(0,info.scope->curvatureRadius(),0),
				   info.scope->curvatureRadius(), 
				   Particle::IP_LATEST, false /* true */);
  if(!good)
    {
      info.status = TS_MISSED_REFLECTOR_SPHERE;
      return 0;
    }

  // Assume mirrors on hexagonal grid - use VVV routines to find which mirror
  double tx = ray.Position().r.x / info.scope->facetSpacing();
  double tz = ray.Position().r.z / info.scope->facetSpacing();
  double costheta = cos(info.scope->reflectorRotation());
  double sintheta = sin(info.scope->reflectorRotation());
  info.reflec_x = tx*costheta - tz*sintheta;  // Align with hex grid in dir of
  info.reflec_z = tz*costheta + tx*sintheta;  // Vec3D(0,-reflectorRotation,0)
  info.reflec_dx = info.reflec_x;
  info.reflec_dz = info.reflec_z;
  if(info.scope->mirrorParity())info.reflec_dx = -info.reflec_dx; // Reverse parity if required
  xy_to_nh(&info.reflec_dx,&info.reflec_dz,&info.mirror_hexid);
  if(info.scope->mirrorParity())info.reflec_dx = -info.reflec_dx; // Reverse parity if required
  
  // Find mirror (if there is a real mirror at that site)
  info.mirror = info.scope->mirrorByHexID(info.mirror_hexid);
  if(info.mirror==0)
    {
      info.status = TS_NO_MIRROR;
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  if(info.mirror->removed())
    {
      info.status = TS_MIRROR_REMOVED;
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  // Propagate to intersection with the mirror sphere
  double mirror_radius = info.mirror->focalLength()*2.0;
  Vec3D mirror_center = 
    info.mirror->pos() + info.mirror->align()*mirror_radius;

  good = ray.PropagateFreeToSphere(mirror_center, mirror_radius,
				   Particle::IP_LATEST, true);
  if(!good)
    {
      info.status = TS_MISSED_MIRROR_SPHERE;
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  // Convert to mirror coordinates
  info.mirror->reflectorToMirror(ray);
#ifdef DEBUG_DIRECTION
  std::cerr << "B: " << ray.Momenta().r/ray.Momenta().r0 << std::endl;
#endif

  // **************************************************************************
  // ******************** RAY IS NOW IN MIRROR COORDINATES ********************
  // **************************************************************************

  info.mirror_x = ray.Position().r.x;
  info.mirror_y = ray.Position().r.y;
  info.mirror_z = ray.Position().r.z;
  
  // Check if ray impacted beyond the edge of this mirror - in which case we
  // should check if it hits neighbour -- this is likely to be unimportant so
  // just throw the ray away
  static const double cos60 = 1.0/2.0;
  static const double sin60 = sqrt(3.0)/2.0;
  double edge = info.scope->facetSize()/2.0;
  double x_0 = info.mirror_x;
  double x_pos60 = cos60*info.mirror_x - sin60*info.mirror_z;
  double x_neg60 = cos60*info.mirror_x + sin60*info.mirror_z;

  if((x_0>edge)||(x_0<-edge)||(x_pos60>edge)||(x_pos60<-edge)||
     (x_neg60>edge)||(x_neg60<-edge))
    {
      info.status = TS_MISSED_MIRROR_EDGE;
      info.mirror->mirrorToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  // Check to see if photon is absorbed at the mirror. Would be faster to 
  // do this check before the edge check, but this way gets info.status 
  // correct .. is this important ?
  if(fRNG.Uniform() > info.mirror->degradingFactor())
    {
      info.status = TS_ABSORBED_AT_MIRROR;
      info.mirror->mirrorToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  // Find normal at the point of intersection
  info.mirror_normal = Vec3D(0,mirror_radius,0) - ray.Position().r;
  info.mirror_normal /= info.mirror_normal.Norm();

  // Scatter the normal to account for the spot size ot the focal length of the
  // radius. The spot size is given as the DIAMETER at the focal distance.
  // Must divide by 2.0 (for reflection)
  info.mirror_normal_dispersion = 
    info.mirror->spotSize()/2.0/info.mirror->focalLength();

  info.mirror_scattered = info.mirror_normal;
  info.mirror_scattered.ScatterDirection(info.mirror_normal_dispersion,fRNG);
 
  // Reflect ray
  ray.Reflect(info.mirror_scattered);
  
  // Back to reflector coordinates
  info.mirror->mirrorToReflector(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "C: " << ray.Momenta().r/ray.Momenta().r0 << std::endl;
#endif

  // **************************************************************************
  // ****************** RAY IS NOW IN REFLECTOR COORDINATES *******************
  // **************************************************************************

  // Translate to focal plane coordinates
  info.scope->reflectorToFocalPlane(ray);

#ifdef DEBUG_DIRECTION
  std::cerr << "D: " << ray.Momenta().r/ray.Momenta().r0 << std::endl;
#endif

  // **************************************************************************
  // ***************** RAY IS NOW IN FOCAL PLANE COORDINATES ******************
  // **************************************************************************

  // Propagate back to camera plane
  good = ray.PropagateFreeToPlane(Vec3D(0,1,0),0,false);
  if(!good)
    {
      info.status = TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE;
      info.mirror->mirrorToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  info.fplane_x = ray.Position().r.x / info.scope->pixelSpacing();
  info.fplane_z = ray.Position().r.z / info.scope->pixelSpacing();
  info.fplane_dx = info.fplane_x;
  info.fplane_dz = info.fplane_z;
  info.fplane_t = ray.Position().r0 / Constants::cgs_c;
  if(info.scope->pixelParity())info.fplane_dx = -info.fplane_dx;
  xy_to_nh(&info.fplane_dx,&info.fplane_dz,&info.pixel_hexid);
  if(info.scope->pixelParity())info.fplane_dx = -info.fplane_dx;

  // Find pixel (if there is a real pixel at that site)
  info.pixel = info.scope->pixelByHexID(info.pixel_hexid);
  if(info.pixel==0)
    {
      info.status = TS_NO_PIXEL;
      info.scope->focalPlaneToReflector(ray);
      info.scope->reflectorToGlobal(ray);
      return 0;
    }

  info.pixel_dist = 
    sqrt(info.fplane_dx*info.fplane_dx + info.fplane_dz*info.fplane_dz) * 
    info.scope->pixelSpacing();

  info.concentrator_hit = 
    (info.pixel_dist > info.scope->cathodeDiameter()/2.0);
  if(info.concentrator_hit)
    {
      if(fRNG.Uniform() > info.scope->concentratorSurvivalProb())
	{
	  info.status = TS_ABSORBED_AT_CONCENTRATOR;
	  return 0;
	}
    }

  // Translate to reflector coordinates
  info.scope->focalPlaneToReflector(ray);

  // **************************************************************************
  // ****************** RAY IS NOW IN REFLECTOR COORDINATES *******************
  // **************************************************************************

  info.pixel_dist = Vec3D(ray.Position().r - info.pixel->pos() -
			  info.scope->focalPlanePosition()).Norm();

  // Transform back to global
  info.scope->reflectorToGlobal(ray);

  // **************************************************************************
  // ******************** RAY IS NOW IN GLOBAL COORDINATES ********************
  // **************************************************************************
    
  info.status = TS_PE_GENERATED;
  return info.pixel;
}

bool VSORayTracer::beam(Physics::Particle& photon,
		     const Physics::Vec3D& origin,
		     const Physics::Vec3D& direction, 
		     double beam_start, double beam_stop, 
		     double beam_radius_in, double beam_radius_out,
		     double beam_angle_lo, double beam_angle_hi,
		     double lambda_nm)
{
  double d_norm = direction.Norm();
  if(d_norm == 0)return false;
  Physics::Vec3D d_hat = direction/d_norm;

  Vec3D tangent_a;
  if((d_hat.x<=d_hat.y)&&(d_hat.x<=d_hat.z))tangent_a = d_hat^Vec3D(1,0,0);
  else if(d_hat.y<=d_hat.z)tangent_a = d_hat^Vec3D(0,1,0);
  else tangent_a = d_hat^Vec3D(0,0,1);
  tangent_a /= tangent_a.Norm();
  
  Vec3D tangent_b(direction^tangent_a);
  tangent_b /= tangent_b.Norm();

  // CHOOSE PHOTON EMISSION POINT
  Physics::Vec3D emission_point(origin);
  emission_point += d_hat*beam_start;

  // SAMPLE PHOTON EMISSION POINT FROM BEAM LENGTH
  if(beam_stop != beam_start)
    {
      emission_point += d_hat*((beam_stop-beam_start)*fRNG.Uniform());
    }

  // SAMPLE PHOTON EMISSION POINT FROM BEAM AREA
  if((beam_radius_in != 0)||(beam_radius_out != 0))
    {
      double theta = 2.0*M_PI*fRNG.Uniform();
      double rho2 = beam_radius_in*beam_radius_in;
      if(beam_radius_in != beam_radius_out)
	rho2 += (beam_radius_out*beam_radius_out-rho2) * fRNG.Uniform();
      double rho = sqrt(rho2);
      emission_point += tangent_a*rho*cos(theta) + tangent_b*rho*sin(theta);
    }

  // SAMPLE PHOTON EMISSION DIRECTION
  double phi = fRNG.Uniform() * 2.0*Constants::num_pi;
  double costheta = cos(beam_angle_lo);
  if(beam_angle_lo != beam_angle_hi)
    {
      costheta += fRNG.Uniform() * (cos(beam_angle_hi)-costheta);
    }
  double theta = acos(costheta);

  if(theta != 0)
    {
#if 1
      d_hat = d_hat*costheta + (tangent_a*cos(phi)+tangent_b*sin(phi))*sin(theta);
#else  
      Vec3D axis = tangent_a*cos(phi) + tangent_b*sin(phi);
      d_hat.Rotate(axis*theta);
#endif
    }

  // MAKE PHOTON
  double E = Constants::cgs_h * Constants::cgs_c / ( lambda_nm * 1e-7 );
  photon = Particle(Physics::Vec4D(0,emission_point),d_hat,E,0,0);
  return true;
}

bool VSORayTracer::laserBeam(Physics::Particle& photon,
			  const Physics::Vec3D& origin,
			  const Physics::Vec3D& direction,
			  double d0, double sampling_radius,
			  double lambda_nm)
{
  return beam(photon, origin, direction, d0, d0, 0, sampling_radius, 0, 0,
	      lambda_nm);
}

bool VSORayTracer::fanBeam(Physics::Particle& photon,
			const Physics::Vec3D& origin, 
			const Physics::Vec3D& direction, 
			double half_angle_spread, 
			double lambda_nm)
{
  return beam(photon, origin, direction, 0, 0, 0, 0, 0, half_angle_spread,
	      lambda_nm);
}

bool VSORayTracer::muonBeam(Physics::Particle& photon,
			 const Physics::Vec3D& origin, 
			 const Physics::Vec3D& direction, 
			 double muon_travel_distance, double opening_angle, 
			 double lambda_nm)
{
  return beam(photon, origin, direction, 0, muon_travel_distance, 0, 0, 
	      opening_angle, opening_angle, lambda_nm);
}

