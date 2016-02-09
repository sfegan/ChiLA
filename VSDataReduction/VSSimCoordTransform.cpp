//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimCoordTransform.hpp

    Perform coordinate transforms on simulation events.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSSimCoordTransform.cpp,v 3.1 2008/09/24 20:01:16 matthew Exp $

*/

#include <VSSimCoordTransform.hpp>

using namespace VERITAS;
using namespace SEphem;

std::pair<double,double> VSSimCoordTransform::s_default_src_radec = 
std::pair<double,double>(0,0);
double VSSimCoordTransform::s_default_wobble_phi_deg = 0;

VSSimCoordTransform::
VSSimCoordTransform(const std::pair<double,double>& src_radec,
		    double wobble_phi_deg):
  m_src_radec(M_PI/2.-src_radec.second,src_radec.first),
  m_obs_radec(),
  m_wobble_phi_rad(wobble_phi_deg*M_PI/180.)
{

}

VSSimCoordTransform::~VSSimCoordTransform()
{

}

void VSSimCoordTransform::transform(SphericalCoords& radec,
				    uint32_t corsika_particle_id,
				    const SphericalCoords& azzn, 
				    const SphericalCoords& obs_azzn)
{
  if(corsika_particle_id == 1)
    {
      SphericalCoords obs_radec = obs_azzn;
      obs_radec.rotate(0,-azzn.thetaRad(),-azzn.phiRad());
      double obs_phi   = obs_radec.phiRad();
      radec.rotate(0,-azzn.thetaRad(),-azzn.phiRad());
      radec.rotate(m_src_radec.phiRad(),m_src_radec.thetaRad(),
		   M_PI-obs_phi+m_wobble_phi_rad);
    }
  else if(corsika_particle_id > 1)
    {
      radec.rotate(0,-obs_azzn.thetaRad(),-obs_azzn.phiRad());
      radec.rotate(m_obs_radec.phiRad(),m_obs_radec.thetaRad(),0);
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSSimCoordTransform::configure(VSOptions& options, 
				    const std::string& profile, 
				    const std::string& opt_prefix)
{
  options.findWithValue(OPTNAME(opt_prefix,"sim_src_radec"),
			s_default_src_radec,
			"Set the RA/DEC coordinate in radians of the "
			"centroid for simulated gamma-ray events.",
			OPTNAME(opt_prefix,"common"));

  options.findWithValue(OPTNAME(opt_prefix,"sim_wobble_deg"),
			s_default_wobble_phi_deg,
			"Set the wobble direction in degrees for simulation "
			"events.",
			OPTNAME(opt_prefix,"common"));
}
