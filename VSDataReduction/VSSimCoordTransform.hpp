//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimCoordTransform.hpp

  Perform coordinate transforms on simulation events.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSSimCoordTransform.hpp,v 3.1 2008/09/24 20:01:16 matthew Exp $

*/

#ifndef VSSIMCOORDTRANSFORM_HPP
#define VSSIMCOORDTRANSFORM_HPP

#include <SphericalCoords.h>
#include <VSOptions.hpp>

namespace VERITAS
{
  class VSSimCoordTransform
  {
  public:
    VSSimCoordTransform(const std::pair<double,double>& src_radec =
			s_default_src_radec,
			double wobble_phi_deg = s_default_wobble_phi_deg);
    virtual ~VSSimCoordTransform();

    virtual void transform(SEphem::SphericalCoords& radec,
			   unsigned corsika_particle_id,
			   const SEphem::SphericalCoords& azzn, 
			   const SEphem::SphericalCoords& obs_azzn);

    virtual void setObsRADec(const SEphem::SphericalCoords& obs_radec)
    {
      m_obs_radec = obs_radec;
    }

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");

  private:
    SEphem::SphericalCoords m_src_radec;
    SEphem::SphericalCoords m_obs_radec;
    double                  m_wobble_phi_rad;

    // Default options --------------------------------------------------------
    static std::pair< double, double > s_default_src_radec;
    static double                      s_default_wobble_phi_deg;
  };
 
} // namespace VERITAS


#endif // VSSIMCOORDTRANSFORM_HPP
