//-*-mode:c++; mode:font-lock;-*-

/*! \file VSStarCatalog.hpp

  Catalog of bright stars.

  Default has list of stars from SKY2000 Catalog, Version 4 (Myers+ 2002)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/27/2007

  $Id: VSStarCatalog.hpp,v 3.5 2008/06/27 21:53:09 matthew Exp $

*/

#ifndef VSSTARCATALOG_HPP
#define VSSTARCATALOG_HPP

#include<vector>
#include<memory>
#include<SphericalCoords.h>

namespace VERITAS
{

  class VSStarCatalog
  {
  public:
    
    struct Star
    {
      Star(): radec_j2000(), v_mag() { }
      Star(const SEphem::SphericalCoords& r, double m): 
	radec_j2000(r), v_mag(m) { }
      SEphem::SphericalCoords radec_j2000;
      double v_mag;
    };

    static void getStars(std::vector<Star>& stars, 
			 const SEphem::SphericalCoords& radec_center,
			 const SEphem::Angle max_distance,
			 double max_mag = 4.5);
    
    static VSStarCatalog* getInstance();

  private:
    VSStarCatalog() { }
    VSStarCatalog(const VSStarCatalog& o);
    VSStarCatalog& operator=(const VSStarCatalog& o);

    static std::auto_ptr<VSStarCatalog> s_instance;
  };

}

#endif // ifndef VSSTARCATALOG_HPP
