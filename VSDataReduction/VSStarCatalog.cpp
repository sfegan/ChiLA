//-*-mode:c++; mode:font-lock;-*-

/*! \file VSStarCatalogchpp

  Catalog of bright stars.

  Default has list of stars from SKY2000 Catalog, Version 4 (Myers+ 2002)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/27/2007

  $Id: VSStarCatalog.cpp,v 3.6 2008/03/13 18:23:18 matthew Exp $

*/

#include <cmath>
#include "VSStarCatalog.hpp"

using namespace VERITAS;
using namespace SEphem;

#define NUMEL(x) (sizeof(x)/sizeof(*x))

struct RawStar
{
  double ra_j2000;
  double dec_j2000;
  double v_mag;
};

static RawStar raw_catalog[] = {
#include "star_catalog_radec_j2000_Vmag_sorted.h"
};

std::auto_ptr<VSStarCatalog> VSStarCatalog::s_instance;

VSStarCatalog* VSStarCatalog::getInstance()
{
  if(!s_instance.get())s_instance.reset(new VSStarCatalog);
  return s_instance.get();
}

void VSStarCatalog::
getStars(std::vector<Star>& stars, const SphericalCoords& radec_center,
	 const Angle max_distance, double max_mag)
{
  stars.clear();
  const unsigned nstar = NUMEL(raw_catalog);
  unsigned istar = 0;
  while(istar<nstar && raw_catalog[istar].v_mag <= max_mag)
    {
      SphericalCoords radec(Angle::frCoAngle(raw_catalog[istar].dec_j2000),
			    raw_catalog[istar].ra_j2000);
      if(radec_center.separation(radec) <= max_distance)
	stars.push_back(Star(radec, raw_catalog[istar].v_mag));
      istar++;
    }
}
