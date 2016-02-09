/*! \file VSNSpaceDataSource.hpp

  Utility class for accessing data from a file.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/08/2007

*/

#include<VSNSpace.hpp>
#include<VSOctaveH5Reader.hpp>
#include<VSDatumElementExtractor.hpp>

#ifndef VSNSPACEDATASOURCE_HPP
#define VSNSPACEDATASOURCE_HPP

namespace VERITAS 
{
  class VSNSpaceDataSource
  {
  public:
    virtual ~VSNSpaceDataSource();
    virtual bool getData(std::vector<VSNSpace::Point>& points, 
			 double& theta0, double& theta1) = 0;

  };
}

#endif // VSNSPACEDATASOURCE_HPP
