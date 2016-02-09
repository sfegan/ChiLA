//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDataReductionsConstants.hpp

  Types needed in more than one place

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/28/2007

*/

#ifndef VSDATAREDUCTIONTYPES_HPP
#define VSDATAREDUCTIONTYPES_HPP

#include<VSDataConverter.hpp>

namespace VERITAS
{

  typedef triple<double, double, double> pos_type;
  typedef std::pair<std::string, std::string> target_type;

}

#endif // not defined VSDATAREDUCTIONTYPES_HPP
