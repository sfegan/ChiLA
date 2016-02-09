//-*-mode:c++; mode:font-lock;-*-

/*! \file VSRunInfoData.hpp

  Data extracted from run

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/13/2006

  $Id: VSRunInfoData.hpp,v 3.1 2007/08/21 19:58:47 sfegan Exp $

*/

#ifndef VSRUNINFODATA_HPP
#define VSRUNINFODATA_HPP

#include <VBFRunInfo.hpp>

namespace VERITAS
{
  typedef VBFRunInfo::Data VSRunInfoData;
  typedef VBFRunInfo::SimData VSSimInfoData;
}

#endif // VSRUNINFODATA_HPP
