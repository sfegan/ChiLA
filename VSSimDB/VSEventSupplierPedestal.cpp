//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventSupplierPedestal.cpp

  Generate a pedestal event.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       08/25/2006
  \note
*/

#include <cassert>
#include <vector>
#include <algorithm>

#include "VSEventSupplierPedestal.hpp"

using namespace VERITAS;

VSEventSupplierPedestal::
VSEventSupplierPedestal(unsigned nevents,
			unsigned nscope, unsigned nchan)
  : VSEventSupplier(), fNEvent(nevents),
    fNScope(nscope), fNChannel(nchan), fIEvent()
{
  // nothing to see here
}

VSEventSupplierPedestal::~VSEventSupplierPedestal()
{
  // nothing to see here
}

std::vector<VSEventSupplier::SimParam> VSEventSupplierPedestal::getSimParam()
{
  std::vector<SimParam> tp_vec;
  return tp_vec;
}

bool VSEventSupplierPedestal::getNextEvent(Event& e)
{

  if((fNEvent >0)&&(fIEvent >= fNEvent))return false;
  
  e.fEventID                 = fIEvent+1;
  e.fTargetZenithRad         = 0;
  e.fTargetAzimuthRad        = 0;
  e.fPrimarySpeciesCORSIKA   = 0;
  e.fPrimaryEnergyTeV        = 0;
  e.fPrimaryZenithRad        = 0;
  e.fPrimaryAzimuthRad       = 0;
  e.fPrimaryCoreEastM        = 0;
  e.fPrimaryCoreNorthM       = 0;
  e.fPrimaryCoreUpASLM       = 0;
  e.fSamplingRadiusM         = 0;
  e.fTableIndex              = 0;
  e.fEventComplete           = 0;

  fIEvent++;

  return true;
}

uint32_t VSEventSupplierPedestal::getNumEvents()
{
  return fNEvent;
}
