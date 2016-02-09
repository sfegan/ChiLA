//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFSimplePeds.cpp

  Calculate pedestal information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFChargeIntegration.cpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#include "VBFChargeIntegration.hpp"

using namespace VERITAS;

// ============================================================================
// VBFChargeIntegration
// ============================================================================

VBFChargeIntegration::~VBFChargeIntegration()
{
  // nothing to see here
}

void VBFChargeIntegration::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		const std::vector<uint32_t>& samples,
		const std::vector<uint32_t>& integrated)
{
  uint32_t my_charge=0;
  unsigned nsample = samples.size(); 
  if((nsample>fSampleN)&&(fSampleN>0))nsample=fSampleN;

  for(unsigned isample=fSample0; isample<nsample; isample++)
    my_charge += samples[isample];
  if(lo_gain)my_charge*=6;
  visitCharge(user_data,channel_num,lo_gain,my_charge,fSample0,nsample);
}

void VBFChargeIntegration::
visitCharge(void* user_data,
	    uint32_t                   channel_num,
	    bool                       lo_gain,
	    uint32_t                   charge, 
	    unsigned                   sample0,
	    unsigned                   samplen)
{
  // nothing to see here
}
