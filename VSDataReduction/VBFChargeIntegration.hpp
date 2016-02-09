//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFChargeIntegration.cpp

  Calculate pedestal information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFChargeIntegration.hpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#ifndef VBFCHARGEINTEGRATION_HPP
#define VBFCHARGEINTEGRATION_HPP

#include<vector>

#include "VSSimpleStat.hpp"
#include "VSSimpleVBF.hpp"

// ============================================================================
// VBFChargeIntegration
// ============================================================================

namespace VERITAS
{

  class VBFChargeIntegration: public VSSimpleVBFVisitor
  {
  public:
    VBFChargeIntegration(unsigned sample_0, unsigned sample_N): 
      VSSimpleVBFVisitor(),
      fSample0(sample_0), fSampleN(sample_N) { }
    virtual ~VBFChargeIntegration();
    
    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 const std::vector<uint32_t>& samples,
				 const std::vector<uint32_t>& integrated);

    virtual void visitCharge(void* user_data,
			     uint32_t                   channel_num,
			     bool                       lo_gain,
			     uint32_t                   charge, 
			     unsigned                   sample0,
			     unsigned                   samplen);

  protected:
    unsigned fSample0;
    unsigned fSampleN;
  };
  
}

#endif // VBFCHARGEINTEGRATION_HPP
