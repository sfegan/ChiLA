//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFSimplePeds.cpp

  Calculate pedestal information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFSimplePeds.hpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#ifndef VBFSIMPLEPEDS_HPP
#define VBFSIMPLEPEDS_HPP

#include<vector>

#include "VSSimpleStat.hpp"
#include "VSSimplePedData.hpp"
#include "VBFChargeIntegration.hpp"

// ============================================================================
// VBFSimplePeds
// ============================================================================

namespace VERITAS
{

  class VBFSimplePeds: public VSSimpleVBFVisitor
  {
  public:
    VBFSimplePeds(unsigned sample_separation = 3, 
		  unsigned print_frequency = 0): 
      VSSimpleVBFVisitor(), 
      m_sample_separation(sample_separation), m_scope_num(), m_data(),
      m_print_frequency(print_frequency)
    { /* nothing to see here */ }
    virtual ~VBFSimplePeds();
    
    virtual void visitArrayTrigger(bool& veto_array_event, void* user_data,
				   uint32_t             event_num,
				   const VEventType&    event_type,
				   uint32_t             trigger_mask,
				   uint32_t             flags,
				   const VSTime&        raw_time,
				   uint32_t             at_flags,
				   uint32_t             config_mask,
				   uint32_t             num_telescopes,
				   uint32_t             num_trigger_telescopes,
				   uint32_t             run_number,
				   const uint32_t*      ten_mhz_clocks,
				   const uint32_t*      cal_count,
				   const uint32_t*      ped_count,
				   const VArrayTrigger* trigger);

    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);
        
    virtual void visitScopeEvent(bool& veto_scope_event, void* user_data,
				 uint32_t               event_num,
				 uint32_t               telescope_num, 
				 const VEventType&      event_type,
				 uint32_t               trigger_mask,
				 uint32_t               flags,
				 const VSTime&          raw_time,
				 uint32_t               num_samples,
				 uint32_t               num_channels_saved,
				 uint32_t               num_channels_total,
				 uint32_t               num_clock_trigger,
				 const VEvent*          event);
    
    void getSimplePedData(VSSimplePedData& data);
    
  protected:
    unsigned m_sample_separation;
    unsigned m_scope_num;
    std::vector<std::vector<std::vector<VSSimpleStat2<double> > > > m_data;
    unsigned m_print_frequency;
  };
  
}

#endif // VBFSIMPLEPEDS_HPP
