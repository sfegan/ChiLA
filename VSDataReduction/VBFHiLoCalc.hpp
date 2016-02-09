//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFHiLoCalc.hpp
     Class designed for use with logain.cpp...

  \author     Timothy C. Arlen            \n
              UCLA                        \n
	      arlen@astro.ucla.edu        \n

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/18/2005

  $Id: VBFHiLoCalc.hpp,v 3.1 2009/12/22 19:43:30 matthew Exp $

*/

#ifndef VBFHILOCALC_HPP
#define VBFHILOCALC_HPP

#include<vector>

// Not sure if all these are needed??
#include "VBFSimplePeds.hpp"
#include "VSSimpleHist.hpp"
#include "VBFLaserCalc.hpp"
#include "VSHiLoData.hpp"
#include "VSFileUtility.hpp"
#include "VSChannelMap.hpp"

// ============================================================================
// VBFHiLoCalc
// ============================================================================

// (Every new ChiLA class goes into the VERITAS namespace)
namespace VERITAS
{
  class VBFHiLoCalc: public VSSimpleVBFVisitor
  {
  public:
    struct ChanData
    {
      ChanData(): ped_stat(), nevent() { }
      
      // Statistics counters -------------------------------------------------
      VSSimpleStat1<double>               ped_stat;     // Amplitude
      unsigned                            nevent;
    };

    struct ScopeData
    {
      ScopeData(unsigned nchan=0): chan(nchan,ChanData()), nevent() { }

      // Channel data
      std::vector<ChanData>               chan;

      // Statistics counters --------------------------------------------------
      unsigned                            nevent;
    };

    typedef std::vector<ScopeData*> ArrayData;

    VBFHiLoCalc(int sample_0, unsigned sample_N);
    virtual ~VBFHiLoCalc();
    
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
    
    virtual void visitChannel(bool& veto_channel, void* user_data,
			      uint32_t                  channel_num, 
			      bool                      hit, 
			      bool                      trigger);
    
    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);
    
    virtual void leaveScopeEvent(bool veto_scope_event, void* user_data);

    virtual void leaveArrayEvent(bool veto_array_event, void* user_data);
    
    void getData(VSHiLoData& data, unsigned nevent_min = 10) const;

    ArrayData               scope;
    
  private:
    VBFHiLoCalc(VBFLaserCalc&);
    VBFHiLoCalc& operator= (const VBFLaserCalc&);

    // Settings
    int                     m_sample_0;
    unsigned                m_sample_N;

    // State
    unsigned                m_runno;
    unsigned                m_scope_id;
  };

}

#endif // VBFHILOCALC_HPP
