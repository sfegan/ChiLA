//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFTimingCalc.cpp

  Calculate timing information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFTOffsetCalc.hpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#ifndef VBFTIMINGCALC_HPP
#define VBFTIMINGCALC_HPP

#include<vector>

#include "VSSimpleStat.hpp"
#include "VSSimpleHist.hpp"
#include "VSTimingCalc.hpp"
#include "VSSimpleVBF.hpp"
#include "VSSimplePedData.hpp"

// ============================================================================
// VBFTOffsetCalc
// ============================================================================

namespace VERITAS
{

  class VBFTOffsetCalc: public VSSimpleVBFVisitor
  {
  public:
    
    struct ChanData
    {
      ChanData(double hist_res):
	stat_time(), stat_signal(), 
	global_suppressed(), global_ped(), global_crate(),
	ev_hit(), ev_time(), ev_signal(), hist_time(hist_res) { }
      
      VSSimpleStat2<double>               stat_time;
      VSSimpleStat2<double>               stat_signal;
      
      bool                                global_suppressed;
      double                              global_ped;
      unsigned                            global_crate;
      
      bool                                ev_hit;
      double                              ev_time;
      double                              ev_signal;
      
      VSSimpleHist<double>                hist_time;
    };
    
    struct ScopeData
    {
      ScopeData(double _hist_res,
		const std::vector<unsigned>& _boards_per_crate,
		const std::vector<unsigned>& _l2_pulse_channel,
		unsigned nchan=0):
	chan(nchan,ChanData(_hist_res)),
	stat_mean_nflash(), stat_mean_signal(), stat_mean_crate_time(),
	max_nsample(0), hist_res(_hist_res),
	boards_per_crate(_boards_per_crate), 
	l2_pulse_channel(_l2_pulse_channel)
      { }
      std::vector<ChanData>               chan;
      VSSimpleStat2<double>               stat_mean_nflash;
      VSSimpleStat2<double>               stat_mean_signal;  
      std::vector<VSSimpleStat2<double> > stat_mean_crate_time;
      unsigned                            max_nsample;
      double                              hist_res;
      std::vector<unsigned>               boards_per_crate;
      std::vector<unsigned>               l2_pulse_channel;
    };

    typedef std::vector<ScopeData*> ArrayData;

    VBFTOffsetCalc(VSTimingCalc* calc,
		   double threshold_dc, unsigned threshold_nchan,
		   const std::vector<std::vector<unsigned> >& boards_per_crate,
		   const std::vector<std::vector<unsigned> >& l2_pulse_channel,
		   const VSSimplePedData& peds,
		   bool fill_hist = false, bool include_lo_gain = false,
		   double hist_res = 0.5);
    
    virtual ~VBFTOffsetCalc();
    
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
    
    ArrayData               m_data;
    
  private:
    VBFTOffsetCalc(VBFTOffsetCalc&);
    VBFTOffsetCalc& operator= (const VBFTOffsetCalc&);

    VSTimingCalc*           m_calc;
    unsigned                m_scope_id;
    double                  m_threshold_dc;
    unsigned                m_threshold_nchan;
    bool                    m_fill_hist;
    bool                    m_include_lo_gain;
    unsigned                m_nflash_chan;
  };
  

}

#endif // VBFTIMINGCALC_HPP
