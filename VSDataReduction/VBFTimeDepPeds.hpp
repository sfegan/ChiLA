//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFTimeDepPeds.cpp

  Calculate time dependent pedestal information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VBFTimeDepPeds.hpp,v 3.3 2007/06/13 23:30:28 sfegan Exp $

*/

#ifndef VBFTIMEDEPPEDS_HPP
#define VBFTIMEDEPPEDS_HPP

#include <vector>
#include <list>

#include <VSSimpleStat.hpp>
#include <VSTimeDepPedData.hpp>
#include <VSSimpleVBF.hpp>
#include <VBFRunInfo.hpp>
#include <VSTimeDepSupData.hpp>
#include <VSTimeDepPedData.hpp>

namespace VERITAS
{

  class VBFTimeDepPeds: public VSSimpleVBFVisitor
  {
  public:
    VBFTimeDepPeds(const std::vector<VBFRunInfo::Slice>& slice,
		   const VSTimeDepSupData* suppress,
		   unsigned window_min = 1, unsigned sample_separation = 3);
    virtual ~VBFTimeDepPeds();
    
    virtual void visitArrayEvent(bool& veto_array_event, void* user_data,
				 uint32_t               num_scope_events,
				 bool                   has_array_trigger,
				 bool                   has_event_number,
				 uint32_t               event_num,
 				 bool                   has_good_event_time,
				 const VSTime&          best_event_time,
				 EventType              event_type,
				 uint32_t               l2_trigger_mask,
				 const VArrayEvent*     array_event);

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

    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);

    void getData(VSTimeDepPedData& data) const;
    
  protected:

    struct Chan
    {
      Chan(): ped_event_count(), window_stat() { }
      unsigned ped_event_count;
      std::vector<VSSimpleStat2<double> > window_stat;
    };

    typedef std::vector<Chan> Scope;

    struct Slice
    {
      Slice()
	: event_num_lo(), event_num_hi(), event_time_lo(), event_time_hi(),
	  ped_event_count(), scope() { /* nothing to see here */ }
      unsigned event_num_lo;
      unsigned event_num_hi;
      VSTime event_time_lo;
      VSTime event_time_hi;
      unsigned ped_event_count;
      std::vector<Scope> scope;
    };

    typedef VSSimpleStat2<double> ChanPed;
    typedef std::vector<ChanPed> ScopePed;

    typedef VSSimpleHist<int32_t> ChanHist;
    typedef std::vector<ChanHist> ScopeHist;

    // SETTINGS ---------------------------------------------------------------
    unsigned                   m_window_min;
    unsigned                   m_sample_separation;
    const VSTimeDepSupData*    m_suppress;

    // DATA -------------------------------------------------------------------
    std::vector<Slice>         m_slice;
    std::vector<ScopePed>      m_hi_ped;
    std::vector<ScopeHist>     m_min_window_hist;

    // INTERMEDIATE -----------------------------------------------------------
    unsigned                   m_islice;
    unsigned                   m_suppress_islice;
    unsigned                   m_suppress_iword;
    unsigned                   m_suppress_mask;
    unsigned                   m_scope_num;
    std::vector<unsigned>      m_nchan;
    std::vector<unsigned>      m_nsamples;

    // TEMPORARY DATA ---------------------------------------------------------
    bool                       m_is_pedestal;

  private:
    VBFTimeDepPeds(const VBFTimeDepPeds&);
    VBFTimeDepPeds& operator= (const VBFTimeDepPeds&);
  };  

}

#endif // VBFTIMEDEPPEDS_HPP
