/*! \file VBFSimpleRunInfoExtractor.hpp

  Extract simple run info from the VBF file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/25/2006

  $Id: VBFSimpleRunInfoExtractor.hpp,v 3.3 2007/12/03 22:14:35 sfegan Exp $

*/

#ifndef VBFSIMPLERUNINFOEXTRACTOR_HPP
#define VBFSIMPLERUNINFOEXTRACTOR_HPP

#include <VSTime.hpp>
#include <VSSimpleVBF.hpp>

namespace VERITAS
{
  // ==========================================================================
  // SIMPLE RUN INFO EXTRACTOR
  // ==========================================================================

  class VBFSimpleRunInfoExtractor: public VSSimpleVBFVisitor
  {
  public:
    VBFSimpleRunInfoExtractor();
    virtual ~VBFSimpleRunInfoExtractor();

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
#warning What do I do with clock trigger data ?
				 const VEvent*          event);

    bool hasRunNumber() const { return m_run_number_set; }
    unsigned runNumber() const { return m_run_number; }

    bool hasApproximateStartTime() const 
    { return m_approximate_start_time_set; }
    const VSTime& approximateStartTime() const 
    { return m_approximate_start_time; }

    bool hasSubarrayTelescopeList() const 
    { return m_subarray_telescope_list_set; }
    const std::vector<unsigned> subarrayTelescopeList() const
    { return m_subarray_telescope_list; }

  private:
    unsigned              m_components_found;
    unsigned              m_events_found;
    bool                  m_run_number_set;
    unsigned              m_run_number;
    bool                  m_approximate_start_time_set;
    VSTime                m_approximate_start_time;
    bool                  m_subarray_telescope_list_set;
    std::vector<unsigned> m_subarray_telescope_list;
  };

}

#endif // defined VBFSIMPLERUNINFOEXTRACTOR_HPP
