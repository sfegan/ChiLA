/*! \file VBFSimpleRunInfoExtractor.cpp

  Extract simple run info from the VBF file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/25/2006

  $Id: VBFSimpleRunInfoExtractor.cpp,v 3.4 2008/06/27 21:53:09 matthew Exp $

*/

#include <algorithm>

#include <VBFSimpleRunInfoExtractor.hpp>

using namespace VERITAS;

VBFSimpleRunInfoExtractor::VBFSimpleRunInfoExtractor(): 
  VSSimpleVBFVisitor(), m_components_found(0), m_events_found(0),
  m_run_number_set(false), m_run_number(0),
  m_approximate_start_time_set(false), m_approximate_start_time(),
  m_subarray_telescope_list_set(false), m_subarray_telescope_list()
{ 
  // nothing to see here
}

VBFSimpleRunInfoExtractor::~VBFSimpleRunInfoExtractor()
{
  // nothing to see here
}

void VBFSimpleRunInfoExtractor::
visitArrayEvent(bool& veto_array_event, void* user_data,
		uint32_t               num_scope_events,
		bool                   has_array_trigger,
		bool                   has_event_number,
		uint32_t               event_num,
		bool                   has_good_event_time,
		const VSTime&          best_event_time,
		EventType              event_type,
		uint32_t               l2_trigger_mask,
		const VArrayEvent*     array_event)
{
  if((!m_approximate_start_time_set)&&(has_good_event_time))
    {
      m_components_found++;
      m_approximate_start_time_set = true;
      m_approximate_start_time = best_event_time;
    }

  if((m_components_found == 3)||(m_events_found == 100))
    {
      if(!m_subarray_telescope_list.empty())
	m_subarray_telescope_list_set = true;
      fDispatchStop->stopProcessingFile();
      veto_array_event = true;
    }

  m_events_found++;
}

void VBFSimpleRunInfoExtractor::
visitArrayTrigger(bool& veto_array_event, void* user_data,
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
		  const VArrayTrigger* trigger)
{
  if((!m_run_number_set)&&(run_number>0))
    {
      m_components_found++;
      m_run_number_set = true;
      m_run_number = run_number;
    }

  if(!m_subarray_telescope_list_set)
    {
      m_components_found++;
      m_subarray_telescope_list_set = true;
      m_subarray_telescope_list.clear();
      for(unsigned iscope=0;iscope<sizeof(config_mask)*8;iscope++)
	if(config_mask & (1<<iscope))
	  m_subarray_telescope_list.push_back(iscope);
    }

  if(m_components_found == 3)
    {
      fDispatchStop->stopProcessingFile();
    }
}

void VBFSimpleRunInfoExtractor::
visitScopeEvent(bool& veto_scope_event, void* user_data,
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
		const VEvent*          event)
{
  if(!m_subarray_telescope_list_set)
    {
      bool this_scope_set = false;
      for(unsigned iscope=0;iscope<m_subarray_telescope_list.size();iscope++)
	if(m_subarray_telescope_list[iscope] == telescope_num)
	  this_scope_set = true;

      if(!this_scope_set)
	{
	  m_subarray_telescope_list.push_back(telescope_num);
	  std::sort(m_subarray_telescope_list.begin(),
		    m_subarray_telescope_list.end());
	}
    }

  veto_scope_event = true;  
}
