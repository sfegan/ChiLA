#include <VSRTraceVisitor.hpp>

using namespace VERITAS;

VSRTraceVisitor::VSRTraceVisitor(std::ostream& stream, 
				 unsigned scope,
				 unsigned channel):
  VSSimpleVBFVisitor(), m_stream(stream), m_channels() 
{ 
  m_channels[scope].insert(channel);
}

VSRTraceVisitor::
VSRTraceVisitor(std::ostream& stream, 
		const std::vector< std::pair< unsigned,unsigned > >& channels):
  VSSimpleVBFVisitor(), m_stream(stream), m_channels() 
{ 
  for(std::vector< std::pair< unsigned,unsigned > >::const_iterator itr = 
	channels.begin(); itr != channels.end(); itr++)
    m_channels[itr->first].insert(itr->second);
}

VSRTraceVisitor::~VSRTraceVisitor()
{
  // nothing to see here
}

void VSRTraceVisitor::
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
		const VEvent*          event)
{
  m_current_scope = telescope_num;
  m_scope_event_type = event_type.trigger;

  
}

void VSRTraceVisitor::visitHitChannel(void* user_data,
				      uint32_t               channel_num,
				      uint32_t               charge, 
				      uint32_t               pedestal,
				      bool                   lo_gain,
				      unsigned               nsample,
				      const uint32_t*        samples,
				      const uint32_t*        integrated)
{
  if(m_current_scope == 0)
    {
      std::cout << " " << channel_num;
      
      for(unsigned isample=0;isample<nsample;isample++)
	{
	  std::cout << ' ' << samples[isample];
	}
      std::cout << std::endl;
    }

  if(m_scope_event_type == VEventType::PED_TRIGGER) return;
  if(m_channels.find(m_current_scope) == m_channels.end()) return;
  else if(m_channels[m_current_scope].find(channel_num) == 
	  m_channels[m_current_scope].end()) return;


  if(m_trace[m_current_scope][channel_num].size() < nsample)      
    m_trace[m_current_scope][channel_num].resize(nsample);

  m_stream 
    << "HITCHAN:  " << channel_num << ' ' << charge << ' ' 
    << pedestal << ' ' << lo_gain << ' ' << m_scope_event_type;

  m_nevent[m_current_scope][channel_num]++;

  for(unsigned isample=0;isample<nsample;isample++)
    {
      m_stream << ' ' << samples[isample];
      m_trace[m_current_scope][channel_num][isample] += samples[isample];
    }

  m_stream << ' ' << integrated[nsample-1] << '\n';
}

std::vector< double > VSRTraceVisitor::getAverageTrace(unsigned scope, 
						       unsigned channel)
{
  vsassert(m_channels.find(scope) != m_channels.end());
  vsassert(m_channels[scope].find(channel) != m_channels[scope].end());

  std::vector< double > avg_trace;

  const unsigned nsample = m_trace[scope][channel].size();
  if(avg_trace.size() < m_trace[scope][channel].size())
    avg_trace.resize(m_trace[scope][channel].size());

  for(unsigned isample = 0; isample < nsample; isample++)
    {
      avg_trace[isample] = 
	m_trace[scope][channel][isample]/m_nevent[scope][channel];
    }
  
  return avg_trace;
}
