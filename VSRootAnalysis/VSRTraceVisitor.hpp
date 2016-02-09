#ifndef VSRTRACEVISITOR_HPP
#define VSRTRACEVISITOR_HPP

#include <VSSimpleVBF.hpp>

namespace VERITAS
{
  class VSRTraceVisitor: public VSSimpleVBFVisitor
  {
  public:
    VSRTraceVisitor(std::ostream& stream, 
		    unsigned scope = 0,
		    unsigned channel = 0);
    VSRTraceVisitor(std::ostream& stream, 
		    const std::vector< std::pair< unsigned, unsigned > >&
		    channels);
    virtual ~VSRTraceVisitor();

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

    std::vector< double > getAverageTrace(unsigned scope, unsigned channel);

  private:
    std::ostream& m_stream;
    std::map< unsigned, std::set< unsigned > > m_channels;
    unsigned                                   m_current_scope;
    VEventType::TriggerType                    m_scope_event_type;

    std::map< unsigned, std::map< unsigned, std::vector< double > > > m_trace;
    std::map< unsigned, std::map< unsigned, double > > m_nevent;
  };

} // namespace VERITAS

#endif // VSRTRACEVISITOR_HPP
