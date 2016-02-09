//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFPrintFrequency.hpp

  Write the packet count to a stream every so often

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/13/2005

  $Id: VBFPrintFrequency.hpp,v 3.1 2007/04/20 23:33:15 sfegan Exp $

*/

#ifndef VBFPRINTFREQUENCY_HPP
#define VBFPRINTFREQUENCY_HPP

#include <iostream>

#include <VSSimpleVBF.hpp>

namespace VERITAS
{

  class VBFPrintFrequency: public VSSimpleVBFVisitor
  {
  public:
    VBFPrintFrequency(std::ostream& stream, unsigned frequency=1000):
      VSSimpleVBFVisitor(), m_stream(stream), m_frequency(frequency) { }
    virtual ~VBFPrintFrequency();

    virtual void visitPacket(bool& veto_packet, void*& user_data,
			     uint32_t                   seq_packet_number,
			     uint32_t                   vbf_packet_number,
			     bool                       has_array_event,
			     bool                       has_sim_header,
			     bool                       has_sim_event,
			     uint32_t                   num_overflow_datum,
			     const VBFPacket*           packet);

  private:
    std::ostream& m_stream;
    unsigned m_frequency;
  };

}

#endif // not defined VBFPRINTFREQUENCY_HPP
