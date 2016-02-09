//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFPrintFrequency.cpp

  Write the packet count to a stream every so often

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/13/2005

  $Id: VBFPrintFrequency.cpp,v 3.1 2007/04/20 23:33:15 sfegan Exp $

*/

#include <VBFPrintFrequency.hpp>

using namespace VERITAS;

VBFPrintFrequency::~VBFPrintFrequency()
{
  // nothing to see here
}

void VBFPrintFrequency::
visitPacket(bool& veto_packet, void*& user_data,
	    uint32_t                   seq_packet_number,
	    uint32_t                   vbf_packet_number,
	    bool                       has_array_event,
	    bool                       has_sim_header,
	    bool                       has_sim_event,
	    uint32_t                   num_overflow_datum,
	    const VBFPacket*           packet)
{
  if((seq_packet_number)&&(seq_packet_number % m_frequency == 0))
    m_stream << seq_packet_number << std::endl;
  veto_packet=true;
}
