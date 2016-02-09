//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFSimplePeds.cpp

  Calculate pedestal information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFSimplePeds.cpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#include "VBFSimplePeds.hpp"

using namespace VERITAS;

// ============================================================================
// VBFSimplePeds
// ============================================================================

VBFSimplePeds::~VBFSimplePeds()
{
  // nothing to see here
}

void VBFSimplePeds::
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
  if(m_print_frequency && (event_num % m_print_frequency == 0))
    std::cerr << event_num << std::endl;
  
  if(event_type.trigger != VEventType::PED_TRIGGER)veto_array_event=true;
}

void VBFSimplePeds::
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
  m_scope_num = telescope_num;
  if(m_scope_num >= m_data.size())
    m_data.resize(m_scope_num+1);
  if(num_channels_total >= m_data[m_scope_num].size())
    m_data[m_scope_num].resize(num_channels_total);
}

void VBFSimplePeds::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		unsigned               nsample,
		const uint32_t*        samples,
		const uint32_t*        integrated)
{
  if(lo_gain)return;

  if(m_data[m_scope_num][channel_num].size() < nsample)
    m_data[m_scope_num][channel_num].resize(nsample);
  
  for(unsigned iwindow=0;iwindow<nsample;iwindow++)
    {
      m_data[m_scope_num][channel_num][iwindow].
	accumulate(integrated[iwindow]);
      unsigned isample=iwindow+1;
      while(isample+iwindow+m_sample_separation < nsample)
	{
	  isample += m_sample_separation;
	  unsigned Q = integrated[isample+iwindow]-integrated[isample-1];
	  m_data[m_scope_num][channel_num][iwindow].accumulate(Q);
	  isample += iwindow+1;
	}
    }
}

void VBFSimplePeds::getSimplePedData(VSSimplePedData& data)
{
  unsigned nscope = m_data.size();
  data.m_data.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nchan = m_data[iscope].size();
      data.m_data[iscope].resize(nchan);
      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  unsigned nsample = m_data[iscope][ichan].size();
	  if(m_data[iscope][ichan][nsample-1].count())
	    {
	      const double ped = 
		m_data[iscope][ichan][nsample-1].mean()/double(nsample);
	      data.m_data[iscope][ichan].ped = ped;
	      data.m_data[iscope][ichan].dev.clear();
	      data.m_data[iscope][ichan].dev.resize(nsample,0);
	      for(unsigned isample=0;isample<nsample;isample++)
		{
		  // Should use the pedestal calculated from the
		  // largest window size, and appropriately scaled, as
		  // the baseline for the calculation of the RMS. This
		  // is done using the chi2() function of the
		  // SimpleStat2 rather than the dev(), which
		  // implicitly uses the sampled mean.
		  const double isampleped = ped*double(isample+1);
		  const double chi2 = 
		    m_data[iscope][ichan][isample].chi2(isampleped);
		  data.m_data[iscope][ichan].dev[isample] = sqrt(chi2);
		}
	      data.m_data[iscope][ichan].suppressed = false;
	    }
	  else
	    {
	      data.m_data[iscope][ichan].ped = 0;
	      data.m_data[iscope][ichan].dev.clear();
	      data.m_data[iscope][ichan].dev.resize(nsample,0);
	      data.m_data[iscope][ichan].suppressed = true;
	    }
	}
    }
}
