//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFTimeDepPeds.cpp

  Calculate time dependent pedestal information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VBFTimeDepPeds.cpp,v 3.7 2007/12/04 18:05:03 sfegan Exp $

*/

#include "VBFTimeDepPeds.hpp"

using namespace VERITAS;

// ============================================================================
// VBFTimeDepPeds
// ============================================================================

VBFTimeDepPeds::
VBFTimeDepPeds(const std::vector<VBFRunInfo::Slice>& slice,
	       const VSTimeDepSupData* suppress,
	       unsigned window_min, unsigned sample_separation)
  : VSSimpleVBFVisitor(),
    m_window_min(window_min), m_sample_separation(sample_separation),
    m_suppress(suppress), m_slice(slice.size()), m_hi_ped(),
    m_min_window_hist(0,ScopeHist(0,1)),
    m_islice(), m_suppress_islice(), m_suppress_iword(), m_suppress_mask(),
    m_scope_num(), m_nchan(), m_nsamples(), m_is_pedestal()
{
  unsigned nslice = slice.size();
  for(unsigned islice=0;islice<nslice;islice++)
    {
      m_slice[islice].event_num_lo  = slice[islice].event_num_lo;
      m_slice[islice].event_num_hi  = slice[islice].event_num_hi;
      m_slice[islice].event_time_lo = slice[islice].event_time_lo;
      m_slice[islice].event_time_hi = slice[islice].event_time_hi;
    }

  if(m_suppress)
    m_suppress->getFastSliceFindInfo(m_suppress_islice,
				     m_suppress_iword, m_suppress_mask);
}

VBFTimeDepPeds::~VBFTimeDepPeds()
{
  // nothing to see here
}

void VBFTimeDepPeds::
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
  m_is_pedestal = false;

  bool lo = event_num>=m_slice[m_islice].event_num_lo;
  bool hi = event_num<m_slice[m_islice].event_num_hi;

  if(!lo)
    {
      while((m_islice)&&(event_num<m_slice[--m_islice].event_num_lo));
      vsassert((m_islice)||(event_num>=m_slice[0].event_num_lo));
    }
  else if(!hi)
    {
      const unsigned nslice = m_slice.size();
      while((m_islice<nslice-1)
	    &&(event_num>=m_slice[++m_islice].event_num_hi));
      vsassert((m_islice<nslice-1)||(event_num<m_slice[nslice-1].event_num_hi));
    }

  if(m_suppress)
    {
      unsigned suppressed_islice = m_suppress->getSliceByEventNum(event_num);
      if(suppressed_islice != m_suppress_islice)
	{
	  m_suppress_islice = suppressed_islice;
	  m_suppress->getFastSliceFindInfo(m_suppress_islice,
					   m_suppress_iword, m_suppress_mask);
	}
    }

  m_slice[m_islice].ped_event_count++;
}

void VBFTimeDepPeds::
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
  if(event_type.trigger != VEventType::PED_TRIGGER)
    {
      veto_array_event=true;
      return;
    }

  m_is_pedestal = true;
}

void VBFTimeDepPeds::
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
  if((m_is_pedestal == false)
     && (event_type.trigger != VEventType::PED_TRIGGER))
    {
      veto_scope_event=true;
      return;
    }

  m_scope_num = telescope_num;

  if(m_scope_num >= m_nchan.size())
    {
      m_nchan.resize(m_scope_num+1);
      m_nsamples.resize(m_scope_num+1);
      m_hi_ped.resize(m_scope_num+1);
      m_min_window_hist.resize(m_scope_num+1);
      for(std::vector<Slice>::iterator islice = m_slice.begin();
	  islice!=m_slice.end(); islice++)
	islice->scope.resize(m_scope_num+1);
    }

  if(num_channels_total > m_nchan[m_scope_num])
    {
      m_nchan[m_scope_num]=num_channels_total;
      m_hi_ped[m_scope_num].resize(num_channels_total);
      m_min_window_hist[m_scope_num].resize(num_channels_total,1);
      for(std::vector<Slice>::iterator islice = m_slice.begin();
	  islice!=m_slice.end(); islice++)
	islice->scope[m_scope_num].resize(num_channels_total);
    }

  if(num_samples > m_nsamples[m_scope_num])
    {
      m_nsamples[m_scope_num]=num_samples;
      for(std::vector<Slice>::iterator islice = m_slice.begin();
	  islice!=m_slice.end(); islice++)
	for(Scope::iterator ichan=islice->scope[m_scope_num].begin();
	    ichan!=islice->scope[m_scope_num].end();ichan++)
	  ichan->window_stat.resize(num_samples-m_window_min+1);
    }
}

void VBFTimeDepPeds::
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

  // The seemingly pointless divide and subsequent multiply allow
  // different events to have different sample sizes
  double pedQ = double(integrated[nsample-1])/double(nsample);
  m_hi_ped[m_scope_num][channel_num].accumulate(pedQ,nsample);

  Chan& chan(m_slice[m_islice].scope[m_scope_num][channel_num]);

  if((m_suppress)
     &&(m_suppress->isSuppressedFast(m_suppress_iword, m_suppress_mask,
				     m_scope_num, channel_num)))return;
  
  chan.ped_event_count++;

  unsigned iwindow=m_window_min;
  if(iwindow<=nsample)
    {
      chan.window_stat[iwindow-m_window_min].accumulate(integrated[iwindow-1]);
      m_min_window_hist[m_scope_num][channel_num].
	accumulate(integrated[iwindow-1]);
      unsigned isample=iwindow;
      while(isample+iwindow+m_sample_separation-1 < nsample)
	{
	  isample += m_sample_separation;
	  unsigned Q = integrated[isample-1+iwindow]-integrated[isample-1];
	  chan.window_stat[iwindow-m_window_min].accumulate(Q);
	  m_min_window_hist[m_scope_num][channel_num].accumulate(Q);
	  isample += iwindow;
	}
      iwindow++;
    }

  while(iwindow<=nsample)
    {
      chan.window_stat[iwindow-m_window_min].accumulate(integrated[iwindow-1]);
      unsigned isample=iwindow;
      while(isample+iwindow+m_sample_separation-1 < nsample)
	{
	  isample += m_sample_separation;
	  unsigned Q = integrated[isample-1+iwindow]-integrated[isample-1];
	  chan.window_stat[iwindow-m_window_min].accumulate(Q);
	  isample += iwindow;
	}
      iwindow++;
    }
}

void VBFTimeDepPeds::getData(VSTimeDepPedData& data) const
{
  data.clear();
  data.m_window_min = m_window_min;

  const unsigned nscope = m_hi_ped.size();
  data.m_hi_ped.resize(nscope);
  data.m_lo_ped.resize(nscope);
  data.m_hi_dev.resize(nscope);
  data.m_min_window_hist.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nchan = m_hi_ped[iscope].size();
      data.m_hi_ped[iscope].resize(nchan);
      data.m_lo_ped[iscope].resize(nchan);
      data.m_hi_dev[iscope].resize(nchan);
      data.m_min_window_hist[iscope].resize(nchan,1);
      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  double hi_ped = -1;
	  double hi_dev = -1;
	  if(m_hi_ped[iscope][ichan].count())
	    {
	      hi_ped = m_hi_ped[iscope][ichan].mean();
	      double N = 
		double(m_hi_ped[iscope][ichan].count()/m_nsamples[iscope]);
	      double meansq = 
		m_hi_ped[iscope][ichan].sumsq()/N*m_nsamples[iscope];
	      double mean   =
		m_hi_ped[iscope][ichan].sum()/N;
	      hi_dev = sqrt(meansq - mean*mean);
	    }
	  data.m_hi_ped[iscope][ichan] = hi_ped;
	  data.m_lo_ped[iscope][ichan] = hi_ped;
	  data.m_hi_dev[iscope][ichan] = hi_dev;
	  data.m_min_window_hist[iscope][ichan] = 
	    m_min_window_hist[iscope][ichan];
	}    
    }

  const unsigned nslice = m_slice.size();
  data.m_slice.resize(nslice);
  for(unsigned islice=0;islice<nslice;islice++)
    {
      data.m_slice[islice].event_num_lo = m_slice[islice].event_num_lo;
      data.m_slice[islice].event_num_hi = m_slice[islice].event_num_hi;
      data.m_slice[islice].event_time_lo = m_slice[islice].event_time_lo;
      data.m_slice[islice].event_time_hi = m_slice[islice].event_time_hi;
      data.m_slice[islice].ped_event_count = m_slice[m_islice].ped_event_count;

      data.m_slice[islice].scope.resize(nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  unsigned nchan = m_slice[islice].scope[iscope].size();
	  if(nchan==0)continue;
	  unsigned nwindow = 
	    m_slice[islice].scope[iscope][0].window_stat.size();

	  data.m_slice[islice].scope[iscope].resize(nchan);
	  for(unsigned ichan=0;ichan<nchan;ichan++)
	    {
	      VSTimeDepPedData::Chan& 
		o_chan(data.m_slice[islice].scope[iscope][ichan]);
	      const Chan& i_chan(m_slice[islice].scope[iscope][ichan]);
	      const double ped = data.m_hi_ped[iscope][ichan];
	      o_chan.ped_event_count = i_chan.ped_event_count;
	      o_chan.window_dev.resize(nwindow);
	      for(unsigned iwindow=0;iwindow<nwindow;iwindow++)
		if(i_chan.window_stat[iwindow].count())
		  {
		    double wped = ped*double(m_window_min+iwindow);
		    o_chan.window_dev[iwindow] =
		      sqrt(i_chan.window_stat[iwindow].chi2(wped));
		  }
		else
		  {
		    o_chan.window_dev[iwindow] = 0;
		  }
	    }
	}
    }
}
