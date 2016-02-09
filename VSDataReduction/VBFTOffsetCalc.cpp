//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFTOffsetCalc.cpp

  Calculate timing information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFTOffsetCalc.cpp,v 3.2 2008/02/23 22:41:34 sfegan Exp $

*/

#include<vsassert>

#include"VBFTOffsetCalc.hpp"

using namespace VERITAS;

VBFTOffsetCalc::
VBFTOffsetCalc(VSTimingCalc* calc,
	       double threshold_dc, unsigned threshold_nchan,
	       const std::vector<std::vector<unsigned> >& boards_per_crate,
	       const std::vector<std::vector<unsigned> >& l2_pulse_channel,
	       const VSSimplePedData& peds, 
	       bool fill_hist, bool include_lo_gain, double hist_res)
  : m_data(), m_calc(calc),  m_scope_id(0),
    m_threshold_dc(threshold_dc), m_threshold_nchan(threshold_nchan),
    m_fill_hist(fill_hist), m_include_lo_gain(include_lo_gain),
    m_nflash_chan()
{ 
  unsigned nscope = peds.m_data.size();
  vsassert(boards_per_crate.size() >= nscope);
  vsassert(l2_pulse_channel.size() >= nscope);

  m_data.resize(nscope);

  for(unsigned iscope=0; iscope<nscope; iscope++)
    if(peds.m_data[iscope].size() != 0)
      {
	ScopeData* data = 
	  new ScopeData(hist_res,
			boards_per_crate[iscope], l2_pulse_channel[iscope]);
	m_data[iscope] = data;
	unsigned nchan = peds.m_data[iscope].size();
	data->chan.resize(nchan,data->hist_res);
	unsigned crate = 0;
	unsigned board = 0;
	for(unsigned ichan=0;ichan<nchan; ichan++)
	  {
	    data->chan[ichan].global_ped = 
	      peds.m_data[iscope][ichan].ped;
	    data->chan[ichan].global_suppressed = 
	      peds.m_data[iscope][ichan].suppressed;
	    data->chan[ichan].global_crate=crate;
	    if(ichan%10 == 9)board++;
	    if((crate<=data->boards_per_crate.size())
	       &&(board==data->boards_per_crate[crate]))board=0,crate++;
	  }
	data->stat_mean_crate_time.resize(data->boards_per_crate.size());
	for(unsigned icrate=0;icrate<data->l2_pulse_channel.size();icrate++)
	  if(data->l2_pulse_channel[icrate]<nchan)
	    data->chan[data->l2_pulse_channel[icrate]].global_suppressed =
	      true;
      }
}

VBFTOffsetCalc::~VBFTOffsetCalc()
{
  // nothing to see here
}

void VBFTOffsetCalc::
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
  if(event_type.trigger == VEventType::PED_TRIGGER)veto_array_event=true;
}

void VBFTOffsetCalc::
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
  m_scope_id = telescope_num;
  vsassert(m_data[m_scope_id] != 0);
  vsassert(m_data[m_scope_id]->chan.size() <= num_channels_total);
  m_nflash_chan=0;
}

void VBFTOffsetCalc::
visitChannel(bool& veto_channel, void* user_data,
	     uint32_t                  channel_num, 
	     bool                      hit, 
	     bool                      trigger)
{
  m_data[m_scope_id]->chan[channel_num].ev_hit = hit;
}

void VBFTOffsetCalc::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		unsigned               nsample,
		const uint32_t*        samples,
		const uint32_t*        integrated)
{
  ScopeData* data = m_data[m_scope_id];
  double signal;
  double time;
  double ped = data->chan[channel_num].global_ped;
  if(nsample > data->max_nsample)data->max_nsample=nsample;
  if((!m_calc->calc(signal,time,lo_gain,6.0,nsample,samples,integrated,ped))||
     (lo_gain && !m_include_lo_gain))
    {
      data->chan[channel_num].ev_hit = false;
    }
  else
    {
      data->chan[channel_num].ev_signal = signal;
      data->chan[channel_num].ev_time = time;

      // ----------------------------------------------------------------------
      // Add the channel to the tally if it exceeds threshold
      // ----------------------------------------------------------------------
      if((signal>=m_threshold_dc)
	 &&(!data->chan[channel_num].global_suppressed))
	m_nflash_chan++;
    }
}

void VBFTOffsetCalc::leaveScopeEvent(bool veto_scope_event, void* user_data)
{
  if(veto_scope_event)return;
  ScopeData* data = m_data[m_scope_id];
  data->stat_mean_nflash.accumulate(m_nflash_chan);
  if(m_nflash_chan>=m_threshold_nchan)
    {
      unsigned ncrate = data->boards_per_crate.size();
      VSSimpleStat2<double> stat_ev_signal;  
      std::vector<VSSimpleStat2<double> > stat_ev_time(ncrate);

      unsigned nchan = data->chan.size();
      for(unsigned ichan=0;ichan<nchan; ichan++)
	{
	  double ev_signal = data->chan[ichan].ev_signal;
	  double ev_time = data->chan[ichan].ev_time;
	  unsigned crate = data->chan[ichan].global_crate;

	  if(data->chan[ichan].ev_hit)
	    {
	      if(ev_signal>=m_threshold_dc)
		{
		  data->chan[ichan].stat_time.accumulate(ev_time);
		  if(m_fill_hist)
		    data->chan[ichan].hist_time.accumulate(ev_time);
		}
	      data->chan[ichan].stat_signal.accumulate(ev_signal);

	      if(!data->chan[ichan].global_suppressed)
		{
		  if(ev_signal>=m_threshold_dc)
		    stat_ev_time[crate].accumulate(ev_time);
		  stat_ev_signal.accumulate(ev_signal);
		}
	    }
	}
      if(stat_ev_signal.count())
	data->stat_mean_signal.accumulate(stat_ev_signal.mean());
      for(unsigned icrate=0; icrate<ncrate; icrate++)
	if(stat_ev_time[icrate].count())
	  data->stat_mean_crate_time[icrate].
	    accumulate(stat_ev_time[icrate].mean());
    }
}
