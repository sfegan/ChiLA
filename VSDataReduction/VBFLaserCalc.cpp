//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFLaserCalc.cpp

  Calculate timing information

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/05/2006

  $Id: VBFLaserCalc.cpp,v 3.10 2009/10/14 22:03:30 matthew Exp $

*/

#include<vsassert>

#include"VBFLaserCalc.hpp"

using namespace VERITAS;

VBFLaserCalc::
VBFLaserCalc(VSTimingCalc* calc,
	     double threshold_dc, unsigned threshold_nchan,
	     unsigned threshold_nscope, unsigned threshold_nflash, 
	     double ped_suppress_lo, double ped_suppress_hi,
	     double singlepe_dev,
	     const std::vector<std::vector<unsigned> >& boards_per_crate,
	     const std::vector<std::vector<unsigned> >& l2_pulse_channel,
	     const VSAnalysisStage1Data& stage1_data, 
	     bool fill_hist, bool include_lo_gain)
  : m_data(), m_calc(calc), m_scope_id(0), m_nflash_scope(),
    m_nscope(),
    m_threshold_dc(threshold_dc), m_threshold_nchan(threshold_nchan),
    m_threshold_nscope(threshold_nscope), 
    m_threshold_nflash(threshold_nflash), 
    m_singlepe_dev(singlepe_dev),
    m_fill_hist(fill_hist), m_include_lo_gain(include_lo_gain),
    m_stage1_data(stage1_data)
{ 
  const VSTimeDepPedData& peds = m_stage1_data.pedestals;
  const VSRunInfoData& run_info = m_stage1_data.run_info;
  unsigned nscope = peds.nscope();
  vsassert(boards_per_crate.size() >= nscope);
  vsassert(l2_pulse_channel.size() >= nscope);

  for(unsigned iscope = 0; iscope < nscope; iscope++)
    if(run_info.nchan[iscope] > 0) m_nscope++;
  if(m_threshold_nscope == 0) m_threshold_nscope = m_nscope;

  const double hist_res_time = 0.5;
  const double hist_res_signal = 5;

  // --------------------------------------------------------------------------
  // Find the median hi-gain ped dev for each telescope
  // --------------------------------------------------------------------------
  std::vector<double> median_dev(nscope);

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      std::vector<double> ped_dev;
      for(unsigned ichan=0;ichan<peds.nchan(iscope);ichan++)
	ped_dev.push_back(peds.hiDev(iscope,ichan));
      median_dev[iscope] = median(ped_dev);	   
    }

  m_data.resize(nscope);

  for(unsigned iscope=0; iscope<nscope; iscope++)
    if(peds.nchan(iscope) != 0)
      {
	double locut = median_dev[iscope]*ped_suppress_lo;
	double hicut = median_dev[iscope]*ped_suppress_hi;
	ScopeData* data = 
	  new ScopeData(hist_res_time,hist_res_signal,
			boards_per_crate[iscope], l2_pulse_channel[iscope]);
	m_data[iscope] = data;
	unsigned nchan = peds.nchan(iscope);
	data->chan.resize(nchan,ChanData(data->hist_res_time,
					 data->hist_res_signal));
	unsigned crate = 0;
	unsigned board = 0;
	for(unsigned ichan=0;ichan<nchan; ichan++)
	  {
	    data->chan[ichan].global_ped = peds.hiPed(iscope,ichan);
	    // ----------------------------------------------------------------
	    // Define as 'excluded' those channels that fall outside
	    // of the acceptable ped dev range.  These channels will
	    // be left out of any calculations of mean/median
	    // telescope quantities.
	    // ----------------------------------------------------------------
	    if(peds.hiDev(iscope,ichan) < locut ||
	       peds.hiDev(iscope,ichan) > hicut)
	      data->chan[ichan].global_excluded = true;
	    else
	      data->chan[ichan].global_excluded = false;
	    data->chan[ichan].global_crate=crate;
	    if(ichan%10 == 9)board++;
	    if((crate<=data->boards_per_crate.size())
	       &&(board==data->boards_per_crate[crate]))board=0,crate++;
	  }
	data->mean_crate_time_stat.resize(data->boards_per_crate.size());
	for(unsigned icrate=0;icrate<data->l2_pulse_channel.size();icrate++)
	  if(data->l2_pulse_channel[icrate]<nchan)
	    data->chan[data->l2_pulse_channel[icrate]].global_excluded = true;
      }
}

VBFLaserCalc::~VBFLaserCalc()
{
  // nothing to see here
}

void VBFLaserCalc::
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
  m_nflash_scope=0;
}

void VBFLaserCalc::
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
  vsassert(m_scope_id < m_data.size());
  vsassert(m_data[m_scope_id] != 0);
  vsassert(m_data[m_scope_id]->chan.size() <= num_channels_total);
  m_data[m_scope_id]->ev_nchan_flash = 0;
  m_data[m_scope_id]->ev_nchan_logain = 0;
  m_data[m_scope_id]->ev_nchan_higain = 0;
}

void VBFLaserCalc::
visitChannel(bool& veto_channel, void* user_data,
	     uint32_t                  channel_num, 
	     bool                      hit, 
	     bool                      trigger)
{
  m_data[m_scope_id]->chan[channel_num].ev_hit = hit;
}

void VBFLaserCalc::
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

  if(lo_gain)
    data->ev_nchan_logain++;
  else
    data->ev_nchan_higain++;

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
	 &&(!data->chan[channel_num].global_excluded))
	{
	  data->ev_nchan_flash++;
	  data->chan[channel_num].nflash++;
	}
      else if(signal>=m_threshold_dc)
	data->chan[channel_num].nflash++;
    }
}

void VBFLaserCalc::leaveScopeEvent(bool veto_scope_event, void* user_data)
{
  if(veto_scope_event)return;
  ScopeData* data = m_data[m_scope_id];
  data->mean_nflash_stat.accumulate(data->ev_nchan_flash);
  data->nchan_flash_hist.accumulate(data->ev_nchan_flash);
  data->nchan_logain_hist.accumulate(data->ev_nchan_logain);
  data->nchan_higain_hist.accumulate(data->ev_nchan_higain);
  if(data->ev_nchan_flash >= m_threshold_nchan) m_nflash_scope++;
}

void VBFLaserCalc::leaveArrayEvent(bool veto_array_event, void* user_data)
{
  if(veto_array_event)return;
  // --------------------------------------------------------------------------
  // Evaluate whether this event qualifies as a laser flash
  // --------------------------------------------------------------------------
  if(m_nflash_scope < m_threshold_nscope) return;
    
  const unsigned nscope = m_data.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      ScopeData* data = m_data[iscope];
      if(data == 0) continue;
      if(data->ev_nchan_flash<m_threshold_nchan) continue;

      unsigned ncrate = data->boards_per_crate.size();
      VSSimpleStat2<double> ev_signal_stat;  
      std::vector<VSSimpleStat2<double> > ev_time_stat(ncrate);

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
		  data->chan[ichan].time_stat.accumulate(ev_time);
		  if(m_fill_hist)
		    {
		      data->chan[ichan].time_hist.accumulate(ev_time);
		      data->chan[ichan].signal_hist.accumulate(ev_signal);
		    }
		}
	      data->chan[ichan].signal_stat.accumulate(ev_signal);

	      if(!data->chan[ichan].global_excluded)
		{
		  if(ev_signal>=m_threshold_dc)
		    ev_time_stat[crate].accumulate(ev_time);
		  ev_signal_stat.accumulate(ev_signal);
		}
	    }
	}
      if(ev_signal_stat.count())
	{
	  data->mean_signal_stat.accumulate(ev_signal_stat.mean());
	  data->signal_hist.accumulate(ev_signal_stat.mean());
	}
      for(unsigned icrate=0; icrate<ncrate; icrate++)
	if(ev_time_stat[icrate].count())
	  data->mean_crate_time_stat[icrate].
	    accumulate(ev_time_stat[icrate].mean());

      // ----------------------------------------------------------------------
      // Loop over all pixels and scale the laser pulse amplitude in
      // each pixel by the average for all unsuppressed pixels
      // ----------------------------------------------------------------------
      for(unsigned ichan=0;ichan<nchan; ichan++)
	{
	  double ev_signal = data->chan[ichan].ev_signal;	  
	  if(data->chan[ichan].ev_hit && ev_signal>=m_threshold_dc)
	    {
	      data->chan[ichan].
		signal_corr_stat.accumulate(ev_signal/
					    ev_signal_stat.mean());
	      data->chan[ichan].
		signal_corr_hist.accumulate(ev_signal/
					    ev_signal_stat.mean());
	    }
	}

    }
}

void VBFLaserCalc::getData(VSLaserData& data) const
{
  data.clear();
  const VSTimeDepPedData& peds = m_stage1_data.pedestals;
  // --------------------------------------------------------------------------
  // Constant to correct for the variance of the single PE amplitude
  // --------------------------------------------------------------------------
  const double singlepe_corr = 1 + pow(m_singlepe_dev,2);

  const unsigned nscope = m_data.size();
  const std::vector<unsigned> nsample = m_stage1_data.run_info.nsample;
  vsassert(m_stage1_data.run_info.nchan.size() == m_data.size());
  // --------------------------------------------------------------------------
  // Find the median/mean laser pulse amplitude, absolute gain, and #
  // PEs for each telescope
  // --------------------------------------------------------------------------
  std::vector<double> median_signal(nscope);
  std::vector<double> median_npe(nscope);
  std::vector<double> median_absgain(nscope);
  std::vector< VSSimpleStat2<double> > npe_stat(nscope);
  std::vector< VSSimpleStat2<double> > absgain_stat(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      if(m_data[iscope] == 0) continue;

      VBFLaserCalc::ScopeData* scope_data = m_data[iscope];
      unsigned nchan = scope_data->chan.size();

      std::vector< double > siglist;
      std::vector< double > npelist;
      std::vector< double > absgainlist;
      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  // ------------------------------------------------------------------
	  // If the channel has fewer than the required number of
	  // flashes set its state to 'uncalibrated'.  Exclude from
	  // the calculation of median/mean telescope parameters both
	  // 'uncalibrated' and 'excluded' channels.
	  // ------------------------------------------------------------------
	  if(scope_data->chan[ichan].nflash < m_threshold_nflash) 
	    {
	      scope_data->chan[ichan].global_uncalibrated = true;
	      continue;
	    } 
	  else if(scope_data->chan[ichan].global_excluded)
	    continue;

	  siglist.push_back(scope_data->chan[ichan].signal_stat.mean());
	  
	  double signal_corr_mean =
	    scope_data->chan[ichan].signal_corr_stat.mean()*
	    scope_data->mean_signal_stat.mean();
	  double signal_corr_dev =
	    scope_data->chan[ichan].signal_corr_stat.dev()*
	    scope_data->mean_signal_stat.mean();
	  double ped_dev = peds.dev(iscope,ichan,nsample[iscope]);
	  double laser_var = pow(signal_corr_dev,2) - pow(ped_dev,2);
	  double npe = singlepe_corr*pow(signal_corr_mean,2)/laser_var;
	  double absgain = laser_var/(singlepe_corr*signal_corr_mean);
	      
	  npelist.push_back(npe);
	  absgainlist.push_back(absgain);
	  npe_stat[iscope].accumulate(npe);
	  absgain_stat[iscope].accumulate(absgain);	  
	}
      
      median_signal[iscope] = median(siglist);
      median_npe[iscope] = median(npelist);
      median_absgain[iscope] = median(absgainlist);
    }
  
  // --------------------------------------------------------------------------
  // Fill run data
  // --------------------------------------------------------------------------
  data.m_runno = m_stage1_data.run_info.run_number;
  data.m_threshold_nchan = m_threshold_nchan;
  data.m_threshold_dc = m_threshold_dc;
  data.m_singlepe_dev = m_singlepe_dev;
  
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      if(m_data[iscope] == 0)
	continue;

      VBFLaserCalc::ScopeData* scope_data = m_data[iscope];
      // ----------------------------------------------------------------------
      // Skip filling this telescope if there are too few laser flashes
      // ----------------------------------------------------------------------
      if(scope_data->mean_signal_stat.count()<m_threshold_nflash)
	{
	  std::cerr 
	    << "T" << iscope+1 << " has only " 
	    << scope_data->mean_signal_stat.count() 
	    << " flashes. Need "
	    << m_threshold_nflash << "." << std::endl
	    << " Num Chan Mean:    " 
	    << scope_data->nchan_flash_hist.mean() << std::endl
	    << " Num LoGain Mean:  " 
	    << scope_data->nchan_logain_hist.mean() << std::endl
	    << " Num HiGain Mean:  " 
	    << scope_data->nchan_higain_hist.mean() << std::endl
	    << " Amplitude Mean:   " 
	    << scope_data->mean_signal_stat.mean() << std::endl 
	    << " Amplitude Dev:    " 
	    << scope_data->mean_signal_stat.dev() << std::endl 
	    << std::endl;
	  continue;
	}

      if(iscope >= data.scope.size())data.scope.resize(iscope+1);

      // ----------------------------------------------------------------------
      // Fill telescope data
      // ----------------------------------------------------------------------
      unsigned nchan = scope_data->chan.size();
      data.scope[iscope].runno = m_stage1_data.run_info.run_number;
      data.scope[iscope].nchan = nchan;
      data.scope[iscope].nevents = scope_data->mean_nflash_stat.count();
      data.scope[iscope].nflash = scope_data->mean_signal_stat.count();
      data.scope[iscope].absgain_mean = absgain_stat[iscope].mean();
      data.scope[iscope].absgain_median = median_absgain[iscope];
      data.scope[iscope].npe_mean = npe_stat[iscope].mean();
      data.scope[iscope].npe_median = median_npe[iscope];
      data.scope[iscope].nchan_mean = scope_data->mean_nflash_stat.mean();
      data.scope[iscope].nchan_dev = scope_data->mean_nflash_stat.dev();
      data.scope[iscope].signal_mean = scope_data->mean_signal_stat.mean();
      data.scope[iscope].signal_dev = scope_data->mean_signal_stat.dev();
      data.scope[iscope].signal_hist = scope_data->signal_hist;
      data.scope[iscope].nchan_flash_hist = scope_data->nchan_flash_hist;
      data.scope[iscope].nchan_logain_hist = scope_data->nchan_logain_hist;
      data.scope[iscope].nchan_higain_hist = scope_data->nchan_higain_hist;

      // ----------------------------------------------------------------------
      // Fill crate data
      // ----------------------------------------------------------------------
      unsigned ncrate = scope_data->l2_pulse_channel.size();
      data.scope[iscope].crate.resize(ncrate);
      for(unsigned icrate=0; icrate<ncrate; icrate++)
	{
	  unsigned l2chan = nchan+1;
	  if(icrate < scope_data->l2_pulse_channel.size())
	    l2chan = scope_data->l2_pulse_channel[icrate];
	  double l2time = ((l2chan<nchan)?
			   scope_data->chan[l2chan].time_stat.mean():0.0);

	  data.scope[iscope].crate[icrate].l2chan = l2chan;
	  data.scope[iscope].crate[icrate].l2time = l2time;
	  data.scope[iscope].crate[icrate].cratetime_mean =	  
	    scope_data->mean_crate_time_stat[icrate].mean()-l2time;
	  data.scope[iscope].crate[icrate].cratetime_dev =	  
	    scope_data->mean_crate_time_stat[icrate].dev();
	}

      // ----------------------------------------------------------------------
      // Fill channel data
      // ----------------------------------------------------------------------
      for(unsigned ichan=0; ichan<nchan; ichan++)
	{	 
	  if(ichan >= data.scope[iscope].chan.size())
	    data.scope[iscope].chan.resize(ichan+1);

	  ChanData& chan_data = scope_data->chan[ichan];
	  VSLaserData::ChanData& chan = data.scope[iscope].chan[ichan];

	  unsigned icrate = chan_data.global_crate;
	  unsigned l2chan = data.scope[iscope].crate[icrate].l2chan;
	  double l2time = ((l2chan<nchan)?
			   scope_data->chan[l2chan].time_stat.mean():0.0);

	  double signal_corr_mean = chan_data.signal_corr_stat.mean()*
	    scope_data->mean_signal_stat.mean();
	  double signal_corr_mean_err = chan_data.signal_corr_stat.dev()*
	    scope_data->mean_signal_stat.mean()/
	    sqrt(chan_data.signal_corr_stat.count());
	  double signal_corr_dev = chan_data.signal_corr_stat.dev()*
	    scope_data->mean_signal_stat.mean();
	  double signal_corr_dev_err = chan_data.signal_corr_stat.dev()*
	    scope_data->mean_signal_stat.mean()/
	    sqrt(2*chan_data.signal_corr_stat.count());

	  chan.nflash = chan_data.signal_stat.count();
	  chan.l2chan = l2chan;
	  chan.crate = icrate;
	  chan.uncalibrated = chan_data.global_uncalibrated;
	  chan.excluded = chan_data.global_excluded;

	  if(!chan_data.global_uncalibrated)
	    {
	      chan.chantime = chan_data.time_stat.mean() - l2time;
	      chan.cratetime =
		scope_data->mean_crate_time_stat[icrate].mean() - l2time;
	      chan.l2time = l2time;
	      chan.gain = chan_data.signal_stat.mean()/median_signal[iscope];
	      chan.gain_err = chan_data.signal_stat.dev()/
		(chan_data.signal_stat.count()*median_signal[iscope]);
	      chan.signal_mean = chan_data.signal_stat.mean();
	      chan.signal_dev = chan_data.signal_stat.dev();
	      chan.signal_corr_mean = signal_corr_mean;
	      chan.signal_corr_dev = signal_corr_dev;

	      double ped_dev = peds.dev(iscope,ichan,nsample[iscope]);
	      double laser_var = std::pow(signal_corr_dev,2) - pow(ped_dev,2);
	      double laser_var_err = 
		std::pow(signal_corr_dev,2)*
		2*signal_corr_dev_err/signal_corr_dev;

	      chan.absgain = laser_var/(singlepe_corr*signal_corr_mean);
	      chan.absgain_err = 
		chan.absgain*
		sqrt(std::pow(laser_var_err/laser_var,2) +
		     std::pow(signal_corr_mean_err/signal_corr_mean,2));
	      
	      chan.eff = singlepe_corr*pow(signal_corr_mean,2)/
		(laser_var*median_npe[iscope]);
	    }
	  else
	    {
	      chan.chantime = 0;
	      chan.cratetime = 0;
	      chan.l2time = 0;
	      chan.gain = 0;
	      chan.gain_err = 0;
	      chan.signal_mean = 0;
	      chan.signal_dev = 0;
	      chan.signal_corr_mean = 0;
	      chan.signal_corr_dev = 0;
	      chan.absgain = 0;
	      chan.absgain_err = 0;
	      chan.eff = 0;
	      chan.eff_err = 0;
	    }

	  chan.signal_hist = chan_data.signal_hist;
	}
    }
}
