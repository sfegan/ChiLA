//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFHiLoCalc.cpp
    Class designed for use with logain.cpp...

  \author     Timothy C. Arlen            \n
              UCLA                        \n
	      arlen@astro.ucla.edu        \n

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  
  \date       07/23/2008

  $Id: VBFHiLoCalc.cpp,v 3.1 2009/12/22 19:43:30 matthew Exp $

*/


// NOTE: Commit changes to this and Makefile and comment top when done!!


#include "VBFHiLoCalc.hpp"

using namespace VERITAS;

// ============================================================================
// VBFHiLoCalc
// ============================================================================

VBFHiLoCalc::VBFHiLoCalc(int sample_0, unsigned sample_N):
  scope(), m_sample_0(sample_0), m_sample_N(sample_N),
  m_runno(), m_scope_id()
{
  // nothing to see here
}

VBFHiLoCalc::~VBFHiLoCalc()
{
  // nothing to see here
}

void VBFHiLoCalc::
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
  m_runno = run_number;
}

void VBFHiLoCalc::
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
  if(m_scope_id >= scope.size())scope.resize(m_scope_id+1);
  if(scope[m_scope_id]==0)scope[m_scope_id]=new ScopeData(num_channels_total);
  if(num_channels_total > scope[m_scope_id]->chan.size())
    scope[m_scope_id]->chan.resize(num_channels_total);
  scope[m_scope_id]->nevent++;
}
    
void VBFHiLoCalc::
visitChannel(bool& veto_channel, void* user_data,
	     uint32_t                  channel_num, 
	     bool                      hit, 
	     bool                      trigger)
{

}
    
void VBFHiLoCalc::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		unsigned               nsample,
		const uint32_t*        samples,
		const uint32_t*        integrated)
{
  if(lo_gain)
    {
      unsigned sample0 = 0;
      unsigned sampleN = m_sample_N;
      if(m_sample_0 > 0)sample0 = unsigned(m_sample_0);
      else if(m_sample_0 < 0 && unsigned(abs(m_sample_0)) < nsample)
	sample0 = nsample - unsigned(abs(m_sample_0));
      if(sampleN > nsample)sampleN=nsample;
      if(sample0 + sampleN > nsample)sample0 = nsample-sampleN;
      if(sampleN == 0)sampleN = nsample-sample0;

      if(sampleN)
	{
	  ChanData& cd(scope[m_scope_id]->chan[channel_num]);
	  unsigned sum = integrated[sample0+sampleN-1];
	  if(sample0)sum -= integrated[sample0-1];
	  cd.ped_stat.accumulate(double(sum)/double(sampleN),sampleN);
	  cd.nevent++;
	}
    }
}
    
void VBFHiLoCalc::leaveScopeEvent(bool veto_scope_event, void* user_data)
{

}

void VBFHiLoCalc::leaveArrayEvent(bool veto_array_event, void* user_data)
{

}
    
void VBFHiLoCalc::getData(VSHiLoData& data, unsigned nevent_min) const
{
  const unsigned nscope = scope.size();

  data.clear();
  data.runno                = m_runno;
  data.has_hi_lo_gain_ratio = false;
  data.has_lo_gain_ped      = true;
  data.has_lo_gain_switch   = false;
  data.scope                .resize(nscope);

  for(unsigned iscope=0; iscope<nscope; iscope++)
    if(scope[iscope] && scope[iscope]->nevent>=nevent_min)
      {
	const unsigned nchan = scope[iscope]->chan.size();
	//const unsigned nevent = scope[iscope]->nevent;
	data.scope[iscope].chan.resize(nchan);
	for(unsigned ichan=0; ichan<nchan; ichan++)
	  {
	    const ChanData& icd(scope[iscope]->chan[ichan]);
	    VSHiLoData::ChanData& ocd(data.scope[iscope].chan[ichan]);

	    ocd.hi_lo_gain_ratio     = 6.0;
	    ocd.lo_gain_ped          = 0.0;
	    ocd.lo_gain_switch_amp   = 0.0;
	    ocd.has_hi_lo_gain_ratio = false;
	    ocd.has_lo_gain_ped      = false;
	    ocd.has_lo_gain_switch   = false;

	    if(icd.nevent>=nevent_min)
	      {
		ocd.lo_gain_ped      = icd.ped_stat.mean();
		ocd.has_lo_gain_ped  = true;
	      }
	  }
      }
}
