//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFTimeDepNSB.cpp

  Calculate time dependent NSB, can be used to flag suppressed channels or
  estimate tracking performance

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VBFTimeDepNSB.cpp,v 3.7 2008/02/29 21:16:34 matthew Exp $

*/

#include <vsassert>

#include <VBFTimeDepNSB.hpp>

using namespace VERITAS;

VBFTimeDepNSB::
VBFTimeDepNSB(unsigned window_start, unsigned window_width,
	      unsigned events_per_slice, unsigned nevent_cut)
  : VSSimpleVBFVisitor(),
    m_window_start(window_start), m_window_width(window_width),
    m_events_per_slice(events_per_slice), m_nevent_cut(nevent_cut),
    m_slice(), m_nchan(), m_ped(),
    m_is_ped_event(false), m_event_num(), m_scope_num(), 
    m_cached_slice(), m_cached_scope(), m_got_event(false),
    m_lo_time(), m_hi_time()
{
  // nothing to see here
}
    
VBFTimeDepNSB::~VBFTimeDepNSB()
{
  for(std::vector<Slice>::iterator islice=m_slice.begin();
      islice!=m_slice.end(); islice++)
    delete islice->scope;
}
    
void VBFTimeDepNSB::
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
  if(!has_event_number)
    {
      veto_array_event = true;
      return;
    }

  if((!m_got_event)||(best_event_time>m_hi_time))
    m_hi_time=best_event_time, m_got_event=true;
  if((!m_got_event)||(best_event_time<m_lo_time))
    m_lo_time=best_event_time, m_got_event=true;
  
  m_event_num = event_num;
  unsigned islice = m_event_num/m_events_per_slice;
  if(islice >= m_slice.size())m_slice.resize(islice+1);
  if(!m_slice[islice].scope)
    {
      m_slice[islice].scope=new Array;
      resize(m_slice[islice].scope);
    }

  if((has_good_event_time)&&(best_event_time.isGood())
     &&((!m_slice[islice].slice_time.isGood())
	||(best_event_time<m_slice[islice].slice_time)))
    m_slice[islice].slice_time = best_event_time;

  m_cached_slice=m_slice[islice].scope;
}

void VBFTimeDepNSB::
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
  if((m_window_width==0)||(m_window_start+m_window_width>num_samples))
    {
      veto_scope_event=true;
      return;
    }

  m_is_ped_event = (event_type.trigger == VEventType::PED_TRIGGER);
  m_scope_num = telescope_num;
  if(m_scope_num >= m_nchan.size())
    {
      m_nchan.resize(m_scope_num+1);
      for(std::vector<Slice>::iterator islice=m_slice.begin();
	  islice!=m_slice.end(); islice++)
       if(islice->scope)resize(islice->scope);
    }
  if(num_channels_total > m_nchan[m_scope_num])
    {
      m_nchan[m_scope_num]=num_channels_total;
      for(std::vector<Slice>::iterator islice=m_slice.begin();
	  islice!=m_slice.end(); islice++)
	if(islice->scope)resize(islice->scope);
    }
  m_cached_scope=&(*m_cached_slice)[m_scope_num];
}

void VBFTimeDepNSB::
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
  if(m_is_ped_event)
    {
      double q = double(integrated[nsample-1])/double(nsample);
      m_ped[m_scope_num][channel_num].accumulate(q,nsample);
    };
  unsigned Q = integrated[m_window_start+m_window_width-1];
  if(m_window_start)Q -= integrated[m_window_start-1];
  (*m_cached_scope)[channel_num].accumulate(double(Q));
}

void VBFTimeDepNSB::
suppressFromPedDev(VSTimeDepSupData& sup, 
		   const std::vector<double>& sup_dev) const
{
  const unsigned nscope = m_nchan.size();
  vsassert(sup_dev.size() >= nscope);

  const unsigned nslice = m_slice.size();
  sup.resize(m_events_per_slice, nslice, m_nchan);
  for(unsigned islice=0;islice<nslice;islice++)
    if(m_slice[islice].scope)
      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  const double cut_dev = sup_dev[iscope]*sup_dev[iscope];
	  for(unsigned ichan=0;ichan<m_nchan[iscope];ichan++)
	    {
	      double ped = m_ped[iscope][ichan].mean()*double(m_window_width);
	      if(((*m_slice[islice].scope)[iscope][ichan].chi2(ped) < cut_dev)
		 &&((*m_slice[islice].scope)[iscope][ichan].count()
		    >= m_nevent_cut))
		sup.setSuppressed(islice, iscope, ichan);
	    }
	}
}

void VBFTimeDepNSB::getData(VBFTimeDepNSB::Data& data) const
{
  data.clear();

  const unsigned nscope = m_nchan.size();
  data.scope.resize(nscope);

  const unsigned nslice = m_slice.size();
  data.slice_time.resize(nslice);
  for(unsigned islice=0;islice<nslice;islice++)
    data.slice_time[islice] = m_slice[islice].slice_time;

  for(unsigned iscope=0;iscope<m_nchan.size();iscope++)
    if(m_nchan[iscope])
      {
	std::vector<double> scope_dev;
	scope_dev.reserve(m_nchan[iscope]);

	data.scope[iscope].slice_median_dev.resize(nslice);

	for(unsigned islice=0;islice<nslice;islice++)
	  if(m_slice[islice].scope)
	    {
	      for(unsigned ichan=0;ichan<m_nchan[iscope];ichan++)
		{
		  const double ped =
		    m_ped[iscope][ichan].mean()*double(m_window_width);
		  const double dev =
		    sqrt((*m_slice[islice].scope)[iscope][ichan].chi2(ped));
		  if((*m_slice[islice].scope)[iscope][ichan].count()
		     >= m_nevent_cut)
		    {
		      data.scope[iscope].all_dev.accumulate(dev);
		      scope_dev.push_back(dev);
		    }
		}
	      if(!scope_dev.empty())
		data.scope[iscope].slice_median_dev[islice] = 
		  median(scope_dev);
	    }

	data.scope[iscope].suppress_dev_suggested = 0.0;

	const int zbin = data.scope[iscope].all_dev.firstBin();
	const int nbin = data.scope[iscope].all_dev.onePastLastBin(); 
	if(nbin-zbin > 5)
	  {
	    int pbin=zbin;
	    unsigned pcount=0;
	    for(int ibin=zbin;ibin<nbin;ibin++)
	      if(pcount < data.scope[iscope].all_dev.count(ibin))
		pcount=data.scope[iscope].all_dev.count(ibin), pbin=ibin;

	    while((pbin>zbin)
		  &&(data.scope[iscope].all_dev.count(pbin)*2 > pcount))pbin--;

	    while(pbin>zbin)
	      {
		unsigned count = data.scope[iscope].all_dev.count(pbin);
#if 0
		std::cout << pbin << ' ' 
			  << data.scope[iscope].all_dev.val(pbin)
			  << ' ' << count << ' '
			  << data.scope[iscope].all_dev.count(pbin-2) << ' ' 
			  << data.scope[iscope].all_dev.count(pbin-1) << ' ' 
			  << data.scope[iscope].all_dev.count(pbin+1) << ' ' 
			  << data.scope[iscope].all_dev.count(pbin+2) << '\n';
#endif
		if((count <= data.scope[iscope].all_dev.count(pbin-2))
		   &&(count <= data.scope[iscope].all_dev.count(pbin-1))
		   &&(count <= data.scope[iscope].all_dev.count(pbin+1))
		   &&(count <= data.scope[iscope].all_dev.count(pbin+2)))
		  break;

		pbin--;
	      }

	    data.scope[iscope].suppress_dev_suggested = 
	      data.scope[iscope].all_dev.val(pbin);
	  }
	else
	  {
	    data.scope[iscope].suppress_dev_suggested = 
	      data.scope[iscope].all_dev.val(0);
	  }
      }
}

void VBFTimeDepNSB::resize(VBFTimeDepNSB::Array* slice)
{
  const unsigned nscope = m_nchan.size();
  if(nscope >= slice->size())slice->resize(nscope);
  if(nscope >= m_ped.size())m_ped.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nchan = m_nchan[iscope];
      if(nchan >= (*slice)[iscope].size())(*slice)[iscope].resize(nchan);
      if(nchan >= m_ped[iscope].size())m_ped[iscope].resize(nchan);
    }
}

void VBFTimeDepNSB::Data::clear()
{
  scope.clear();
}

void VBFTimeDepNSB::Data::load(VSOctaveH5ReaderStruct* reader)
{
  clear();
  if(!reader->isCellVector("scope"))
    {
      const unsigned nscope = 4;
      scope.resize(nscope);
      for(unsigned iscope=0;iscope<scope.size();iscope++)
	{
	  scope[iscope].all_dev.load(reader->readStruct("all_dev"));
	  reader->readScalar("suppress_dev_suggested",
			     scope[iscope].suppress_dev_suggested);
	  reader->readScalar("suppress_dev",scope[iscope].suppress_dev);
	}
    }
  else
    {
      if(reader->isValid("slice_time_mjd"))
	{
	  std::vector<uint32_t> t_mjd;
	  std::vector<uint64_t> t_dns;
	  reader->readVector("slice_time_mjd",t_mjd);
	  reader->readVector("slice_time_dns",t_dns);
	  unsigned nslice = t_mjd.size();
	  slice_time.resize(nslice);
	  for(unsigned islice=0;islice<nslice;islice++)
	    slice_time[islice].
	      setFromMJDIntAndNS(t_mjd[islice],t_dns[islice]);
	}

      VSOctaveH5ReaderCellVector* c = reader->readCellVector("scope");
      vsassert(c);
      unsigned nscope = c->dimensions();
      scope.resize(nscope);
      for(unsigned iscope=0;iscope<scope.size();iscope++)
	if(!c->isEmpty(iscope))
	  {
	    VSOctaveH5ReaderStruct* s = c->readStruct(iscope);
	    scope[iscope].all_dev.load(s->readStruct("all_dev"));
	    s->readScalar("suppress_dev_suggested",
			  scope[iscope].suppress_dev_suggested);
	    s->readScalar("suppress_dev",scope[iscope].suppress_dev);
	    if(s->isValid("slice_median_dev"))
	      s->readVector("slice_median_dev",scope[iscope].slice_median_dev);
	    delete s;
	  }
      delete c;
    }
}

void VBFTimeDepNSB::Data::save(VSOctaveH5WriterStruct* writer) const
{
  unsigned nslice = slice_time.size();
  if(nslice)
    {
      std::vector<uint32_t> t_mjd(nslice);
      std::vector<uint64_t> t_dns(nslice);
      for(unsigned islice=0;islice<nslice;islice++)
	t_mjd[islice] = slice_time[islice].getMJDInt(),
	  t_dns[islice] = slice_time[islice].getDayNS();
      writer->writeVector("slice_time_mjd",t_mjd);
      writer->writeVector("slice_time_dns",t_dns);
    }

  const unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* c = writer->writeCellVector("scope",nscope);
  for(unsigned iscope=0;iscope<scope.size();iscope++)
    if(!scope[iscope].all_dev.empty())
      {
	VSOctaveH5WriterStruct* s = c->writeStruct(iscope);
	scope[iscope].all_dev.save(s->writeStruct("all_dev"));
	s->writeScalar("suppress_dev_suggested",
		       scope[iscope].suppress_dev_suggested);
	s->writeScalar("suppress_dev",scope[iscope].suppress_dev);
	if(nslice)
	  s->writeVector("slice_median_dev",scope[iscope].slice_median_dev);
	delete s;
      }
  delete c;
}

