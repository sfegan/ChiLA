//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTimeDepPedData.hpp
  TimeDep pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/10/2006

  $Id: VSTimeDepPedData.cpp,v 3.5 2008/02/29 21:16:34 matthew Exp $

*/

#include<VSSimpleStat.hpp>
#include<VSTimeDepPedData.hpp>

using namespace VERITAS;

void VSTimeDepPedData::clear()
{
  m_window_min = 0;
  m_hi_ped.clear();
  m_lo_ped.clear();
  m_hi_dev.clear();
  m_min_window_hist.clear();
  m_slice.clear();
}

double VSTimeDepPedData::
medianDev(unsigned iscope, unsigned ichan, unsigned iwindow,
	  unsigned min_event_count) const
{
  const unsigned nslice = m_slice.size();
  std::vector<std::pair<bool, double> > alldev;
  unsigned ngood=0;
  for(unsigned islice=0;islice<nslice;islice++)
    {
      std::pair<bool, double> onedev;
      onedev.first = (min_event_count==0)
	||(eventCount(islice, iscope, ichan)>min_event_count);
      onedev.second = dev(islice, iscope, ichan, iwindow);
      alldev.push_back(onedev);
      if(onedev.first)ngood++;
    }
  if(ngood)return median(alldev);
  else return 0;
}

void VSTimeDepPedData::load(VSOctaveH5ReaderStruct* reader)
{
  reader->readScalar("window_min",m_window_min);

  VSOctaveH5ReaderCellVector* hiped = reader->readCellVector("hi_ped");
  vsassert(hiped);
  unsigned hi_nscope = hiped->dimensions();
  m_hi_ped.resize(hi_nscope);
  for(unsigned iscope=0;iscope<hi_nscope;iscope++)
    if(!hiped->isEmpty(iscope))hiped->readVector(iscope,m_hi_ped[iscope]);
  delete hiped;

  VSOctaveH5ReaderCellVector* loped = reader->readCellVector("lo_ped");
  vsassert(loped);
  unsigned lo_nscope = loped->dimensions();
  m_lo_ped.resize(lo_nscope);
  for(unsigned iscope=0;iscope<lo_nscope;iscope++)
    if(!loped->isEmpty(iscope))loped->readVector(iscope,m_lo_ped[iscope]);
  delete loped;

  VSOctaveH5ReaderCellVector* hidev = reader->readCellVector("hi_dev");
  vsassert(hidev);
  hi_nscope = hidev->dimensions();
  m_hi_dev.resize(hi_nscope);
  for(unsigned iscope=0;iscope<hi_nscope;iscope++)
    if(!hidev->isEmpty(iscope))hidev->readVector(iscope,m_hi_dev[iscope]);
  delete hidev;

  if(reader->isValid("min_window_hist"))
    {
      VSOctaveH5ReaderCellVector* mwsc_hist = 
	reader->readCellVector("min_window_hist");
      vsassert(mwsc_hist);
      unsigned mw_nscope = mwsc_hist->dimensions();
      m_min_window_hist.resize(mw_nscope);
      for(unsigned iscope=0;iscope<mw_nscope;iscope++)
	if(!mwsc_hist->isEmpty(iscope))
	  {
	    VSOctaveH5ReaderCellVector* mwch_hist = 
	      mwsc_hist->readCellVector(iscope);
	    unsigned mw_nchan = mwch_hist->dimensions();
	    m_min_window_hist[iscope].resize(mw_nchan,1);
	    for(unsigned ichan=0;ichan<mw_nchan;ichan++)
	      {
		VSOctaveH5ReaderStruct* s = mwch_hist->readStruct(ichan);
		m_min_window_hist[iscope][ichan].load(s);
		delete s;
	      }
	      
	    delete mwch_hist;
	  }
      delete mwsc_hist;
    }

  VSOctaveH5ReaderCellVector* slice = reader->readCellVector("slice");
  vsassert(slice);
  unsigned nslice = slice->dimensions();
  m_slice.resize(nslice);
  for(unsigned islice=0;islice<nslice;islice++)
    {
      VSOctaveH5ReaderStruct* s = slice->readStruct(islice);
      s->readScalar("event_num_lo",m_slice[islice].event_num_lo);
      s->readScalar("event_num_hi",m_slice[islice].event_num_hi);
      H5READ_VATIME(s,"event_time_lo",m_slice[islice].event_time_lo);
      H5READ_VATIME(s,"event_time_hi",m_slice[islice].event_time_hi);
      s->readScalar("ped_event_count",m_slice[islice].ped_event_count);

      VSOctaveH5ReaderCellVector* ccnt =
	s->readCellVector("channel_ped_event_count");
      vsassert(ccnt);

      VSOctaveH5ReaderCellVector* cdev =	
	s->readCellVector("channel_dev");
      vsassert(cdev);

      unsigned slice_nscope = cdev->dimensions();
      m_slice[islice].scope.resize(slice_nscope);
      for(unsigned iscope=0;iscope<slice_nscope;iscope++)
	if(!cdev->isEmpty(iscope))
	  {
	    unsigned nchan = 0;
	    unsigned nwindow = 0;
	    cdev->dimensions(iscope,nwindow,nchan);	  
	    m_slice[islice].scope[iscope].resize(nchan);

	    std::vector<unsigned> _ped_event_count;
	    ccnt->readVector(iscope,_ped_event_count);
	    for(unsigned ichan=0;ichan<nchan;ichan++)
	      m_slice[islice].scope[iscope][ichan].ped_event_count =
		_ped_event_count[ichan];
	    
	    if(nwindow)
	      {
		double* data = new double[nchan*nwindow];
		std::vector<unsigned> _ped_event_count;
		cdev->readMatrix(iscope,data);
		for(unsigned ichan=0;ichan<nchan;ichan++)
		  {
		    m_slice[islice].scope[iscope][ichan].
		      window_dev.resize(nwindow);
		    for(unsigned iwindow=0;iwindow<nwindow;iwindow++)
		      m_slice[islice].scope[iscope][ichan].
			window_dev[iwindow] = data[iwindow*nchan+ichan];
		  }
		delete[] data;
	      }
	  }
      delete cdev;
      delete ccnt;
      delete s;
    }
  delete slice;
}

void VSTimeDepPedData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeScalar("window_min",m_window_min);

  const unsigned hi_nscope = m_hi_ped.size();
  VSOctaveH5WriterCellVector* hiped = 
    writer->writeCellVector("hi_ped",hi_nscope);
  for(unsigned iscope=0;iscope<hi_nscope;iscope++)
    if(m_hi_ped[iscope].size())
      hiped->writeVector(iscope,m_hi_ped[iscope]);
  delete hiped;

  const unsigned lo_nscope = m_lo_ped.size();
  VSOctaveH5WriterCellVector* loped = 
    writer->writeCellVector("lo_ped",lo_nscope);
  for(unsigned iscope=0;iscope<lo_nscope;iscope++)
    if(m_lo_ped[iscope].size())
      loped->writeVector(iscope,m_lo_ped[iscope]);
  delete loped;

  const unsigned h2_nscope = m_hi_dev.size();
  VSOctaveH5WriterCellVector* hidev = 
    writer->writeCellVector("hi_dev",h2_nscope);
  for(unsigned iscope=0;iscope<h2_nscope;iscope++)
    if(m_hi_dev[iscope].size())
      hidev->writeVector(iscope,m_hi_dev[iscope]);
  delete hidev;

  unsigned mw_nscope = m_min_window_hist.size();
  VSOctaveH5WriterCellVector* mwsc_hist = 
    writer->writeCellVector("min_window_hist",mw_nscope);
  for(unsigned iscope=0;iscope<mw_nscope;iscope++)
    if(m_min_window_hist[iscope].size())
      {
	unsigned mw_nchan = m_min_window_hist[iscope].size();
	VSOctaveH5WriterCellVector* mwch_hist = 
	  mwsc_hist->writeCellVector(iscope,mw_nchan);
	for(unsigned ichan=0;ichan<mw_nchan;ichan++)
	  {
	    VSOctaveH5WriterStruct* s = mwch_hist->writeStruct(ichan);
	    m_min_window_hist[iscope][ichan].save(s);
	    delete s;
	  }
	delete mwch_hist;
      }
  delete mwsc_hist;

  const unsigned nslice = m_slice.size();
  VSOctaveH5WriterCellVector* slice = writer->writeCellVector("slice",nslice);
  for(unsigned islice=0;islice<nslice;islice++)
    {
      VSOctaveH5WriterStruct* s = slice->writeStruct(islice);
      s->writeScalar("event_num_lo",m_slice[islice].event_num_lo);
      s->writeScalar("event_num_hi",m_slice[islice].event_num_hi);
      H5WRITE_VATIME(s,"event_time_lo",m_slice[islice].event_time_lo);
      H5WRITE_VATIME(s,"event_time_hi",m_slice[islice].event_time_hi);
      s->writeScalar("ped_event_count",m_slice[islice].ped_event_count);
      const unsigned slice_nscope = m_slice[islice].scope.size();
      VSOctaveH5WriterCellVector* ccnt = 
	s->writeCellVector("channel_ped_event_count",slice_nscope);
      VSOctaveH5WriterCellVector* cdev = 
	s->writeCellVector("channel_dev",slice_nscope);
      for(unsigned iscope=0;iscope<slice_nscope;iscope++)
	{
	  unsigned nchan = m_slice[islice].scope[iscope].size();
	  if(nchan)
	    {
	      unsigned nwindow = 
		m_slice[islice].scope[iscope][0].window_dev.size();
	      std::vector<unsigned> _ped_event_count(nchan);
	      double* data = new double[nchan*nwindow];
	      for(unsigned ichan=0;ichan<nchan;ichan++)
		{
		  _ped_event_count[ichan] = 
		    m_slice[islice].scope[iscope][ichan].ped_event_count;
		  for(unsigned iwindow=0;iwindow<nwindow;iwindow++)
		    data[iwindow*nchan+ichan] = 
		      m_slice[islice].scope[iscope][ichan].window_dev[iwindow];
		}
	      ccnt->writeVector(iscope,_ped_event_count);
	      cdev->writeMatrix(iscope,nwindow,nchan,data);
	      delete[] data;
	    }
	}
      delete cdev;
      delete ccnt;
      delete s;
    }
  delete slice;
}

void VSTimeDepPedData::
suppressScopeByDev(const unsigned iscope, 
		   const double cut_lo, const double cut_hi,
		   VSTimeDepSupData& sup, unsigned min_event_count)
{
  const unsigned nslice = m_slice.size();
  for(unsigned islice=0;islice<nslice;islice++)
    {
      const unsigned ev_slice_lo = 
	sup.getSliceByEventNum(sliceEventNumLo(islice));
      const unsigned ev_slice_hi = 
	sup.getSliceByEventNum(sliceEventNumHi(islice)-1)+1;
      const unsigned nchan = m_hi_ped[iscope].size();
      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  const unsigned iwindow = max_window(islice, iscope, ichan);
	  const unsigned count = eventCount(islice,iscope,ichan);
	  const double onedev = dev(islice, iscope, ichan, iwindow);
	  if(((min_event_count!=0)&&(count<=min_event_count))
	     ||(onedev<cut_lo)||(onedev>cut_hi))
	    for(unsigned ev_islice=ev_slice_lo;ev_islice<ev_slice_hi;
		ev_islice++)
	      sup.setSuppressed(ev_islice,iscope,ichan);
	}
    }
}

void VSTimeDepPedData::
suppressByMedianFraction(const double lo, const double hi,
			 VSTimeDepSupData& sup,
			 std::vector<std::pair<double,double> >& cuts, 
			 unsigned min_event_count)
{
  const unsigned nscope = m_hi_ped.size();
  cuts.clear();
  cuts.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      std::vector<std::pair<bool, double> > alldev;
      getAllDev(iscope,alldev,min_event_count);

      const double median_dev = median(alldev);
      const double cut_lo = median_dev*lo;
      const double cut_hi = median_dev*hi;

      cuts[iscope].first = cut_lo;
      cuts[iscope].second = cut_hi;

      suppressScopeByDev(iscope, cut_lo, cut_hi, sup, min_event_count);
    }
}

void VSTimeDepPedData::
suppressByContainmentFraction(const double fraction, const double scale,
			      VSTimeDepSupData& sup,
			      std::vector<std::pair<double,double> >& cuts, 
			      unsigned min_event_count)
{
  const unsigned nscope = m_hi_ped.size();
  cuts.clear();
  cuts.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      std::vector<std::pair<bool, double> > alldev;
      getAllDev(iscope,alldev,min_event_count);

      double cut_lo;
      double cut_hi;
      containmentInterval(alldev,fraction,scale,cut_lo,cut_hi);

      cuts[iscope].first = cut_lo;
      cuts[iscope].second = cut_hi;

      suppressScopeByDev(iscope, cut_lo, cut_hi, sup, min_event_count);
    }
}

void VSTimeDepPedData::
getAllDev(unsigned iscope, std::vector<std::pair<bool, double> >& alldev,
	  unsigned min_event_count) const
{
  alldev.clear();
  const unsigned nchan = m_hi_ped[iscope].size();
  const unsigned nslice = m_slice.size();
  for(unsigned islice=0;islice<nslice;islice++)
    for(unsigned ichan=0;ichan<nchan;ichan++)
      {
	const unsigned iwindow = max_window(islice, iscope, ichan);
	std::pair<bool, double> onedev;
	onedev.first = (min_event_count==0)
	  ||(eventCount(islice, iscope, ichan)>min_event_count);
	onedev.second = dev(islice, iscope, ichan, iwindow);
	alldev.push_back(onedev);
      }
}
