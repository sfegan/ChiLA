//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTimeDepPedData.hpp
  TimeDep pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/10/2006

  $Id: VSTimeDepPedData.hpp,v 3.5 2008/02/29 21:16:34 matthew Exp $

*/

#ifndef VSTIMEDEPPEDDATA_HPP
#define VSTIMEDEPPEDDATA_HPP

#include<vector>
#include<map>

#include<VSOctaveIO.hpp>
#include<VSTime.hpp>
#include<VSTimeDepSupData.hpp>
#include <VSSimpleHist.hpp>

namespace VERITAS
{

  // ==========================================================================
  // TIME DEPENDENT PEDESTAL DATA
  // ==========================================================================

  class VBFTimeDepPeds;
  
  class VSTimeDepPedData
  {
  public:
    VSTimeDepPedData(): 
      m_window_min(), m_hi_ped(), m_lo_ped(), m_hi_dev(), 
      m_min_window_hist(0,std::vector<VSSimpleHist<int32_t> >(0,1)), 
      m_slice() { }

    void clear();
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    unsigned nslice() const { return m_slice.size(); }
    unsigned nscope() const { return m_hi_ped.size(); }
    unsigned nchan(unsigned iscope) const { return m_hi_ped[iscope].size(); }

    unsigned max_window(unsigned islice, unsigned iscope, 
			unsigned ichan) const
    {
      return m_slice[islice].scope[iscope][ichan].
	window_dev.size()+m_window_min-1;
    }
    unsigned min_window() const { return m_window_min; }

    double hiPed(unsigned iscope, unsigned ichan) const
    { return m_hi_ped[iscope][ichan]; }
    double loPed(unsigned iscope, unsigned ichan) const
    { return m_lo_ped[iscope][ichan]; }
    double hiDev(unsigned iscope, unsigned ichan) const
    { return dev(iscope,ichan,max_window(0,iscope,ichan)); }

    double dev(unsigned islice, unsigned iscope, 
	       unsigned ichan, unsigned iwindow) const
    { 
      return m_slice[islice].scope[iscope][ichan]
	.window_dev[iwindow-m_window_min]; 
    }

    double dev(unsigned iscope, unsigned ichan,
	       unsigned iwindow) const
    { 
      double dev_avg = 0;
      for(unsigned islice = 0; islice < nslice(); islice++)
	dev_avg+=dev(islice,iscope,ichan,iwindow);
      return dev_avg/(double)nslice();
    }

    double medianDev(unsigned iscope, unsigned ichan, unsigned iwindow,
		     unsigned min_event_count = 0) const;

    unsigned eventCount(unsigned islice, unsigned iscope, 
			unsigned ichan) const
    {
      return m_slice[islice].scope[iscope][ichan].ped_event_count;
    }

    inline unsigned getSliceByEventNum(unsigned event_num,
				       unsigned islice_hint = 0) const;
    inline unsigned getSliceByTime(const VSTime& event_time,
				   unsigned islice_hint = 0) const;

    unsigned sliceEventNumLo(unsigned islice) const
    { return m_slice[islice].event_num_lo; }
    unsigned sliceEventNumHi(unsigned islice) const
    { return m_slice[islice].event_num_hi; }
    unsigned sliceEventTimeLo(unsigned islice) const
    { return m_slice[islice].event_time_lo; }
    unsigned sliceEventTimeHi(unsigned islice) const
    { return m_slice[islice].event_time_hi; }

    void suppressScopeByDev(const unsigned iscope,
			    const double cut_lo, const double cut_hi,
			    VSTimeDepSupData& sup, 
			    unsigned min_event_count = 0);

    void suppressByMedianFraction(const double lo, const double hi,
				  VSTimeDepSupData& sup, 
				 std::vector<std::pair<double,double> >& cuts, 
				  unsigned min_event_count = 0);

    void suppressByContainmentFraction(const double fraction,
				       const double scale,
				       VSTimeDepSupData& sup, 
				 std::vector<std::pair<double,double> >& cuts, 
				       unsigned min_event_count = 0);

  private:
    
    void getAllDev(unsigned iscope, 
		   std::vector<std::pair<bool, double> >& alldev,
		   unsigned min_event_count) const;

    friend class VBFTimeDepPeds;

    struct Chan
    {
      Chan(): ped_event_count(), window_dev() { }
      unsigned ped_event_count;
      std::vector<double> window_dev;
    };

    typedef std::vector<Chan> Scope;

    struct Slice
    {
      Slice()
	: event_num_lo(), event_num_hi(), event_time_lo(), event_time_hi(),
	  ped_event_count(), scope() { /* nothing to see here */ }
      unsigned event_num_lo;
      unsigned event_num_hi;
      VSTime event_time_lo;
      VSTime event_time_hi;
      unsigned ped_event_count;
      std::vector<Scope> scope;
    };

    unsigned m_window_min;
    std::vector<std::vector<double> > m_hi_ped;
    std::vector<std::vector<double> > m_lo_ped;
    std::vector<std::vector<double> > m_hi_dev;
    std::vector<std::vector<VSSimpleHist<int32_t> > > m_min_window_hist;
    std::vector<Slice> m_slice;
  };

  inline unsigned VSTimeDepPedData::
  getSliceByEventNum(unsigned event_num, unsigned islice_hint) const
  {
    unsigned nslice=m_slice.size();
    if(islice_hint<nslice)
      {
	bool gt_lo = ((islice_hint==0)
		      ||(event_num>=m_slice[islice_hint].event_num_lo));
	bool lt_hi = ((islice_hint==nslice-1)
		      ||(event_num<m_slice[islice_hint].event_num_hi));
	if(gt_lo)
	  {
	    if(lt_hi)return islice_hint;
	    else if((islice_hint+1<nslice)
		    &&((islice_hint==nslice-2)
		       ||(event_num<m_slice[islice_hint+1].event_num_hi)))
	      return islice_hint+1;
	  }
	else
	  if((islice_hint)
	     &&((islice_hint==1)||
		(event_num>=m_slice[islice_hint-1].event_num_lo)))
	    return islice_hint-1;
      }

    for(unsigned islice=0;islice<nslice;islice++)
      if(event_num<m_slice[islice].event_num_hi)return islice;
    return nslice-1;
  }

  inline unsigned VSTimeDepPedData::
  getSliceByTime(const VSTime& event_time, unsigned islice_hint) const
  {
    unsigned nslice=m_slice.size();
    if(islice_hint<nslice)
      {
	bool gt_lo = ((islice_hint==0)
		      ||(event_time>=m_slice[islice_hint].event_time_lo));
	bool lt_hi = ((islice_hint==nslice-1)
		      ||(event_time<m_slice[islice_hint].event_time_hi));
	if(gt_lo)
	  {
	    if(lt_hi)return islice_hint;
	    else if((islice_hint+1<nslice)
		    &&((islice_hint==nslice-2)
		       ||(event_time<m_slice[islice_hint+1].event_time_hi)))
	      return islice_hint+1;
	  }
	else
	  if((islice_hint)
	     &&((islice_hint==1)||
		(event_time>=m_slice[islice_hint-1].event_time_lo)))
	    return islice_hint-1;
      }

    for(unsigned islice=0;islice<nslice;islice++)
      if(event_time<m_slice[islice].event_time_hi)return islice;
    return nslice-1;
  }

};

#endif // VSTIMEDEPPEDDATA_HPP
