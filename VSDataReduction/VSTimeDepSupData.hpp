//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTimeDepSupData.hpp

  Time dependant suppressed data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VSTimeDepSupData.hpp,v 3.0 2007/04/12 17:25:53 sfegan Exp $

*/

#ifndef VSTIMEDEPSUPDATA_HPP
#define VSTIMEDEPSUPDATA_HPP

#include<vector>

#include<VSRunInfoData.hpp>
#include<VSOctaveIO.hpp>

namespace VERITAS
{

  // ==========================================================================
  // TIME DEPENDENT SUPPRESSED DATA
  // ==========================================================================

  class VSTimeDepSupData
  {
  public:

    VSTimeDepSupData(unsigned events_per_slice=0)
      : m_events_per_slice(events_per_slice), 
	m_nslice(), m_scope() { /* nothing to see here */ }
    
    bool isSuppressedFast(const unsigned iword, const unsigned mask,
			  const unsigned iscope, const unsigned ichan) const
    {
      return m_scope[iscope][ichan][iword]&mask; 
    }

    bool isSuppressed(const unsigned islice, const unsigned iscope, 
		      const unsigned ichan) const
    { 
      const unsigned iword = islice>>5;
      const unsigned ibit  = islice&0x1F; 
      const unsigned mask  = 1<<ibit;
      return m_scope[iscope][ichan][iword]&mask; 
    }
    
    void setSuppressed(const unsigned islice, const unsigned iscope, 
		       const unsigned ichan, const bool suppressed = true)
    {
      const unsigned iword = islice>>5;
      const unsigned ibit  = islice&0x1F; 
      const unsigned mask  = 1<<ibit;
      if(suppressed)m_scope[iscope][ichan][iword] |= mask; 
      else m_scope[iscope][ichan][iword] &= ~mask;
    }

    void getFastSliceFindInfo(const unsigned islice,
			      unsigned& iword, unsigned& mask) const
    {
      iword                = islice>>5;
      const unsigned ibit  = islice&0x1F; 
      mask                 = 1<<ibit;
    }

    unsigned getSliceByEventNum(const unsigned event_num) const 
    {
      const unsigned islice = event_num/m_events_per_slice; 
      if(islice >= m_nslice)return m_nslice-1;
      return islice;
    }

    void resize(const unsigned events_per_slice,
		const unsigned nslice, const std::vector<unsigned>& nchans)
    {
      if(m_events_per_slice != events_per_slice)
	{
	  m_events_per_slice = events_per_slice;
	  m_scope.clear();
	}
      //m_nslice = nevents?((nevents-1)/events_per_slice)+1:0;
      m_nslice = nslice;
      unsigned nword = nslice?((nslice-1)>>5)+1:0;
      unsigned nscope = nchans.size();
      m_scope.resize(nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  unsigned nchan = nchans[iscope];
	  m_scope[iscope].resize(nchan);
	  for(unsigned ichan=0;ichan<nchan;ichan++)
	    m_scope[iscope][ichan].resize(nword);
	}
    }

    unsigned nslice() const { return m_nslice; }
    unsigned nscope() const { return m_scope.size(); }
    unsigned nchan(unsigned iscope) const { return m_scope[iscope].size(); }

    unsigned getSuppressedCount(unsigned iscope, unsigned ichan) const;

    bool isAlwaysSuppressed(unsigned iscope, unsigned ichan) const
    {
      return getSuppressedCount(iscope, ichan) == nslice();
    }

    void suppressFromIMon(const std::vector<VBFRunInfo::Slice>& slice);

    void clear();
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

  private:

    typedef std::vector<uint32_t> Chan;
    typedef std::vector<Chan> Scope;
    typedef std::vector<Scope> Array;

    unsigned m_events_per_slice;
    unsigned m_nslice;
    Array m_scope;
  };

} // namespace VERITAS

#endif // not defined VSTIMEDEPSUPDATA_HPP
