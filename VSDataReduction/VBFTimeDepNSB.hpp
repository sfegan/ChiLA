//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFTimeDepNSB.hpp

  Calculate time dependent NSB, can be used to flag suppressed channels or
  estimate tracking performance

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VBFTimeDepNSB.hpp,v 3.1 2007/06/13 23:30:28 sfegan Exp $

*/

#ifndef VBFTIMEDEPNSB_HPP
#define VBFTIMEDEPNSB_HPP

#include <VSSimpleStat.hpp>
#include <VSSimpleHist.hpp>
#include <VSOctaveIO.hpp>
#include <VSSimpleVBF.hpp>
#include <VSTimeDepSupData.hpp>

namespace VERITAS
{

  class VBFTimeDepNSB: public VSSimpleVBFVisitor
  {
  public:
    VBFTimeDepNSB(unsigned window_start = 0, unsigned window_width = 2,
		  unsigned events_per_slice = 1000, unsigned nevent_cut = 400);

    virtual ~VBFTimeDepNSB();
    
    virtual void visitArrayEvent(bool& veto_array_event, void* user_data,
				 uint32_t               num_scope_events,
				 bool                   has_array_trigger,
				 bool                   has_event_number,
				 uint32_t               event_num,
 				 bool                   has_good_event_time,
				 const VSTime&          best_event_time,
				 EventType              event_type,
				 uint32_t               l2_trigger_mask,
				 const VArrayEvent*     array_event);

    virtual void visitScopeEvent(bool& veto_scope_event, void* user_data,
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
				 const VEvent*          event);

    virtual void visitHitChannel(void* user_data,
				 uint32_t               channel_num,
				 uint32_t               charge, 
				 uint32_t               pedestal,
				 bool                   lo_gain,
				 unsigned               nsample,
				 const uint32_t*        samples,
				 const uint32_t*        integrated);

    class Data
    {
    public:
      Data(): scope(), slice_time()
      { /* nothing to see here */ }

      class Scope
      {
      public:
	Scope(): all_dev(0.05), suppress_dev_suggested(), suppress_dev(),
		 slice_median_dev()
	{ /* nothing to see here */ }
	VSSimpleHist<double>    all_dev;
	double                  suppress_dev_suggested;
	double                  suppress_dev;
	std::vector<double>     slice_median_dev;
      };

      std::vector<Scope>        scope;
      std::vector<VATime>       slice_time;

      void clear();
      void load(VSOctaveH5ReaderStruct* reader);
      void save(VSOctaveH5WriterStruct* writer) const;
    };

    void getData(Data& data) const;
    void suppressFromPedDev(VSTimeDepSupData& sup, 
			    const std::vector<double>& sup_dev) const;
    unsigned nslice() const { return m_slice.size(); }

  private:
    VBFTimeDepNSB(const VBFTimeDepNSB&);
    VBFTimeDepNSB& operator= (const VBFTimeDepNSB&);

    typedef VSSimpleStat2<double> Channel;
    typedef std::vector<Channel> Scope;
    typedef std::vector<Scope> Array;

    struct Slice
    {
      Slice(): scope(), slice_time() { /* nothing to see here */ }
      Array* scope;
      VATime slice_time;
    };
    
    typedef VSSimpleStat1<double> ChannelPed;
    typedef std::vector<ChannelPed> ScopePed;
    typedef std::vector<ScopePed> ArrayPed;

    void resize(Array* slice);

    unsigned                 m_window_start;
    unsigned                 m_window_width;
    unsigned                 m_events_per_slice;
    unsigned                 m_nevent_cut;
    std::vector<Slice>       m_slice;
    std::vector<unsigned>    m_nchan;
    ArrayPed                 m_ped;
    bool                     m_is_ped_event;
    unsigned                 m_event_num;
    unsigned                 m_scope_num;
    Array*                   m_cached_slice;
    Scope*                   m_cached_scope;
    bool                     m_got_event;
    VSTime                   m_lo_time;
    VSTime                   m_hi_time;
  };

  typedef VBFTimeDepNSB::Data VSTimeDepNSBData;

}

#endif // not defined VBFTIMEDEPNSB_HPP
