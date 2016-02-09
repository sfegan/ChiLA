//-*-mode:c++; mode:font-lock;-*-

/*! \file VSHiLoData.hpp
  Simple pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/18/2005

  $Id: VSHiLoData.hpp,v 3.2 2008/02/23 23:28:35 sfegan Exp $

*/

#ifndef VSHILODATA_HPP
#define VSHILODATA_HPP

#include<string>
#include<vector>

#include<VSOctaveIO.hpp>
#include<VSSimpleHist.hpp>

namespace VERITAS
{

  class VSHiLoData
  {
  public:

    struct ChanData
    {
      ChanData(): 
	hi_lo_gain_ratio(), lo_gain_ped(), lo_gain_switch_amp(), 
	has_hi_lo_gain_ratio(), has_lo_gain_ped(), has_lo_gain_switch()
      { }

      // Data -----------------------------------------------------------------
      double   hi_lo_gain_ratio;   // Ratio of High:Low gain
      double   lo_gain_ped;        // Low gain pedestal value
      double   lo_gain_switch_amp; // Switch-over amplitude (high gain channel)
      bool     has_hi_lo_gain_ratio;
      bool     has_lo_gain_ped;
      bool     has_lo_gain_switch;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,ChanData,hi_lo_gain_ratio);
	H5_ADDMEMBER(c,ChanData,lo_gain_ped);
	H5_ADDMEMBER(c,ChanData,lo_gain_switch_amp);
	H5_ADDMEMBER(c,ChanData,has_hi_lo_gain_ratio);
	H5_ADDMEMBER(c,ChanData,has_lo_gain_ped);
	H5_ADDMEMBER(c,ChanData,has_lo_gain_switch);
      }
    };

    struct ScopeData
    {
      ScopeData(unsigned _nchan = 0): chan(_nchan) {}

      // Vectors of structures ------------------------------------------------
      std::vector<ChanData>               chan;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
      }
    };

    // Data -------------------------------------------------------------------
    unsigned runno;
    bool     has_hi_lo_gain_ratio;
    bool     has_lo_gain_ped;
    bool     has_lo_gain_switch;

    // Vectors of structures --------------------------------------------------
    std::vector<ScopeData> scope;

    VSHiLoData(): 
      runno(), has_hi_lo_gain_ratio(), has_lo_gain_ped(), has_lo_gain_switch(),
      scope() { }

    bool load(const std::string& filename);

    void clear();
    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSHiLoData,runno);
      H5_ADDMEMBER(c,VSHiLoData,has_hi_lo_gain_ratio);
      H5_ADDMEMBER(c,VSHiLoData,has_lo_gain_ped);
      H5_ADDMEMBER(c,VSHiLoData,has_lo_gain_switch);
    }

  private:
    VSHiLoData(const VSHiLoData&);
    VSHiLoData& operator= (const VSHiLoData&);
  };

};

#endif // VSHILODATA_HPP
