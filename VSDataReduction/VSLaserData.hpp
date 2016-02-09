//-*-mode:c++; mode:font-lock;-*-

/*! \file VSLaserData.hpp
  Simple pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/18/2005

  $Id: VSLaserData.hpp,v 3.8 2009/11/24 01:11:45 matthew Exp $

*/

#ifndef VSLASERDATA_HPP
#define VSLASERDATA_HPP

#include<string>
#include<vector>

#include<VSOctaveIO.hpp>
#include<VSSimpleHist.hpp>

namespace VERITAS
{

  class VSLaserData
  {
  public:

    struct CrateData
    {
      CrateData():
	l2chan(), l2time(), cratetime_mean(), cratetime_dev() {}

      // Data -----------------------------------------------------------------
      unsigned l2chan;         // L2 channel 
      double   l2time;         // L2 arrival time 
      double   cratetime_mean; // Crate arrival time - L2 time
      double   cratetime_dev;  // Crate arrival time dev

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,CrateData,l2chan);
	H5_ADDMEMBER(c,CrateData,l2time);
	H5_ADDMEMBER(c,CrateData,cratetime_mean);
	H5_ADDMEMBER(c,CrateData,cratetime_dev);
      }
    };

    struct ChanData
    {
      ChanData(): 
	nflash(), l2chan(), crate(),
	uncalibrated(), excluded(), chantime(), cratetime(), l2time(), 
	gain(), gain_err(), absgain(), absgain_err(), eff(), eff_err(),
	signal_mean(), signal_dev(), signal_corr_mean(), signal_corr_dev(),
	signal_hist(0) { }

      // Data -----------------------------------------------------------------
      unsigned nflash;            // Number of flashes recorded
      unsigned l2chan;            // L2 channel
      unsigned crate;             // Crate number
      bool     uncalibrated;      // Has less than required # of laser flashes
      bool     excluded;          // Excluded from mean/median calculations
      double   chantime;          // Channel arrival time - L2 time
      double   cratetime;         // Crate arrival time - L2 time
      double   l2time;            // L2 arrival time for crate
      double   gain;              // Total gain relative to median
      double   gain_err;          
      double   absgain;           // Absolute single PE gain in DC 
      double   absgain_err;       
      double   eff;               // Efficiency relative to median
      double   eff_err;               
      double   signal_mean;       // Laser amplitude mean
      double   signal_dev;        // Laser amplitude dev
      double   signal_corr_mean;  // Corrected laser amplitude mean 
      double   signal_corr_dev;   // Corrected laser amplitude dev
     
      // Histograms -----------------------------------------------------------
      VSSimpleHist<double> signal_hist; // Histogram of signal_mean

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,ChanData,nflash);
	H5_ADDMEMBER(c,ChanData,l2chan);
	H5_ADDMEMBER(c,ChanData,crate);
	H5_ADDMEMBER(c,ChanData,uncalibrated);
	H5_ADDMEMBER(c,ChanData,excluded);
	H5_ADDMEMBER(c,ChanData,chantime);
	H5_ADDMEMBER(c,ChanData,cratetime);
	H5_ADDMEMBER(c,ChanData,l2time);
	H5_ADDMEMBER(c,ChanData,gain);
	H5_ADDMEMBER(c,ChanData,gain_err);
	H5_ADDMEMBER(c,ChanData,absgain);
	H5_ADDMEMBER(c,ChanData,absgain_err);
	H5_ADDMEMBER(c,ChanData,eff);
	H5_ADDMEMBER(c,ChanData,eff_err);
	H5_ADDMEMBER(c,ChanData,signal_mean);
	H5_ADDMEMBER(c,ChanData,signal_dev);
	H5_ADDMEMBER(c,ChanData,signal_corr_mean);
	H5_ADDMEMBER(c,ChanData,signal_corr_dev);
      }
    };

    struct ScopeData
    {
      ScopeData(): 
	runno(), nchan(0), nevents(0), nflash(0),
	absgain_mean(), absgain_median(), npe_mean(), npe_median(),
	nchan_mean(), nchan_dev(), signal_mean(), signal_dev(),
	signal_hist(0), nchan_flash_hist(0),
	nchan_logain_hist(0), nchan_higain_hist(0),
	chan(), crate() {}

      // Data -----------------------------------------------------------------
      unsigned runno;
      unsigned nchan;             // Number channels in this telescope
      unsigned nevents;           // Number events in laser run
      unsigned nflash;            // Number events passing laser criteria
      double   absgain_mean;      // Mean absolute gain
      double   absgain_median;    // Median absolute gain
      double   npe_mean;          // Mean number of PEs per pixel
      double   npe_median;        // Median number of PEs per pixel
      double   nchan_mean;        // Number of channels in each flash mean
      double   nchan_dev;         // Number of channels in each flash dev
      double   signal_mean;       // Average telescope laser amplitude mean 
      double   signal_dev;        // Average telescope laser amplitude dev 

      // Histograms -----------------------------------------------------------
      VSSimpleHist<double>   signal_hist;
      VSSimpleHist<unsigned> nchan_flash_hist;   
      VSSimpleHist<unsigned> nchan_logain_hist;
      VSSimpleHist<unsigned> nchan_higain_hist;

      // Vectors of structures ------------------------------------------------
      std::vector<ChanData>               chan;
      std::vector<CrateData>              crate;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,ScopeData,runno);
	H5_ADDMEMBER(c,ScopeData,nchan);
	H5_ADDMEMBER(c,ScopeData,nevents);
	H5_ADDMEMBER(c,ScopeData,nflash);
	H5_ADDMEMBER(c,ScopeData,absgain_mean);
	H5_ADDMEMBER(c,ScopeData,absgain_median);
	H5_ADDMEMBER(c,ScopeData,npe_mean);
	H5_ADDMEMBER(c,ScopeData,npe_median);
	H5_ADDMEMBER(c,ScopeData,nchan_mean);
	H5_ADDMEMBER(c,ScopeData,nchan_dev);
	H5_ADDMEMBER(c,ScopeData,signal_mean);
	H5_ADDMEMBER(c,ScopeData,signal_dev);
      }
    };

    // Data -------------------------------------------------------------------
    unsigned m_runno;

    // Settings ---------------------------------------------------------------
    unsigned m_threshold_nchan;
    double   m_threshold_dc;
    double   m_singlepe_dev;

    // Vectors of structures --------------------------------------------------
    std::vector<ScopeData> scope;

    VSLaserData(): m_runno(), m_threshold_nchan(), m_threshold_dc(), 
		   m_singlepe_dev(), scope() 
    { /* nothing to see here */ }

    bool load(const std::string& filename);
    void suppress(const double lo, const double hi);

    void clear();
    void load(VSOctaveH5ReaderStruct* reader);
    bool load(unsigned iscope,VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSLaserData,m_runno);
      H5_ADDMEMBER(c,VSLaserData,m_threshold_nchan);
      H5_ADDMEMBER(c,VSLaserData,m_threshold_dc);
      H5_ADDMEMBER(c,VSLaserData,m_singlepe_dev);
    }

  private:
    VSLaserData(const VSLaserData&);
    VSLaserData& operator= (const VSLaserData&);
  };

};

#endif // VSLASERDATA_HPP
