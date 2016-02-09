//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMergedCalibrationData.hpp

  Time independent merged calibration data (laser, ped, suppression)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/15/2006

  $Id: VSMergedCalibrationData.hpp,v 3.8 2008/11/12 04:52:23 matthew Exp $

*/

#ifndef MERGEDCALIBRATIONDATA_HPP
#define MERGEDCALIBRATIONDATA_HPP

#include<vector>
#include<VSOctaveIO.hpp>
#include<VSSimpleHist.hpp>

namespace VERITAS
{

  class VSChannelSuppressedCount
  {
  public:
    VSChannelSuppressedCount():
      none(), ped(), pad(), laser(), ped_pad(), ped_laser(), pad_laser(), 
      all(), any() 
    { /* nothing to see here */ }

    unsigned                                     none;
    unsigned                                     ped;
    unsigned                                     pad;
    unsigned                                     laser;
    unsigned                                     ped_pad;
    unsigned                                     ped_laser;
    unsigned                                     pad_laser;
    unsigned                                     all;
    unsigned                                     any;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSChannelSuppressedCount,none);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,ped);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,pad);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,laser);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,ped_pad);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,ped_laser);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,pad_laser);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,all);
      H5_ADDMEMBER(c,VSChannelSuppressedCount,any);
    }
  };

  class VSChannelMergedCalibrationData
  {
  public:
    static const unsigned SR_LASER              = 0x00000001;
    static const unsigned SR_PED                = 0x00000002;
    static const unsigned SR_PAD_PED            = 0x00000004;
    static const unsigned SR_DEMAND             = 0x00000008;
    static const unsigned SR_NO_PMT             = 0x00000010;
    static const unsigned SR_ILLEGAL_PEDVAL     = 0x00000020;
    static const unsigned SR_LO_GAIN_PEDVAL     = 0x00000040;

    VSChannelMergedCalibrationData(unsigned nchan=0):
      ped_hi(), ped_lo(), median_dev(), hi_lo_gain_ratio(6.0), gain(1.0),
      has_pmt(true), suppress_all_events(true), suppress_lo_gain(false), 
      suppress_nslice(), suppress_reason(), 
      crate(5), l2chan(nchan), l2time(), chantime(), cratetime(), 
      median_current(), median_l1_rate(), 
      pmt_multiplier_gain(1.0), pmt_efficiency(1.0),
      pmt_coord_E_deg(), pmt_coord_U_deg()
    { /* nothing to see here */ }

    double                                       ped_hi;
    double                                       ped_lo;	      
    double                                       median_dev;	      
    double                                       hi_lo_gain_ratio;
    double                                       gain;
    bool                                         has_pmt;	      
    bool                                         suppress_all_events;
    bool                                         suppress_lo_gain;
    unsigned                                     suppress_nslice;     
    unsigned                                     suppress_reason;     
    unsigned                                     crate;
    unsigned                                     l2chan;	      
    double                                       l2time;	      
    double                                       chantime;	      
    double                                       cratetime;	      
    double                                       median_current;      
    double                                       median_l1_rate;          
    double                                       pmt_multiplier_gain;
    double                                       pmt_efficiency;
    double                                       pmt_coord_E_deg;
    double                                       pmt_coord_U_deg;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,ped_hi);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,ped_lo);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,median_dev);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,hi_lo_gain_ratio);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,gain);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,has_pmt);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,suppress_all_events);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,suppress_lo_gain);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,suppress_nslice);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,suppress_reason);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,crate);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,l2chan);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,l2time);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,chantime);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,cratetime);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,median_current);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,median_l1_rate);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,pmt_multiplier_gain);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,pmt_efficiency);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,pmt_coord_E_deg);
      H5_ADDMEMBER(c,VSChannelMergedCalibrationData,pmt_coord_U_deg);
    }
  };

  class VSScopeMergedCalibrationBase
  {
  public:
    VSScopeMergedCalibrationBase(unsigned nchan = 0):
      nchan(nchan), suppress(true), gain(1.0), 
      hi_sample_zero(), lo_sample_zero(), 
      position_E_m(), position_N_m(), position_U_m(), 
      median_median_dev(), scaled_dev(),
      suppressed_count(), median_dev_hist(0.1)
    { /* nothing to see here */ }

    unsigned                                     nchan;
    bool                                         suppress;
    double                                       gain;
    unsigned                                     hi_sample_zero;
    unsigned                                     lo_sample_zero;
    double                                       position_E_m;
    double                                       position_N_m;
    double                                       position_U_m;
    double                                       median_median_dev;
    double                                       scaled_dev;
    VSChannelSuppressedCount                     suppressed_count;

    VSSimpleHist<double>                         median_dev_hist;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,nchan);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,suppress);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,gain);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,hi_sample_zero);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,lo_sample_zero);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,position_E_m);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,position_N_m);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,position_U_m);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,median_median_dev);
      H5_ADDMEMBER(c,VSScopeMergedCalibrationBase,scaled_dev);
    }

    void load(VSOctaveH5ReaderStruct* s);
    void save(VSOctaveH5WriterStruct* s) const;
  };

  class VSArrayMergedCalibrationBase
  {
  public:
    VSArrayMergedCalibrationBase():
      mean_scaled_dev(), rms_scaled_dev()
    { /* nothing to see here */ }

    double                                       mean_scaled_dev;
    double                                       rms_scaled_dev;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSArrayMergedCalibrationBase,mean_scaled_dev);
      H5_ADDMEMBER(c,VSArrayMergedCalibrationBase,rms_scaled_dev);
    }
  };

  class VSScopeMergedCalibrationData: public VSScopeMergedCalibrationBase
  {
  public:
    VSScopeMergedCalibrationData(unsigned nchan = 0):
      VSScopeMergedCalibrationBase(nchan), channel(nchan,nchan)
    { /* nothing to see here */ }
      
    std::vector<VSChannelMergedCalibrationData>  channel;
  };

  class VSArrayMergedCalibrationData: public VSArrayMergedCalibrationBase
  {
  public:
    VSArrayMergedCalibrationData(const std::vector<unsigned>& nchan =
				 std::vector<unsigned>()):
      VSArrayMergedCalibrationBase(), scope()
    { /* nothing to see here */ }

    std::vector<VSScopeMergedCalibrationData>    scope;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSArrayMergedCalibrationBase,mean_scaled_dev);
    }

    void clear() { scope.clear(); }
    void load(VSOctaveH5ReaderStruct* s);
    void save(VSOctaveH5WriterStruct* s) const;    
  };

}

#endif // not defined MERGEDCALIBRATIONDATA_HPP
