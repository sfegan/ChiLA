//-*-mode:c++; mode:font-lock;-*-

/*! \file VQMCArrayData.hpp

  Data classes for VQMultiCamera code. Classes to hold channel/camera/array
  data and limits.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/30/2005

  \see VQMultiCamera
*/

#ifndef VQMCARRAYDATA_HPP
#define VQMCARRAYDATA_HPP

namespace VERITAS
{

  class VQMCChannelData
  {
  public:
    VQMCChannelData(bool no_data = false):
      fValue(), fType(), fSuppress(), fHighlight(), fNoDraw(), fNoData(no_data)
    { /* nothing to see here */ }
    double    fValue;
    unsigned  fType;
    bool      fSuppress;
    bool      fHighlight;
    bool      fNoDraw;
    bool      fNoData;
    inline bool operator== (const VQMCChannelData& o) const;
    bool operator!= (const VQMCChannelData& o) const { return !(*this==o); }
  };

  class VQMCCameraData
  {
  public:
    VQMCCameraData(unsigned nchannels=0, 
		   const VQMCChannelData& _default = VQMCChannelData()):
      fChannelData(nchannels, _default) { /* nothing to see here */ }
    std::vector<VQMCChannelData> fChannelData;
  };

  class VQMCArrayData
  {
  public:
    VQMCArrayData(unsigned ncameras=0, unsigned nchannels=0,
		  const VQMCChannelData& _default = VQMCChannelData()): 
      fCameraData(ncameras, VQMCCameraData(nchannels, _default)) 
    { /* nothing to see here */ }
    std::vector<VQMCCameraData> fCameraData;
  };

  class VQMCCameraDataLimits
  {
  public:
    VQMCCameraDataLimits(): 
      fLimitLo(), fLimitHi() { /* nothing to see here */ }
    VQMCCameraDataLimits(double val_lo, double val_hi): 
      fLimitLo(), fLimitHi() 
    { 
      fLimitLo.fValue = val_lo;
      fLimitHi.fValue = val_hi;
    }
    VQMCChannelData fLimitLo;
    VQMCChannelData fLimitHi;
    inline void startLimitComputation(const VQMCChannelData& data);
    inline void extendLimits(const VQMCChannelData& data);
    bool operator== (const VQMCCameraDataLimits& o) const
    { return (fLimitLo==o.fLimitLo)&&(fLimitHi==o.fLimitHi); }
    bool operator!= (const VQMCCameraDataLimits& o) const 
    { return !(*this==o); }
  private:
    bool fFoundData;
  };

  class VQMCArrayDataLimits
  {
  public:
    VQMCArrayDataLimits(unsigned ncameras=0): 
      fArrayLimits(ncameras) { /* nothing to see here */ }
    std::vector<VQMCCameraDataLimits> fArrayLimits;
  };

  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //
  // FUNCTION DEFINITIONS
  //
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------

  inline bool VQMCChannelData::operator== (const VQMCChannelData& o) const
  {
    return 
      (fValue == o.fValue)
      &&(fType == o.fType)
      &&(fSuppress == o.fSuppress)
      &&(fHighlight == o.fHighlight)
      &&(fNoData == o.fNoData)
      &&(fNoDraw == o.fNoDraw);
  }

  inline void VQMCCameraDataLimits::
  startLimitComputation(const VQMCChannelData& data)
  {
    fFoundData = false;
    fLimitLo.fValue = 0;
    fLimitHi.fValue = 0;
    fLimitLo.fType = 0;
    fLimitHi.fType = 0;
    extendLimits(data);
  }
  
  inline void VQMCCameraDataLimits::extendLimits(const VQMCChannelData& data)
  {
    if(data.fNoData)return;
    if(fFoundData)
      {
	if(fLimitHi.fValue < data.fValue)fLimitHi.fValue = data.fValue;
	if(fLimitHi.fType < data.fType)fLimitHi.fType = data.fType;
	if(fLimitLo.fValue > data.fValue)fLimitLo.fValue = data.fValue;
	if(fLimitLo.fType > data.fType)fLimitLo.fType = data.fType;
      }
    else
      {
	fFoundData = true;
	fLimitLo.fValue = data.fValue;
	fLimitHi.fValue = data.fValue;
	fLimitLo.fType = data.fType;
	fLimitHi.fType = data.fType;
      }
  }

}

#endif // VQMCARRAYDATA_HPP
