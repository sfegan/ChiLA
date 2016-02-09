//-*-mode:c++; mode:font-lock;-*-

/*! \file VSimpleCameraLayout.hpp

  Simple class to contain tube positions. This class should have very little
  in it so that it can be adapted to the widest possible camera classes that
  exist out there, for example those from the OAWG and from VSOptics.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/30/2005
*/

#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>

#ifndef VSIMPLECAMERALAYOUT_HPP
#define VSIMPLECAMERALAYOUT_HPP

#ifndef NO_WHIPPLE_CAMERAS
#define USE_WHIPPLE_CAMERAS
#endif

namespace VERITAS
{

  // --------------------------------------------------------------------------
  // VSimpleCameraChannel
  // --------------------------------------------------------------------------

  class VSimpleCameraChannel
  {
  public:
    VSimpleCameraChannel(double x, double y, double r, 
			 const std::string& name):
      fXDeg(x), fYDeg(y), fRDeg(r), fName(name) { }
    VSimpleCameraChannel(double x, double y, double r, unsigned id):
      fXDeg(x), fYDeg(y), fRDeg(r), fName() 
    { char buffer[20]; buffer[0]='\0'; sprintf(buffer,"%u",id); fName=buffer; }
    virtual ~VSimpleCameraChannel();
    
    double distanceSquared(double x, double y) const
    { return (x-fXDeg)*(x-fXDeg)+(y-fYDeg)*(y-fYDeg); }
    bool isInside(double x, double y) const
    { return distanceSquared(x,y)<=fRDeg*fRDeg; }

    double xDeg() const { return fXDeg; }
    double yDeg() const { return fYDeg; }
    double rDeg() const { return fRDeg; }
    const std::string& name() const { return fName; }

  private:
    double fXDeg;
    double fYDeg;
    double fRDeg;
    std::string fName;
  };

  // --------------------------------------------------------------------------
  // VSimpleCameraLayout
  // --------------------------------------------------------------------------

  class VSimpleCameraLayout
  {
  public:
    typedef VSimpleCameraChannel                            channel_type;
    typedef std::vector<channel_type>                       channel_list;
    typedef channel_list::const_iterator                    channel_iterator;
    typedef channel_list::size_type                         channel_size_type;

    VSimpleCameraLayout(): 
      fChannels(), fExtentRight(), fExtentLeft(), fExtentTop(),
      fExtentBottom(), fExtentRadial() { }
    virtual ~VSimpleCameraLayout();

    void addChannel(const channel_type& channel);

    channel_size_type numChannels() const { return fChannels.size(); }
    channel_iterator beginChannels() const { return fChannels.begin(); }
    channel_iterator endChannels() const { return fChannels.end(); }
    const channel_type& channel(channel_size_type i) const
    { return fChannels[i]; }
    const channel_list& channels() const { return fChannels; }

    double extentRight() const { return fExtentRight; }
    double extentLeft() const { return fExtentLeft; }
    double extentTop() const { return fExtentTop; }
    double extentBottom() const { return fExtentBottom; }
    double extentRadial() const { return fExtentRadial; }

  private:
    channel_list fChannels;
    double fExtentRight;
    double fExtentLeft;
    double fExtentTop;
    double fExtentBottom;
    double fExtentRadial;
  };

  // --------------------------------------------------------------------------
  // VSimpleCameraFactory
  // --------------------------------------------------------------------------
  
  class VSimpleCameraFactory
  {
  public:
    virtual ~VSimpleCameraFactory();
#ifdef USE_WHIPPLE_CAMERAS
    enum WhippleCamera { WC109, WC151, WC331, WC490, VC499, VC499BIGTUBES,
			 VC500, VC500BIGTUBES };
    VSimpleCameraLayout* getWhippleCamera(WhippleCamera wc) const;
#endif
    VSimpleCameraLayout* getHexCamera(unsigned nrings, double r) const;
    VSimpleCameraLayout* getRoundedHexCamera(double camera_r, double r) const;
    
    static VSimpleCameraFactory* getInstance();
  protected:
    VSimpleCameraFactory() { /* nothing to see here */ }
  private:
    static std::auto_ptr<VSimpleCameraFactory> sInstance;
  };

};

#endif // VSIMPLECAMERALAYOUT_HPP
