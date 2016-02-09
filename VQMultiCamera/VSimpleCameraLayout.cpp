//-*-mode:c++; mode:font-lock;-*-

/*! \file VSimpleCameraLayout.cpp

  Simple class to contain tube positions. This class should have very little
  in it so that it can be adapted to the widest possible camera classes that
  exist out there, for example those from the OAWG and from VSOptics.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/30/2005
*/

#include <cmath>

#include <xytohex.h>

#include "VSimpleCameraLayout.hpp"

using namespace VERITAS;

// ----------------------------------------------------------------------------
// VSimpleCameraChannel
// ----------------------------------------------------------------------------

VSimpleCameraChannel::~VSimpleCameraChannel()
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// VSimpleCameraLayout
// ----------------------------------------------------------------------------

VSimpleCameraLayout::~VSimpleCameraLayout()
{
  // nothing to see here
}

void VSimpleCameraLayout::addChannel(const channel_type& channel)
{
  float l = channel.xDeg()-channel.rDeg();
  float r = channel.xDeg()+channel.rDeg();
  float b = channel.yDeg()-channel.rDeg();
  float t = channel.yDeg()+channel.rDeg();
  float radius = sqrt(channel.distanceSquared(0,0))+channel.rDeg();
  fChannels.push_back(channel);
  if(l<fExtentLeft)fExtentLeft=l;
  if(r>fExtentRight)fExtentRight=r;
  if(b<fExtentBottom)fExtentBottom=b;
  if(t>fExtentTop)fExtentTop=t;
  if(radius>fExtentRadial)fExtentRadial=radius;
}

// ----------------------------------------------------------------------------
// VSimpleCameraFactory
// ----------------------------------------------------------------------------

std::auto_ptr<VSimpleCameraFactory> VSimpleCameraFactory::sInstance;

VSimpleCameraFactory::~VSimpleCameraFactory()
{
  // nothing to see here
}

VSimpleCameraFactory* VSimpleCameraFactory::getInstance()
{
  if(!sInstance.get())sInstance.reset(new VSimpleCameraFactory);
  return sInstance.get();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Functions supporting creation of standard Whipple (VERITAS) cameras
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#ifdef USE_WHIPPLE_CAMERAS

#include "WhippleCams.h"

#define NUMOF(x) (sizeof((x))/sizeof((*x)))

VSimpleCameraLayout* VSimpleCameraFactory::
getWhippleCamera(VSimpleCameraFactory::WhippleCamera wc) const
{
  float* x;
  float* y;
  float* r;
  unsigned n;

  switch(wc)
    {
    case WC109:
      x = WC109Xcoord;
      y = WC109Ycoord;
      r = WC109Radius;
      n = NUMOF(WC109Xcoord);
      break;
    case WC151:
      x = WC151Xcoord;
      y = WC151Ycoord;
      r = WC151Radius;
      n = NUMOF(WC151Xcoord);
      break;
    case WC331:
      x = WC331Xcoord;
      y = WC331Ycoord;
      r = WC331Radius;
      n = NUMOF(WC331Xcoord);
      break;
    case WC490:
      x = WC490Xcoord;
      y = WC490Ycoord;
      r = WC490Radius;
      n = NUMOF(WC490Xcoord);
      break;
    case VC499:
    case VC500:
      x = VC499Xcoord;
      y = VC499Ycoord;
      r = VC499Radius;
      n = NUMOF(VC499Xcoord);
      break;
    case VC499BIGTUBES:
    case VC500BIGTUBES:
      x = VC499Xcoord;
      y = VC499Ycoord;
      r = VC499BigRadius;
      n = NUMOF(VC499Xcoord);
      break;
    default:
      return 0;
      break;
    }

  VSimpleCameraLayout* cam = new VSimpleCameraLayout;
  
  for(unsigned ichan=0;ichan<n;ichan++)
    {
      VSimpleCameraChannel chan(x[ichan],y[ichan],r[ichan],ichan+1);
      cam->addChannel(chan);
    }

  unsigned N = n;
  if((wc==VC500)||(wc==VC500BIGTUBES))N=500;
  
  if(N>n)
    {
      float x0 = x[0];
      float y0 = y[0];
      for(unsigned ichan=1;ichan<n;ichan++)
	{
	  if(x[ichan]<x0)x0=x[ichan];
	  if(y[ichan]>y0)y0=y[ichan];
	}
      float delta = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
     
      for(unsigned ichan=n;ichan<N;ichan++)
	{
	  VSimpleCameraChannel chan(x0, y0-delta*(ichan-n),r[0],ichan+1);
	  cam->addChannel(chan);
	}
    }

  return cam;
}

#endif // USE_WHIPPLE_CAMERAS

VSimpleCameraLayout* VSimpleCameraFactory::
getHexCamera(unsigned nrings, double r) const
{
  VSimpleCameraLayout* cam = new VSimpleCameraLayout;

  unsigned nchan = 3*nrings*(nrings+1)+1;
  for(unsigned ichan=0;ichan<nchan;ichan++)
    {
      int n = ichan+1;
      double x = 0;
      double y = 0;
      nh_to_xy(&n, &x, &y);
      x*=2.0*r;
      y*=2.0*r;
      VSimpleCameraChannel chan(x,y,r,ichan+1);
      cam->addChannel(chan);
    }

  return cam;
}

VSimpleCameraLayout*  VSimpleCameraFactory::
getRoundedHexCamera(double camera_r, double r) const
{
  VSimpleCameraLayout* cam = new VSimpleCameraLayout;

  // Add the center channel
  VSimpleCameraChannel chan(0,0,r,1);
  cam->addChannel(chan);

  // Keep adding hexagonal rings until no more tubes are added to the camera

  unsigned iring = 1;
  unsigned ichan = 1;
  unsigned nring_chan_added;
  do
    {
      nring_chan_added = 0;
      unsigned nring_chan = 6*iring;
      for(unsigned iring_chan=0;iring_chan<nring_chan;iring_chan++)
	{
	  int n = ichan+1;
	  double x = 0;
	  double y = 0;
	  nh_to_xy(&n, &x, &y);
	  x*=2.0*r;
	  y*=2.0*r;
	  if(((x*x+y*y)<=camera_r*camera_r))
	    {
	      VSimpleCameraChannel chan(x,y,r,ichan+1);
	      cam->addChannel(chan);
	      nring_chan_added++;
	    }
	  ichan++;
	}
      iring++;
    }while(nring_chan_added);

  return cam;
}
