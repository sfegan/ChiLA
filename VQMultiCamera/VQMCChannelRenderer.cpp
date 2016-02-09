//-*-mode:c++; mode:font-lock;-*-

/*! \file VQMCChannelRenderer.cpp

  QT4 class which does the acrtual drawing of the channels in the camera

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/30/2005

  \see VQMultiChannel
*/

#include "VQMCChannelRenderer.hpp"

using namespace VERITAS;

// ----------------------------------------------------------------------------
// VQMCChannelRenderer
// ----------------------------------------------------------------------------

VQMCChannelRenderer::~VQMCChannelRenderer()
{
  // nothing to see here
}

void VQMCChannelRenderer::drawAnnotations(QPainter& painter)
{
  // nothing to see here
}

// ----------------------------------------------------------------------------
// VQMCSimpleChannelRendererBase
// ----------------------------------------------------------------------------

VQMCSimpleChannelRendererBase::~VQMCSimpleChannelRendererBase()
{
  // nothing to see here
}

void 
VQMCSimpleChannelRendererBase::blendColors(qreal x, QColor& c, const QColor& b)
{
  if(x<0)x=0;
  if(x>1)x=1;
  qreal omx=1-x;
  QColor temp;
  temp.setRgbF(x*c.redF()   + omx*b.redF(),
	       x*c.greenF() + omx*b.greenF(),
	       x*c.blueF()  + omx*b.blueF());
  c=temp;
}

void 
VQMCSimpleChannelRendererBase::
blendColors2(qreal x, QColor& c, const QColor& b)
{
  blendColors(x,c,b);
}
