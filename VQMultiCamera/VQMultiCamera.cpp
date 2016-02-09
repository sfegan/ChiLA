//-*-mode:c++; mode:font-lock;-*-

/*! \file VQMultiCamera.hpp

  QT4 camera viewing/drawing classes. These classes support drawing of
  multiple cameras in one widget. Camera drawing is done by the class
  VQMultiCameraRenderer with the help of VQMCChannelRenderer. The widget
  class which can be used in Qt GUIs is VQMultiCamera.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/30/2005
*/

#include <cmath>
#include <cfloat>
#include <algorithm>

#include <Qt/QtGui>

#include <VSTestTimer.hpp>

#include "VQMultiCamera.hpp"

using namespace VERITAS;

// ----------------------------------------------------------------------------
// UTILITY FUNCTIONS
// ----------------------------------------------------------------------------

template<typename T>
static inline void MAXIMIZE(T& x, const T& y)
{
  if(x<y)x=y;
}

template<typename T>
static inline void MINIMIZE(T& x, const T& y)
{
  if(y<x)x=y;
}

template<typename T>
static inline void UPDATE_AND_NOTE_CHANGE(T& x, const T& y, bool& change)
{
  if(x!=y)
    {
      x=y;
      change=true;
    }
}

template<typename T>
static inline void RESIZE_AND_NOTE_CHANGE(T& x, typename T::size_type count, 
					  bool& change)
{
  if(x.size() != count)
    {
      x.resize(count);
      change=true;
    }
}

// ----------------------------------------------------------------------------
// VQMultiCameraRenderer
// ----------------------------------------------------------------------------

VQMultiCameraRenderer::
VQMultiCameraRenderer(VQMCChannelRenderer* channel_renderer,
		      const camera_list& cameras):
  fChannelRenderer(channel_renderer), fCameras(cameras),
  fDemandGlobalRotationRad(),
  fDemandGlobalCameraRotationRad(),
  fDemandCameraRotationRad(),
  fDemandCameraOffsetXDeg(),
  fDemandCameraOffsetYDeg(),
  fDemandRadialExtentDeg(),
  fDemandZoom(1.0),
  fDemandPanX(0.0),
  fDemandPanY(0.0),
  fDemandUseCommonLimits(),
  fDemandUseCameraLimits(),
  fDemandCameraLimits(),
  fDemandCRM(CRM_ALL_VALUES_BY_CAMERA_AND_CHANNEL),
  fDemandClipToLimits(false),
  fDemandLabelChannels(true),
  fDemandFontFamily("helvetica"),
  fDemandFontScale(0.9),
  fCachedExtentL(),
  fCachedExtentR(),
  fCachedExtentD(),
  fCachedExtentU(),
  fCachedGlobalRotationCos(),
  fCachedGlobalRotationSin(),
  fCachedCameraRotationCos(),
  fCachedCameraRotationSin(),
  fCachedChannelCenterX(),
  fCachedChannelCenterY(),
  fCachedChannelScale(),
  fCachedChannelRotation(),
  fCachedHasCameraLimits(),
  fCachedCameraLimits(),
  fCachedFontPointSize(),
  fCachedPainterScale(),
  fCachedLabelLengthToRadiusRatio(),
  fCachedRenderOrder(),
  fCachedRenderCount()
{
  calculateCachedParameters();
}

VQMultiCameraRenderer::~VQMultiCameraRenderer()
{
  // nothing to see here
}

bool VQMultiCameraRenderer::calculateCachedParameters()
{
  bool change=false;

  unsigned ncam = fCameras.size();
  RESIZE_AND_NOTE_CHANGE(fDemandCameraRotationRad,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fDemandChannelRotationRad,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fDemandCameraOffsetXDeg,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fDemandCameraOffsetYDeg,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fDemandUseCameraLimits,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fDemandCameraLimits.fArrayLimits,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedCameraRotationCos,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedCameraRotationSin,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedChannelRotationCos,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedChannelRotationSin,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedChannelCenterX,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedChannelCenterY,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedChannelScale,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedChannelRotation,ncam,change);
  RESIZE_AND_NOTE_CHANGE(fCachedCameraLimits.fArrayLimits,ncam,change);

  unsigned total_nchan = 0;
  for(unsigned icam=0;icam!=ncam;icam++)
    {
      unsigned nchan = 0;
      if(fCameras[icam])nchan=fCameras[icam]->numChannels();
      RESIZE_AND_NOTE_CHANGE(fCachedChannelCenterX[icam],nchan,change);
      RESIZE_AND_NOTE_CHANGE(fCachedChannelCenterY[icam],nchan,change);
      RESIZE_AND_NOTE_CHANGE(fCachedChannelScale[icam],nchan,change);
      total_nchan+=nchan;
    }

  RESIZE_AND_NOTE_CHANGE(fCachedRenderOrder,total_nchan,change);

  double extent_l = -1;
  double extent_r =  1;
  double extent_d = -1;
  double extent_u =  1;

  if(fDemandRadialExtentDeg>0)
    {
      extent_l=fDemandRadialExtentDeg;
      extent_r=fDemandRadialExtentDeg;
      extent_d=fDemandRadialExtentDeg;
      extent_u=fDemandRadialExtentDeg;
    }
  else
    {
      bool initialized = false;
      for(unsigned icam=0;icam!=ncam;icam++)
	if(fCameras[icam])
	  {
	    double l = 
	      fDemandCameraOffsetXDeg[icam]+fCameras[icam]->extentLeft();
	    double r = 
	      fDemandCameraOffsetXDeg[icam]+fCameras[icam]->extentRight();
	    double d = 
	      fDemandCameraOffsetYDeg[icam]+fCameras[icam]->extentBottom();
	    double u = 
	      fDemandCameraOffsetYDeg[icam]+fCameras[icam]->extentTop();

	    if(initialized)
	      {
		MINIMIZE(extent_l,l);
		MAXIMIZE(extent_r,r);
		MINIMIZE(extent_d,d);
		MAXIMIZE(extent_u,u);
	      }
	    else
	      {
		extent_l = l;
		extent_r = r;
		extent_d = d;
		extent_u = u;
		initialized = true;
	      }
	}
      if(initialized)
	{
	  extent_l -= (extent_r-extent_l)*0.02;
	  extent_r += (extent_r-extent_l)*0.02;
	  extent_d -= (extent_u-extent_d)*0.02;
	  extent_u += (extent_u-extent_d)*0.02;
	}
    }

  UPDATE_AND_NOTE_CHANGE(fCachedExtentL, extent_l, change);
  UPDATE_AND_NOTE_CHANGE(fCachedExtentR, extent_r, change);
  UPDATE_AND_NOTE_CHANGE(fCachedExtentD, extent_d, change);
  UPDATE_AND_NOTE_CHANGE(fCachedExtentU, extent_u, change);

  UPDATE_AND_NOTE_CHANGE(fCachedGlobalRotationCos, 
			 cos(fDemandGlobalRotationRad), change);
  UPDATE_AND_NOTE_CHANGE(fCachedGlobalRotationSin, 
			 sin(fDemandGlobalRotationRad), change);

  for(unsigned icam=0;icam!=ncam;icam++)
    {
      param_type rot = 
	fDemandGlobalRotationRad
	+fDemandGlobalCameraRotationRad
	+fDemandCameraRotationRad[icam];

      UPDATE_AND_NOTE_CHANGE(fCachedCameraRotationCos[icam], cos(rot), change);
      UPDATE_AND_NOTE_CHANGE(fCachedCameraRotationSin[icam], sin(rot), change);

      rot += fDemandChannelRotationRad[icam];
      UPDATE_AND_NOTE_CHANGE(fCachedChannelRotationCos[icam], 
			     cos(rot), change);
      UPDATE_AND_NOTE_CHANGE(fCachedChannelRotationSin[icam], 
			     sin(rot), change);
    }

  double max_label_length_to_radius_ratio = 0;
  
  for(unsigned icam=0;icam!=ncam;icam++)
    if(fCameras[icam])
      {
	const VSimpleCameraLayout* cam = fCameras[icam];

	qreal camera_origin_x = 
	  fDemandCameraOffsetXDeg[icam]*fCachedGlobalRotationCos 
	  - fDemandCameraOffsetYDeg[icam]*fCachedGlobalRotationSin;
	qreal camera_origin_y = 
	  fDemandCameraOffsetYDeg[icam]*fCachedGlobalRotationCos 
	  + fDemandCameraOffsetXDeg[icam]*fCachedGlobalRotationSin;

	const qreal camera_rotation_deg = getActualCameraRotationDeg(icam);
	const qreal camera_rotation_cos = fCachedCameraRotationCos[icam];
	const qreal camera_rotation_sin = fCachedCameraRotationSin[icam];
	const bool has_camera_rotation = 
	  fabs(camera_rotation_deg)>DBL_EPSILON;

	const qreal channel_rotation_deg = getActualChannelRotationDeg(icam);
 	fCachedChannelRotation[icam]        = channel_rotation_deg;
 	
	unsigned nchan = fCameras[icam]->numChannels();
	for(unsigned ichan=0;ichan!=nchan;ichan++)
	  {
	    const VSimpleCameraChannel* chan = &cam->channel(ichan);

	    qreal channel_origin_x;
	    qreal channel_origin_y;
	    qreal channel_scale = chan->rDeg()/1.0;
	    
	    if(has_camera_rotation)
	      {
		channel_origin_x = 
		  chan->xDeg()*camera_rotation_cos
		  - chan->yDeg()*camera_rotation_sin;
		channel_origin_y = 
		  chan->yDeg()*camera_rotation_cos
		  + chan->xDeg()*camera_rotation_sin;
	      }
	    else
	      {
		channel_origin_x = chan->xDeg();
		channel_origin_y = chan->yDeg();
	      }

	    fCachedChannelCenterX[icam][ichan]  = 
	      camera_origin_x+channel_origin_x;
	    fCachedChannelCenterY[icam][ichan]  = 
	      camera_origin_y+channel_origin_y;
	    fCachedChannelScale[icam][ichan]    =
	      channel_scale;

	    const double label_length_to_radius_ratio = 
	      double(chan->name().length())/(2.0*chan->rDeg());
	    MAXIMIZE(max_label_length_to_radius_ratio, 
		     label_length_to_radius_ratio);
	  }
      }

  UPDATE_AND_NOTE_CHANGE(fCachedLabelLengthToRadiusRatio,
			 max_label_length_to_radius_ratio, change);

  return change;
}

bool VQMultiCameraRenderer::setChannelRenderer(VQMCChannelRenderer* renderer)
{
  fChannelRenderer=renderer;
  configurationChanged();
  return true;
}

bool VQMultiCameraRenderer::setCameras(const camera_list& cameras)
{
  setNumCameras(cameras.size());
  for(unsigned icam=0; icam!=cameras.size(); icam++)
    setCamera(icam, cameras[icam]);
  configurationChanged();
  return true;
}

bool VQMultiCameraRenderer::setNumCameras(camera_size_type ncam)
{
  if(ncam != fCameras.size())
    {
      fCameras.resize(ncam);
      calculateCachedParameters();
      fCachedHasCameraLimits=false;
      configurationChanged();
      return true;
    }
  else return false;
}

bool VQMultiCameraRenderer::
setCamera(camera_size_type icam, camera_type camera)
{
  fCameras[icam]=camera;
  calculateCachedParameters();
  fCachedHasCameraLimits=false;
  configurationChanged();
  return true;
}

bool VQMultiCameraRenderer::setAutomaticExtent()
{
  fDemandRadialExtentDeg=0;
  bool change = calculateCachedParameters();
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::setDemandRadialExtentDeg(param_type extent)
{
  fDemandRadialExtentDeg=extent;
  bool change = calculateCachedParameters();
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::setGlobalRotationDeg(param_type angle)
{
  fDemandGlobalRotationRad=angle/180.0*M_PI;
  return calculateCachedParameters();  
}

bool VQMultiCameraRenderer::setGlobalCameraRotation(param_type angle)
{
  fDemandGlobalCameraRotationRad=angle/180.0*M_PI;
  bool change = calculateCachedParameters();
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::
setCameraRotationDeg(camera_size_type icam, param_type angle)
{
  fDemandCameraRotationRad[icam]=angle/180.0*M_PI;
  bool change = calculateCachedParameters();
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::
setChannelRotationDeg(camera_size_type icam, param_type angle)
{
  fDemandChannelRotationRad[icam]=angle/180.0*M_PI;
  bool change = calculateCachedParameters();
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::
setCameraOffsetXDeg(camera_size_type icam, param_type offset)
{
  bool change=false;
  UPDATE_AND_NOTE_CHANGE(fDemandCameraOffsetXDeg[icam], offset, change);
  change |= calculateCachedParameters();
  return change;
}

bool VQMultiCameraRenderer::
setCameraOffsetYDeg(camera_size_type icam, param_type offset)
{
  bool change=false;
  UPDATE_AND_NOTE_CHANGE(fDemandCameraOffsetYDeg[icam], offset, change);
  change |= calculateCachedParameters();
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::setZoom(param_type zoom)
{
  if(zoom<0)zoom=0;
  bool change=false;
  UPDATE_AND_NOTE_CHANGE(fDemandZoom, zoom, change);
  if(change)configurationChanged();
  return change;  
}

bool VQMultiCameraRenderer::setPanX(param_type pan_x)
{
  if(pan_x<-1.0)pan_x=-1.0;
  else if(pan_x>1.0)pan_x=1.0;
  bool change=false;
  UPDATE_AND_NOTE_CHANGE(fDemandPanX, pan_x, change);
  if(change)configurationChanged();
  return change;  
}

bool VQMultiCameraRenderer::setPanY(param_type pan_y)
{
  if(pan_y<-1.0)pan_y=-1.0;
  else if(pan_y>1.0)pan_y=1.0;
  bool change=false;
  UPDATE_AND_NOTE_CHANGE(fDemandPanY, pan_y, change);
  if(change)configurationChanged();
  return change;  
}

bool VQMultiCameraRenderer::setUseCommonLimits(bool use_common_limits)
{
  bool change=false;
  UPDATE_AND_NOTE_CHANGE(fDemandUseCommonLimits, use_common_limits, change);
  if(change)fCachedHasCameraLimits=false;
  if(change)configurationChanged();
  return change;  
}

bool VQMultiCameraRenderer::clearCameraLimits(camera_size_type icam)
{
  bool change=false;
  if(fDemandUseCameraLimits[icam])
    {
      fDemandUseCameraLimits[icam]=false;
      fCachedHasCameraLimits=false;
      configurationChanged();
      change=true;
    }
  return change;  
}

bool VQMultiCameraRenderer::
setCameraLimits(camera_size_type icam, 
		const VQMCCameraDataLimits& camera_limits)
{
  bool change=false;
  if(!fDemandUseCameraLimits[icam])
    {
      fDemandUseCameraLimits[icam]=true;
      change=true;
    }

  UPDATE_AND_NOTE_CHANGE(fDemandCameraLimits.fArrayLimits[icam], 
			 camera_limits, change);

  if(change)fCachedHasCameraLimits=false;
  if(change)configurationChanged();
  return change;  
}

bool VQMultiCameraRenderer::clearAllCameraLimits()
{
  bool change = false;
  for(unsigned icam=0;icam!=fCameras.size();icam++)
    change |= clearCameraLimits(icam);
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::
setAllCameraLimits(const VQMCCameraDataLimits& camera_limits)
{
  bool change = false;
  for(unsigned icam=0;icam!=fCameras.size();icam++)
    change |= setCameraLimits(icam, camera_limits);
  if(change)configurationChanged();
  return change;
}

bool VQMultiCameraRenderer::setChannelRenderMode(ChannelRenderMode mode)
{
  if(fDemandCRM != mode)
    {
      fDemandCRM = mode;
      fCachedHasCameraLimits=false;
      configurationChanged();
      return true;
    }
  else return false;
}

bool VQMultiCameraRenderer::setClipToLimits(bool clip)
{
  if(clip != fDemandClipToLimits)
    {
      fDemandClipToLimits = clip;
      configurationChanged();
      return true;
    }
  else return false;
}

bool VQMultiCameraRenderer::setLabelChannels(bool label)
{
  if(fDemandLabelChannels != label)
    {
      fDemandLabelChannels = label;
      return true;
    }
  else return false;
}

bool VQMultiCameraRenderer::setFontFamily(const std::string& family)
{
  if(fDemandFontFamily != family)
    {
      fDemandFontFamily = family;
      fCachedPainterScale = 0;
      return true;
    }
  else return false;
}

bool VQMultiCameraRenderer::setFontScale(double scale)
{
  if(fDemandFontScale != scale)
    {
      fDemandFontScale = scale;
      fCachedPainterScale = 0;
      return true;
    }
  else return false;
}

void VQMultiCameraRenderer::
computeAndCacheLimits(const VQMCArrayData& array_data)
{
  doCacheLimits(array_data);
  fCachedHasCameraLimits=true;
}

void VQMultiCameraRenderer::clearCachedLimits()
{
  fCachedHasCameraLimits=false;
}

void VQMultiCameraRenderer::
renderCamera(QPainter& painter, 
	     const VQMCArrayData& array_data, bool draw_fast)
{
  doRender(painter,array_data,draw_fast,false);
}

void VQMultiCameraRenderer::
renderCamera(QPaintDevice* pd, 
	     const VQMCArrayData& array_data, bool draw_fast)
{
  QPainter painter(pd);
  configurePainter(painter,draw_fast);
  doRender(painter,array_data,draw_fast,true);
}

void VQMultiCameraRenderer::
renderCamera(QPainter& painter, 
	     const VQMCCameraData& camera_data, bool draw_fast)
{
  VQMCArrayData array_data(1);
  array_data.fCameraData[0]=camera_data;
  doRender(painter,array_data,draw_fast,false);
}

void VQMultiCameraRenderer::
renderCamera(QPaintDevice* pd, 
	     const VQMCCameraData& camera_data, bool draw_fast)
{
  VQMCArrayData array_data(1);
  array_data.fCameraData[0]=camera_data;
  renderCamera(pd,array_data,draw_fast);
}

void VQMultiCameraRenderer::doCacheLimits(const VQMCArrayData& array_data)
{
#ifdef TESTTIME
  VSTestTimer L("L");
  L.start();
#endif

  unsigned ncam = fCameras.size();
  for(unsigned icam=0;icam!=ncam;icam++)
    {
      if(!fDemandUseCameraLimits[icam])
	{
	  bool has_camera_data = (array_data.fCameraData.size()>icam);
	  
	  if(has_camera_data)
	    {
	      const VQMCCameraData& camera_data(array_data.fCameraData[icam]);
	      
	      unsigned nchan = camera_data.fChannelData.size();
	      if(!fCameras[icam])
		nchan = 0;
	      else if(fCameras[icam]->numChannels() > nchan)
		nchan = fCameras[icam]->numChannels();
	      
	      if(nchan > 0)
		{
		  fCachedCameraLimits.fArrayLimits[icam].
		    startLimitComputation(camera_data.fChannelData[0]);
		  
		  for(unsigned ichan=1;ichan!=nchan;ichan++)
		    fCachedCameraLimits.fArrayLimits[icam].
		      extendLimits(camera_data.fChannelData[ichan]);
		}
	      else
		{
		  fCachedCameraLimits.fArrayLimits[icam] 
		    = VQMCCameraDataLimits();
		}
	    }
	  else
	    {
	      fCachedCameraLimits.fArrayLimits[icam] 
		= VQMCCameraDataLimits();
	    }
	}
      else
	{
	  fCachedCameraLimits.fArrayLimits[icam] 
	    = fDemandCameraLimits.fArrayLimits[icam];
	}
    }

  if(fDemandUseCommonLimits)
    {
      VQMCCameraDataLimits common_limits;
      common_limits.
	startLimitComputation(fCachedCameraLimits.fArrayLimits[0].fLimitLo);
      common_limits.
	extendLimits(fCachedCameraLimits.fArrayLimits[0].fLimitHi);

      for(unsigned icam=1;icam!=fCameras.size();icam++)
	{
	  common_limits.
	    extendLimits(fCachedCameraLimits.fArrayLimits[icam].fLimitLo);
	  common_limits.
	    extendLimits(fCachedCameraLimits.fArrayLimits[icam].fLimitHi);
	}

      for(unsigned icam=0;icam!=fCameras.size();icam++)
	fCachedCameraLimits.fArrayLimits[icam] = common_limits;
    }

  fCachedRenderCount=0;
  if(fDemandCRM == CRM_MAX_VALUE_PER_CHANNEL)
    {
      for(unsigned icam=0;icam!=ncam;icam++)
	{
	  const VSimpleCameraLayout* cam = fCameras[icam];
	  if(!cam)continue;
      
	  bool has_camera_data = (array_data.fCameraData.size()>icam);
	  
	  unsigned nchan = cam->numChannels();
	  for(unsigned ichan=0;ichan!=nchan;ichan++)
	    {
	      if(fCachedRenderCount==ichan)
		{
		  fCachedRenderOrder[fCachedRenderCount].fICamera = icam;
		  fCachedRenderOrder[fCachedRenderCount].fIChannel = ichan;
		  fCachedRenderOrder[fCachedRenderCount].fValue = 0;
		  fCachedRenderOrder[fCachedRenderCount].fNoData = true;
		  fCachedRenderCount++;
		}
	      
	      if((has_camera_data)&&
		 (array_data.fCameraData[icam].fChannelData.size()>ichan)&&
		 ((fCachedRenderOrder[ichan].fNoData)
		  ||(fCachedRenderOrder[ichan].fValue<
		     array_data.fCameraData[icam].fChannelData[ichan].fValue)))
		{
		  fCachedRenderOrder[ichan].fICamera = icam;
		  fCachedRenderOrder[ichan].fIChannel = ichan;
		  fCachedRenderOrder[ichan].fValue =
		    array_data.fCameraData[icam].fChannelData[ichan].fValue;
		  fCachedRenderOrder[ichan].fNoData = false;
		}
	    }
	}
    }
  else
    {
      for(unsigned icam=0;icam!=ncam;icam++)
	{
	  const VSimpleCameraLayout* cam = fCameras[icam];
	  if(!cam)continue;
	  
	  bool has_camera_data = (array_data.fCameraData.size()>icam);
	  
	  unsigned nchan = cam->numChannels();
	  for(unsigned ichan=0;ichan!=nchan;ichan++)
	    {
	      fCachedRenderOrder[fCachedRenderCount].fICamera = icam;
	      fCachedRenderOrder[fCachedRenderCount].fIChannel = ichan;
	      
	      if((has_camera_data)&&
		 (array_data.fCameraData[icam].fChannelData.size()>ichan))
		{
		  fCachedRenderOrder[fCachedRenderCount].fValue =
		    array_data.fCameraData[icam].fChannelData[ichan].fValue;
		  fCachedRenderOrder[fCachedRenderCount].fNoData = false;
		}
	      else
		{
		  fCachedRenderOrder[fCachedRenderCount].fValue = 0;
		  fCachedRenderOrder[fCachedRenderCount].fNoData = true;
		}
	      
	      fCachedRenderCount++;
	    }
	}
    }

  if(fDemandCRM == CRM_ALL_VALUES_SORTED_ASC)
    std::sort(fCachedRenderOrder.begin(), fCachedRenderOrder.end(),
	      VQMCChannelSortOrderGreaterCompare());
  else if(fDemandCRM == CRM_ALL_VALUES_SORTED_DEC)
    std::sort(fCachedRenderOrder.begin(), fCachedRenderOrder.end(),
	      VQMCChannelSortOrderSmallerCompare());
  
#ifdef TESTTIME
  L.stop();
  std::cout << L << std::endl;
#endif
}

void VQMultiCameraRenderer::
doCacheFontSize(QPainter& painter, const QMatrix& matrix)
{
  double scale;
  if(matrix.m12() == 0)scale=matrix.m11();
  else scale=sqrt(matrix.m11()*matrix.m11()+matrix.m12()*matrix.m12());

  if(scale != fCachedPainterScale)
    {
      fCachedPainterScale = scale;
      fCachedFontPointSize = 
	unsigned(round(scale/fCachedLabelLengthToRadiusRatio
		       /(painter.device()->logicalDpiX()/72)));
      bool first = true;
      bool increase = true;
      do
	{
	  QFont font(fDemandFontFamily.c_str(),fCachedFontPointSize);
	  QFontMetrics metrics(font, painter.device());
	  const double pix20 = metrics.width("88888888888888888888");
	  const double pix = pix20 / 20.0 / fDemandFontScale;
	  const double deg = pix/scale;
	  const double factor = fCachedLabelLengthToRadiusRatio*deg;
#if 0
	  std::cerr << fCachedFontPointSize << ' ' << pix << ' ' << deg 
		    << ' ' << factor << std::endl;
#endif
	  if((factor<1)&&((first)||(increase)))
	    fCachedFontPointSize++,first=false,increase=true;
	  else if((factor>1)&&((first)||(!increase)))
	    fCachedFontPointSize--,first=false,increase=false;
	  else
	    {
	      fCachedFontPointSize = font.pointSize();
	      if((!first)&&(increase))fCachedFontPointSize--;
	      else if((!first)&&(!increase))fCachedFontPointSize++;
	      break;
	    }
	}while(1);
    }
}

void VQMultiCameraRenderer::
configureMatrix(QMatrix& matrix, int width, int height) const
{
  matrix.translate(width/2, height/2);

  qreal painterwidth = qreal(width);
  qreal painterheight = qreal(height);

  qreal scale = 
    qMin(painterwidth/(fCachedExtentR-fCachedExtentL),
	 painterheight/(fCachedExtentU-fCachedExtentD))*fDemandZoom;
  
  qreal xc = (fCachedExtentR+fCachedExtentL)/2.0 
    + fabs(fCachedExtentR-fCachedExtentL-painterwidth/scale)/2.0*fDemandPanX;
  
  qreal yc = (fCachedExtentU+fCachedExtentD)/2.0
    + fabs(fCachedExtentU-fCachedExtentD-painterheight/scale)/2.0*fDemandPanY;
  
  matrix.scale(scale,-scale);
  matrix.translate(-xc,-yc);  
}

void VQMultiCameraRenderer::
configurePainter(QPainter& painter, bool draw_fast)
{
  if(fCameras.size()==0)return;

  QMatrix matrix(painter.matrix());
  configureMatrix(matrix,painter.device()->width(),painter.device()->height());
  painter.setMatrix(matrix);
  
  if(!draw_fast)
    {
      painter.setRenderHint(QPainter::Antialiasing);
      painter.setRenderHint(QPainter::TextAntialiasing);
      painter.setRenderHint(QPainter::SmoothPixmapTransform);
    }
}

void VQMultiCameraRenderer::
doRender(QPainter& painter, const VQMCArrayData& array_data, bool draw_fast,
	 bool i_configured_painter)
{
  if(fCameras.size()==0)return;

#ifdef TESTTIME
  VSTestTimer A("A");
  VSTestTimer B("B");
  VSTestTimer C("C");
  A.start();
#endif
  
  if(!fCachedHasCameraLimits)doCacheLimits(array_data);

  const QMatrix saved_matrix = painter.matrix();

  unsigned ntotal = fCachedRenderCount;
  for(unsigned itotal=0;itotal!=ntotal;itotal++)
    {
      unsigned icam = fCachedRenderOrder[itotal].fICamera;
      unsigned ichan = fCachedRenderOrder[itotal].fIChannel;
      bool has_data = !fCachedRenderOrder[itotal].fNoData;

      const VQMCCameraDataLimits& limit(getCachedCameraLimits(icam));

      qreal channel_rotation = fCachedChannelRotation[icam];
      bool has_channel_rotation = (fabs(channel_rotation) > DBL_EPSILON);

      static const VQMCChannelData no_data(true);
      const VQMCChannelData* data = 
	has_data?&array_data.fCameraData[icam].fChannelData[ichan]:&no_data;
      
      if(!data->fNoDraw)
	{
	  qreal x = fCachedChannelCenterX[icam][ichan];
	  qreal y = fCachedChannelCenterY[icam][ichan];
	  qreal s = fCachedChannelScale[icam][ichan];
	  
#ifdef TESTTIME
	  B.start();
#endif
	  QMatrix matrix(saved_matrix);
	  matrix.translate(x,y);
	  matrix.scale(s,s);
	  if(has_channel_rotation)matrix.rotate(channel_rotation);
	  painter.setMatrix(matrix);
#ifdef TESTTIME
	  B.stop();

	  C.start();
#endif
	  fChannelRenderer->renderChannel(icam, fCameras.size(),
					  *data, limit, fDemandClipToLimits,
					  painter, draw_fast);
#ifdef TESTTIME
	  C.stop();
#endif	  

	}
    }
  
  if(fDemandLabelChannels)
    {
      doCacheFontSize(painter, saved_matrix);
      QFont saved_font = painter.font();
      QFont font(fDemandFontFamily.c_str(),fCachedFontPointSize);
      painter.setFont(font);

      QMatrix matrix(1,0,0,1,saved_matrix.dx(),saved_matrix.dy());
      painter.setMatrix(matrix);

      for(unsigned itotal=0;itotal!=ntotal;itotal++)
	{
	  unsigned icam = fCachedRenderOrder[itotal].fICamera;
	  unsigned ichan = fCachedRenderOrder[itotal].fIChannel;
	  bool has_data = !fCachedRenderOrder[itotal].fNoData;

	  if((has_data)&&
	     (!array_data.fCameraData[icam].fChannelData[ichan].fNoDraw))
	    {
	      qreal X = fCachedChannelCenterX[icam][ichan];
	      qreal Y = fCachedChannelCenterY[icam][ichan];

	      qreal x = saved_matrix.m11()*X + saved_matrix.m21()*Y;
	      qreal y = saved_matrix.m12()*X + saved_matrix.m22()*Y;

	      painter.drawText(QRectF(x,y,0,0),
			       Qt::AlignCenter|Qt::TextDontClip, 
			       fCameras[icam]->channel(ichan).name().c_str());
	    }
	}

      painter.setFont(saved_font);
    }

  painter.setMatrix(saved_matrix);

#ifdef TESTTIME
  A.stop();
  std::cout << A << std::endl;
  std::cout << B << std::endl;
  std::cout << C << std::endl;
#endif
}

void VQMultiCameraRenderer::
findHitChannels(const QRect& bounding_rect, const QPoint& pos,
		std::vector<HitList>& resultsAndHint) const
{
  QMatrix saved_matrix;
  saved_matrix.translate(-bounding_rect.x(),-bounding_rect.y());
  configureMatrix(saved_matrix,bounding_rect.width(),bounding_rect.height());

  qreal widget_x = pos.x();
  qreal widget_y = pos.y();

  qreal fov_x;
  qreal fov_y;
  saved_matrix.inverted().map(widget_x, widget_y, &fov_x, &fov_y);

  unsigned ncam = fCameras.size();
  if(resultsAndHint.size() != ncam)resultsAndHint=std::vector<HitList>(ncam);
  for(unsigned icam=0;icam!=ncam;icam++)
    {
      const VSimpleCameraLayout* cam = fCameras[icam];
      if(!cam){ resultsAndHint[icam].fCameraHit=false; continue; }

      qreal camera_x = fov_x - fDemandCameraOffsetXDeg[icam];
      qreal camera_y = fov_y - fDemandCameraOffsetYDeg[icam];

      if(camera_x*camera_x+camera_y*camera_y > 
	 cam->extentRadial()*cam->extentRadial()*1.05*1.05)
	{
	  resultsAndHint[icam].fCameraHit=false;
	  resultsAndHint[icam].fChannel=0;
	  continue;
	}
      
      qreal channel_rotation = fCachedChannelRotation[icam];
      bool has_channel_rotation = (fabs(channel_rotation) > DBL_EPSILON);
      
      unsigned nchan = cam->numChannels();
      unsigned ichan;
      bool hint;
      if((resultsAndHint[icam].fCameraHit)
	 &&(resultsAndHint[icam].fChannel<nchan))
	ichan=resultsAndHint[icam].fChannel, hint=true;
      else
	ichan=0, hint=false;
      resultsAndHint[icam].fCameraHit=false;
      resultsAndHint[icam].fChannel=0;

      while(ichan!=nchan)
	{
	  qreal x = fCachedChannelCenterX[icam][ichan];
	  qreal y = fCachedChannelCenterY[icam][ichan];
	  qreal s = fCachedChannelScale[icam][ichan];

	  QMatrix matrix(saved_matrix);
	  matrix.translate(x,y);
	  matrix.scale(s,s);
	  if(has_channel_rotation)matrix.rotate(channel_rotation);
	  
	  qreal channel_x;
	  qreal channel_y;

	  matrix.inverted().map(widget_x,widget_y,&channel_x,&channel_y);

	  if(fChannelRenderer->isChannelHit(icam,ncam,channel_x,channel_y))
	    {
	      resultsAndHint[icam].fCameraHit=true;
	      resultsAndHint[icam].fChannel=ichan;
	      break;
	    }

	  if(hint)ichan=0,hint=false;
	  else ichan++;
	}
    }
}

void VQMultiCameraRenderer::configurationChanged()
{
  // nothing to see here
}


// ----------------------------------------------------------------------------
// VQMultiCamera
// ----------------------------------------------------------------------------

VQMultiCamera::~VQMultiCamera()
{
  // nothing to see here
}

bool VQMultiCamera::setChannelRenderer(VQMCChannelRenderer* renderer)
{
  bool change = 
    VQMultiCameraRenderer::setChannelRenderer(renderer);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setCameras(const camera_list& cameras)
{
  bool change = 
    VQMultiCameraRenderer::setCameras(cameras);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setNumCameras(camera_size_type ncam)
{
  bool change = 
    VQMultiCameraRenderer::setNumCameras(ncam);
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::setCamera(camera_size_type icam, camera_type camera)
{
  bool change = 
    VQMultiCameraRenderer::setCamera(icam, camera);
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::setGlobalRotationDeg(param_type angle)
{
  bool change = 
    VQMultiCameraRenderer::setGlobalRotationDeg(angle);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setGlobalCameraRotation(param_type angle)
{
  bool change = 
    VQMultiCameraRenderer::setGlobalCameraRotation(angle);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::
setCameraRotationDeg(camera_size_type icam, param_type angle)
{
  bool change = 
    VQMultiCameraRenderer::setCameraRotationDeg(icam, angle);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::
setChannelRotationDeg(camera_size_type icam, param_type angle)
{
  bool change = 
    VQMultiCameraRenderer::setChannelRotationDeg(icam, angle);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::
setCameraOffsetXDeg(camera_size_type icam, param_type offset)
{
  bool change = 
    VQMultiCameraRenderer::setCameraOffsetXDeg(icam, offset);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::
setCameraOffsetYDeg(camera_size_type icam, param_type offset)
{
  bool change = 
    VQMultiCameraRenderer::setCameraOffsetYDeg(icam, offset);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setAutomaticExtent()
{
  bool change = 
    VQMultiCameraRenderer::setAutomaticExtent();
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setDemandRadialExtentDeg(param_type extent)
{
  bool change = 
    VQMultiCameraRenderer::setDemandRadialExtentDeg(extent);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setZoom(param_type zoom)
{
  bool change = 
    VQMultiCameraRenderer::setZoom(zoom);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setPanX(param_type pan_x)
{
  bool change = 
    VQMultiCameraRenderer::setPanX(pan_x);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setPanY(param_type pan_y)
{
  bool change = 
    VQMultiCameraRenderer::setPanY(pan_y);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setUseCommonLimits(bool use_common_limits)
{
  bool change = 
    VQMultiCameraRenderer::setUseCommonLimits(use_common_limits);
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::clearCameraLimits(camera_size_type icam)
{
  bool change = 
    VQMultiCameraRenderer::clearCameraLimits(icam);
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::
setCameraLimits(camera_size_type icam, 
		const VQMCCameraDataLimits& camera_limits)
{
  bool change = 
    VQMultiCameraRenderer::setCameraLimits(icam, camera_limits);
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::clearAllCameraLimits()
{
  bool change = 
    VQMultiCameraRenderer::clearAllCameraLimits();
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::
setAllCameraLimits(const VQMCCameraDataLimits& camera_limits)
{
  bool change = 
    VQMultiCameraRenderer::setAllCameraLimits(camera_limits);
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::setArrayData(const VQMCArrayData& array_data)
{
  bool change = true;
  fArrayData=array_data;
  computeAndCacheLimits(fArrayData);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setCameraData(camera_size_type icam, 
				  const VQMCCameraData& camera_data)
{
  bool change = true;
  fArrayData.fCameraData[icam]=camera_data;
  computeAndCacheLimits(fArrayData);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setChannelRenderMode(ChannelRenderMode mode)
{
  bool change = 
    VQMultiCameraRenderer::setChannelRenderMode(mode);
  if(change){ computeAndCacheLimits(fArrayData); scheduleUpdate(); }
  return change;
}

bool VQMultiCamera::setClipToLimits(bool clip)
{
  bool change = 
    VQMultiCameraRenderer::setClipToLimits(clip);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setLabelChannels(bool label)
{
  bool change = 
    VQMultiCameraRenderer::setLabelChannels(label);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setFontFamily(const std::string& family)
{
  bool change = 
    VQMultiCameraRenderer::setFontFamily(family);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setFontScale(double scale)
{
  bool change = 
    VQMultiCameraRenderer::setFontScale(scale);
  if(change)scheduleUpdate();
  return change;
}

bool VQMultiCamera::setFastRender(bool fast_render)
{
  bool change = fFastRender != fast_render;
  if(change) { fFastRender = fast_render; scheduleUpdate(); }
  return change;
}

void VQMultiCamera::paintEvent(QPaintEvent* e)
{
  renderCamera(this, fArrayData, fFastRender);
}

void VQMultiCamera::mousePressEvent ( QMouseEvent * e )
{
  std::vector<HitList> hit;
  findHitChannels(rect(), e->pos(), hit);
  for(unsigned icam=0;icam!=hit.size();icam++)
    if(hit[icam].fCameraHit)
      std::cout << icam+1 << ' ' << hit[icam].fChannel+1 << std::endl;
  e->accept();
}

