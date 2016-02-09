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

  \see VQMCArrayData
  \see VQMCChannelRenderer
*/

#include <cmath>
#include <vector>
#include <iostream>

#include <Qt/QtGui>

#include "VSimpleCameraLayout.hpp"
#include "VQMCArrayData.hpp"
#include "VQMCChannelRenderer.hpp"

#ifndef VQMULTICAMERA_HPP
#define VQMULTICAMERA_HPP

namespace VERITAS
{

  class VQMultiCameraRenderer
  {
  public:
    typedef VSimpleCameraLayout*                          camera_type;
    typedef std::vector<camera_type>                      camera_list;
    typedef camera_list::size_type                        camera_size_type;
    typedef double                                        param_type;
    typedef std::vector<param_type>                       param_list;

    class HitList
    {
    public:
      HitList(): fCameraHit(false), fChannel() { }
      HitList(unsigned channel): fCameraHit(true), fChannel(channel) { }
      bool       fCameraHit;
      unsigned   fChannel;
    };
    
    enum ChannelRenderMode { 
      CRM_ALL_VALUES_BY_CAMERA_AND_CHANNEL, 
      CRM_ALL_VALUES_SORTED_ASC, 
      CRM_ALL_VALUES_SORTED_DEC, 
      CRM_MAX_VALUE_PER_CHANNEL };
    
    VQMultiCameraRenderer(VQMCChannelRenderer* channel_renderer,
			  const camera_list& cameras);
    virtual ~VQMultiCameraRenderer();
    
    // GETTERS

    VQMCChannelRenderer* getChannelRenderer() const 
    { return fChannelRenderer; }

    camera_type getCameraType(camera_size_type icam) const 
    { return fCameras[icam]; }

    param_type getGlobalRotationDeg() const 
    { return fDemandGlobalRotationRad/M_PI*180.0; }
    param_type getGlobalCameraRotationDeg() const
    { return fDemandGlobalCameraRotationRad/M_PI*180.0; }
    param_type getCameraRotationDeg(camera_size_type icam) const
    { return fDemandCameraRotationRad[icam]/M_PI*180.0; }
    param_type getChannelRotationDeg(camera_size_type icam) const
    { return fDemandChannelRotationRad[icam]/M_PI*180.0; }

    param_type getCameraOffsetXDeg(camera_size_type icam) const
    { return fDemandCameraOffsetXDeg[icam]; }
    param_type getCameraOffsetYDeg(camera_size_type icam) const
    { return fDemandCameraOffsetYDeg[icam]; }

    param_type getDemandRadialExtentDeg() const 
    { return fDemandRadialExtentDeg; }
    param_type getZoom() const { return fDemandZoom; }
    param_type getPanX() const { return fDemandPanX; }
    param_type getPanY() const { return fDemandPanY; }
    
    param_type getExtentLDeg() const { return fCachedExtentL; }
    param_type getExtentRDeg() const { return fCachedExtentR; }
    param_type getExtentDDeg() const { return fCachedExtentD; }
    param_type getExtentUDeg() const { return fCachedExtentU; }

    inline param_type getActualCameraRotationDeg(camera_size_type icam) const;
    inline param_type getActualChannelRotationDeg(camera_size_type icam) const;
    inline param_type getActualCameraOffsetXDeg(camera_size_type icam) const;
    inline param_type getActualCameraOffsetYDeg(camera_size_type icam) const;

    bool getUseCommonLimits() const { return fDemandUseCommonLimits; }
    bool getUseCameraLimits(camera_size_type icam) const
    { return fDemandUseCameraLimits[icam]; }
    const VQMCCameraDataLimits& getCameraLimits(camera_size_type icam) const
    { return fDemandCameraLimits.fArrayLimits[icam]; }

    const VQMCCameraDataLimits& getCachedCameraLimits(camera_size_type icam) 
      const { return fCachedCameraLimits.fArrayLimits[icam]; }

    ChannelRenderMode getChannelRenderMode() const { return fDemandCRM; }
    bool getClipToLimits() const { return fDemandClipToLimits; }

    bool getLabelChannels() const { return fDemandLabelChannels; }
    const std::string& getFontFamily() const { return fDemandFontFamily; }
    param_type getFontScale() const { return fDemandFontScale; }

    // SETTERS

    bool setChannelRenderer(VQMCChannelRenderer* renderer);
    bool setCameras(const camera_list& cameras);
    bool setNumCameras(camera_size_type ncam);
    bool setCamera(camera_size_type icam, camera_type camera);

    bool setGlobalRotationDeg(param_type angle);
    bool setGlobalCameraRotation(param_type angle);
    bool setCameraRotationDeg(camera_size_type icam, param_type angle);
    bool setChannelRotationDeg(camera_size_type icam, param_type angle);

    bool setCameraOffsetXDeg(camera_size_type icam, param_type offset);
    bool setCameraOffsetYDeg(camera_size_type icam, param_type offset);

    bool setAutomaticExtent();
    bool setDemandRadialExtentDeg(param_type extent);

    bool setZoom(param_type zoom);
    bool setPanX(param_type pan_x);
    bool setPanY(param_type pan_y);

    bool setUseCommonLimits(bool use_common_limits);
    bool clearCameraLimits(camera_size_type icam);
    bool setCameraLimits(camera_size_type icam, 
			 const VQMCCameraDataLimits& camera_limits);
    bool clearAllCameraLimits();
    bool setAllCameraLimits(const VQMCCameraDataLimits& camera_limits);

    void computeAndCacheLimits(const VQMCArrayData& array_data);
    void clearCachedLimits();

    bool setChannelRenderMode(ChannelRenderMode mode);

    bool setClipToLimits(bool clip);
    
    bool setLabelChannels(bool label);
    bool setFontFamily(const std::string& family);
    bool setFontScale(double scale);

    // RENDER FUNCTIONS

    void configurePainter(QPainter& painter, bool draw_fast = false);

    void renderCamera(QPainter& painter, const VQMCArrayData& array_data,
		      bool draw_fast = false);
    void renderCamera(QPaintDevice* pd, const VQMCArrayData& array_data,
		      bool draw_fast = false);

    void renderCamera(QPainter& painter, const VQMCCameraData& camera_data,
		      bool draw_fast = false);
    void renderCamera(QPaintDevice* pd, const VQMCCameraData& camera_data,
		      bool draw_fast = false);

    // CHANNEL FIND FUNCTIONS

    void findHitChannels(const QRect& bounding_rect, const QPoint& pos,
			 std::vector<HitList>& resultsAndHint) const;

  protected:
    typedef std::vector<param_list>                       array_param_list;
    
    class VQMCChannelSortOrder
    {
    public:
      VQMCChannelSortOrder(): fValue(), fNoData(), fICamera(), fIChannel() { }
      double     fValue;
      bool       fNoData;
      unsigned   fICamera;
      unsigned   fIChannel;
      bool operator < (const VQMCChannelSortOrder& o) const 
      { return fNoData?!o.fNoData:o.fValue<fValue; }
      bool operator > (const VQMCChannelSortOrder& o) const 
      { return fNoData?!o.fNoData:o.fValue>fValue; }
    };

    class VQMCChannelSortOrderGreaterCompare
    {
    public:
      int operator()(const VQMCChannelSortOrder& x1,
		     const VQMCChannelSortOrder& x2)
      { return x1>x2; }
    };

    class VQMCChannelSortOrderSmallerCompare
    {
    public:
      int operator()(const VQMCChannelSortOrder& x1,
		     const VQMCChannelSortOrder& x2)
      { return x1<x2; }
    };

    void configureMatrix(QMatrix& matrix, int width, int height) const;

    virtual void configurationChanged();

    void doCacheFontSize(QPainter& painter, const QMatrix& matrix);

    void doCacheLimits(const VQMCArrayData& array_data);

    void doRender(QPainter& painter, const VQMCArrayData& array_data,
		  bool draw_fast, bool i_configured_painter);

    bool calculateCachedParameters();
    
    VQMCChannelRenderer*  fChannelRenderer;

    camera_list           fCameras;

    // Settable parameters which control how cameras appear

    param_type            fDemandGlobalRotationRad;
    param_type            fDemandGlobalCameraRotationRad;
    param_list            fDemandCameraRotationRad;
    param_list            fDemandChannelRotationRad;
    param_list            fDemandCameraOffsetXDeg;
    param_list            fDemandCameraOffsetYDeg;
    param_type            fDemandRadialExtentDeg;
    param_type            fDemandZoom;
    param_type            fDemandPanX;
    param_type            fDemandPanY;
    bool                  fDemandUseCommonLimits;
    std::vector<bool>     fDemandUseCameraLimits;
    VQMCArrayDataLimits   fDemandCameraLimits;
    ChannelRenderMode     fDemandCRM;
    bool                  fDemandClipToLimits;
    bool                  fDemandLabelChannels;
    std::string           fDemandFontFamily;
    param_type            fDemandFontScale;

    // Cached parameters to speed rendering

    param_type            fCachedExtentL;
    param_type            fCachedExtentR;
    param_type            fCachedExtentD;
    param_type            fCachedExtentU;
    param_type            fCachedGlobalRotationCos;
    param_type            fCachedGlobalRotationSin;
    param_list            fCachedCameraRotationCos;
    param_list            fCachedCameraRotationSin;
    param_list            fCachedChannelRotationCos;
    param_list            fCachedChannelRotationSin;
    array_param_list      fCachedChannelCenterX;
    array_param_list      fCachedChannelCenterY;
    array_param_list      fCachedChannelScale;
    param_list            fCachedChannelRotation;
    bool                  fCachedHasCameraLimits;
    VQMCArrayDataLimits   fCachedCameraLimits;
    unsigned              fCachedFontPointSize;
    param_type            fCachedPainterScale;
    param_type            fCachedLabelLengthToRadiusRatio;

    std::vector<VQMCChannelSortOrder> fCachedRenderOrder;
    unsigned                          fCachedRenderCount;
  };

  class VQMultiCamera: public QWidget, public VQMultiCameraRenderer
  {
    Q_OBJECT
  public:
    VQMultiCamera(VQMCChannelRenderer* channel_renderer,
		  const camera_list& cameras, QWidget* parent = 0):
      QWidget(parent),
      VQMultiCameraRenderer(channel_renderer, cameras),
      fArrayData(), fFastRender(false) { }
    virtual ~VQMultiCamera();

    // OVERRIDE VQMultiCameraRenderer SETTER FUNCTIONS

    bool setChannelRenderer(VQMCChannelRenderer* renderer);
    bool setCameras(const camera_list& cameras);
    bool setNumCameras(camera_size_type ncam);
    bool setCamera(camera_size_type icam, camera_type camera);

    bool setGlobalRotationDeg(param_type angle);
    bool setGlobalCameraRotation(param_type angle);
    bool setCameraRotationDeg(camera_size_type icam, param_type angle);
    bool setChannelRotationDeg(camera_size_type icam, param_type angle);

    bool setCameraOffsetXDeg(camera_size_type icam, param_type offset);
    bool setCameraOffsetYDeg(camera_size_type icam, param_type offset);

    bool setAutomaticExtent();
    bool setDemandRadialExtentDeg(param_type extent);

    bool setZoom(param_type zoom);
    bool setPanX(param_type pan_x);
    bool setPanY(param_type pan_y);

    bool setUseCommonLimits(bool use_common_limits);
    bool clearCameraLimits(camera_size_type icam);
    bool setCameraLimits(camera_size_type icam, 
			 const VQMCCameraDataLimits& camera_limits);
    bool clearAllCameraLimits();
    bool setAllCameraLimits(const VQMCCameraDataLimits& camera_limits);

    bool setArrayData(const VQMCArrayData& array_data);
    bool setCameraData(camera_size_type icam, 
		       const VQMCCameraData& camera_data);
 
    bool setChannelRenderMode(ChannelRenderMode mode);

    bool setClipToLimits(bool clip);

    bool setLabelChannels(bool label);
    bool setFontFamily(const std::string& family);
    bool setFontScale(double scale);

    bool setFastRender(bool fast_render);

    // QOBJECT INTERFACE

    void scheduleUpdate() { update(); }

    virtual void paintEvent(QPaintEvent* e);
    virtual void mousePressEvent ( QMouseEvent * e );

  protected:
    VQMCArrayData      fArrayData;
    bool               fFastRender;
  };  

  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //
  // FUNCTION DEFINITIONS
  //
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------

  inline VQMultiCameraRenderer::param_type VQMultiCameraRenderer::
  getActualCameraRotationDeg(VQMultiCameraRenderer::camera_size_type icam)
    const
  {
    param_type rot = 
      fDemandGlobalRotationRad
      +fDemandGlobalCameraRotationRad
      +fDemandCameraRotationRad[icam];
    rot = fmod(fmod(rot/M_PI*180+360,360)+360,360);
    if(rot>=180)rot-=360;
    return rot;
  }

  inline VQMultiCameraRenderer::param_type VQMultiCameraRenderer::
  getActualChannelRotationDeg(VQMultiCameraRenderer::camera_size_type icam)
    const
  {
    param_type rot = 
      fDemandGlobalRotationRad
      +fDemandGlobalCameraRotationRad
      +fDemandCameraRotationRad[icam]
      +fDemandChannelRotationRad[icam];
    rot = fmod(fmod(rot/M_PI*180+360,360)+360,360);
    if(rot>=180)rot-=360;
    return rot;
  }

  inline VQMultiCameraRenderer::param_type VQMultiCameraRenderer::
  getActualCameraOffsetXDeg(camera_size_type icam) const
  {
    return 
      fDemandCameraOffsetXDeg[icam]*fCachedGlobalRotationCos
      -fDemandCameraOffsetYDeg[icam]*fCachedGlobalRotationSin;
  }

  inline VQMultiCameraRenderer::param_type VQMultiCameraRenderer::
  getActualCameraOffsetYDeg(camera_size_type icam) const
  {
    return 
      fDemandCameraOffsetYDeg[icam]*fCachedGlobalRotationCos
      +fDemandCameraOffsetXDeg[icam]*fCachedGlobalRotationSin;
  }
  
}

#endif // VQMULTICAMERA_HPP
