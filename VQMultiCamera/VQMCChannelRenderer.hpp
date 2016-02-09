//-*-mode:c++; mode:font-lock;-*-

/*! \file VQMCChannelRenderer.hpp

  QT4 class which does the acrtual drawing of the channels in the camera

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/30/2005

  \see VQMultiChannel
*/

#include <cmath>
#include <vector>
#include <iostream>

#include <Qt/QtGui>

#include "VQMCArrayData.hpp"

#ifndef VQMCCHANNELRENDERER_HPP
#define VQMCCHANNELRENDERER_HPP

namespace VERITAS
{

  class VQMCChannelRenderer
  {
  public:
    virtual ~VQMCChannelRenderer();
    virtual QColor renderChannel(unsigned camera, unsigned ncameras, 
				 const VQMCChannelData& data, 
				 const VQMCCameraDataLimits& limits, bool clip,
				 QPainter& painter, bool draw_fast) = 0;
    virtual bool isChannelHit(unsigned camera, unsigned ncameras,
			      qreal x, qreal y) = 0;
    virtual void drawAnnotations(QPainter& painter);
  };
  
  class VQMCSimpleChannelRendererBase: public VQMCChannelRenderer
  {
  public:
    enum ColorMode { CM_BY_CAMERA, CM_BY_TYPE };
    enum ValueMode { VM_RADIUS, VM_AREA, VM_ALPHA, VM_COLOR };
    
    virtual ~VQMCSimpleChannelRendererBase();

    inline bool getChannelBaseColor(QColor& color,
				    unsigned camera,
				    const VQMCChannelData& data, 
				    const VQMCCameraDataLimits& limits);
    
    virtual QColor renderChannel(unsigned camera, unsigned ncameras, 
				 const VQMCChannelData& data, 
				 const VQMCCameraDataLimits& limits, bool clip,
				 QPainter& painter, bool draw_fast) = 0;
    virtual bool isChannelHit(unsigned camera, unsigned ncameras,
			      qreal x, qreal y) = 0;
    
  protected:
    VQMCSimpleChannelRendererBase(const std::vector<QColor>& colors, 
				  ColorMode cmode, ValueMode vmode,
				  bool draw_outline, 
				  const QColor& outline_color):
      fColorMode(cmode), fValueMode(vmode), fColors(colors),
      fDrawOutline(draw_outline), fOutlineColor(outline_color) 
    { /* nothing to see here */ }
    
    void blendColors(qreal x, QColor& c, const QColor& b);
    void blendColors2(qreal x, QColor& c, const QColor& b);
    
    ColorMode             fColorMode;
    ValueMode             fValueMode;
    std::vector<QColor>   fColors;
    bool                  fDrawOutline;
    QColor                fOutlineColor;
  };
  
  template<typename ShapeDraw>
  class VQMCSimpleChannelRenderer: public VQMCSimpleChannelRendererBase
  {
  public:
    VQMCSimpleChannelRenderer(const std::vector<QColor>& colors, 
			      ColorMode cmode = CM_BY_CAMERA, 
			      ValueMode vmode = VM_AREA,
			      bool draw_outline = true,
			      const QColor& outline_color = Qt::black):
      VQMCSimpleChannelRendererBase(colors,cmode,vmode,
				    draw_outline,outline_color)
    { /* nothing to see here */ }
    virtual ~VQMCSimpleChannelRenderer();
    virtual QColor renderChannel(unsigned camera, unsigned ncameras, 
				 const VQMCChannelData& data, 
				 const VQMCCameraDataLimits& limits, bool clip,
				 QPainter& painter, bool draw_fast);
    virtual bool isChannelHit(unsigned camera, unsigned ncameras,
			      qreal x, qreal y);
  };
  
  class VQMCNullDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast)
    { /* nothing to see here */ }
    static bool isHit(unsigned icam, unsigned ncam,qreal x, qreal y)
    { return false; }
  };

  // I have optimized VQMCCircleDraw to be quite fast when "draw_fast" is
  // enabled, however it will not produce accurate results when the
  // coordinate system is sheared, which should never happen unless
  // you are renderering onto a painter that you have independently set
  // to have a sheared world matrix... now why would you do that ??
  
  class VQMCCircleDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { 
      if(draw_fast)
	{
	  const QMatrix& matrix(painter.matrix());
	  // Determinant is product of eigenvalues
	  const int r = int(floor(sqrt(fabs(matrix.det()))));
	  const int x = int(round(matrix.dx()))-r;
	  const int y = int(round(matrix.dy()))-r;
	  const int d = 2*r;
	  painter.setMatrixEnabled(false);
	  painter.drawEllipse(QRect(x,y,d,d));
	  painter.setMatrixEnabled(true);	  
	}
      else
	{
	  painter.drawEllipse(QRectF(-1,-1,2,2));
	} 
    }
    static bool isHit(unsigned icam, unsigned ncam,
		      qreal x, qreal y)
    { return x*x+y*y<=1; }
  };
  
  class VQMCSmallSquareDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { painter.drawRect(QRectF(-1/M_SQRT2,-1/M_SQRT2,2/M_SQRT2,2/M_SQRT2)); }
    static bool isHit(unsigned icam, unsigned ncam,qreal x, qreal y)
    { return (fabs(x)<=1/M_SQRT2)&&(fabs(y)<=1/M_SQRT2); }
  };

  class VQMCSquareDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { painter.drawRect(QRectF(-1,-1,2,2)); }
    static bool isHit(unsigned icam, unsigned ncam,qreal x, qreal y)
    { return (fabs(x)<=1)&&(fabs(y)<=1); }
  };
  
#ifndef M_SQRT3
#define M_SQRT3 1.73205080756887729352
#endif
  
  class VQMCSmallHexDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { static const QPointF hex[6] = { QPointF(0.0,1.0), 
				      QPointF(M_SQRT3/2.0, 0.5),
				      QPointF(M_SQRT3/2.0, -0.5),
				      QPointF(0.0,-1.0), 
				      QPointF(-M_SQRT3/2.0, -0.5),
				      QPointF(-M_SQRT3/2.0, 0.5) };
      painter.drawConvexPolygon(hex,sizeof(hex)/sizeof(*hex)); }
    static bool isHit(unsigned icam, unsigned ncam,qreal x, qreal y)
    { 
      static const double cos60 = 1.0/2.0;
      static const double sin60 = M_SQRT3/2.0;
      double x_pos60 = cos60*x - sin60*y;
      double x_neg60 = cos60*x + sin60*y;
      return 
	(fabs(x)<=M_SQRT3/2.0)
	&&(fabs(x_pos60)<=M_SQRT3/2.0)
	&&(fabs(x_neg60)<=M_SQRT3/2.0); 
    }
  };
  
  class VQMCHexDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { static const QPointF hex[6] = 
	{ QPointF( 0.0,  2.0/M_SQRT3), 
	  QPointF( 1.0,  1.0/M_SQRT3),
	  QPointF( 1.0, -1.0/M_SQRT3),
	  QPointF( 0.0, -2.0/M_SQRT3), 
	  QPointF(-1.0, -1.0/M_SQRT3),
	  QPointF(-1.0,  1.0/M_SQRT3) };
      painter.drawConvexPolygon(hex,sizeof(hex)/sizeof(*hex)); }
    static bool isHit(unsigned icam, unsigned ncam,qreal x, qreal y)
    { 
      static const double cos60 = 1.0/2.0;
      static const double sin60 = M_SQRT3/2.0;
      double x_pos60 = cos60*x - sin60*y;
      double x_neg60 = cos60*x + sin60*y;
      return (fabs(x)<=1)&&(fabs(x_pos60)<=1)&&(fabs(x_neg60)<=1);
    }
  };
  
  class VQMCStarDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { 
      static const QPointF star[12] = 
	{ QPointF(  0.0,      2.0/M_SQRT3),   // OUTER
	  QPointF(  1.0/4.0,  M_SQRT3/4.0 ),  // INNER
	  QPointF(  1.0,      1.0/M_SQRT3),   // OUTER
	  QPointF(  1.0/2.0,  0.0 ),          // INNER
	  QPointF(  1.0,     -1.0/M_SQRT3),   // OUTER
	  QPointF(  1.0/4.0, -M_SQRT3/4.0 ),  // INNER
	  QPointF(  0.0,     -2.0/M_SQRT3),   // OUTER
	  QPointF( -1.0/4.0, -M_SQRT3/4.0 ),  // INNER
	  QPointF( -1.0,     -1.0/M_SQRT3 ),  // OUTER
	  QPointF( -1.0/2.0,  0.0 ),          // INNER
	  QPointF( -1.0,      1.0/M_SQRT3 ),  // OUTER
	  QPointF( -1.0/4.0,  M_SQRT3/4.0 )   // INNER
	};              
      painter.drawConvexPolygon(star,sizeof(star)/sizeof(*star)); 
    }
    static bool isHit(unsigned icam, unsigned ncam,qreal x, qreal y)
    { 
      // Star has 12-fold symmetry which is used to reduce the problem
      // to one of whether the point is inside a triangle. The first
      // two fabs() apply 4-fold x,y symmetry to map the point to the
      // first quadrant. If this mapped point does not have theta<30
      // it is rotated clockwise by 60 degrees and a third fabs() is
      // applied (a 3-fold symmetry operation). Draw a diagram and all
      // will become clear! :-)

      static const double cos60 = 1.0/2.0;
      static const double sin60 = M_SQRT3/2.0;
      x = fabs(x);
      y = fabs(y);
      if(y>x/M_SQRT3)
	{
	  qreal tx=x; 
	  qreal ty=y;
	  x = cos60*tx + sin60*ty;
	  y = fabs(cos60*ty - sin60*tx);
	}
      return (x<=1)&&(y>=(x-0.5)*(2.0/M_SQRT3));
    }
  };

#if 0
  class VQMCShamrockDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { 
      static const QPointF shamrock[] = 
	{ QPointF(  0.0,      2.0/M_SQRT3),   // OUTER
	  QPointF(  1.0/4.0,  M_SQRT3/4.0 ),  // INNER
	  QPointF(  1.0,      1.0/M_SQRT3),   // OUTER
	  QPointF(  1.0/2.0,  0.0 ),          // INNER
	  QPointF(  1.0,     -1.0/M_SQRT3),   // OUTER
	  QPointF(  1.0/4.0, -M_SQRT3/4.0 ),  // INNER
	  QPointF(  0.0,     -2.0/M_SQRT3),   // OUTER
	  QPointF( -1.0/4.0, -M_SQRT3/4.0 ),  // INNER
	  QPointF( -1.0,     -1.0/M_SQRT3 ),  // OUTER
	  QPointF( -1.0/2.0,  0.0 ),          // INNER
	  QPointF( -1.0,      1.0/M_SQRT3 ),  // OUTER
	  QPointF( -1.0/4.0,  M_SQRT3/4.0 )   // INNER
	};              
      painter.drawConvexPolygon(star,sizeof(shamrock)/sizeof(*shamrock)); 
    }
    static bool isHit(unsigned icam, unsigned ncam,qreal x, qreal y)
    { 
      return false;
    }
  };  
#endif

  class VQMCPieDraw
  {
  public:
    static void draw(unsigned icam, unsigned ncam,
		     QPainter& painter, bool draw_fast) 
    { 
      int a0 = int(round(double(icam)/double(ncam)*5760));
      int a1 = int(round(double(icam+1)/double(ncam)*5760));
      if(painter.brush().style() == Qt::NoBrush)
	{
	  painter.drawArc(QRectF(-1,-1,2,2),a0,a1-a0);
	}
      else if(painter.pen().color() != painter.brush().color())
	{
	  QPen pen = painter.pen();
	  painter.setPen(QPen(painter.brush().color()));
	  painter.drawPie(QRectF(-1,-1,2,2),a0,a1-a0);
	  painter.setPen(pen);
	  painter.drawArc(QRectF(-1,-1,2,2),a0,a1-a0);
	}
      else
	{
	  painter.drawPie(QRectF(-1,-1,2,2),a0,a1-a0);
	}
    }
    static bool isHit(unsigned icam, unsigned ncam,
		      qreal x, qreal y)
    { if(x*x+y*y>1)return false;
      double theta = atan2(-y,x)/(2.0*M_PI);
      if(theta<0)theta+=1.0;
      double tcam = theta*double(ncam);
      return((tcam>=double(icam))&&(tcam<double(icam+1))); }
  };

  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //
  // (TEMPLATE) FUNCTION DEFINITIONS
  //
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  
  inline bool VQMCSimpleChannelRendererBase::
  getChannelBaseColor(QColor& color,
		      unsigned camera,
		      const VQMCChannelData& data, 
		      const VQMCCameraDataLimits& limits)
  {
    if(fColors.size() == 0)return false;

    if(fValueMode == VM_COLOR)
      {
	double y = (data.fValue - limits.fLimitLo.fValue)/
	  (limits.fLimitHi.fValue - limits.fLimitLo.fValue);
	if(y>1)y=1;
	else if(y<0)y=0;
	double z = y*double(fColors.size()-1);
	unsigned index = unsigned(floor(z));
	double d = fmod(z,1);
	color = fColors[index];
	if(index != fColors.size()-1)
	  blendColors2(1-d,color,fColors[index+1]);
      }
    else
      {
	unsigned color_index = 0;
	switch(fColorMode)
	  {
	  case CM_BY_CAMERA:
	    color_index = camera;
	    break;
	  case CM_BY_TYPE:
	    color_index = data.fType;
	    break;
	  }    
	color_index = color_index%fColors.size();
	color = fColors[color_index];
      }

    return true;
  }
  
  template<typename ShapeDraw>
  VQMCSimpleChannelRenderer<ShapeDraw>::~VQMCSimpleChannelRenderer()
  {
    // nothing to see here
  }
  
  template<typename ShapeDraw>
  QColor VQMCSimpleChannelRenderer<ShapeDraw>::
  renderChannel(unsigned camera, unsigned ncameras, 
		const VQMCChannelData& data, 
		const VQMCCameraDataLimits& limits,
		bool clip, QPainter& painter, bool draw_fast)
  {
    QColor background_color = painter.background().color();
    QColor return_color = background_color;
    if(data.fNoDraw)return return_color;
    
    if((!data.fNoData)&&(!data.fSuppress))
      {
	double lim = qMax(fabs(limits.fLimitLo.fValue),
			  fabs(limits.fLimitHi.fValue));
	double x = data.fValue/lim;
	
	if(clip)
	  {
	    if(x>1)x=1;
	    else if(x<-1)x=-1;
	  }

	double radius = 0;

	QColor color = background_color;
	bool found_color = getChannelBaseColor(color, camera, data, limits);
	if((draw_fast)&&(color.alphaF()<1.0))
	  blendColors(1.0-color.alphaF(),color,background_color);
	QBrush brush(color,Qt::SolidPattern);
	
	switch(fValueMode)
	  {
	  case VM_RADIUS:
	    radius = fabs(x);
	    if(x<0)brush.setStyle(Qt::NoBrush);
	    else return_color = color;
	    break;
	  case VM_AREA:
	    radius = sqrt(fabs(x));
	    if(x<0)brush.setStyle(Qt::NoBrush);
	    else return_color = color;
	    break;
	  case VM_ALPHA:
	    if(x>=1)
	      {
		return_color = color;
		radius = 1;
	      }
	    else if(x>0)
	      {
		if(draw_fast)blendColors(x,color,background_color);
		else color.setAlphaF(color.alphaF()*x);
		brush.setColor(color);
		return_color = color;
		radius = 1;
	      }
	    else
	      {
		radius = 0;
	      }
	    break;
	  case VM_COLOR:
	    radius = 1;
	    return_color = color;
	    break;
	  }

	if((found_color)&&(radius>0))
	  {
	    QPen saved_pen = painter.pen();
	    QBrush saved_brush = painter.brush();
	    
	    QPen pen(saved_pen);
	    pen.setWidth(0);
	    if((radius<=1.0)||(!fDrawOutline))pen.setColor(color);
	    else pen.setColor(fOutlineColor);
	    painter.setPen(pen);
	    painter.setBrush(brush);

	    painter.scale(radius,radius);
	    ShapeDraw::draw(camera,ncameras,painter,draw_fast);
	    
	    if((fDrawOutline)&&(radius<=1.0))
	      {
		const qreal rescale = 1.0/radius;
		painter.scale(rescale,rescale);
		pen.setColor(fOutlineColor);
		painter.setPen(pen);
		painter.setBrush(QBrush());
		ShapeDraw::draw(camera,ncameras,painter,draw_fast);
	      }

	    painter.setBrush(saved_brush);
	    painter.setPen(saved_pen);
	  }
	else if(fDrawOutline)
	  {
	    QPen saved_pen = painter.pen();
	    QPen pen = painter.pen();
	    pen.setWidth(0);
	    pen.setColor(fOutlineColor);
	    painter.setPen(pen);
	    ShapeDraw::draw(camera,ncameras,painter,draw_fast);
	    painter.setPen(saved_pen);
	  }
      }
    else if(fDrawOutline)
      {
	QPen saved_pen = painter.pen();
	QPen pen = painter.pen();
	pen.setWidth(0);
	pen.setColor(fOutlineColor);
	painter.setPen(pen);
	ShapeDraw::draw(camera,ncameras,painter,draw_fast);
	painter.setPen(saved_pen);
      }

    return return_color;
  }

  template<typename ShapeDraw>
  bool VQMCSimpleChannelRenderer<ShapeDraw>::
  isChannelHit(unsigned camera, unsigned ncameras, qreal x, qreal y)
  {
    return ShapeDraw::isHit(camera,ncameras,x,y);
  }

} // namespace VERITAS

#endif // VQMCCHANNELRENDERER_HPP
