//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRCANVASELEMENT_HPP
#define VSRCANVASELEMENT_HPP

#include <iostream>
#include <vector>

#include <TLegend.h>

class VSRCanvasElement
{
public:

  struct Options
  {
    Options(): 
      loZ(), hiZ(),
      lineWidth(1), 
      lineColor(1), lineStyle(1), 
      fillColor(),
      markerColor(1), markerStyle(21),
      drawOptions(), legendStyle("lp"),
      colorPalette(),
      showErrors(true), showChannelLabels(false)
    { }

    double                     loZ;
    double                     hiZ;
    unsigned                   lineWidth;
    unsigned                   lineColor;
    unsigned                   lineStyle;
    unsigned                   fillColor;
    unsigned                   markerColor;
    unsigned                   markerStyle;
    std::string                drawOptions;
    std::string                legendStyle;
    std::string                colorPalette;
    bool                       showErrors;
    bool                       showChannelLabels;
  };

  VSRCanvasElement(const Options& options = s_default_options);
  virtual ~VSRCanvasElement() {}

  virtual void dump() { }
  virtual void draw() = 0;
  virtual VSRCanvasElement* clone() const = 0;

  virtual void addLegendEntry(TLegend *leg, 
			      const std::string& label) const = 0;
  
  virtual double getLoX() const = 0;
  virtual double getHiX() const = 0;
  virtual double getLoY() const = 0;
  virtual double getHiY() const = 0;
  virtual double getPosLoX() const = 0;
  virtual double getPosHiX() const = 0;
  virtual double getPosLoY() const = 0;
  virtual double getPosHiY() const = 0;

  virtual bool is1D() const = 0;
  virtual bool is2D() const = 0;

  virtual void multiply(double x) { }

  virtual VSRCanvasElement& operator*= (double x)
  {
    return *this;
  }

  
  const Options& getOptions() { return m_options; }

  void setShowErrors(bool showErrors) 
  { 
    m_options.showErrors = false; 
  }

  void setDrawOptions(const std::string& drawOptions)
  {
    m_options.drawOptions = drawOptions;
  }

  void setLoZ(double loZ) { m_options.loZ = loZ; }
  void setHiZ(double hiZ) { m_options.hiZ = hiZ; }

  void setLineWidth(unsigned lineWidth) { m_options.lineWidth = lineWidth; }
  void setMarkerStyle(unsigned markerStyle) 
  { m_options.markerStyle = markerStyle; }
  void setFillColor(unsigned fillColor) { m_options.fillColor = fillColor; }
  void setOptions(const Options& options) { m_options = options; }

protected:

  Options m_options;

  static Options s_default_options;
};

#endif // VSRCANVASELEMENT_HPP
