//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRCANVAS_HPP
#define VSRCANVAS_HPP

#include <vector>
#include <map>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSOptions.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TPad.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram1D.hpp"
#include "VSRHistogram2D.hpp"
#include "VSRGraph1D.hpp"
#include "VSRHistogramCamera.hpp"
#include "VSRText.hpp"
#include "VSRCanvasElement.hpp"

class VSRCanvasPad
{
public:
  struct Options
  {
    Options();

    unsigned     lineWidth;
    std::vector< unsigned > markerStyle;
    std::vector< unsigned > lineStyle;
    
    bool         logX;
    bool         logY;
    bool         logZ;
    bool         drawGrid;
    bool         showErrors;
    bool         showChannelLabels;

    std::string  drawOptions;
    std::string  legendStyle;
    std::string  colorPalette;

    std::vector< std::string > binLabels;
  };

  VSRCanvasPad(const Options& opt = s_default_options);
  ~VSRCanvasPad();

  void draw();
  void add(VSRCanvasElement* element, const std::string& label = "");


  unsigned nElements() const { return m_elements.size(); }

  VSRCanvasElement* getElement(unsigned iel);

  std::string getTitle() { return m_title; }

  double getLoX() const { return m_loX; }
  double getHiX() const { return m_hiX; }
  double getLoY() const { return m_loY; }
  double getHiY() const { return m_hiY; }
  double getLoZ() const { return m_loZ; }
  double getHiZ() const { return m_hiZ; }

  double getPosLoX() const { return m_posLoX; }
  double getPosHiX() const { return m_posHiX; }
  double getPosLoY() const { return m_posLoY; }
  double getPosHiY() const { return m_posHiY; }

  bool getLogX() { return m_options.logX; }
  bool getLogY() { return m_options.logY; }

  bool is1D() { return m_is1D; }
  bool is2D() { return m_is2D; }
  
  const std::vector< std::string >& getBinLabels() { return m_binLabels; }

  // Setters ------------------------------------------------------------------
  void setTitle(const std::string& title) { m_title = title; }

  void setLoX(double loX) { m_loX = loX; }
  void setHiX(double hiX) { m_hiX = hiX; }
  void setLoY(double loY) { m_loY = loY; }
  void setHiY(double hiY) { m_hiY = hiY; }
  void setLoZ(double loZ) { m_loZ = loZ; }
  void setHiZ(double hiZ) { m_hiZ = hiZ; }

  void setLogX(bool logX = true) { m_options.logX = logX; }
  void setLinX(bool linX = true) { m_options.logX = !linX; }
  void setLogY(bool logY = true) { m_options.logY = logY; }
  void setLinY(bool linY = true) { m_options.logY = !linY; }
  void setGrid(bool s) { m_options.drawGrid = s; }

  void setBinLabels(const std::vector< std::string >& labels)
  {
    m_binLabels = labels;
  }

  // Static Methods -----------------------------------------------------------
  static void configure(VERITAS::VSOptions& options);
  static void setDefaultDrawOptions(const std::string& draw_options)
  { s_default_options.drawOptions = draw_options; }
  static void setDefaultShowErrors(bool show_errors)
  { s_default_options.showErrors = show_errors; }
  static void setDefaultLineWidth(unsigned line_width)
  { s_default_options.lineWidth = line_width; }
  static void setDefaultLegendStyle(const std::string& legend_style)
  { s_default_options.legendStyle = legend_style; }

private:

  double                           m_loX;
  double                           m_hiX;
  double                           m_loY;
  double                           m_hiY;
  double                           m_loZ;
  double                           m_hiZ;
  double                           m_posLoX;
  double                           m_posHiX;
  double                           m_posLoY;
  double                           m_posHiY;
  
  bool                             m_is1D;
  bool                             m_is2D;

  std::string                      m_title;
  std::vector< VSRCanvasElement* > m_elements;
  std::vector< std::string >       m_labels;
  std::vector< std::string >       m_binLabels;
  TLegend*                         m_legend;

  static Options s_default_options;
  Options        m_options;
};

class VSRCanvas 
{
public:
  class Options
  {
  public:
    Options();

    std::pair< unsigned, unsigned > canvas_size;
  };

  VSRCanvas(const std::string& title,
	    const Options& opt = s_default_options);

  ~VSRCanvas();

  void add(unsigned ipad, VSRCanvasElement* element,
	   const std::string& label = "");

  void setTitle(unsigned ipad, const std::string& title);
  void setDimensions(unsigned nxpix, unsigned nypix)
  {
    m_npixelsX = nxpix;
    m_npixelsY = nypix;
  }

  void setLimitsX(unsigned ipad, double lo, double hi);
  void setLimitsY(unsigned ipad, double lo, double hi);
  void setLimitsZ(unsigned ipad, double lo, double hi);

  void setBinLabels(unsigned ipad, const std::vector< std::string >& labels);

  void setLogX(unsigned ipad, bool logX = true);
  void setLinX(unsigned ipad, bool linX = true);
  void setLogY(unsigned ipad, bool logY = true);
  void setLinY(unsigned ipad, bool linY = true);
  void setGrid(unsigned ipad, bool grid = true);

  void draw();    
  VSRCanvasElement* getElement(unsigned ipad, unsigned iel);

  static void configure(VERITAS::VSOptions& options);

private:

  unsigned    m_npixelsX;
  unsigned    m_npixelsY;
  std::string m_title;

  std::map< unsigned, VSRCanvasPad* > m_pads;

  TCanvas* m_canvas;

  static Options s_default_options;
  Options        m_options;

};


#endif // VSRCANVAS_HPP
