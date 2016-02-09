#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cassert>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <WhippleCams.h>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TEllipse.h>
#include <TPaletteAxis.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TBox.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TMath.h>
#include <THStack.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TPaveLabel.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRCanvas.hpp"

using namespace VERITAS;

// ============================================================================
// VSRCanvasPad
// ============================================================================

VSRCanvasPad::Options VSRCanvasPad::s_default_options;

VSRCanvasPad::Options::Options():
  lineWidth(1),
  logX(false),
  logY(false),
  logZ(false),
  drawGrid(true),
  showErrors(true),
  showChannelLabels(false),
  drawOptions(),
  legendStyle("lp"),
  colorPalette("default")
{

}

VSRCanvasPad::VSRCanvasPad(const Options& opt): 
  m_loX(), m_hiX(), m_loY(), m_hiY(), m_loZ(), m_hiZ(),
  m_posLoX(), m_posHiX(), m_posLoY(), m_posHiY(), m_is1D(), m_is2D(),
  m_title(), m_elements(), m_options(opt)
{
  m_legend = new TLegend(0.6,0.6,0.8,0.8);
}

VSRCanvasPad::~VSRCanvasPad()
{
  const unsigned nelements = m_elements.size();
  for(unsigned ielement = 0; ielement < nelements; ielement++)
    delete m_elements[ielement];

  delete m_legend;
}

void VSRCanvasPad::draw()
{  
  bool draw_legend = false;
  for(unsigned ilabel = 0; ilabel < m_labels.size(); ilabel++)
    if(!m_labels[ilabel].empty()) draw_legend = true;

  const unsigned nelements = m_elements.size();
  for(unsigned ielement = 0; ielement < nelements; ielement++)
    {
      if(!m_labels[ielement].empty())
	m_elements[ielement]->addLegendEntry(m_legend,m_labels[ielement]);

      m_elements[ielement]->setLoZ(m_loZ);
      m_elements[ielement]->setHiZ(m_hiZ);

      m_elements[ielement]->draw();
    }

  if(m_options.logX) gPad->SetLogx();
  if(m_options.logY) gPad->SetLogy();
  if(m_options.logZ) gPad->SetLogz();
  if(m_options.drawGrid) gPad->SetGrid();
  if(draw_legend) m_legend->Draw();
}

void VSRCanvasPad::add(VSRCanvasElement* element, const std::string& label)
{
  if(m_elements.size() == 0)
    {
      m_is2D = element->is2D();
      m_is1D = element->is1D();

      m_loX = element->getLoX();
      m_hiX = element->getHiX();
      m_loY = element->getLoY();
      m_hiY = element->getHiY();

      if(element->getPosLoX()) m_posLoX = element->getPosLoX();
      if(element->getPosHiX()) m_posHiX = element->getPosHiX();
      if(element->getPosLoY()) m_posLoY = element->getPosLoY();
      if(element->getPosHiY()) m_posHiY = element->getPosHiY();
    }
  else 
    {
      //      vsassert(element->is2D() == m_is2D);

      m_loX = std::min(m_loX,element->getLoX());
      m_hiX = std::max(m_hiX,element->getHiX());
      m_loY = std::min(m_loY,element->getLoY());
      m_hiY = std::max(m_hiY,element->getHiY());

      if(element->getLoX()>0 && m_posLoX) 
	m_posLoX = std::min(m_posLoX,element->getLoX());
      else if(element->getLoX()>0) m_posLoX = element->getLoX();

      if(element->getHiX()>0 && m_posHiX) 
	m_posHiX = std::max(m_posHiX,element->getHiX());
      else if(element->getHiX()>0)m_posHiX = element->getHiX();

      if(element->getLoY()>0 && m_posLoY) 
	m_posLoY = std::min(m_posLoY,element->getLoY());
      else if(element->getLoY()>0) m_posLoY = element->getLoY();

      if(element->getHiY()>0 && m_posHiY) 
	m_posHiY = std::max(m_posHiY,element->getHiY());
      else if(element->getHiY()>0)m_posHiY = element->getHiY();
    }      

  if(element->is1D())
    {
      m_loY -= fabs(m_loY)*0.2;
      m_hiY += fabs(m_hiY)*0.2;
      m_posLoY *= 0.8;
      m_posHiY *= 1.2;
    }

  const unsigned colors[12] = {1,2,3,4,6,7,1,2,3,4,6,7};
  const unsigned markerStyles[12] = {20,21,22,23,24,25,20,21,22,23,24,25};
  const unsigned lineStyles[12] = {1,1,1,1,1,1,2,2,2,2,2,2};

  unsigned ielement = m_elements.size();

  VSRCanvasElement::Options options = element->getOptions();

  options.lineWidth = m_options.lineWidth;
  options.lineColor = colors[ielement];

  if(m_options.lineStyle.size() > ielement)
    options.lineStyle = m_options.lineStyle[ielement];
  else
    options.lineStyle = lineStyles[ielement];

  options.markerColor = colors[ielement];

  if(m_options.markerStyle.size() > ielement)
    options.markerStyle = m_options.markerStyle[ielement];
  else
    options.markerStyle = markerStyles[ielement];
  options.showErrors = m_options.showErrors;
  options.drawOptions = m_options.drawOptions;
  options.colorPalette = m_options.colorPalette;
  options.showChannelLabels = m_options.showChannelLabels;

  VSRCanvasElement* e = element->clone();
  e->setOptions(options);
  m_elements.push_back(e);
  m_labels.push_back(label);
}

VSRCanvasElement* VSRCanvasPad::getElement(unsigned iel)
{
  vsassert(iel < m_elements.size());
  return m_elements[iel];
}

void VSRCanvasPad::configure(VSOptions& options)
{
  options.findBoolValue("show_channel_labels", 
			s_default_options.showChannelLabels, false,
			"Enable channel labels for all camera plots.");

  options.findWithValue("line_width", s_default_options.lineWidth,
			"Line width.");

  options.findWithValue("marker_style", s_default_options.markerStyle,
			"Marker style.");

  options.findWithValue("line_style", s_default_options.lineStyle,
			"Line style.");

  options.findBoolValue("logx", s_default_options.logX, true,
                        "Plot x axis with log scale.");

  options.findBoolValue("logy", s_default_options.logY, true,
                        "Plot y axis with log scale.");

  options.findBoolValue("logz", s_default_options.logZ, true,
                        "Plot z axis with log scale.");

  options.findWithValue("draw_options", s_default_options.drawOptions,
			"Options for drawing.");

  options.findBoolValue("errors", s_default_options.showErrors, true,
			"Show error bars.");

  options.findWithValue("palette", s_default_options.colorPalette, 
			"Choose color palette to use for plotting 2D "
			"when plotting 2D distributions.");
}

// ============================================================================
// VSRCanvas
// ============================================================================

VSRCanvas::Options VSRCanvas::s_default_options;

VSRCanvas::Options::Options():
  canvas_size(800,600)
{

}

VSRCanvas::VSRCanvas(const std::string& title, const Options& opt):
  m_npixelsX(opt.canvas_size.first), 
  m_npixelsY(opt.canvas_size.second), 
  m_title(title), m_pads(),
  m_canvas()
{

}

VSRCanvas::~VSRCanvas() 
{
  const unsigned npads = m_pads.size();
  for(unsigned ipad = 0; ipad < npads; ipad++)
    delete m_pads[ipad];

  delete m_canvas;
}

void VSRCanvas::add(unsigned ipad, VSRCanvasElement* element,
		    const std::string& label)
{
  if(m_pads.find(ipad) == m_pads.end())
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->add(element,label);
}

void VSRCanvas::setTitle(unsigned ipad, const std::string& title)
{
  vsassert(m_pads.find(ipad) != m_pads.end());
  m_pads[ipad]->setTitle(title);
}

void VSRCanvas::setLimitsX(unsigned ipad, double lo, double hi)
{
  if(m_pads.find(ipad) == m_pads.end()) 
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->setLoX(lo);
  m_pads[ipad]->setHiX(hi);
}

void VSRCanvas::setLimitsY(unsigned ipad, double lo, double hi)
{
  if(m_pads.find(ipad) == m_pads.end()) 
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->setLoY(lo);
  m_pads[ipad]->setHiY(hi);
}

void VSRCanvas::setLimitsZ(unsigned ipad, double lo, double hi)
{
  vsassert(m_pads.find(ipad) != m_pads.end());
  m_pads[ipad]->setLoZ(lo);
  m_pads[ipad]->setHiZ(hi);
}

void VSRCanvas::setBinLabels(unsigned ipad,
			     const std::vector< std::string >& labels)
{
  vsassert(m_pads.find(ipad) != m_pads.end());
  m_pads[ipad]->setBinLabels(labels);
}

void VSRCanvas::setLogX(unsigned ipad, bool logX)
{
  if(m_pads.find(ipad) == m_pads.end()) 
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->setLogX(logX);
}

void VSRCanvas::setLinX(unsigned ipad, bool linX)
{
  if(m_pads.find(ipad) == m_pads.end()) 
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->setLinX(linX);
}

void VSRCanvas::setLogY(unsigned ipad, bool logY)
{
  if(m_pads.find(ipad) == m_pads.end()) 
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->setLogY(logY);
}

void VSRCanvas::setLinY(unsigned ipad, bool linY)
{
  if(m_pads.find(ipad) == m_pads.end()) 
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->setLinY(linY);
}

void VSRCanvas::setGrid(unsigned ipad, bool grid)
{
  if(m_pads.find(ipad) == m_pads.end()) 
    m_pads[ipad] = new VSRCanvasPad;

  m_pads[ipad]->setGrid(grid);
}

void VSRCanvas::draw()
{
  m_canvas = new TCanvas(m_title.c_str(),m_title.c_str(),
			 m_npixelsX,m_npixelsY);
      
  if(m_pads.size() == 2)
    m_canvas->Divide(2,1);
  else if(m_pads.size() > 2 && m_pads.size() <= 4)
    m_canvas->Divide(2,2);
  else if(m_pads.size() > 4 && m_pads.size() <= 9)
    m_canvas->Divide(3,3);
  else if(m_pads.size() > 9 && m_pads.size() <= 16)
    m_canvas->Divide(4,4);
  else if(m_pads.size() > 16 && m_pads.size() <= 25)
    m_canvas->Divide(5,5);
  else if(m_pads.size() > 25 && m_pads.size() <= 36)
    m_canvas->Divide(6,6);
  else if(m_pads.size() > 36 && m_pads.size() <= 64)
    m_canvas->Divide(8,8);     
//     } 
//   else
//     m_canvas->Divide(npadx,npady);
  
  for(std::map< unsigned, VSRCanvasPad* >::iterator itr = m_pads.begin();
      itr != m_pads.end(); ++itr) 
    {
      VSRCanvasPad& pad = *itr->second;

      double lox = pad.getLoX();
      double hix = pad.getHiX();
      double loy = pad.getLoY();
      double hiy = pad.getHiY();

      std::cout << loy << " " << hiy << std::endl;

      if(pad.getLogY() && !pad.is2D())
	{
	  loy = pad.getPosLoY();
	  hiy = pad.getPosHiY();
	}
      else if(!pad.is2D())
	{
	  loy = pad.getLoY();
	  hiy = pad.getHiY();
	}

      if(loy == hiy) 
	{
	  hiy += 1.0;
	  loy -= 1.0;
	} else if (lox == hix) continue;
      
      TVirtualPad* vpad = m_canvas->cd(itr->first+1);

      gStyle->SetOptStat(0);

      std::vector< std::string > binLabels = pad.getBinLabels();

      if(binLabels.size())
	{
	  unsigned nlabel = binLabels.size();

	  TH1F *h = new TH1F("","",nlabel,lox,hix);

	  for(unsigned ilabel = 0; ilabel < nlabel; ilabel++)
	    h->GetXaxis()->SetBinLabel(ilabel+1,binLabels[ilabel].c_str());

	  h->Draw();

	  h->SetMinimum(loy);
	  h->SetMaximum(hiy);

	  if(!pad.getTitle().empty())
	    h->SetTitle(pad.getTitle().c_str());  
	}
      else if(pad.is1D() || pad.is2D())
	{
	  TH1F* h = vpad->DrawFrame(lox,loy,hix,hiy);       
	  if(!pad.getTitle().empty())
	    h->SetTitle(pad.getTitle().c_str());  
	}
      else
	{
	  TPaveText* title = new TPaveText(0.05,0.9,0.35,0.95);

	  title->AddText(pad.getTitle().c_str());
	  title->SetBorderSize(1);
	  title->Draw();
	}

      pad.draw();
    }

  m_canvas->Draw();
  //    gPad->SetFixedAspectRatio();
}

VSRCanvasElement* VSRCanvas::getElement(unsigned ipad, unsigned iel) 
{
  vsassert(m_pads.find(ipad) != m_pads.end());
  return m_pads[ipad]->getElement(iel);
}

void VSRCanvas::configure(VSOptions& options)
{
  VSRCanvasPad::configure(options);

  options.findWithValue("canvas_size", s_default_options.canvas_size, 
			"");
}

