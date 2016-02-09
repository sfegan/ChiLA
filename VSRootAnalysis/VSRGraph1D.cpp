#include <iostream>
#include <algorithm>
#include <sstream>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <WhippleCams.h>
#include <VSSimpleHist.hpp>

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

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRGraph1D.hpp"

using namespace VERITAS;

// ============================================================================
// Constructors
// ============================================================================

VSRGraph1D::VSRGraph1D():
  VSRCanvasElement(),
  m_loX(), m_hiX(), m_loY(), m_hiY(), 
  m_posLoX(), m_posHiX(), m_posLoY(), m_posHiY(), 
  m_nPoints()
{
  m_graph.reset(new TGraphErrors);
}

VSRGraph1D::~VSRGraph1D() {

}

void VSRGraph1D::draw()
{
  if(!m_nPoints) return;

  m_graph->SetMarkerStyle(m_options.markerStyle);
  m_graph->SetMarkerColor(m_options.markerColor);
  m_graph->SetLineStyle(m_options.lineStyle);
  m_graph->SetLineColor(m_options.lineColor);
  m_graph->SetLineWidth(m_options.lineWidth);

  std::string drawOptions = "P";

  if(!m_options.drawOptions.empty()) drawOptions = m_options.drawOptions;

  if(!m_options.showErrors) drawOptions += "X";

  m_graph->Draw(drawOptions.c_str());
}

VSRGraph1D* VSRGraph1D::clone() const
{
  VSRGraph1D* hist = new VSRGraph1D(*this);
  return hist;
}

VSRGraph1D& VSRGraph1D::operator*= (double x)
{
  const unsigned np = getN();
  for(unsigned ip = 0; ip < np; ip++)
    {
      m_graph->GetY()[ip] *= x;
      m_graph->GetEY()[ip] *= x;
    }

  m_hiY *= x;
  m_loY *= x;
  m_posHiY *= x;
  m_posLoY *= x;

  return *this;
}

double VSRGraph1D::getX(unsigned ipoint) const
{
  return m_graph->GetX()[ipoint];
}

double VSRGraph1D::getErrX(unsigned ipoint) const
{
  return m_graph->GetErrorX(ipoint);
}

double VSRGraph1D::getY(unsigned ipoint) const
{
  return m_graph->GetY()[ipoint];
}

double VSRGraph1D::getErrY(unsigned ipoint) const
{
  return m_graph->GetErrorY(ipoint);
}

void VSRGraph1D::set(double x, double y) 
{
  if(!m_nPoints || x > m_hiX) m_hiX = x;
  if(!m_nPoints || x < m_loX) m_loX = x;
  if(!m_nPoints || y > m_hiY) m_hiY = y;
  if(!m_nPoints || y < m_loY) m_loY = y;

  if(!m_posHiX || (x && x > m_posHiX)) m_posHiX = x;
  if(!m_posLoX || (x && x < m_posLoX)) m_posLoX = x;
  if(!m_posHiY || (y && y > m_posHiY)) m_posHiY = y;
  if(!m_posLoY || (y && y < m_posLoY)) m_posLoY = y;

  m_graph->SetPoint(m_nPoints,x,y);
  m_graph->SetPointError(m_nPoints,0,0);

  m_nPoints++;
}

void VSRGraph1D::set(double x, double y, double yerr) 
{
  if(!m_nPoints || x > m_hiX) m_hiX = x;
  if(!m_nPoints || x < m_loX) m_loX = x;
  if(!m_nPoints || y + yerr > m_hiY) m_hiY = y+yerr;
  if(!m_nPoints || y - yerr < m_loY) m_loY = y-yerr;

  if(!m_posHiX || (x && x > m_posHiX)) m_posHiX = x;
  if(!m_posLoX || (x && x < m_posLoX)) m_posLoX = x;
  if(!m_posHiY || (y && y > m_posHiY)) m_posHiY = y;
  if(!m_posLoY || (y && y < m_posLoY)) m_posLoY = y;

  m_graph->SetPoint(m_nPoints,x,y);
  m_graph->SetPointError(m_nPoints,0,yerr);

  m_nPoints++;
}

void VSRGraph1D::print()
{
  for(int n = 0; n < m_graph->GetN(); n++)
    {
      std::cout << m_graph->GetX()[n] << " "
		<< m_graph->GetY()[n] << std::endl;
    }
}
void VSRGraph1D::dump()
{
  print();
}


void VSRGraph1D::clear()
{
  m_nPoints = 0;
  m_loX = 0;
  m_hiX = 0;
  m_loY = 0;
  m_hiY = 0;
  m_graph->Clear();
}

// ============================================================================
// Copy Constructor
// ============================================================================
VSRGraph1D::VSRGraph1D(const VSRGraph1D& graph) 
{
  *this = graph;
  VSRCanvasElement::operator= (graph);
}

// ============================================================================
// Assignment Operator
// ============================================================================
VSRGraph1D& VSRGraph1D::operator=(const VSRGraph1D &graph) 
{
  m_graph.reset( (TGraphErrors*)graph.m_graph->Clone() );
  m_loX = graph.m_loX;
  m_hiX = graph.m_hiX;
  m_loY = graph.m_loY;
  m_hiY = graph.m_hiY;
  m_nPoints = graph.m_nPoints;
  return *this;
}

VSRGraph1D& VSRGraph1D::operator/=(const VSRGraph1D & graph)
{
  std::vector< double > x;
  std::vector< double > y;
  std::vector< double > yerr;

  const unsigned np = getN();
  for(unsigned ip = 0; ip < np; ip++)
    {
      if(ip >= graph.getN()) break;

      double x1 = getX(ip);
      double x2 = graph.getX(ip);

      //if(fabs(x1-x2) > DBL_EPSILON) continue;

      double y1 = getY(ip);
      double y1_err = getErrY(ip);
      double y2 = graph.getY(ip);
      double y2_err = graph.getErrY(ip);

      if(y2 == 0) 
	{
	  x.push_back(x1);
	  y.push_back(1);
	  yerr.push_back(0);
	}
      else
	{
	  x.push_back(x1);
	  y.push_back(y1/y2);
	  yerr.push_back(y1/y2*sqrt(std::pow(y1_err/y1,2) + 
				    std::pow(y2_err/y2,2)));
	}
    }

  clear();

  for(unsigned ip = 0; ip < np; ip++)
    set(x[ip],y[ip],yerr[ip]);

  return *this;
}
