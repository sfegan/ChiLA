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
#include <TColor.h>
#include <TROOT.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram2D.hpp"
#include "VSRPalette.hpp"

using namespace VERITAS;

// ============================================================================
// Factory Methods
// ============================================================================

VSRHistogram2D* VSRHistogram2D::create(VERITAS::VSNSpace& nspace)
{
  vsassert(nspace.space().ndim >= 2);

  VERITAS::VSNSpace nspace_proj = nspace;
  std::set< unsigned > dims;

  dims.insert(0);
  dims.insert(1);
    
  nspace_proj.project(dims);
    
  VSRHistogram2D* hist = 
    new VSRHistogram2D(nspace_proj.space().axes[0].nbin,
		       nspace_proj.space().axes[0].lo_bound,
		       nspace_proj.space().axes[0].hi_bound,
		       nspace_proj.space().axes[1].nbin,
		       nspace_proj.space().axes[1].lo_bound,
		       nspace_proj.space().axes[1].hi_bound);
    
  VERITAS::VSNSpace::Cell c(2);
    
  for(c.i[0]=0; c.i[0]<nspace.space().axes[0].nbin; c.i[0]++)
    {
      for(c.i[1]=0; c.i[1]<nspace.space().axes[1].nbin; c.i[1]++)
	{	      
	  VERITAS::VSNSpace::Weight w = -1; // initialize any value
	  nspace_proj.getWeight(c,w);
	    
	  if(std::isfinite(w))
	    hist->setBinValue(c.i[0],c.i[1],w);	
	}	
    }
    
  return hist;
}

VSRHistogram2D* VSRHistogram2D::create(TMatrixT<double>& m, 
				       double xmin, double xmax, 
				       double ymin, double ymax)
{
  const unsigned nrow = m.GetNrows();
  const unsigned ncol = m.GetNrows();
  
  VSRHistogram2D* hist;

  if(xmin != xmax)
    hist = new VSRHistogram2D(nrow, xmin, xmax, ncol, ymin, ymax);
  else
    hist = new VSRHistogram2D(nrow, 0, (double)nrow, ncol, 0, (double)ncol);

  for(unsigned irow = 0; irow < nrow; irow++)
    {
      for(unsigned icol = 0; icol < ncol; icol++)
	{
	  hist->setBinValue(irow,icol,m[irow][icol]);
	}
    }

  return hist;
}

// ============================================================================
// Constructors
// ============================================================================

VSRHistogram2D::VSRHistogram2D():
  m_title(), m_nBinsX(), m_loX(),m_hiX(),
  m_nBinsY(),m_loY(),m_hiY()
{
  m_hist.reset(new TH2F("","",m_nBinsX,m_loX,m_hiX,m_nBinsY,m_loY,m_hiY));
  m_hist->SetDirectory(0);
  m_hist->SetContour(50);
}

VSRHistogram2D::VSRHistogram2D(unsigned nBinsX, double loX, double hiX,
			       unsigned nBinsY, double loY, double hiY,
			       const std::string& title):
  m_title(title),
  m_nBinsX(nBinsX),m_loX(loX),m_hiX(hiX),
  m_nBinsY(nBinsY),m_loY(loY),m_hiY(hiY)
{
  m_binWidthX = (m_hiX-m_loX)/(double)m_nBinsX;
  m_binWidthY = (m_hiY-m_loY)/(double)m_nBinsY;

  m_hist.reset(new TH2F(m_title.c_str(),m_title.c_str(),
			m_nBinsX,m_loX,m_hiX,m_nBinsY,m_loY,m_hiY));
  m_hist->SetDirectory(0);
  m_hist->SetContour(50);
}

VSRHistogram2D::~VSRHistogram2D() 
{

}

void setPalette(const std::vector<int>& colors)
{
  
}

void VSRHistogram2D::draw()
{
  if(m_options.loZ != m_options.hiZ)
    {
      m_hist->SetMinimum(m_options.loZ);
      m_hist->SetMaximum(m_options.hiZ);
    } 

  double min = m_hist->GetMinimum();
  double max = m_hist->GetMaximum();

  VSRPalette::setPalette(m_options.colorPalette,min,max);

  std::string drawOptions = "SAME";

  if(m_options.drawOptions.empty()) drawOptions += "COLZ";
  else drawOptions += m_options.drawOptions;

  m_hist->Draw(drawOptions.c_str());
}

VSRHistogram2D* VSRHistogram2D::clone() const
{
  VSRHistogram2D* hist = new VSRHistogram2D(*this);
  return hist;
}

VSRHistogram2D::VSRHistogram2D(const VSRHistogram2D& hist) 
{
  *this = hist;
}

TH2F* VSRHistogram2D::getHistogram() { return m_hist.get(); }

VSRHistogram2D& VSRHistogram2D::operator=(const VSRHistogram2D & hist) 
{
  m_hist.reset( (TH2F*)hist.m_hist->Clone() );
  m_hist->SetDirectory(0);
  m_title = hist.m_title;
  m_nBinsX = hist.m_nBinsX;
  m_loX = hist.m_loX;
  m_hiX = hist.m_hiX;
  m_binWidthX = hist.m_binWidthX;
  m_nBinsY = hist.m_nBinsY;
  m_loY = hist.m_loY;
  m_hiY = hist.m_hiY;
  m_binWidthY = hist.m_binWidthY;
  return *this;
}

void VSRHistogram2D::fill(double x, double y, double w) 
{
  if(!std::isfinite(x) || !std::isfinite(y)) return;
  else if( x < m_loX || x > m_hiX || y < m_loY || y > m_hiY )
    return;

  m_hist->Fill(x,y,w);
}

void VSRHistogram2D::setBinValue(unsigned ibinx, unsigned ibiny, double z)
{
  m_hist->SetBinContent(ibinx+1,ibiny+1,z);
}

void VSRHistogram2D::setBinError(unsigned ibinx, unsigned ibiny, double err)
{
  m_hist->SetBinError(ibinx+1,ibiny+1,err);
}

unsigned VSRHistogram2D::getIndexX(double x) const
{
  return (unsigned)((x-m_loX)/m_binWidthX);
}

unsigned VSRHistogram2D::getIndexY(double y) const
{
  return (unsigned)((y-m_loY)/m_binWidthY);
}

void VSRHistogram2D::fit(const std::string& fname)
{
  m_hist->Fit(fname.c_str(),"NR","E");
}

void VSRHistogram2D::normalizeMax()
{
  double max = m_hist->GetMaximum();

  for(unsigned ix = 0; ix < m_nBinsX; ix++) 
    for(unsigned iy = 0; iy < m_nBinsY; iy++) 
      setBinValue(ix,iy,get(ix,iy)/max);
}

VSRHistogram1D* VSRHistogram2D::projectX(unsigned lobin, unsigned hibin)
{
  if(lobin == hibin) lobin = 0, hibin = m_nBinsY;

  VSRHistogram1D* h = new VSRHistogram1D(m_nBinsX,m_loX,m_hiX);
  for(unsigned iy = lobin; iy <hibin; iy++)
    for(unsigned ix = 0; ix < m_nBinsX; ix++)
      h->fillBin(ix,get(ix,iy));

  return h;
}

VSRHistogram1D* VSRHistogram2D::projectY(unsigned lobin, unsigned hibin)
{
  if(lobin == hibin) lobin = 0, hibin = m_nBinsX;

  VSRHistogram1D* h = new VSRHistogram1D(m_nBinsY,m_loY,m_hiY);
  for(unsigned ix = lobin; ix <hibin; ix++)
    for(unsigned iy = 0; iy < m_nBinsY; iy++)
      h->fillBin(iy,get(ix,iy));

  return h;
}

void VSRHistogram2D::smooth(unsigned n)
{
  std::vector< std::vector< double > > tmp(m_nBinsX);

  for(unsigned ix = 0; ix < m_nBinsX; ix++)
    {
      for(unsigned iy = 0; iy < m_nBinsY; iy++)
	{
	  tmp[ix].push_back(sum(ix,iy,n));
	}
    }

  for(unsigned ix = 0; ix < m_nBinsX; ix++)
    for(unsigned iy = 0; iy < m_nBinsY; iy++)
      setBinValue(ix,iy,tmp[ix][iy]);
}

double VSRHistogram2D::sum(unsigned _ix, unsigned _iy, unsigned n)
{
  unsigned xlo = std::max((int)(_ix-n),0);
  unsigned xhi = std::min((int)(_ix+n),(int)m_nBinsX);
  unsigned ylo = std::max((int)(_iy-n),0);
  unsigned yhi = std::min((int)(_iy+n),(int)m_nBinsY);

  double sum = 0;

  for(unsigned ix = xlo; ix < xhi; ix++)
    for(unsigned iy = ylo; iy < yhi; iy++)
      sum+=get(ix,iy);

  return sum;
}

void VSRHistogram2D::multiply(double x)
{
  for(unsigned ix = 0; ix < m_nBinsX; ix++)
    for(unsigned iy = 0; iy < m_nBinsY; iy++)
      {
	setBinValue(ix,iy,x*get(ix,iy));
	setBinError(ix,iy,x*getError(ix,iy));
      }
}

VSRHistogram2D VSRHistogram2D::operator-(const VSRHistogram2D & hist) const 
{
  return VSRHistogram2D(*this) -= hist;
}

VSRHistogram2D& VSRHistogram2D::operator-=(const VSRHistogram2D & hist)
{
  vsassert(getNBinsX() == hist.getNBinsX() &&
	   getNBinsY() == hist.getNBinsY());

  for(unsigned ix = 0; ix < getNBinsX(); ix++)
    for(unsigned iy = 0; iy < getNBinsY(); iy++)
      {
	double x = get(ix,iy);
	double xerr = getError(ix,iy);
	double y = hist.get(ix,iy);
	double yerr = hist.getError(ix,iy);
	setBinValue(ix,iy,x-y);
	setBinError(ix,iy,xerr+yerr);
      }
  return *this;
}

VSRHistogram2D VSRHistogram2D::operator/(const VSRHistogram2D & hist) const 
{
  return VSRHistogram2D(*this) /= hist;
}

VSRHistogram2D& VSRHistogram2D::operator/=(const VSRHistogram2D & hist)
{
  vsassert(getNBinsX() == hist.getNBinsX() &&
	   getNBinsY() == hist.getNBinsY());

  for(unsigned ix = 0; ix < getNBinsX(); ix++)
    for(unsigned iy = 0; iy < getNBinsY(); iy++)
      {
	if(hist.get(ix,iy) == 0 || get(ix,iy) == 0)
	  {
	    setBinValue(ix,iy,0);
	    setBinError(ix,iy,0);
	  }
	else
	  {
	    double x = get(ix,iy);
	    double xerr = getError(ix,iy);
	    double y = hist.get(ix,iy);
	    double yerr = hist.getError(ix,iy);
	    double err = x/y*sqrt(pow(xerr/x,2) + pow(yerr/y,2));
	    setBinValue(ix,iy,x/y);
	    setBinError(ix,iy,err);
	  }
      }
  return *this;
}
