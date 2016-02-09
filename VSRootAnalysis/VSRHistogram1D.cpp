//-*-mode:c++; mode:font-lock;-*-

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
#include <TVectorT.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram1D.hpp"


using namespace std;
using namespace VERITAS;

// ============================================================================
// Factory Methods
// ============================================================================

VSRHistogram1D* VSRHistogram1D::create(TVectorT<double>& v, 
				       double min, double max)
{
  const unsigned nrow = v.GetNrows();
  
  VSRHistogram1D* hist;
  
  if(min != max)
    hist = new VSRHistogram1D(nrow, min, max);
  else
    hist = new VSRHistogram1D(nrow, 0, (double)nrow);
  
  for(unsigned irow = 0; irow < nrow; irow++)
    hist->set(irow,v[irow]);
  
  return hist;
}

VSRHistogram1D* VSRHistogram1D::create(TVectorT<double>& v, 
				       TVectorT<double>& var, 
				       double min, double max)
{
  const unsigned nrow = v.GetNrows();
        
  VSRHistogram1D* hist;

  if(min != max)
    hist = new VSRHistogram1D(nrow, min, max);
  else
    hist = new VSRHistogram1D(nrow, 0, (double)nrow);

  for(unsigned irow = 0; irow < nrow; irow++)
    {
      hist->set(irow,v[irow]);
      hist->setError(irow,sqrt(var[irow]));
    }

  return hist;
}

VSRHistogram1D* VSRHistogram1D::create(VERITAS::VSNSpace& nspace)
{
  VERITAS::VSNSpace nspace_proj = nspace;
  std::set< unsigned > dims;

  dims.insert(0);
    
  nspace_proj.project(dims);
    
  VSRHistogram1D* hist = 
    new VSRHistogram1D(nspace_proj.space().axes[0].nbin,
		       nspace_proj.space().axes[0].lo_bound,
		       nspace_proj.space().axes[0].hi_bound);
    
  VERITAS::VSNSpace::Cell c(1);
    
  for(c.i[0]=0; c.i[0]<nspace.space().axes[0].nbin; c.i[0]++)
    {
      VERITAS::VSNSpace::Weight w = -1; // initialize any value
      nspace_proj.getWeight(c,w);
      
      if(std::isfinite(w))
	hist->set(c.i[0],w);	
    }
    
  return hist;
}

// ============================================================================
// Constructors
// ============================================================================

VSRHistogram1D::VSRHistogram1D():
  m_title(), m_nBins(), m_loX(), m_hiX(), m_binWidth()
{
  m_hist.reset(new TH1F("","",100,0.0,1.0));
  m_hist->SetDirectory(0);
}

VSRHistogram1D::VSRHistogram1D(unsigned nbins, double lo, double hi,
			       const std::string& title):
  m_title(title), m_nBins(nbins), m_loX(lo), m_hiX(hi), m_binWidth()
{
  m_binWidth = (m_hiX-m_loX)/(double)m_nBins;
  m_hist.reset(new TH1F(title.c_str(),title.c_str(),nbins,lo,hi));
  m_hist->SetDirectory(0);
}

VSRHistogram1D::~VSRHistogram1D() 
{

}

void VSRHistogram1D::dump()
{
  for(unsigned ibin = 0; ibin < getNBins(); ibin++)
    {
      std::cout << std::setw(20) << getX(ibin) 
		<< std::setw(20) << get(ibin) 
		<< std::setw(20) << getError(ibin) 
		<< std::endl;
    }
}

void VSRHistogram1D::draw()
{
  m_hist->SetMarkerStyle(m_options.markerStyle);
  m_hist->SetMarkerColor(m_options.markerColor);
  m_hist->SetLineStyle(m_options.lineStyle);
  m_hist->SetLineColor(m_options.lineColor);
  m_hist->SetLineWidth(m_options.lineWidth);
  m_hist->SetFillColor(m_options.fillColor);

  std::string drawOptions = "SAME";  // Forces not to redraw axes!

  if(m_options.showErrors) drawOptions += "E";
  else drawOptions += "HIST";

  drawOptions += m_options.drawOptions;
  m_hist->Draw(drawOptions.c_str());
}

VSRHistogram1D* VSRHistogram1D::clone() const
{
  VSRHistogram1D* hist = new VSRHistogram1D(*this);
  return hist;
}

double VSRHistogram1D::getIntegral()
{

  double bin_size = (m_hiX-m_loX)/(double)m_nBins;

  double sum = 0;
  for(unsigned ibin = 0; ibin < getNBins(); ibin++)
    sum += this->get(ibin)*bin_size;

  return sum;
}

TH1F* VSRHistogram1D::getHistogram() { return m_hist.get(); }

void VSRHistogram1D::initialize(unsigned n, double lo, double hi) 
{
  m_hist.reset(new TH1F(m_title.c_str(),
			m_title.c_str(),
			n,lo,hi));
  m_hist->SetDirectory(0);
}

VSRHistogram1D::VSRHistogram1D(const VSRHistogram1D& hist) 
{
  *this = hist;
}

VSRHistogram1D& VSRHistogram1D::operator=(const VSRHistogram1D & hist) 
{
  m_hist.reset( (TH1F*)hist.m_hist->Clone() );
  m_hist->SetDirectory(0);
  m_title = hist.m_title;
  m_nBins = hist.m_nBins;
  m_loX = hist.m_loX;
  m_hiX = hist.m_hiX;
  return *this;
}

void VSRHistogram1D::subtract(VSRHistogram1D* hist, double weight) 
{
  m_hist->Sumw2();
  hist->getHistogram()->Sumw2();

  m_hist->Add(hist->getHistogram(),-weight);
}


VSRHistogram1D VSRHistogram1D::operator-(const VSRHistogram1D & hist) const 
{
  return VSRHistogram1D(*this) -= hist;
}

VSRHistogram1D& VSRHistogram1D::operator-=(const VSRHistogram1D & hist)
{
  m_hist->Add(hist.m_hist.get(),-1.0);
  return *this;
}

VSRHistogram1D VSRHistogram1D::operator/(const VSRHistogram1D & hist) const 
{
  return VSRHistogram1D(*this) /= hist;
}

VSRHistogram1D& VSRHistogram1D::operator/=(const VSRHistogram1D & hist)
{
  for(unsigned ibin = 0; ibin < getNBins(); ibin++)
    {
      if(hist.get(ibin) == 0 || get(ibin) == 0)
	{
	  set(ibin,0);
	  setError(ibin,0);
	}
      else
	{
	  double x = get(ibin);
	  double xerr = getError(ibin);
	  double y = hist.get(ibin);
	  double yerr = hist.getError(ibin);

	  double err = x/y*sqrt(pow(xerr/x,2) + pow(yerr/y,2));
	  set(ibin,x/y);
	  setError(ibin,err);
	}
    }

  return *this;
}

void VSRHistogram1D::fill(double x, double w) 
{
  m_hist->Fill(x,w);
}

void VSRHistogram1D::fillBin(unsigned ibin, double w) 
{
  set(ibin,w+get(ibin));
}

void VSRHistogram1D::rebin(unsigned ngroup)
{
  m_hist->Rebin(ngroup);
}

void VSRHistogram1D::set(unsigned ibin, double x) 
{
  m_hist->SetBinContent(ibin+1,x);
}

void VSRHistogram1D::setError(unsigned ibin, double err) 
{
  m_hist->SetBinError(ibin+1,err);
}

double VSRHistogram1D::getX(unsigned ibin) const
{
  return m_hist->GetBinCenter(ibin+1);
}

double VSRHistogram1D::get(unsigned ibin) const
{
  return m_hist->GetBinContent(ibin+1);
}

double VSRHistogram1D::getError(unsigned ibin) const
{
  return m_hist->GetBinError(ibin+1);
}

double VSRHistogram1D::getSum() 
{
  double sum = 0;

  for(int i = 1; i <= m_hist->GetNbinsX(); i++) 
    sum += m_hist->GetBinContent(i);

  return sum;
}

double VSRHistogram1D::getRMS() 
{
  return m_hist->GetRMS();
}

void VSRHistogram1D::fit(string fname) 
{

  m_hist->Fit(fname.c_str(),"N","E");

}

void VSRHistogram1D::fit(string fname, double lo, double hi) 
{

  m_hist->Fit(fname.c_str(),"N","E",lo,hi);

}

void VSRHistogram1D::normalize() 
{
  double sum = getSum();

  if(sum == 0)
    return;

  for(int i = 1; i <= m_hist->GetNbinsX(); i++) 
    {
      double content = m_hist->GetBinContent(i);;
      double err = m_hist->GetBinError(i);
      m_hist->SetBinContent(i,content/sum);
      m_hist->SetBinError(i,err/sum);      
    }
}

void VSRHistogram1D::normalizeMax() 
{
  double max = m_hist->GetMaximum();

  for(int i = 1; i <= m_hist->GetNbinsX(); i++) 
    {
      set(i-1,get(i-1)/max);
      setError(i-1,getError(i-1)/max);
    }
}

void VSRHistogram1D::transform(const std::string& transform) 
{
  for(int i = m_hist->GetNbinsX(); i >= 1; i--) 
    {
      double content = m_hist->GetBinContent(i);;
      double err = m_hist->GetBinError(i);
      
      if(transform == "sq")
	{
	  content = std::pow(content,2);
	  err = 2*err;
	}
      
      m_hist->SetBinContent(i,content);
      m_hist->SetBinError(i,err);      
    }

}

void VSRHistogram1D::cumulative(const std::string& lr) 
{

  double sum = 0;
  double sum_var = 0;

  if(lr == "r")
    {
      for(int i = m_hist->GetNbinsX(); i >= 1; i--) 
	{
	  double content = m_hist->GetBinContent(i);;
	  double err = m_hist->GetBinError(i);
	  
	  sum += content;
	  sum_var += err*err;
	  
	  m_hist->SetBinContent(i,sum);
	  m_hist->SetBinError(i,sqrt(sum_var));      
	}
    }
  else
    {
      for(int i = 1; i <= m_hist->GetNbinsX(); i++) 
	{
	  double content = m_hist->GetBinContent(i);;
	  double err = m_hist->GetBinError(i);
	  
	  sum += content;
	  sum_var += err*err;
	  
	  m_hist->SetBinContent(i,sum);
	  m_hist->SetBinError(i,sqrt(sum_var));      
	}
    }
}

void VSRHistogram1D::multiply(double x)
{
  for(int i = 1; i <= m_hist->GetNbinsX(); i++) 
    {
      double content = m_hist->GetBinContent(i);;
      double err = m_hist->GetBinError(i);
      m_hist->SetBinContent(i,content*x);
      m_hist->SetBinError(i,err*x);      
    }
}

void VSRHistogram1D::print()
{
  for(int ibin = 1; ibin <= m_hist->GetNbinsX(); ibin++) 
    {
      double content = m_hist->GetBinContent(ibin);;
      double err = m_hist->GetBinError(ibin);

      double x = m_hist->GetBinCenter(ibin);

      std::cout << x << " " << content << " " << err << std::endl;
    }
}
