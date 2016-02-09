//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRHISTOGRAM2D_HPP
#define VSRHISTOGRAM2D_HPP

#include <vector>
#include <memory>

#include <VSSimple2DHist.hpp>
#include <VSNSpace.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TPad.h>
#include <TH2F.h>
#include <TMatrixT.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram.hpp"
#include "VSRHistogram1D.hpp"

class VSRHistogram2D : public VSRHistogram 
{
public:
  friend class VSRH5Loader;

  VSRHistogram2D();
  VSRHistogram2D(unsigned nbinsx,double lox,double hix,
		 unsigned nbinsy,double loy,double hiy,
		 const std::string& title = "");
  ~VSRHistogram2D();

  // --------------------------------------------------------------------------
  // Virtual Methods
  // --------------------------------------------------------------------------
  virtual void draw();
  virtual VSRHistogram2D* clone() const;

  virtual void addLegendEntry(TLegend *leg, const std::string& label) const
  {
    leg->AddEntry(m_hist.get(),label.c_str(),m_options.legendStyle.c_str());
  }
  
  virtual double getLoX() const { return m_loX; }
  virtual double getHiX() const { return m_hiX; }
  virtual double getLoY() const { return m_loY; }
  virtual double getHiY() const { return m_hiY; }
  virtual double getPosLoX() const { return m_loX; }
  virtual double getPosHiX() const { return m_hiX; }
  virtual double getPosLoY() const { return m_loY; }
  virtual double getPosHiY() const { return m_hiY; }

  virtual bool is1D() const { return false; }
  virtual bool is2D() const { return true; }

  void fill(double x, double y, double w = 1);

  std::string getTitle() { return m_title; }

  TH2F* getHistogram();

  unsigned getIndexX(double x) const;
  unsigned getIndexY(double y) const;

  double getValueX(unsigned ix) const
  { 
    return m_loX + (ix+0.5)*m_binWidthX;
  }

  double getValueY(unsigned iy) const
  { 
    return m_loY + (iy+0.5)*m_binWidthY;
  }

  unsigned getNBins() const { return m_nBinsX*m_nBinsY; }
  unsigned getNBinsX() const { return m_nBinsX; }
  unsigned getNBinsY() const { return m_nBinsY; }

  double get(unsigned ix, unsigned iy) const 
  { return m_hist->GetBinContent(ix+1,iy+1); }

  double getError(unsigned ix, unsigned iy) const 
  { return m_hist->GetBinError(ix+1,iy+1); }

  void setBinValue(unsigned ibinx, unsigned ibiny, double z);
  void setBinError(unsigned ibinx, unsigned ibiny, double err);

  void fit(const std::string& fname);

  void normalize() { }
  void normalizeMax();

  VSRHistogram1D* projectX(unsigned lobin = 0, unsigned hibin = 0);
  VSRHistogram1D* projectY(unsigned lobin = 0, unsigned hibin = 0);

  void smooth(unsigned n);

  void setRange(double xlo, double xhi, double ylo, double yhi)
  {
    m_hist->GetXaxis()->SetRangeUser(xlo,xhi);
    m_hist->GetYaxis()->SetRangeUser(ylo,yhi);
  }

  void multiply(double x);

  // --------------------------------------------------------------------------
  // Factory methods
  //---------------------------------------------------------------------------
  static VSRHistogram2D* create(TMatrixT<double>& m, 
				double xmin = 0, double xmax = 0, 
				double ymin = 0, double ymax = 0);
  static VSRHistogram2D* create(VERITAS::VSNSpace& nspace);

  template< typename T1, typename T2 >
  static VSRHistogram2D* create(VERITAS::VSSimple2DHist<T1,T1,T2>* hist);

  // --------------------------------------------------------------------------
  // Copy Constructor
  //---------------------------------------------------------------------------
  VSRHistogram2D(const VSRHistogram2D& hist);

  // --------------------------------------------------------------------------
  // Assignment Operator
  //---------------------------------------------------------------------------
  VSRHistogram2D& operator=(const VSRHistogram2D & hist);

  // --------------------------------------------------------------------------
  // Arithmetic Operators
  //---------------------------------------------------------------------------
  VSRHistogram2D operator-(const VSRHistogram2D & hist) const;
  VSRHistogram2D& operator-=(const VSRHistogram2D & hist);

  VSRHistogram2D operator/(const VSRHistogram2D & hist) const;
  VSRHistogram2D& operator/=(const VSRHistogram2D & hist);

private:

  double sum(unsigned ix, unsigned iy, unsigned n);

  std::string m_title;
  unsigned m_nBinsX;
  double m_loX;
  double m_hiX;
  double m_binWidthX;
  unsigned m_nBinsY;
  double m_loY;
  double m_hiY;
  double m_binWidthY;

  std::auto_ptr<TH2F> m_hist;
};

template< typename T1, typename T2 >
VSRHistogram2D* VSRHistogram2D::create(VERITAS::VSSimple2DHist<T1,T1,T2>* hist)
{
  VSRHistogram2D* output_hist = new VSRHistogram2D;

  output_hist->m_nBinsX = hist->nXBins();
  output_hist->m_loX = hist->xLoLimit();
  output_hist->m_hiX = hist->xHiLimit();
  output_hist->m_nBinsY = hist->nYBins();
  output_hist->m_loY = hist->yLoLimit();
  output_hist->m_hiY = hist->yHiLimit();

  output_hist->m_hist.reset(new TH2F("","",
				     output_hist->m_nBinsX,
				     output_hist->m_loX,
				     output_hist->m_hiX,
				     output_hist->m_nBinsY,
				     output_hist->m_loY,
				     output_hist->m_hiY));
  output_hist->m_hist->SetDirectory(0);

  for(typename VERITAS::VSSimple2DHist<T1,T1,T2>::iterator itr = 
	hist->begin(); itr != hist->end(); ++itr)
    {
      double count = 0;
	
      if(std::isfinite(itr->count()))
	count = (double)itr->count();
	
      unsigned xbin = hist->indexToXBin(itr->bin());
      unsigned ybin = hist->indexToYBin(itr->bin());
      output_hist->m_hist->SetBinContent(xbin+1,ybin+1,count);
    }
    
  return output_hist;
}


#endif // VSRHISTOGRAM2D_HPP
