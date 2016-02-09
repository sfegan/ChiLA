//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRHISTOGRAM1D_HPP
#define VSRHISTOGRAM1D_HPP

#include <vector>
#include <memory>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSSimpleHist.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSNSpace.hpp>
#include <VSSimpleStat.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram.hpp"

class VSRHistogram1D : public VSRHistogram 
{
public:
  friend class VSRH5Loader;

  VSRHistogram1D();
  VSRHistogram1D(unsigned nbins, double lo, double hi,
		 const std::string& title = "");

  ~VSRHistogram1D();

  // Virtual Methods ----------------------------------------------------------
  virtual void dump();
  virtual void draw();
  virtual VSRHistogram1D* clone() const;

  virtual void addLegendEntry(TLegend *leg, const std::string& label) const
  {
    leg->AddEntry(m_hist.get(),label.c_str(),
		  m_options.legendStyle.c_str());
  }
  
  virtual double getLoX() const { return m_loX; }
  virtual double getHiX() const { return m_hiX; }
  virtual double getLoY() const { return m_hist->GetMinimum(); }
  virtual double getHiY() const { return m_hist->GetMaximum(); }
  virtual double getPosLoX() const { return m_loX; }
  virtual double getPosHiX() const { return m_hiX; }
  virtual double getPosLoY() const 
  { 
    double min = 0;

    for(unsigned ibin = 0; ibin < getNBins(); ibin++)
      {
	if(get(ibin) && !min) min = get(ibin);
	else if(get(ibin) && min && get(ibin) < min) min = get(ibin);
      }

    return min;
  }
  virtual double getPosHiY() const { return m_hist->GetMaximum(); }

  virtual bool is1D() const { return true; }
  virtual bool is2D() const { return false; }

  void fill(double x, double w = 1);
  void fillBin(unsigned ibin, double w = 1);

  void rebin(unsigned ngroup);

  void initialize(unsigned n, double lo, double hi);

  void set(unsigned ibin, double x);
  void setError(unsigned ibin, double err);

  std::string getTitle() const { return m_title; }

  void subtract(VSRHistogram1D *hist, double weight);

  void normalize();
  void normalizeMax();
  void transform(const std::string& transform) ;
  void cumulative(const std::string& lr);

  void fit(std::string fname);
  void fit(std::string fname, double lo, double hi);

  double getX(unsigned ibin) const;
  double get(unsigned ibin) const;
  double getError(unsigned ibin) const;

  double getSum();
  double getRMS();

  void multiply(double x);

  unsigned getNBins() const { return m_hist->GetNbinsX(); }

  double getIntegral();

  TH1F* getHistogram();

  void print();

  // --------------------------------------------------------------------------
  // Factory methods
  //---------------------------------------------------------------------------
  static VSRHistogram1D* create(TVectorT<double>& v, 
				double min = 0, double max = 0);
  static VSRHistogram1D* create(TVectorT<double>& v, 
				TVectorT<double>& var, 
				double min = 0, double max = 0);

  static VSRHistogram1D* create(VERITAS::VSNSpace& nspace);

  template< typename T1, typename T2, typename T3 >
  static VSRHistogram1D* create(VERITAS::VSLimitedHist<T1,T2,T3>* hist);

  template< typename T1, typename T2 >
  static VSRHistogram1D* create(VERITAS::VSLimitedErrorsHist<T1,T1,T2>* hist);

  // --------------------------------------------------------------------------
  // Copy Constructor
  //---------------------------------------------------------------------------
  VSRHistogram1D(const VSRHistogram1D& hist);

  // --------------------------------------------------------------------------
  // Assignment Operator
  //---------------------------------------------------------------------------
  VSRHistogram1D& operator=(const VSRHistogram1D & hist);

  // --------------------------------------------------------------------------
  // Arithmetic Operator
  //---------------------------------------------------------------------------
  VSRHistogram1D operator-(const VSRHistogram1D & hist) const;
  VSRHistogram1D& operator-=(const VSRHistogram1D & hist);

  VSRHistogram1D operator/(const VSRHistogram1D & hist) const;
  VSRHistogram1D& operator/=(const VSRHistogram1D & hist);

private:

  std::string         m_title;
  unsigned            m_nBins;
  double              m_loX;
  double              m_hiX;
  double              m_binWidth;
  VERITAS::VSSimpleStat2<double> m_stat;
  std::auto_ptr<TH1F> m_hist;
};

// (Bin coord, Type of contents, Bin calculator method)
template< typename T1, typename T2, typename T3 >
VSRHistogram1D* VSRHistogram1D::create(VERITAS::VSLimitedHist<T1,T2,T3>* hist)
{
  VSRHistogram1D* output_hist = new VSRHistogram1D;

  output_hist->m_nBins = 
    lround((hist->hiLimit()-hist->loLimit())/hist->binSize());
  output_hist->m_loX    = hist->loLimit();
  output_hist->m_hiX    = hist->hiLimit();

  int offset = hist->bin(hist->loLimit());

  output_hist->m_hist.reset(new TH1F("","",
				     output_hist->m_nBins,
				     output_hist->m_loX,
				     output_hist->m_hiX) );
  output_hist->m_hist->SetDirectory(0);

  for(typename VERITAS::VSLimitedHist<T1,T2,T3>::iterator itr = 
	hist->begin(); itr != hist->end(); ++itr)
    {
      int bin = itr->bin() - hist->firstBin();
      output_hist->set(bin-offset,(double)itr->count());
    }

  return output_hist;
}

template< typename T1, typename T2 >
VSRHistogram1D* 
VSRHistogram1D::create(VERITAS::VSLimitedErrorsHist<T1,T1,T2>* hist)
{
  VSRHistogram1D* output_hist = new VSRHistogram1D;

  output_hist->m_nBins = 
    lround((hist->hiLimit()-hist->loLimit())/hist->binSize());
  output_hist->m_loX    = hist->loLimit();
  output_hist->m_hiX    = hist->hiLimit();

  int offset = hist->bin(hist->loLimit());

  output_hist->m_hist.reset(new TH1F("","",
				     output_hist->m_nBins,
				     output_hist->m_loX,
				     output_hist->m_hiX) );
  output_hist->m_hist->SetDirectory(0);

  for(typename VERITAS::VSLimitedErrorsHist<T1,T1,T2>::iterator itr = 
	hist->begin(); itr != hist->end(); ++itr)
    {
      

      int bin = itr->bin() - hist->firstBin();

      if(itr->count() > 1E18)
	{
	  output_hist->set(bin-offset,0);
	  output_hist->setError(bin-offset,0);
	}
      else
	{
	  
	  output_hist->set(bin-offset,(double)itr->count());
	  output_hist->setError(bin-offset,(double)itr->err());
	}
    }

  return output_hist;
}


#endif // VSRHISTOGRAM1D_HPP
