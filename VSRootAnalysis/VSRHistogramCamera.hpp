//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRHISTOGRAMCAMERA_HPP
#define VSRHISTOGRAMCAMERA_HPP

#include <vector>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TPad.h>

#include "VSRHistogram.hpp"

class VSRHistogramCamera : public VSRHistogram
{

public:
  VSRHistogramCamera(double lo = 0, double hi = 0,
		     const std::string& title = "");
  ~VSRHistogramCamera();

  virtual void draw();
  virtual VSRHistogramCamera* clone() const;

  virtual void addLegendEntry(TLegend *leg, const std::string& label) const
  {
    return;
  }

  virtual double getLoX() const { return 0; }
  virtual double getHiX() const { return 1; }
  virtual double getLoY() const { return 0; }
  virtual double getHiY() const { return 1; }
  virtual double getPosLoX() const { return 0; }
  virtual double getPosHiX() const { return 1; }
  virtual double getPosLoY() const { return 0; }
  virtual double getPosHiY() const { return 1; }

  virtual bool is1D() const { return false; }
  virtual bool is2D() const { return false; }

  void setPixel(unsigned ichan, double value)
  {
    //   m_channel_values[ichan] = value;
    m_pixelValues[ichan] = value;
    m_isSuppressed[ichan] = false;
  }

  void setSuppressed(unsigned ichan, bool suppressed)
  {
    if(m_isSuppressed.size() <= ichan)
      m_isSuppressed.resize(ichan+1);

    m_isSuppressed[ichan] = suppressed;
  }

  double getMedian(std::vector< double > &pixelValues,
		   std::vector< bool > &isPixelSuppressed);
  double getMean(std::vector< double > &pixelValues,
		 std::vector< bool > &isPixelSuppressed);

  void smooth(double radius);

private:
  double      m_lo;
  double      m_hi;
  std::string m_title;

  std::vector< double > m_pixelValues;
  std::vector< bool >   m_isSuppressed;

  //  std::map< unsigned, double > m_channel_values;

};


#endif // VHISTOGRAMCAMERA_HPP
