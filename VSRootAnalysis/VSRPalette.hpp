//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRPALETTE_HPP
#define VSRPALETTE_HPP

#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>

#include "VSAssert.hpp"

class VSRPalette
{
public:
  VSRPalette(unsigned ncont, unsigned nrgb): 
    m_ncont(ncont), m_nrgb(nrgb), m_colors(nrgb), m_stops(nrgb)
  { }

  virtual ~VSRPalette() { }
  


  unsigned nrgb() const { return m_nrgb; }
  unsigned ncont() const { return m_ncont; }

  int getColor(unsigned icolor) const
  {
    vsassert(icolor < m_nrgb);
    return m_colors[icolor];
  }

  float getStop(unsigned icolor) const
  {
    vsassert(icolor < m_nrgb);
    return m_stops[icolor];
  }

  void setColor(unsigned icolor, int color, float stop) 
  {
    vsassert(icolor < m_nrgb);
    m_colors[icolor] = color; 
    m_stops[icolor] = stop; 
  }

  virtual void setPalette()
  {
    const unsigned nrgb = VSRPalette::nrgb();    
    double red[nrgb];
    double green[nrgb];
    double blue[nrgb];
    double stops[nrgb];
    
    for(unsigned irgb = 0; irgb < nrgb; irgb++)
      {
	TColor *c = gROOT->GetColor(getColor(irgb));
	float rgb[3];
	c->GetRGB(rgb[0],rgb[1],rgb[2]);

	red[irgb]   = rgb[0];
	green[irgb] = rgb[1];
	blue[irgb]  = rgb[2];
	stops[irgb] = getStop(irgb);
      }

    gStyle->SetNumberContours(ncont());
    TColor::CreateGradientColorTable(nrgb, stops, red, green, blue, ncont());
  }

  static void setPalette(const std::string& palette_name, 
			 double lo, double hi);


private:
  unsigned             m_ncont;
  unsigned             m_nrgb;
  std::vector< int >   m_colors;
  std::vector< float > m_stops;
};


class VSRPaletteSkyMap : public VSRPalette
{
public:
  VSRPaletteSkyMap(double lo, double hi): VSRPalette(255,5) 
  {     
    double zp = fabs(lo)/(hi-lo);

    setColor(0,kBlue+4,0.0);
    setColor(1,kBlue,zp);
    setColor(2,kRed,zp + (1-zp)*0.3);
    setColor(3,kOrange,zp + (1-zp)*0.6);
    setColor(4,kYellow,1.0);
  }
};

class VSRPaletteDefault : public VSRPalette
{
public:
  VSRPaletteDefault(double lo, double hi);
};

#endif
