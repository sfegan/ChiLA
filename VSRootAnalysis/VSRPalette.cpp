#include "VSRPalette.hpp"

void VSRPalette::setPalette(const std::string& palette_name, 
			    double lo, double hi)
{
  VSRPalette* palette;
  
  if(palette_name == "skymap")
    palette = new VSRPaletteSkyMap(lo,hi);
  else 
    palette = new VSRPaletteDefault(lo,hi);

  palette->setPalette();
  
  delete palette;
}


VSRPaletteDefault::VSRPaletteDefault(double lo, double hi): VSRPalette(255,6) 
{     
  setColor(0,kViolet,0.0);
  setColor(1,kBlue,0.2);
  setColor(2,kCyan,0.4);
  setColor(3,kGreen,0.6);
  setColor(4,kYellow,0.8);
  setColor(5,kRed,1.0);
}
