#include "VQColormap.hpp"

using namespace VERITAS;

Colormap VQColormapFactory::jet() 
{
  std::vector<QColor> c;
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.5625 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.6250 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.6875 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.7500 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.8125 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.8750 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.9375 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0625 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.1250 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.1875 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.2500 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.3125 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.3750 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.4375 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.5000 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.5625 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.6250 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.6875 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.7500 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.8125 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.8750 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.9375 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 1.0000 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 0.0625 , 1.0000 , 0.9375 ));
  c.push_back(QColor::fromRgbF( 0.1250 , 1.0000 , 0.8750 ));
  c.push_back(QColor::fromRgbF( 0.1875 , 1.0000 , 0.8125 ));
  c.push_back(QColor::fromRgbF( 0.2500 , 1.0000 , 0.7500 ));
  c.push_back(QColor::fromRgbF( 0.3125 , 1.0000 , 0.6875 ));
  c.push_back(QColor::fromRgbF( 0.3750 , 1.0000 , 0.6250 ));
  c.push_back(QColor::fromRgbF( 0.4375 , 1.0000 , 0.5625 ));
  c.push_back(QColor::fromRgbF( 0.5000 , 1.0000 , 0.5000 ));
  c.push_back(QColor::fromRgbF( 0.5625 , 1.0000 , 0.4375 ));
  c.push_back(QColor::fromRgbF( 0.6250 , 1.0000 , 0.3750 ));
  c.push_back(QColor::fromRgbF( 0.6875 , 1.0000 , 0.3125 ));
  c.push_back(QColor::fromRgbF( 0.7500 , 1.0000 , 0.2500 ));
  c.push_back(QColor::fromRgbF( 0.8125 , 1.0000 , 0.1875 ));
  c.push_back(QColor::fromRgbF( 0.8750 , 1.0000 , 0.1250 ));
  c.push_back(QColor::fromRgbF( 0.9375 , 1.0000 , 0.0625 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.9375 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.8750 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.8125 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.7500 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.6875 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.6250 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.5625 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.5000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.4375 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.3750 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.3125 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.2500 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.1875 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.1250 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.0625 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.9375 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.8750 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.8125 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.7500 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.6875 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.6250 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.5625 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.5000 , 0.0000 , 0.0000 ));
  return c;
}

Colormap VQColormapFactory::hot() 
{
  std::vector<QColor> c;
  c.push_back(QColor::fromRgbF( 0.0417 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.0833 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.1250 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.1667 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.2083 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.2500 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.2917 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.3333 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.3750 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.4167 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.4583 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.5000 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.5417 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.5833 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.6250 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.6667 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.7083 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.7500 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.7917 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.8333 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.8750 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.9167 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 0.9583 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.0417 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.0833 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.1250 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.1667 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.2083 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.2500 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.2917 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.3333 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.3750 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.4167 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.4583 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.5000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.5417 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.5833 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.6250 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.6667 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.7083 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.7500 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.7917 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.8333 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.8750 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.9167 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.9583 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.0625 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.1250 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.1875 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.2500 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.3125 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.3750 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.4375 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.5000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.5625 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.6250 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.6875 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.7500 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.8125 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.8750 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.9375 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 1.0000 ));
  return c;
}

Colormap VQColormapFactory::blue_red_yellow() 
{
  std::vector<QColor> c;
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 0.2000 ));
  c.push_back(QColor::fromRgbF( 0.0000 , 0.0000 , 1.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 0.0000 , 0.0000 ));
  c.push_back(QColor::fromRgbF( 1.0000 , 1.0000 , 0.3000 ));
  return c;
}

Colormap VQColormapFactory::hsv(unsigned entries, double sat, double val)
{
  std::vector<QColor> c;
  for(unsigned ientry=0; ientry<entries;ientry++)
    c.push_back(QColor::fromHsvF(double(ientry)/double(entries),sat,val));
  return c;
}

Colormap VQColormapFactory::getColormap(const std::string name)
{
  if(name=="jet")return jet();
  if(name=="hot")return hot();
  if(name=="blue_red_yellow")return blue_red_yellow();
  if(name=="hsv")return hsv();
  return Colormap();
}
