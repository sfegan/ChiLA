//-*-mode:c++; mode:font-lock;-*-

/*! \file VQColormapFactory.hpp

  Class to return a vector of colors that can be used as a colormap

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/30/2005

  \see VQMCArrayData
  \see VQMCChannelRenderer
*/

#include <vector>
#include <string>

#include <Qt/QtGui>

#ifndef VQCOLORMAP_HPP
#define VQCOLORMAP_HPP

namespace VERITAS
{

  typedef std::vector<QColor> Colormap;

  class VQColormapFactory
  {
  public:
    static Colormap getColormap(const std::string name);
    static Colormap jet();
    static Colormap hot();
    static Colormap hsv(unsigned entries=64, double sat=1.0, double val=1.0);
    static Colormap blue_red_yellow();
  };

}

#endif // VQCOLORMAP_HPP
