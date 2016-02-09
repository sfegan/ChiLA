//-*-mode:c++; mode:font-lock;-*-

/*! \file VSUncopyable.hpp

  Header for all uncopyable classes

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       12/62/2006
*/

#ifndef VSUNCOPYABLE_HPP
#define VSUNCOPYABLE_HPP

namespace VERITAS
{

  class VSUncopyable
  {
  protected:
    VSUncopyable() { /* nothing to see here */ }
  private:
    VSUncopyable(const VSUncopyable&);
    VSUncopyable& operator=(const VSUncopyable&);
  };

}

#endif // VSUNCOPYABLE_HPP
