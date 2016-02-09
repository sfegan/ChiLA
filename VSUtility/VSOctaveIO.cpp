//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOctaveIO.cpp

  Class to write to GNU/Octave file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/20/2006
*/

#include<cctype>
#include<vsassert>
#include<stdint.h>
#include<string>
#include<sstream>

#include"VSOctaveIO.hpp"

using namespace VERITAS;

// ============================================================================
// EXCEPTION
// ============================================================================

VSOctaveH5Exception::~VSOctaveH5Exception()
{
  // nothing to see here
}

// ============================================================================
// VSOCTAVEH5COMMON
// ============================================================================


unsigned VSOctaveH5Common::cellPrecision(unsigned nel)
{
  //vsassert(nel>0);
  //nel--;
  unsigned name_precision = 0;
  while(nel)name_precision++, nel/=10 ;
  return name_precision;
}
