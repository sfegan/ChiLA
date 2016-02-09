//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAStatistics.hpp
  Various statistics functions

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#ifndef VSASTATISTICS_HPP
#define VSASTATISTICS_HPP

#include <vector>

namespace VERITAS
{
  namespace VSAStatistics
  {
    double heleneUL(double e, double err, double cl);
    double limaSignificance(unsigned on, unsigned off, double alpha);
    double limaSignificance(std::vector< double >& on, 
			    std::vector< double >& off, 
			    std::vector< double >& alpha);
  }
}

#endif // VSASTATISTICS_HPP
