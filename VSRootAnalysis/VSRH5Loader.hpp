//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRH5LOADER_HPP
#define VSRH5LOADER_HPP

#include <iostream>

#include <TMatrixT.h>
#include <TVectorT.h>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSSimple2DHist.hpp>
#include <VSSimpleHist.hpp>
#include <VSSimpleGraph.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSNSpace.hpp>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram2D.hpp"
#include "VSRHistogram1D.hpp"
#include "VSRGraph1D.hpp"

class VSRH5Loader
{
public:
  static TVectorT<double> createVector(const VERITAS::VSNSpace& nspace);  
  static TMatrixT<double> createMatrix(const VERITAS::VSNSpace& nspace);  
  static VERITAS::VSNSpace createNSpace(const TVectorT<double>& v,
					double lo, double hi);

};

#endif // VSRH5LOADER_HPP
