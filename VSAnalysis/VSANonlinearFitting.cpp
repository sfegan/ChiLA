//-*-mode:c++; mode:font-lock;-*-

/*! \file VSALinearLeastSquares.cpp

  Generalized Linear Least Squares Fitting

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n        

  \version    0.1
  \date       01/20/2009
*/

#include <vector>
#include <VSAAlgebra.hpp>
#include <VSANonlinearFitting.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAMath;
using namespace VERITAS::VSAAlgebra;

#ifdef TEST_MAIN

// g++ -O3 -g -Wall -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DUSEALLOCA -D H5_USE_16_API -fPIC -I/usr/include/mysql -I../VSUtility -I../VSSimDB -I../Physics -I../VSShower -I../VSOptics -I../VSAnalysis -I../VSNSpace -I../VSDataReduction -I../SEphem -I. -I../VSCommon -I/usr/local/veritas/include -DTEST_MAIN -o test VSANonlinearFitting.cpp VSASVD.cpp ../Physics/RandomNumbers_TNG.o VSAMath.o VSAAlgebra.o -L../VSCommon -lpthread -lVSCommon

#include<iostream>
#include<iomanip>
#include<VSAMath.hpp>
#include<RandomNumbers.hpp>
#include<VSACoord.hpp>

using namespace VSACoord;
using namespace VERITAS;

int main(int argc, char**argv)
{
  VSAFunction::Constant<Coord2D> fn;
  VSAMath::Data<Coord2D> data;

  VSAMath::DataPoint<Coord2D> p;

  data.insert(p);

  {
    VSAMath::FitLM< VSAFunction::Constant<Coord2D>,Coord2D > fitter;
    
    fitter.set(data,fn);

  }

}

#endif
