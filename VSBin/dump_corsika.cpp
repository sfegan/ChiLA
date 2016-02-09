//-*-mode:c++; mode:font-lock;-*-

/*! \file dump_corsika.cpp

  Dump CORSIKA output file to console

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    1.0
  \date       03/02/2005
*/

#include <iostream>
#include <vector>
#include <cmath>

#include "VSCORSIKAEvent.hpp"

using namespace VERITAS;

int main(int argc, char**argv)
{
  char* progname = *argv;
  argc--,argv++;

  VSCORSIKAFileDumper visitor(std::cout);
  VSCORSIKAEventDispatcher dispatcher(&visitor);

  while(argc)
    {
      dispatcher.processFile(*argv);
      argc--,argv++;
    }
}
