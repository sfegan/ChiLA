//-*-mode:c++; mode:font-lock;-*-

/*! \file dump_optics.cpp

  Dump array from the database

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \date    12/07/2004
  \version 0.2
  \note
*/

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>

#include <VSOptions.hpp>
#include <VSDBFactory.hpp>

#include "VSOTelescopeArray.hpp"

using namespace VERITAS;

int main(int argc, char ** argv)
{
  VSOptions options(argc,argv);
  VSOArrayParameters::setCanonicalValuesFromOptions(options);
  VSDBFactory::configure(&options);

  bool use_file = false;
  if(options.find("file") != VSOptions::FS_NOT_FOUND)use_file=true;

  char *progname = *argv;
  argv++, argc--;

  // --------------------------------------------------------------------------
  // Get the database name or filename from the command line arguements
  // --------------------------------------------------------------------------
  if(argc<2)
    {
      std::cerr << "Usage: " << progname << " database optics_id" 
		<< std::endl;
      return EXIT_FAILURE;
    }

  std::string database(*argv);
  argc--; argv++;

  unsigned optics_id = 0;
  VSDataConverter::fromString(optics_id, *argv);
  argc--; argv++;

  VSOTelescopeArray array;

  if(use_file)
    {
      std::ifstream stream(database.c_str());
      if(stream)array.readFromShortDump(stream);
    }
  else
    {
      // ----------------------------------------------------------------------
      // Use the factory to create the database
      // ----------------------------------------------------------------------
  
      VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
      db->useDatabase(database);

      // ----------------------------------------------------------------------
      // Retrieve an array configuration
      // ----------------------------------------------------------------------
      
      array.readFromDatabase(db, optics_id);
    }

  // --------------------------------------------------------------------------
  // Dump the array
  // --------------------------------------------------------------------------

  array.dumpShort(std::cout);
}
