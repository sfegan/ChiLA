/*! \file extract_egy_kernel.cpp

  Program to extract energy kernel and effective area curves from stage3
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/25/2007

  $Id: extract_egy_kernel.cpp,v 3.6 2009/05/13 06:04:03 matthew Exp $

*/

#include <iomanip>
#include <fstream>

#include <VSOptions.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSOctaveH5Writer.hpp>
#include <VSFileUtility.hpp>
#include <VSNSpace.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSResultsData.hpp>
#include <VSResultsSimData.hpp>

using namespace VERITAS;
using namespace SEphem;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] stage3_file" 
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  std::vector<std::string> command_line(argc);
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)command_line[iarg]=argv[iarg];

  std::cout << "Command line:";
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
    std::cout << ' ' << argv[iarg];
  std::cout << std::endl;
 
  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  std::string output_file = "egy_kernel.h5";
  options.findWithValue("o", output_file, 
			"Set the output file name."); 

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
        if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  int arg_req = 1;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string stage3_file = *argv;
  argv++,argc--;

  VSFileUtility::expandFilename(stage3_file);
  if(!output_file.empty())VSFileUtility::expandFilename(output_file);

  try
    {
      VSOctaveH5Reader reader(stage3_file);

      if(!reader.isStruct("sim"))
	{
	  std::cerr << progname << ": " << stage3_file
		    << ": no \"sim\" structure" << std::endl;
	  exit(EXIT_SUCCESS);
	}

      VSStage3SimDatum sim;
      sim.load(reader.readStruct("sim"));

      VSOctaveH5Writer* writer = 
	new VSOctaveH5Writer(output_file,true,
			     std::string("Energy kernel from: ")
			     + stage3_file);

      VSScaledParameterLibraryWriter* library =
	new VSScaledParameterLibraryWriter(writer);

      library->write(sim.egy_kernel,
		     VSScaledParameterLibraryWriter::SPS_ENERGY_KERNEL,
		     true, 0, 0, true, 0, 0, true, 0);

      library->write(sim.effarea,
		     VSScaledParameterLibraryWriter::SPS_EFFECTIVE_AREA,
		     true, 0, 0, true, 0, 0, true, 0);

      delete library;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cout << "Caught instance of VSOctaveH5Exception" << std::endl
		<< e.message() << std::endl;
      exit(EXIT_FAILURE);
    }
  catch(const std::exception& x)
    {
      std::cout << "Caught instance of " << x.what() << std::endl;
      exit(EXIT_FAILURE);
    }
  catch(...)
    {
      std::cout << "Caught some exception" << std::endl;
      exit(EXIT_FAILURE);
    }

  H5close();
}
