/*! \file combine_sp_tables.cpp

  Program to combine individual scaled parameter tables made with 
  "make_sp_tables" into a library
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/12/2007

  $Id: combine_sp_tables.cpp,v 3.1 2009/12/22 19:43:30 matthew Exp $

*/

#include <iomanip>
#include <fstream>
#include <sstream>
#include <cerrno>

#include <VSOptions.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSFileUtility.hpp>
#include <VSEventData.hpp>
#include <VSSimulationData.hpp>
#include <VSNSpace.hpp>
#include <VSGenSpace.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <Angle.h>

using namespace VERITAS;
using namespace SEphem;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] list_file" 
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
  stream << std::endl
	 << 
   "The input file (list_file) should contain a list of individual scaled\n"
   "parameter, energy lookup or energy kernel tables, such as those produced\n"
   "by \"make_sp_tables\", \"make_egy_tables\", and \"extract_kernel\",\n"
   "and the zenith angle and azimuth angles that they represent. If either\n"
   "angle is given as -1 the file will be used as default. In practice,\n"
   "only the azimuth angle should be set to default. For example:\n\n"
   "0    -1 base_g_wob05_2008_PL2.5_Zn00.0_Wob0.5_2007A01_NSB0.07_1M_scp.h5\n"
   "13.1 -1 base_g_wob05_2008_PL2.5_Zn13.1_Wob0.5_2007A01_NSB0.07_1M_scp.h5\n"
   "20.2 -1 base_g_wob05_2008_PL2.5_Zn20.2_Wob0.5_2007A01_NSB0.07_1M_scp.h5\n"
   "\n"
   "would combine three scaled parameter tables, corresponding to 0, 13.1\n"
   "and 20.2 degrees zenith angle into one library. Each of the table sets\n"
   "would apply at all azimuth angles."
	 << std::endl;
}

typedef triple<double,double,std::string> ListItem;

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

  std::string output_file = "sp_library.h5";
  options.findWithValue("o", output_file, "Set the output file name."); 

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

  std::string list_file = *argv;
  argv++,argc--;

  VSFileUtility::expandFilename(list_file);
  if(!output_file.empty())VSFileUtility::expandFilename(output_file);

  std::ifstream stream(list_file.c_str());
  if(!stream)
    {
      std::cerr << progname << ": could not open list file: " << list_file
		<< std::endl
		<< progname << ": " << strerror(errno) << std::endl;
      exit(EXIT_FAILURE);
    }

  try
    {
      VSOctaveH5Writer* writer =
	new VSOctaveH5Writer(output_file,true,
			     std::string("Scaled parameter table from: ")
			     + list_file);

      VSScaledParameterLibraryWriter wlibrary(writer);

      VSOctaveH5WriterCompositeVector<ListItem>* liw =
	writer->writeCompositeExpandableVector<ListItem>("list_items");
      
      std::string line;
      std::getline(stream,line);
      while(stream)
	{
	  std::istringstream linestream(line);
	  std::getline(stream,line);
	  double zn;
	  double az;
	  double ped_rms;
	  std::string file;
	  linestream >> zn >> az >> ped_rms >> file;
	  if(file.empty())continue;
	  ListItem li(zn,az,file);
	  liw->append(li);
	  
	  VSOctaveH5Reader* reader = new VSOctaveH5Reader(file);
	  VSScaledParameterLibraryReader rlibrary(reader);
	  
	  rlibrary.copyTables(wlibrary, zn, az, ped_rms);
	}

      delete liw;
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
