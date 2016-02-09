/*! \file extract_egy_kernel.cpp

  Program to extract data required by mat_diagnostics so that file
  fits into memory in wasteful Matlab
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/25/2007

  $Id: extract_diagnostics.cpp,v 3.2 2010/01/01 20:47:54 matthew Exp $

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

using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] stage2_file" 
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  std::string output_file = "diagnostics.h5";
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

  std::string stage2_file = *argv;
  argv++,argc--;

  VSFileUtility::expandFilename(stage2_file);
  if(!output_file.empty())VSFileUtility::expandFilename(output_file);

  try
    {
      VSOctaveH5Reader reader(stage2_file);

      VSOctaveH5Writer* writer = 
	new VSOctaveH5Writer(output_file,true,
			     std::string("Diagnostics from: ")
			     + stage2_file);

      std::vector<std::string> fieldnames = reader.variables();
      unsigned nfield = fieldnames.size();
      for(unsigned ifield=0;ifield<nfield;ifield++)
	if(fieldnames[ifield] != "events")
	  reader.copyTo(fieldnames[ifield], writer);

      VSOctaveH5ReaderStruct* res = reader.readStruct("events");
      VSOctaveH5WriterStruct* wes = writer->writeStruct("events");

      res->copyTo("nscope_image",wes);
      res->copyTo("mean_derotated_fov_x",wes);
      res->copyTo("mean_derotated_fov_y",wes);
      res->copyTo("Rx",wes);
      res->copyTo("Ry",wes);

      unsigned ntel;
      VSOctaveH5ReaderCellVector* rtc = res->readCellVector("t");
      ntel = rtc->dimensions();
      VSOctaveH5WriterCellVector* wtc = wes->writeCellVector("t",ntel);
      for(unsigned itel=0;itel<ntel;itel++)
	{
	  VSOctaveH5ReaderStruct* rts = rtc->readStruct(itel);
	  VSOctaveH5WriterStruct* wts = wtc->writeStruct(itel);
	  rts->copyTo("has_image",wts);
	  rts->copyTo("fp_N",wts);
	  rts->copyTo("fp_xc",wts);
	  rts->copyTo("fp_yc",wts);
	  rts->copyTo("delta2l",wts);
	  rts->copyTo("d1",wts);
	  rts->copyTo("d2",wts);
	}
      
      delete writer;
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
