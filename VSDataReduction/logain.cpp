/*! \file logain.cpp

  Low gain analysis swiss army knife program

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/25/2006

  $Id: logain.cpp,v 3.2 2009/07/27 18:21:37 tarlen Exp $

*/

#include<vector>
#include<string>
#include<iostream>
#include<vsassert>

#include <SphericalCoords.h>

#include "VSOptions.hpp"
#include "VBFHiLoCalc.hpp"

using namespace VERITAS;


// ============================================================================
// MAIN
// ============================================================================

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
	 << " [options] vbf_file" 
	 << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  //signal(SIGSEGV,sigabort);

  try
    {
      std::string progname(*argv);
      std::vector<std::string> command_line(argc);
      for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
	command_line[iarg]=argv[iarg];

      VSOptions options(argc, argv, true);

      bool print_usage = false;
      if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
	print_usage=true;
      if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
	print_usage=true;

      std::string output_file("x?_hilo.h5");
      options.findWithValue("o", output_file, 
			"Set the output file name. If the output file is "
			"empty a name defined by the run number will be used, "
			"for example: x12345_hilo.h5. The first instance of "
			"a question mark in the name will be replaced with "
			"the run number.");

      // ----------------------------------------------------------------------
      // Hi-Lo specific configuration
      // ----------------------------------------------------------------------

      int sample_0 = -10;
      unsigned sample_N = 0;
      unsigned nevent_min = 10;

      options.findWithValue("sample0", sample_0,
			    "Set the pedestal integration window start. A "
			    "negative value specifies that the window should "
			    "start some number of samples from the end.");

      options.findWithValue("sampleN", sample_N,
			    "Set the pedestal integration window size. A "
			    "value of zero means that the window size should "
			    "be automatically chosen based on the size of the "
			    "recorded trace and \"sample_0\".");

      options.findWithValue("nevent_min", nevent_min,
			    "Minimum number of events required for each "
			    "channel.");

      // ----------------------------------------------------------------------
      // Finish option processing and get arguments
      // ----------------------------------------------------------------------

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
	  usage(progname, options, std::cerr);
	  exit(EXIT_SUCCESS);
	}

      if(argc!=1)
	{
	  std::cerr << progname << ": need exactly one argumenent (got " 
		    << argc << ')' << std::endl << std::endl;
	  usage(progname, options, std::cerr);
	  exit(EXIT_FAILURE);
	}

      std::string vbf_file = *argv;
      argv++,argc--;

      // ----------------------------------------------------------------------
      // Run the calculator
      // ----------------------------------------------------------------------

      VBFHiLoCalc hilo_calc(sample_0, sample_N);
      VSSimpleVBFDispatcher dispatcher(&hilo_calc);
      dispatcher.processFile(vbf_file.c_str());      // Runs the analysis!

      // ----------------------------------------------------------------------
      // Get the results
      // ----------------------------------------------------------------------
      
      VSHiLoData hilo_data;
      hilo_calc.getData(hilo_data,nevent_min);
      
      // ----------------------------------------------------------------------
      // Print some stuff
      // ----------------------------------------------------------------------

      const unsigned nscope = hilo_data.scope.size();
      for(unsigned iscope=0; iscope<nscope; iscope++)
	if(!hilo_data.scope[iscope].chan.empty())
	  {
	    const unsigned nchan = hilo_data.scope[iscope].chan.size();
	    unsigned ncal = 0;
	    for(unsigned ichan=0; ichan<nchan; ichan++)
	      if(hilo_data.scope[iscope].chan[ichan].has_lo_gain_ped)
		ncal++;
	    std::cout << 'T' << iscope+1 << ": " << ncal
		      << " low gain peds calculated\n";
	  }

      // ----------------------------------------------------------------------
      // Write data
      // ----------------------------------------------------------------------

      if(!output_file.empty())
	VSFileUtility::replaceQuestionWithNumber(output_file, hilo_data.runno),
	  VSFileUtility::expandFilename(output_file);

      time_t thetime = time(NULL);
      std::string thetimestring(ctime(&thetime));
      std::string comment =
	std::string("# VERITAS Hi/Lo pedetsal data, written by ") + progname +
	std::string(" - ") + thetimestring.substr(0,thetimestring.length()-1);
      
      VSOctaveH5Writer writer(output_file, true, comment);
      VSOctaveH5WriterStruct* hilo_struct = writer.writeStruct("hilo");
      hilo_data.save(hilo_struct);
      delete hilo_struct;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cout << "Caught instance of VSOctaveH5Exception" << std::endl
		<< e.message() << std::endl;
      exit(EXIT_FAILURE);
    }
  catch(const VSAssert& a)
    {
      std::cout << a.what() << std::endl;
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
