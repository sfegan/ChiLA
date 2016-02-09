//-*-mode:c++; mode:font-lock;-*-

/*! \file collect_muon.cpp

  Collect muon analysis entries from muliple runs into one file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/15/2007

  $Id: collect_muon.cpp,v 3.3 2007/08/23 17:36:13 sfegan Exp $

*/

#include <stdlib.h>
#include <signal.h>

#include <VSOptions.hpp>
#include <VSOctaveIO.hpp>
#include <VSMuonAnalysisData.hpp>

using namespace VERITAS;

// ============================================================================
// MAIN
// ============================================================================

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
	 << " [options] stage2_file [stage2_file...]"
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

  std::string output_file = "muons.h5";
  options.findWithValue("o", output_file,
			"Set the ouput file name");
    
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

  if(argc==0)
    {
      std::cerr << "Need at lead one input file" << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_FAILURE);
    }

  VSOctaveH5Writer writer(output_file, true);
  VSOctaveH5WriterCellVector* c = writer.writeCellVector("muon_analysis",argc);

  for(int iarg=0;iarg<argc;iarg++)
    {
      try
	{
	  VSOctaveH5Reader reader(argv[iarg]);
	  if(reader.isStruct("muon_analysis"))
	    {
	      VSOctaveH5ReaderStruct* sr = reader.readStruct("muon_analysis");
	      VSMuonAnalysisData data;
	      data.load(sr);
	      delete sr;
	      
	      VSOctaveH5WriterStruct* sw = c->writeStruct(iarg);
	      data.save(sw);
	      delete sw;
	    }
	  else 
	    {
	      std::cerr << argv[iarg] 
			<< ": does not contain muon analysis data"
			<< std::endl;
	    }
	}
      catch(const VSOctaveH5Exception& e)
	{
	  std::cerr << argv[iarg] 
		    << ": caught exception: " << e.message() << std::endl;
	}
    }

  delete c;
}
