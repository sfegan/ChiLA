/*! \file optimize_msc.cpp

  Calculate optimal theta^2, msc_width and msc_length parameter values
  from a bunch of files
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/06/2007

  $Id: optimize_msc.cpp,v 1.39 2009/11/20 02:50:26 matthew Exp $

*/

#include <iomanip>
#include <fstream>

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
#include <VSScaledParameterCalc.hpp>
#include <VSRunInfoData.hpp>
#include <VSSimEnergyWeightCalc.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSCutOptimizer.hpp>
#include <VSSimpleCutOptimizer.hpp>
#include <VSCutsEvaluator.hpp>
#include <VSLineTokenizer.hpp>

using namespace VERITAS;
using namespace SEphem;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream 
    << "Usage: " << progname 
    << " [options] stage2_file [stage2_file..] [file_list..]" 
    << std::endl << std::endl
    << "Arguments: The argument list can include a list of stage2 "
    << std::endl
    << "files or text files containing lists of stage2 files." << std::endl
    << std::endl << std::endl;

  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  std::vector<std::string> command_line(argc);
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)command_line[iarg]=argv[iarg];

  std::cerr << "Command line:";
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
    std::cerr << ' ' << argv[iarg];
  std::cerr << std::endl;
 
  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  std::string output_file = "opt_tables.h5";
  options.findWithValue("o", output_file, 
			"Set the output file name."); 

  VSCutOptimizerFactory::configure(options);
  VSCutsEvaluator::configure(options);

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
        if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  if(print_usage)
    {
      usage(progname, options, std::cerr);
      exit(EXIT_SUCCESS);
    }

  argv++,argc--;

  int arg_req = 1;
  if(argc < arg_req)
    {
      std::cerr << progname << ": need at least " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_SUCCESS);
    }

  // --------------------------------------------------------------------------
  // Insert the list of files to process into file_list
  // --------------------------------------------------------------------------
  std::list<std::string> file_list;
  for(int i=0;i<argc;i++) file_list.push_back(argv[i]);

  try
    {
      VSOctaveH5Writer* writer = NULL;
      
      if(!output_file.empty())
	writer = new VSOctaveH5Writer(output_file,true);

      VSCutOptimizer* optimizer = 
	VSCutOptimizerFactory::getInstance()->createCutOptimizer();

      VSEventDataDispatcher* dispatcher = 
	new VSEventDataDispatcher(optimizer);
      
      dispatcher->processFiles(file_list);  
      
      optimizer->optimize();

      if(writer) optimizer->save(writer);

      // ----------------------------------------------------------------------
      // Write configuration etc..
      // ----------------------------------------------------------------------
      VSOctaveH5WriterCellVector* wc = 
	writer->writeCellVector("command_line", command_line.size());
      for(unsigned iarg=0;iarg<command_line.size();iarg++)
	wc->writeString(iarg,command_line[iarg]);
      delete wc;
      
      unsigned noptionrecords = options.getOptionRecords().size();
      wc = writer->writeCellVector("configuration", noptionrecords);
      for(unsigned irecord=0;irecord<noptionrecords;irecord++)
	{
	  VSOctaveH5WriterStruct* ws = wc->writeStruct(irecord);
	  ws->writeString("key", options.getOptionRecords()[irecord].key);
	  ws->writeScalar("status", 
			  int(options.getOptionRecords()[irecord].status));
	  ws->writeScalar("val_requested", 
			  options.getOptionRecords()[irecord].val_requested);
	  ws->writeString("val", options.getOptionRecords()[irecord].val);
	  delete ws;
	}
      delete wc;  

      std::vector<std::string> files = dispatcher->getFiles();

      unsigned nfile = files.size();
      wc = writer->writeCellVector("files", nfile);
      for(unsigned ifile = 0; ifile < nfile; ifile++)
	{
	  wc->writeString(ifile,files[ifile]);
	}
      delete wc;

      delete dispatcher;
      delete optimizer;      
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
}
