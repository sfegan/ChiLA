/*! \file stage3.cpp

  Program for calculating rates and significances.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/05/2007

  $Id: stage3.cpp,v 3.17 2010/01/12 20:24:10 matthew Exp $

*/

#include <iomanip>
#include <string>
#include <fstream>

#include <VSOptions.hpp>
#include <VSLineTokenizer.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSFileUtility.hpp>
#include <VSEventData.hpp>
#include <VSRunInfoData.hpp>
#include <VSSimulationData.hpp>
#include <VSSimpleGraph.hpp>
#include <VSDatumElementExtractor.hpp>
#include <VSAnalysisStage3.hpp>

using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream 
    << "Usage: " << progname 
    << " [options] stage2_file [stage2_file..] [file_list..]" 
    << std::endl << std::endl
    << "Arguments: The argument list can include one or more stage2 "
    << std::endl
    << "files and/or single column lists of stage2 files.  " << std::endl
    << std::endl << std::endl;
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

  options.addCatagory("file", "Input and output files");

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  bool compress = true;
  std::string output_file = "";

  VSAnalysisStage3::configure(options);

  options.findWithValue("o",
			output_file,
			"Name of stage3 output file.",
			"file");

  options.findBoolValue("compress", compress, true,
			"Compress data in output files where possible.",
			"file"); 

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
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  if(!output_file.empty()) 
    VSFileUtility::expandFilename(output_file);

  // --------------------------------------------------------------------------
  // Read the data files
  // --------------------------------------------------------------------------
  VSAnalysisStage3 stage3;

  // --------------------------------------------------------------------------
  // Insert the list of files to process into file_list
  // --------------------------------------------------------------------------
  std::list<std::string> file_list;
  for(int i=1;i<argc;i++) file_list.push_back(argv[i]);

  if(!file_list.size())
    {
      std::cerr << "Error: No input files." << std::endl;
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // Write the output HDF5 file
  // --------------------------------------------------------------------------
  if(output_file.empty())
    {
      if(file_list.size() == 1)
	{
	  VSOctaveH5Reader *reader = new VSOctaveH5Reader(*file_list.begin());

	  unsigned run_number;
	  reader->readScalar("stage1.run_info.run_number",run_number);
	  output_file = "x" + VSDataConverter::toString(run_number) + "_s3.h5";

	  delete reader;
	}
      else
	{
	  VSOctaveH5Reader *reader = new VSOctaveH5Reader(*file_list.begin());

	  std::string source_name;
	  reader->readString("observation.name",source_name);
	  output_file = source_name + "_s3.h5";

	  delete reader;
	}
    }
  
  VSOctaveH5Writer::setDefaultCompress(compress);
  VSOctaveH5Writer *writer = new VSOctaveH5Writer(output_file,true);

  stage3.runStage3(file_list,writer);

  // --------------------------------------------------------------------------
  // Write configuration etc..
  // --------------------------------------------------------------------------
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
  delete writer;
}
