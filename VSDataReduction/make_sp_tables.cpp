/*! \file make_sp_tables.cpp

  Program to make scaled parameter tables from stage2 file
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/25/2007

  $Id: make_sp_tables.cpp,v 3.26 2009/10/01 17:53:13 matthew Exp $

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
#include <VSSimEnergyWeightCalc.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <Angle.h>
#include <VSLTLibraryVisitor.hpp>

using namespace VERITAS;
using namespace SEphem;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] stage2_file [stage2_file..]" 
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

  std::string output_file = "sp_tables.h5";
  options.findWithValue("o", output_file, 
			"Set the output file name."); 

  VSScaledParameterLibraryVisitor::configure(options);


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
  if(argc < arg_req)
    {
      std::cerr << progname << ": need at least " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::list<std::string> file_list;
  for(int i=0;i<argc;i++) file_list.push_back(argv[i]);

  try
    {
      VSOctaveH5Writer* writer = new VSOctaveH5Writer(output_file,true);
      
      VSScaledParameterLibraryVisitor* sp_visitor =
	new VSScaledParameterLibraryVisitor;

      VSEventDataDispatcher* dispatcher = 
	new VSEventDataDispatcher(sp_visitor);
      
      dispatcher->processFiles(file_list);  
      
      sp_visitor->createLibrary();
      sp_visitor->save(writer);
      
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

      delete sp_visitor;
      delete dispatcher;
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





