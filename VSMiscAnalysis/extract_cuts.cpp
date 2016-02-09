/*! \file extract_h5struct.cpp

  Program to extract an h5 structure and write it to a new file
  
  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       08/25/2007

  $Id: extract_cuts.cpp,v 1.1 2008/12/04 02:53:37 matthew Exp $

*/

#include <VSAssert.hpp>
#include <VSOptions.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSOctaveH5Writer.hpp>
#include <VSFileUtility.hpp>
#include <VSNSpaceCutsCalc.hpp>
#include <VSSimpleCutsCalc.hpp>

using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] h5_file path" 
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

  std::string output_file = "nspace_cuts.h5";
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
  if(argc < arg_req)
    {
      std::cerr << progname << ": need at least " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  std::string file, path;

  file = *argv;
  argv++,argc--;
  if(argc) path = *argv;
  
  std::cout << file << " " << path << std::endl;
  
  VSOctaveH5Reader reader(file);
  VSOctaveH5ReaderStruct* s = reader.readStruct(path);

  VSOctaveH5ReaderStruct* simple_cuts_struct = 
    s->readStruct("simple_cuts");
  VSSimpleCutsCalc* simple_cuts = new VSSimpleCutsCalc;

  std::cout << simple_cuts_struct << std::endl;

  if(simple_cuts_struct) simple_cuts->load(simple_cuts_struct);
  delete simple_cuts_struct;

  VSOctaveH5ReaderStruct* nspace_cuts_struct = 
    s->readStruct("nspace_cuts");
  VSNSpaceCutsCalc* nspace_cuts = new VSNSpaceCutsCalc;
  if(nspace_cuts_struct) nspace_cuts->load(nspace_cuts_struct);
  delete nspace_cuts_struct;

  delete s;
  
  VSOctaveH5Writer* writer = new VSOctaveH5Writer(output_file,true);
  nspace_cuts->save(writer->writeStruct("nspace_cuts"));
  simple_cuts->save(writer->writeStruct("simple_cuts"));
  delete writer;
}
