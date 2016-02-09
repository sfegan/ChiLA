/*! \file s2opt.cpp

  Program to print out or compare options from stage2 file
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/12/2007

  $Id: s2opt.cpp,v 3.5 2007/12/04 18:05:04 sfegan Exp $

*/

#include <iomanip>
#include <iostream>
#include <map>
#include <string>

#include <VSOptions.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSFileUtility.hpp>

using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] stage2_file [stage2_file]" 
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
  stream << std::endl
	 << 
    "Print list of configuration options and their values from a stage2 file\n"
    "or compare two stage2 files and print only options which differ."
	 << std::endl;
}

typedef std::map<std::string,std::string> OptMap;

OptMap getOptMap(std::string filename)
{
  VSFileUtility::expandFilename(filename);
  VSOctaveH5Reader reader(filename);
  OptMap option;
  VSOctaveH5ReaderCellVector* cv = reader.readCellVector("configuration");
  vsassert(cv);
  unsigned nrecord = cv->dimensions();
  for(unsigned irecord=0;irecord<nrecord;irecord++)
    {
      VSOctaveH5ReaderStruct* rs = cv->readStruct(irecord);
      std::string key;
      std::string val;
      rs->readString("key", key);
      rs->readString("val", val);
      option[key]=val;
    }
  return option;
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

  if(argc == 1)
    {
      OptMap opt = getOptMap(*argv);
      for(OptMap::const_iterator iopt=opt.begin();iopt!=opt.end();iopt++)
	std::cout << std::setw(30) << std::left
		  << iopt->first << ' ' 
		  << iopt->second << '\n';
    }
  else if(argc == 2)
    {
      OptMap opt1 = getOptMap(*argv++);
      OptMap opt2 = getOptMap(*argv);

      std::cout << "Differences in options:" << '\n';
      for(OptMap::const_iterator iopt=opt1.begin();iopt!=opt1.end();iopt++)
	if((opt2.find(iopt->first) != opt2.end())
	   &&(iopt->second != opt2[iopt->first]))
	  std::cout << std::setw(30) << std::left
		    << iopt->first << ' ' 
		    << iopt->second << " <---> "
		    << opt2[iopt->first] << '\n';
      
      std::cout << "\nOptions only in first file:" << '\n';
      for(OptMap::const_iterator iopt=opt1.begin();iopt!=opt1.end();iopt++)
	if(opt2.find(iopt->first) == opt2.end())
	  std::cout << std::setw(30) << std::left
		    << iopt->first << ' ' 
		    << iopt->second << '\n';

      std::cout << "\nOptions only in second file:" << '\n';
      for(OptMap::const_iterator iopt=opt2.begin();iopt!=opt2.end();iopt++)
	if(opt1.find(iopt->first) == opt1.end())
	  std::cout << std::setw(30) << std::left
		    << iopt->first << ' ' 
		    << iopt->second << '\n';
    }
  else
    {
      std::cerr << progname << ": need 1 or 2 arguments, got "
		<< argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }
}
