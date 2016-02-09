/*! \file optimize_msc.cpp

  Test a set of cuts, calculated with optimize_msc or otherwise, on a
  bunch of files
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/06/2007

  $Id: stat.cpp,v 1.1 2007/12/14 23:48:38 sfegan Exp $

*/

#include <iomanip>
#include <fstream>
#include <limits>

#include <VSOptions.hpp>
#include <VSFileUtility.hpp>
#include <VSSimpleStat.hpp>

using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname << " [options] file" << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);

#if 0
  std::vector<std::string> command_line(argc);
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)command_line[iarg]=argv[iarg];

  std::cerr << "Command line:";
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
    std::cerr << ' ' << argv[iarg];
  std::cerr << std::endl;
 #endif

  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  double cut_hi = std::numeric_limits<double>::infinity();
  double cut_lo = -std::numeric_limits<double>::infinity();
  
  options.findWithValue("cut_hi", cut_hi,
			"Set the upper bound of data to accept.");

  options.findWithValue("cut_lo", cut_lo,
			"Set the lower bound of data to accept.");

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

  std::vector<double> data;
  data.reserve(1000);
  VSSimpleStat2<double> stat;

  if(argc==0)
    {
      while(std::cin)
	{
	  double x;
	  if((std::cin >> x)&&(x>=cut_lo)&&(x<=cut_hi))
	    data.push_back(x), stat.accumulate(x);
	}
    }
  while(argc)
    {
      std::string file = *argv;
      argv++,argc--;
      VSFileUtility::expandFilename(file);
      
      std::ifstream stream(file.c_str());
      while(stream)
	{
	  double x;
	  if((std::cin >> x)&&(x>=cut_lo)&&(x<=cut_hi))
	    data.push_back(x), stat.accumulate(x);
	}
    }

  std::cout << stat.count() << ' ' << stat.mean() << ' ' << stat.dev() << ' '
	    << median(data) << '\n';
}
