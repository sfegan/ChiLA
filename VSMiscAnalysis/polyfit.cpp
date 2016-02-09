/*! \file polyfit.cpp

  Polynomial fit to a data set
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       12/04/2007

  $Id: polyfit.cpp,v 1.5 2009/04/01 00:52:36 matthew Exp $

*/

#include<iostream>
#include<cstdlib>
#include<sstream>
#include<fstream>
#include<string>
#include<cerrno>

#include<VSOptions.hpp>
#include<VSAData.hpp>
#include<VSALinearLeastSquares.hpp>

using namespace VERITAS;

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] [filename]" 
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

  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  bool dump_data = false;
  options.findBoolValue("dump_data",dump_data,true,"Dump data.");

  bool quiet = false;
  options.findBoolValue("q",quiet,true,"Be quiet.");

  unsigned poly_order = 0;
  options.findWithValue("n", poly_order,
			"Specify the polynomial order to fit.");
  
  triple<double,double,double> func_domain;
  options.findWithValue("func", func_domain,
			"Output the values of the fitted function over some "
			"range, rather than the polynomial coefficients. "
			"Specify the range as xlo/dx/xhi.");

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

  if(argc > 1)
    {
      std::cerr << progname << ": need at most 1 argument, got " 
		<< argc << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  const char* filename = 0;
  
  if(argc)
    {
      filename = *argv;
      argc--,argv++;
    }
  
  std::istream* stream = &std::cin;  
  std::ifstream* filestream = 0;

  if(filename)stream = filestream = new std::ifstream(filename);
  
  if(!stream->good())
    {
      std::cerr << progname << ": " << filename << ": " << strerror(errno)
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  VSAMath::Data<double> data;

  std::string line;
  std::getline(*stream, line);
  while(*stream)
    {
      std::istringstream linestream(line);
      double x;
      double y;
      double s = 1.0;
      if(linestream >> x >> y)
	{
	  linestream >> s;
	  data.insert(VSAMath::DataPoint<double>(x,y,s));
	}
      std::getline(*stream, line);
    }
  delete filestream;

  if(dump_data)
    {
      for(unsigned idata=0;idata<data.size();idata++)
	std::cout << data[idata].x << ' ' << data[idata].y << ' '
		  << data[idata].sigma << '\n';
      exit(EXIT_FAILURE);
    }

  if(!quiet)
    {
      if(filename)std::cerr << filename;
      else std::cerr << "[stdin]";
      std::cerr << ": loaded " << data.size() << " points\n";
    }

  VSAAlgebra::MatrixND cov;
  VSAAlgebra::VecND p;
  double chi2;

  try
    {
      chi2 = VSAMath::PolyFit::fit(poly_order,data,p,&cov);
    }
  catch(const std::exception& x)
    {
      std::cerr << progname << ": exception: " << x.what() << std::endl;
      exit(EXIT_FAILURE);
    }

  if(!quiet)
    std::cerr<< "chi^2: " << chi2 << " [reduced: " 
	     << chi2/double(data.size()-p.size()) << "]\n";

  if((func_domain.first<func_domain.third)&&(func_domain.second>0.0))
    {
      if(!quiet)
	for(unsigned ip=0;ip<p.size();ip++)
	  std::cerr << "x^" << ip << ": " << p[ip] << " +/- "
		    << std::sqrt(cov(ip,ip)) << '\n';

      double x = func_domain.first;
      while(x<func_domain.third)
	{
	  std::cout << x << ' ' << VSAMath::PolyFit::val(p,x) << '\n';
	  x += func_domain.second;
	}
    }
  else
    {
      for(unsigned ip=0;ip<p.size();ip++)
	std::cout << p[ip] << ' ' << std::sqrt(cov(ip,ip)) << '\n';
    }
}
