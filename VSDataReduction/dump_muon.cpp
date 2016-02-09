//-*-mode:c++; mode:font-lock;-*-

/*! \file dump_muon.cpp

  Dump muon analysis entries to the console

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/15/2007

  $Id: dump_muon.cpp,v 3.3 2007/12/04 18:05:04 sfegan Exp $

*/

#include <stdlib.h>
#include <signal.h>

#include <VSOptions.hpp>
#include <VSOctaveIO.hpp>
#include <VSMuonAnalysisData.hpp>
#include <Angle.h>

using namespace SEphem;
using namespace VERITAS;

// ============================================================================
// MAIN
// ============================================================================

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
	 << " [options] stage2_file|muon_file [stage2_file|muon_file...]"
	 << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

void dumpMuonData(const VSMuonAnalysisData& data)
{
  unsigned nscope = data.scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nmuon = data.scope[iscope].size();
      for(unsigned imuon=0;imuon<nmuon;imuon++)
	{
	  const VSMuonAnalysisDatum& datum(data.scope[iscope][imuon]);
	  std::cout << std::left
		    << std::setw(6) << datum.ievent << ' '
		    << std::right
		    << iscope << ' ' 
		    << std::fixed
		    << std::setw(5) << std::setprecision(2) 
		    << Angle::toDeg(datum.x0) << ' ' 
		    << std::setw(5) << std::setprecision(2) 
		    << Angle::toDeg(datum.y0) << ' '
		    << std::setw(5) << std::setprecision(2) 
		    << Angle::toDeg(datum.r0) << ' '
		    << std::setw(7) << std::setprecision(4)
		    << datum.chi2 << ' ' 

		    << std::setw(3) 
		    << datum.c_nimage << ' '
		    << std::setw(7) << std::setprecision(1)
		    << datum.c_signal << ' '
		    << std::setw(6) << std::setprecision(4)
		    << Angle::toDeg(datum.c_rms) << ' '
		    << std::setw(5) << std::setprecision(2) 
		    << Angle::toDeg(datum.c_cx) << ' ' 
		    << std::setw(5) << std::setprecision(2) 
		    << Angle::toDeg(datum.c_cy) << ' '
		    << std::setw(6) << std::setprecision(4) 
		    << datum.c_xi << ' '
		    << std::setw(7) << std::setprecision(1)
		    << datum.c_U0_r0 << ' '
		    << std::setw(3) 

		    << datum.r_nimage << ' '
		    << std::setw(7) << std::setprecision(1)
		    << datum.r_signal << ' '
		    << std::setw(6) << std::setprecision(4)
		    << Angle::toDeg(datum.r_rms) << ' '
		    << std::setw(5) << std::setprecision(2) 
		    << Angle::toDeg(datum.r_cx) << ' ' 
		    << std::setw(5) << std::setprecision(2) 
		    << Angle::toDeg(datum.r_cy) << ' '
		    << std::setw(6) << std::setprecision(4) 
		    << datum.r_xi << ' '
		    << std::setw(7) << std::setprecision(1)
		    << datum.r_U0_r0 << ' '
		    << std::setw(7) << std::setprecision(1)
		    << datum.r_U0_r0_corr

#ifdef MUON_TEST_NON_UNIFORMITY
		    << ' '
		    << std::setw(7) << std::setprecision(4)
		    << datum.nu_xi << ' '
		    << std::setw(7) << std::setprecision(1)
		    << datum.nu_U0_r0 << ' '
		    << std::setw(7) << std::setprecision(1)
		    << datum.nu_U0_r0_corr
#endif
		    << '\n';
	}
    }
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

  if(argc==0)
    {
      std::cerr << "Need at lead one input file" << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_FAILURE);
    }

  for(int iarg=0;iarg<argc;iarg++)
    {
      try
	{
	  VSOctaveH5Reader reader(argv[iarg]);
	  if(reader.isStruct("muon_analysis"))
	    {
	      VSOctaveH5ReaderStruct* sr = reader.readStruct("muon_analysis");
	      vsassert(sr);
	      VSMuonAnalysisData data;
	      data.load(sr);
	      delete sr;
	      dumpMuonData(data);
	    }
	  else if(reader.isCellVector("muon_analysis"))
	    {
	      VSOctaveH5ReaderCellVector* c = 
		reader.readCellVector("muon_analysis");
	      vsassert(c);
	      unsigned nfile = c->dimensions();
	      for(unsigned ifile=0;ifile<nfile;ifile++)
		{
		  VSOctaveH5ReaderStruct* sr = c->readStruct(ifile);
		  vsassert(sr);
		  VSMuonAnalysisData data;
		  data.load(sr);
		  delete sr;
		  dumpMuonData(data);
		}
	      delete c;
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
}
