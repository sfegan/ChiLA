/*! \file logain.cpp

  Program which flatfields the camera.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/25/2006

  $Id: flatfield.cpp,v 3.2 2008/10/19 02:08:51 matthew Exp $

*/

#include<vector>
#include<string>
#include<iostream>
#include<vsassert>
#include<fstream>

#include <SphericalCoords.h>

#include "VSOptions.hpp"
#include "VBFSimplePeds.hpp"
#include "VSSimpleHist.hpp"
#include "VBFLaserCalc.hpp"
#include "VSHiLoData.hpp"
#include "VSMiscellaneousDBData.hpp"
#include "VSFileUtility.hpp"
#include "VSChannelMap.hpp"

using namespace VERITAS;

// ============================================================================
// MAIN
// ============================================================================

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
	 << " [options] stage1_file" 
	 << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

double get_avg_hv(const std::vector<VSMiscellaneousDBData::OneHVMeasScope>& hv,
		  unsigned ichan)
{
  if(hv.empty()) return 0;

  double sum = 0;
  double n = 0;

  for(std::vector<VSMiscellaneousDBData::OneHVMeasScope>::const_iterator itr =
	hv.begin(); itr != hv.end(); ++itr)
    {
      if(itr->chan[ichan].voltage == 0) return 0;

      sum += itr->chan[ichan].voltage;
      n++;
    }

  return sum/n;
}

struct ChanData
{
  ChanData(): chan(), gain(), absgain(), hv(), uncalibrated() {}
  ChanData(unsigned _chan, double _gain, double _absgain, double _hv,
	   bool _uncalibrated): 
    chan(_chan), gain(_gain), absgain(_absgain), hv(_hv), 
    uncalibrated(_uncalibrated) {}

  unsigned chan;
  double   gain;
  double   absgain;
  double   hv;
  bool     uncalibrated;

  static bool sort_absgain(const ChanData & lhs, const ChanData& rhs)
  {
    return lhs.absgain < rhs.absgain;
  }

  double calc_hv(double _gain,double _gamma)
  {
    return exp((1./_gamma)*log(_gain/gain)+log(hv));
  }
};

int main(int argc, char** argv)
{
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

      // ----------------------------------------------------------------------
      // Flatfield options
      // ----------------------------------------------------------------------

      double gamma = 7.4;
      std::vector<double> scope_absgain(4,5.0);
      double gain_exclude_lo = 0.5;
      double gain_exclude_hi = 1.5;
      std::string pixel_exclude_file("");
      std::vector< std::pair<unsigned,unsigned> > pixel_exclude;

      options.findWithValue("pixel_exclude_file", pixel_exclude_file,
			    "File containing a list of pixels in each "
			    "telescope that should be excluded from "
			    "flat-fielding.");

      options.findWithValue("gain_exclude_lo", gain_exclude_lo,
			    "Set the gain threshold below which "
			    "pixels will be excluded from flatfielding.");

      options.findWithValue("gain_exclude_hi", gain_exclude_hi,
			    "Set the gain threshold above which "
			    "all pixels will be excluded from flatfielding.");

      options.findWithValue("gamma", gamma,
			    "Set the value of the exponent for the gain-"
			    "voltage relation.");

      options.findWithValue("absgain", scope_absgain,
			    "Set the target absolute gain of each telescope "
			    "in DC/PE.");

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

      std::string stage1_file = *argv;
      argv++,argc--;


      VSAnalysisStage1Data stage1;

      if(!VSFileUtility::isFile(stage1_file))
	{
	  std::cerr << "Error: " << stage1_file
		    << " is not a valid file." << std::endl;
	  exit(EXIT_FAILURE);
	}
      else if(VSOctaveH5ReaderStruct::isHDF5(stage1_file))
	{
	  VSOctaveH5Reader *reader = new VSOctaveH5Reader(stage1_file);
	  stage1.load(reader->readStruct("stage1"));
	  delete reader;
	}


      if(VSFileUtility::isFile(pixel_exclude_file))
	{
	  std::ifstream file(pixel_exclude_file.c_str());
	  unsigned scope, pixel;

	  while(file >> scope >> pixel)
	    pixel_exclude.push_back(std::make_pair(scope,pixel));	  
	}

      std::cout << std::setw(6) << "SCOPE"
		<< std::setw(6) << "PIXEL"
		<< std::setw(10) << "GAIN"
		<< std::setw(10) << "ABSGAIN"
		<< std::setw(10) << "OLDHV"
		<< std::setw(10) << "NEWHV"
		<< std::setw(10) << "HVDIFF"
		<< std::endl;

      const unsigned nscope = stage1.run_info.nchan.size();
      for(unsigned iscope = 0; iscope < nscope; iscope++)
	{
	  std::vector<VSMiscellaneousDBData::OneHVMeasScope>& scope_hv = 
	    stage1.misc_db->scope[iscope].hv_status;

	  std::vector< ChanData > sort_data;
	  std::vector< ChanData > all_data;

	  const unsigned nchan = stage1.run_info.nchan[iscope];
	  for(unsigned ichan = 0; ichan < nchan; ichan++)
	    {
	      double avg_hv = get_avg_hv(scope_hv,ichan);
	      double gain = stage1.laser->scope[iscope].chan[ichan].gain;
	      double absgain = 
		stage1.laser->scope[iscope].chan[ichan].absgain;
	      bool uncalibrated = 
		stage1.laser->scope[iscope].chan[ichan].uncalibrated;

	      if(!stage1.laser->scope[iscope].chan[ichan].uncalibrated)
		sort_data.push_back(ChanData(ichan,gain,absgain,avg_hv,
					     uncalibrated));

	      all_data.push_back(ChanData(ichan,gain,absgain,avg_hv,
					  uncalibrated));
	    }

	  std::sort(sort_data.begin(),sort_data.end(),ChanData::sort_absgain);
	  ChanData ref_pixel = *(sort_data.begin()+sort_data.size()/2);

	  double ref_gain = 
	    ref_pixel.gain*scope_absgain[iscope]/ref_pixel.absgain;
	  
	  std::string date = 
	    VSDataConverter::toString(stage1.run_info.first_event_time.
				       getDBTimeStamp()).substr(0,8);

	  std::string hv_filename = 
	    "T" + VSDataConverter::toString(iscope+1) + "_" + date + "_" +
	    VSDataConverter::toString(stage1.run_info.run_number) + ".hv";

	  std::ofstream hv_file(hv_filename.c_str());

	  for(std::vector<ChanData>::iterator itr = all_data.begin();
	      itr != all_data.end(); ++itr)
	    {
	      if(itr->chan+1 == 500) continue; // skip pixel 500

	      bool exclude = false;

	      for(std::vector<std::pair<unsigned,unsigned> >::iterator 
		    exclude_itr = pixel_exclude.begin(); 
		  exclude_itr != pixel_exclude.end(); exclude_itr++)
		{
		  if(iscope+1 == exclude_itr->first && 
		     itr->chan+1 == exclude_itr->second)
		    {
		      exclude = true;
		      break;
		    }
		}

	      int old_hv = lround(itr->hv);
	      int hv = 0;
	      if(itr->uncalibrated || exclude || 
		 itr->gain < gain_exclude_lo ||
		 itr->gain > gain_exclude_hi) 
		hv = lround(itr->hv);
	      else hv = lround(itr->calc_hv(ref_gain,7.0));

	      std::cout 
		<< std::setw(6) << iscope+1 
		<< std::setw(6) << itr->chan+1 
		<< std::setw(10) << itr->gain
		<< std::setw(10) << itr->absgain
		<< std::setw(10) << old_hv
		<< std::setw(10) << hv
		<< std::setw(10) << hv-old_hv
		<< std::endl;

	      hv_file << std::setw(6) << itr->chan+1 
		      << std::setw(10) << hv
		      << std::endl;
	    }
	}     
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
}


