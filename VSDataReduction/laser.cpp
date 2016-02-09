//-*-mode:c++; mode:font-lock;-*-

/**
 * \class laser.cpp
 * \ingroup tools
 * \brief Calculate laser gains and timing offsets
 *
 * 
 * Calculate timing offsets and gains using laser flashes.
 *
 * Original Author: Stephen Fegan
 * $Author: matthew $
 * $Date: 2008/09/30 03:21:38 $
 * $Revision: 3.11 $
 * $Tag$
 *
 * $Id: laser.cpp,v 3.11 2008/09/30 03:21:38 matthew Exp $
 *
 **/

#include<vector>
#include<string>
#include<iostream>
#include<vsassert>

#include <SphericalCoords.h>

#include "VSOptions.hpp"
#include "VBFSimplePeds.hpp"
#include "VSSimpleHist.hpp"
#include "VBFLaserCalc.hpp"
#include "VSLaserData.hpp"
#include "VSAnalysisStage1.hpp"
#include "VSFileUtility.hpp"
#include "VSChannelMap.hpp"

static const char*const VERSION = 
  "$Id: laser.cpp,v 3.11 2008/09/30 03:21:38 matthew Exp $";

static const char*const REVISION = 
  "$Revision: 3.11 $";

static const char* const LOGO =
  "  ____      ________    _ __    ___     ____    \n"
  " / / /     / ____/ /_  (_) /   /   |    \\ \\ \\   \n"
  "/ / /     / /   / __ \\/ / /   / /| |     \\ \\ \\          __ _ _  \n"
  "\\ \\ \\    / /___/ / / / / /___/ ___ |     / / /  |   /\\ (_ |_|_) \n"
  " \\_\\_\\   \\____/_/ /_/_/_____/_/  |_|    /_/_/   |__/--\\__)|_| \\ \n"
  "                                                                 \n"
  "\n"
  "$Id: laser.cpp,v 3.11 2008/09/30 03:21:38 matthew Exp $\n";

using namespace VERITAS;
using namespace SEphem;

// ============================================================================
// MAIN
// ============================================================================

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname << " [options] filename" << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

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

      std::string profile = "veritas";
      options.findWithValue("profile", profile,
			    "Set the default values for a given profile. The "
			    "profiles presently defined are: \"veritas\", "
			    "\"simulations\" and \"whipple\"");

      bool no_db = false;
      bool no_l3 = false;
      bool no_threads = false;
      std::string output_file("");
      double threshold_frac = 0.5;
      double threshold_dc = 40;
      unsigned threshold_nchan = 400;
      unsigned threshold_nscope = 1;
      unsigned threshold_nflash = 50;
      bool integral_timing = false;
      bool fill_hist = true;
      bool include_low_gain = false;
      unsigned sample_separation = 3;
      double ped_suppress_hi = 3.0;
      double ped_suppress_lo = 1/ped_suppress_hi;
      double singlepe_dev = 0.469;

#ifdef NOVDB
      VSDBFactory::configure(&options);
#endif

      VSAnalysisStage1::configure(options,profile,"");

      triple<std::string, std::string, double> 
	array_pos("31d40m29.04s","-110d57m10.08s",1270);

      options.findWithValue("array_pos", array_pos,
			    "Set the reference position for the array. The "
			    "coordinate triple should be specified as "
			    "latitude/longitude/elevation.");

      Angle latitude;
      Angle longitude;
      latitude.setFromDMSString(array_pos.first);
      longitude.setFromDMSString(array_pos.second);
      SphericalCoords earth_position;
      earth_position.setLatLong(latitude,longitude);
      //  double earth_elevation = array_pos.third;

      options.findBoolValue("no_db", no_db, true,
			    "Do not collect information from the database.");
      options.findBoolValue("no_l3", no_l3, true, 
			    "Do not use L3 time.");

#ifndef NOTHREADS
      options.findBoolValue("no_threads", no_threads, true,
			"Disable multi-threaded operation. This may slow the "
			"analysis on machines with multiple processors, "
			"multiple cores, hyperthreading and when reading data "
			"from a slow source (such as NFS).");
#endif

      options.findWithValue("o", output_file, 
			"Set the output file name. If the output file is "
			"empty a name defined by the run number will be used, "
			"for example: x12345_s1.h5.");

      if(!output_file.empty())VSFileUtility::expandFilename(output_file);

      // ----------------------------------------------------------------------
      // Laser-specific configuration
      // ----------------------------------------------------------------------

      options.findBoolValue("integral", integral_timing, true,
		        "Use integral timing calculator.");
      if(integral_timing==true)threshold_frac=0.2;

      options.findWithValue("threshold_frac", threshold_frac,
			"Set the fraction of the pulse at which "
			"the time calculation is done.");

      options.findWithValue("threshold_dc", threshold_dc,
			"Set the per channel laser flash threshold (in DC).");

      options.findWithValue("threshold_nchan", threshold_nchan,
			"Set the threshold number of channels required to "
			"for an event to be considered a laser flash.");

      options.findWithValue("threshold_nscope", threshold_nscope,
			"Set the threshold number of telescopes meeting the "
			"channel number criterion required for an event to be "
			"considered a laser flash.  Setting this to zero is "
			"equivalent to requiring that all of the telescopes "
			"in the array be present.");

      options.findWithValue("threshold_nflash", threshold_nflash,
			"Set the threshold number of flashes required for "
			"a telescope to be considered as calibrated.");

      options.findBoolValue("fill_hist",fill_hist, true,
			"Fill offset histogram.");

      options.findBoolValue("include_low_gain",include_low_gain, true,
			"Include LOW gain channels in statistics.");
  
      options.findWithValue("ped_sample_separation", sample_separation,
			"Specify the number of samples between independent "
			"smaples of the pedestal variance from within the "
			"same trace. This value should be approximately the "
			"single PE pulse width, in samples");

      options.findWithValue("ped_suppress_lo", ped_suppress_lo,
                        "Fraction of median pedestal variance to below which "
                        "channel is suppressed.");

      options.findWithValue("ped_suppress_hi", ped_suppress_hi,
                        "Fraction of median pedestal variance to above which "
                        "channel is suppressed.");

      options.findWithValue("singlepe_dev", singlepe_dev,
                        "Standard deviation of single PE pulse amplitude "
			"as a fraction of the mean amplitude.  Used to "
			"estimate the absolute gain and relative efficiency "
			"of each pixel.");

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

      std::string vbf_file = *argv;
      argv++,argc--;

      // ----------------------------------------------------------------------
      // Calculate stage 1 information
      // ----------------------------------------------------------------------

      VSSimpleVBFDispatcher::catchSignalAndStopProcessingFile(SIGINT);
      VSAnalysisStage1 stage1_calc;
      VSAnalysisStage1Data stage1_data;
      stage1_calc.runStage1(vbf_file, stage1_data, earth_position,
			    no_db, no_l3, no_threads);

      // ----------------------------------------------------------------------
      // Calculate laser data
      // ----------------------------------------------------------------------

      VSChannelMap chmap(stage1_data.run_info.lo_event_time);
      std::vector<std::vector<unsigned> > boards_per_crate(chmap.nscope());
      std::vector<std::vector<unsigned> > l2_pulse_channel(chmap.nscope());
      for(unsigned iscope=0;iscope<boards_per_crate.size();iscope++)
	{
	  unsigned ncrate = chmap.ncratesForScope(iscope);
	  boards_per_crate[iscope].resize(ncrate);
	  l2_pulse_channel[iscope].resize(ncrate);
	  for(unsigned icrate=0;icrate<ncrate;icrate++)
	    {
	      boards_per_crate[iscope][icrate] = 
		chmap.nboardsForCrate(iscope,icrate);
	      l2_pulse_channel[iscope][icrate] = 
		chmap.l2ChannelForCrate(iscope,icrate);
	    }
	}

      VSTimingCalc* calc = 0;
      if(integral_timing)calc = new VSIntegralTimingCalc(threshold_frac);
      else calc = new VSPeakTimingCalc(threshold_frac);

      std::cout << LOGO << std::endl;
      VBFLaserCalc laser_calc(calc, 
			      threshold_dc,threshold_nchan,threshold_nscope,
			      threshold_nflash,
			      ped_suppress_lo, ped_suppress_hi, singlepe_dev,
			      boards_per_crate, l2_pulse_channel,
			      stage1_data, fill_hist, include_low_gain);

      VSSimpleVBFDispatcher dispatcher(&laser_calc);
      dispatcher.processFile(vbf_file.c_str());

      VSLaserData* laser = new VSLaserData;
      laser_calc.getData(*laser);

      unsigned nscope = laser->scope.size();
      for(unsigned iscope = 0; iscope < nscope; iscope++)
	std::cout 
	  << "T" << iscope+1 << " ************************" << std::endl
	  << " Num Flashes:      " 
	  << laser->scope[iscope].nflash << std::endl
	  << " Num Chan Mean:    " 
	  << laser->scope[iscope].nchan_flash_hist.mean() << std::endl
	  << " Num LoGain Mean:  " 
	  << laser->scope[iscope].nchan_logain_hist.mean() << std::endl
	  << " Num HiGain Mean:  " 
	  << laser->scope[iscope].nchan_higain_hist.mean() << std::endl
	  << " Amplitude Mean:   " 
	  << laser->scope[iscope].signal_mean << std::endl 
	  << " Amplitude Dev:    " 
	  << laser->scope[iscope].signal_dev << std::endl 
	  << std::endl;

      // ----------------------------------------------------------------------
      // Write stage 1 data
      // ----------------------------------------------------------------------

      std::ostringstream s1_filename;
      s1_filename << 'x' << stage1_data.run_info.run_number << "_s1.h5";
      if(output_file.empty())output_file=s1_filename.str();
  
      time_t thetime = time(NULL);
      std::string thetimestring(ctime(&thetime));
      std::string comment =
	std::string("# VERITAS analysis stage 1, written by ") + progname +
	std::string(" - ") + thetimestring.substr(0,thetimestring.length()-1);

      VSOctaveH5Writer writer(output_file, true, comment);
      stage1_data.laser = laser;
      VSOctaveH5WriterStruct* stage1_struct = writer.writeStruct("stage1");
      stage1_data.save(stage1_struct);
      delete stage1_struct;

      // ----------------------------------------------------------------------
      // Write configuration etc..
      // ----------------------------------------------------------------------

      VSOctaveH5WriterCellVector* wc = 
	writer.writeCellVector("command_line", command_line.size());
      for(unsigned iarg=0;iarg<command_line.size();iarg++)
	wc->writeString(iarg,command_line[iarg]);
      delete wc;

      unsigned noptionrecords = options.getOptionRecords().size();
      wc = writer.writeCellVector("configuration", noptionrecords);
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

  H5close();
}
