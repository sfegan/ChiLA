//-*-mode:c++; mode:font-lock;-*-

/**
 * \class toffset.cpp
 * \ingroup tools
 * \brief Calculate timing offsets
 *
 * Calculate channel-to-channel timing offsets for any (one) telescope
 * in the file.
 *
 * Original Author: Stephen Fegan
 * $Author: sfegan $
 * $Date: 2007/12/04 18:05:04 $
 * $Revision: 3.2 $
 * $Tag$
 *
 * $Id: toffset.cpp,v 3.2 2007/12/04 18:05:04 sfegan Exp $
 *
 **/

#include<vector>
#include<string>
#include<iostream>
#include<vsassert>

#include "VSOptions.hpp"
#include "VBFSimplePeds.hpp"
#include "VSSimpleHist.hpp"
#include "VBFTOffsetCalc.hpp"

using namespace VERITAS;

#define HIST_RES 0.5

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
  std::string progname(*argv);
  VSOptions options(argc, argv, true);

  std::vector<unsigned> boards_per_crate_0;
  boards_per_crate_0.push_back(13);
  boards_per_crate_0.push_back(12);
  boards_per_crate_0.push_back(13);
  boards_per_crate_0.push_back(12);
  std::vector<unsigned> l2_pulse_channel_0;
  l2_pulse_channel_0.push_back(128);
  l2_pulse_channel_0.push_back(249);
  l2_pulse_channel_0.push_back(259);
  l2_pulse_channel_0.push_back(498);
  std::vector<unsigned> boards_per_crate_1;
  boards_per_crate_1.push_back(13);
  boards_per_crate_1.push_back(12);
  boards_per_crate_1.push_back(13);
  boards_per_crate_1.push_back(12);
  std::vector<unsigned> l2_pulse_channel_1;
  l2_pulse_channel_1.push_back(128);
  l2_pulse_channel_1.push_back(249);
  l2_pulse_channel_1.push_back(259);
  l2_pulse_channel_1.push_back(498);
  std::vector<unsigned> boards_per_crate_2;
  boards_per_crate_2.push_back(13);
  boards_per_crate_2.push_back(12);
  boards_per_crate_2.push_back(13);
  boards_per_crate_2.push_back(12);
  std::vector<unsigned> l2_pulse_channel_2;
  l2_pulse_channel_2.push_back(128);
  l2_pulse_channel_2.push_back(249);
  l2_pulse_channel_2.push_back(259);
  l2_pulse_channel_2.push_back(498);
  std::vector<unsigned> boards_per_crate_3;
  boards_per_crate_3.push_back(13);
  boards_per_crate_3.push_back(12);
  boards_per_crate_3.push_back(13);
  boards_per_crate_3.push_back(12);
  std::vector<unsigned> l2_pulse_channel_3;
  l2_pulse_channel_3.push_back(128);
  l2_pulse_channel_3.push_back(249);
  l2_pulse_channel_3.push_back(259);
  l2_pulse_channel_3.push_back(498);

  double threshold_frac = 0.5;
  double threshold_dc = 40;
  unsigned threshold_nchan = 100;
  unsigned threshold_nflash = 50;
  bool integral_timing = false;
  bool fill_hist = false;
  bool include_low_gain = false;
  std::string ped_file = "";
  unsigned sample_separation = 3;
  double ped_suppress_hi = 2.0;
  double ped_suppress_lo = 1/ped_suppress_hi;

  bool print_usage = false;

  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  options.findWithValue("t0_boards_per_crate", boards_per_crate_0,
			"Number of boards per crate in T1. To be specified as "
			"a comma seperated list.");

  options.findWithValue("t0_l2_channel", l2_pulse_channel_0,
			"ID of channels containing L2 pulse in T1. To be "
			"specified as a comma seperated list of channel IDs "
			"(counting from ZERO). If any channel ID is out of "
			"range then L2 signal compensation is disabled.");

  options.findWithValue("t1_boards_per_crate", boards_per_crate_1,
			"Number of boards per crate in T2. To be specified as "
			"a comma seperated list.");

  options.findWithValue("t1_l2_channel", l2_pulse_channel_1,
			"ID of channels containing L2 pulse in T2. To be "
			"specified as a comma seperated list of channel IDs "
			"(counting from ZERO). If any channel ID is out of "
			"range then L2 signal compensation is disabled.");

  options.findWithValue("t2_boards_per_crate", boards_per_crate_2,
			"Number of boards per crate in T3. To be specified as "
			"a comma seperated list.");

  options.findWithValue("t2_l2_channel", l2_pulse_channel_2,
			"ID of channels containing L2 pulse in T3. To be "
			"specified as a comma seperated list of channel IDs "
			"(counting from ZERO). If any channel ID is out of "
			"range then L2 signal compensation is disabled.");

  options.findWithValue("t3_boards_per_crate", boards_per_crate_3,
			"Number of boards per crate in T4. To be specified as "
			"a comma seperated list.");

  options.findWithValue("t3_l2_channel", l2_pulse_channel_3,
			"ID of channels containing L2 pulse in T4. To be "
			"specified as a comma seperated list of channel IDs "
			"(counting from ZERO). If any channel ID is out of "
			"range then L2 signal compensation is disabled.");
  
  std::vector<std::vector<unsigned> > boards_per_crate;
  boards_per_crate.push_back(boards_per_crate_0);
  boards_per_crate.push_back(boards_per_crate_1);
  boards_per_crate.push_back(boards_per_crate_2);
  boards_per_crate.push_back(boards_per_crate_3);
  std::vector<std::vector<unsigned> > l2_pulse_channel;
  l2_pulse_channel.push_back(l2_pulse_channel_0);
  l2_pulse_channel.push_back(l2_pulse_channel_1);
  l2_pulse_channel.push_back(l2_pulse_channel_2);
  l2_pulse_channel.push_back(l2_pulse_channel_3);
  
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

  options.findWithValue("threshold_nflash", threshold_nflash,
			"Set the threshold number of flashes required for "
			"a telescope to be considered as calibrated.");
  
  options.findBoolValue("fill_hist",fill_hist, true, "Fill offset histogram.");

  options.findBoolValue("include_low_gain",include_low_gain, true,
			"Include LOW gain channels in statistics.");

  options.findWithValue("ped_sample_separation", sample_separation,
			"Specify the number of samples between independent "
			"smaples of the pedestal variance from within the "
			"same trace. This value should be approximately the "
			"single PE pulse width, in samples");

  options.findWithValue("load_peds", ped_file,
                        "Load pedestals information from file rather than "
                        "calculating them.");

  options.findWithValue("ped_suppress_lo", ped_suppress_lo,
                        "Fraction of median pedestal variance to below which "
                        "channel is suppressed.");

  options.findWithValue("ped_suppress_hi", ped_suppress_hi,
                        "Fraction of median pedestal variance to above which "
                        "channel is suppressed.");

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

  std::string filename(*argv);
  argv++,argc--;

  // ==========================================================================
  // CALCULATE PEDESTALS
  // ==========================================================================

  VSSimplePedData ped_data;
  bool ped_loaded = false;

  if(!ped_file.empty())
    {
      ped_loaded = ped_data.load(ped_file);
    }

  if(!ped_loaded)
    {
      std::cerr << "Calculating pedestals..." << std::endl;
      VBFSimplePeds ped_calc(sample_separation);
      try
        {
          VSSimpleVBFDispatcher dispatcher(&ped_calc);
          dispatcher.processFile(filename.c_str());
        }
      catch (const std::exception& e)
        {
          std::cerr << e.what() << std::endl;
          return EXIT_FAILURE;
        }

      ped_calc.getSimplePedData(ped_data);
    }

  ped_data.suppress(ped_suppress_lo, ped_suppress_hi);

  // ==========================================================================
  // CALCULATE LASER DATA
  // ==========================================================================

  VSTimingCalc* calc = 0;
  if(integral_timing)calc = new VSIntegralTimingCalc(threshold_frac);
  else calc = new VSPeakTimingCalc(threshold_frac);

  std::cerr << "Calculating time offsets.." << std::endl;
  VBFTOffsetCalc toffset(calc, threshold_dc,threshold_nchan,
			 boards_per_crate, l2_pulse_channel,
			 ped_data, fill_hist, include_low_gain,
			 HIST_RES);
  try 
    {
      VSSimpleVBFDispatcher dispatcher(&toffset);
      dispatcher.processFile(filename.c_str());
    }
  catch (const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      return EXIT_FAILURE;
    }

  unsigned nscope = toffset.m_data.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(toffset.m_data[iscope] != 0)
      {
	VBFTOffsetCalc::ScopeData* data = toffset.m_data[iscope];

	std::vector<std::vector<std::pair<bool, double> > > median_calc;
	median_calc.resize(boards_per_crate[iscope].size());
	unsigned nchan = data->chan.size();
	
	std::cout << "# number of events     : " 
		  << data->stat_mean_nflash.count() << std::endl;
	
	std::cout << "# number of flashes    : " 
		  << data->stat_mean_signal.count() << std::endl;
	
	std::cout << "# mean nchan per event : " 
		  << data->stat_mean_nflash.mean() << " +/- "
		  << data->stat_mean_nflash.dev() << std::endl;
	
	std::cout << "# mean signal          : " 
		  << data->stat_mean_signal.mean() << " +/- " 
		  << data->stat_mean_signal.dev() << std::endl;

	unsigned npedchan = ped_data.m_data[iscope].size();
	unsigned nsuppressed = 0;
	for(unsigned ichan=0;ichan<npedchan;ichan++)
	  if(ped_data.m_data[iscope][ichan].suppressed)nsuppressed++;

	std::cout << "# num ped suppressed   : " 
		  << nsuppressed << std::endl;

	nsuppressed = 0;
	for(unsigned ichan=0;ichan<npedchan;ichan++)
	  if((!ped_data.m_data[iscope][ichan].suppressed)
	     &&(data->chan[ichan].stat_time.count() < threshold_nflash))
	    {
	      data->chan[ichan].global_suppressed=true;
	      nsuppressed++;
	    }

	std::cout << "# num flash suppressed : " 
		  << nsuppressed << std::endl;
	
	if(data->stat_mean_signal.count()<threshold_nflash)
	  {
	    std::cerr << "T-" << iscope << " has only " 
		      << data->stat_mean_signal.count() << " flashes. Need "
		      << threshold_nflash << std::endl;
	    continue;
	  }

	for(unsigned ichan=0; ichan<nchan; ichan++)
	  {
	    unsigned icrate = data->chan[ichan].global_crate;
	    vsassert(icrate < median_calc.size());
	    std::pair<bool, double> chan;
	    chan.first = !data->chan[ichan].global_suppressed;
	    chan.second = data->chan[ichan].stat_time.mean();
	    median_calc[icrate].push_back(chan);
	  }
	
	std::vector<double> median_time;
	median_time.resize(boards_per_crate[iscope].size());
	for(unsigned icrate=0; icrate<median_time.size(); icrate++)
	  median_time[icrate] = median(median_calc[icrate]);
	
	unsigned ncrate = data->stat_mean_crate_time.size();
	for(unsigned icrate=0; icrate<ncrate; icrate++)
	  std::cout << "# crate " << icrate << " mean time    : " 
		    << data->stat_mean_crate_time[icrate].mean() << " +/- " 
		    << data->stat_mean_crate_time[icrate].dev() 
		    << " (" << median_time[icrate] << ')' << std::endl;

	for(unsigned ichan=0; ichan<nchan; ichan++)
	  {
	    unsigned icrate = data->chan[ichan].global_crate;
	    unsigned l2chan = nchan+1;
	    if(icrate < l2_pulse_channel[iscope].size())
	      l2chan = l2_pulse_channel[iscope][icrate];

	    double dev = ped_data.m_data[iscope][ichan].
	      dev[ped_data.m_data[iscope][ichan].dev.size()-1];

	    std::cout /*  0 */ << iscope << ' '
	              /*  1 */ << ichan << ' '
		      /*  2 */ << ped_data.m_data[iscope][ichan].ped << ' '
	              /*  3 */ << dev << ' '
		      /*  4 */ << data->chan[ichan].global_suppressed << ' '
		      /*  5 */ << data->chan[ichan].stat_time.mean() << ' '
		      /*  6 */ << data->chan[ichan].stat_time.dev() << ' '
		      /*  7 */ << data->chan[ichan].stat_signal.mean() << ' '
		      /*  8 */ << data->chan[ichan].stat_signal.dev() << ' '
		      /*  9 */ << data->stat_mean_crate_time[icrate].mean() << ' '
		      /* 10 */ << l2chan << ' '
		      /* 11 */ << ((l2chan<nchan)?
				   data->chan[l2chan].stat_time.mean():0.0);
	    
	    if(fill_hist)
	      for(double itime=0; itime<double(data->max_nsample);
		  itime+=HIST_RES)
		std::cout << ' ' << data->chan[ichan].hist_time.
		  countForVal(itime);
	    
	    std::cout << std::endl;
	  }
      }
}
