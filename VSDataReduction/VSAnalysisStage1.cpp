//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage1.cpp

  Stage 1 analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/12/2006

  $Id: VSAnalysisStage1.cpp,v 3.18 2009/06/24 19:59:31 matthew Exp $

*/

#include <Angle.h>
#include <SphericalCoords.h>
#include <VSAnalysisStage1.hpp>

#include <VSFileUtility.hpp>
#include <VSSimpleVBF.hpp>
#include <VBFRunInfo.hpp>
#include <VBFPrintFrequency.hpp>
#include <VBFTimeDepNSB.hpp>
#include <VSPointing.hpp>
#include <VBFDumper.hpp>
#include <VBFTimeDepPeds.hpp>

static const char*const VERSION = 
  "$Id: VSAnalysisStage1.cpp,v 3.18 2009/06/24 19:59:31 matthew Exp $";

static const char*const REVISION = 
  "$Revision: 3.18 $";

static const char*const LOGO = 
  "  ____      ________    _ __    ___     ____    \n"
  " / / /     / ____/ /_  (_) /   /   |    \\ \\ \\   \n"
  "/ / /     / /   / __ \\/ / /   / /| |     \\ \\ \\   __               \n"
  "\\ \\ \\    / /___/ / / / / /___/ ___ |     / / /  (__|_ _. _  _  /| \n"
  " \\_\\_\\   \\____/_/ /_/_/_____/_/  |_|    /_/_/   __)|_(_|(_|(/_  | \n"
  "                                                        ._|       \n"
  "\n"
  "$Id: VSAnalysisStage1.cpp,v 3.18 2009/06/24 19:59:31 matthew Exp $\n";

#define vstream std::cout

using namespace VERITAS;
using namespace SEphem;

unsigned VSAnalysisStage1::s_default_print_frequency          = 10000;
unsigned VSAnalysisStage1::s_default_npackets                 = 0;
bool     VSAnalysisStage1::s_default_no_pedestals_in_core     = false;
unsigned VSAnalysisStage1::s_default_nsb_window_start         = 0;
unsigned VSAnalysisStage1::s_default_nsb_window_width         = 2;
unsigned VSAnalysisStage1::s_default_events_per_slice         = 1000;
std::vector<double> 
         VSAnalysisStage1::s_default_nsb_suppress_dev           (4,0);
unsigned VSAnalysisStage1::s_default_nsb_min_events_per_slice = 400;
bool     VSAnalysisStage1::s_default_no_nsb_suppress          = false;
double   VSAnalysisStage1::s_default_ped_slice_time           = 180;
unsigned VSAnalysisStage1::s_default_ped_min_window           = 1;
unsigned VSAnalysisStage1::s_default_ped_sample_separation    = 3;
bool     VSAnalysisStage1::s_default_ped_no_suppress          = false;

VSAnalysisStage1::VSAnalysisStage1()
  : m_print_frequency(s_default_print_frequency),
    m_npackets(s_default_npackets),
    m_no_pedestals_in_core(s_default_no_pedestals_in_core),
    m_nsb_window_start(s_default_nsb_window_start),
    m_nsb_window_width(s_default_nsb_window_width),
    m_events_per_slice(s_default_events_per_slice),
    m_nsb_suppress_dev(s_default_nsb_suppress_dev),
    m_nsb_min_events_per_slice(s_default_nsb_min_events_per_slice),
    m_no_nsb_suppress(s_default_no_nsb_suppress),
    m_ped_slice_time(s_default_ped_slice_time),
    m_ped_min_window(s_default_ped_min_window),
    m_ped_sample_separation(s_default_ped_sample_separation),
    m_ped_no_suppress(s_default_ped_no_suppress)
{
  // nothing to see here
}

VSAnalysisStage1::~VSAnalysisStage1()
{
  // nothing to see here
}

void VSAnalysisStage1::
runStage1(const std::string& vbf_filename, Data& data,
	  const SphericalCoords& earth_position,
	  bool no_db, bool no_l3, bool no_threads, bool no_verbose,
	  bool do_use_overflow)
{
  if(!no_verbose)vstream << LOGO << std::endl;

  VSAnalysisData analysis_data(REVISION,VERSION);

  // --------------------------------------------------------------------------
  // Run the STAGE 1 analysis
  // --------------------------------------------------------------------------

  data.clear();

  // Print frequency writes the event count at a given frequency --------------
  VBFPrintFrequency visitor_printfrequency(std::cout,m_print_frequency);

  // Run info extractor collects lots of useful information -------------------
  VBFRunInfo visitor_runinfo(!m_no_pedestals_in_core,
			     0);

  // Time dependent NSB calculator, used to suppress channels -----------------
  VBFTimeDepNSB visitor_nsb(m_nsb_window_start, m_nsb_window_width,
			    m_events_per_slice, m_nsb_min_events_per_slice);

  // Extract L3 pointing information ------------------------------------------
  VSL3PointingDataSource visitor_l3pointing;

  // Broadcast visitor distributes events to different visitors ---------------
  VBFBroadcastVisitor visitor_master;
  if((!no_verbose)&&(m_print_frequency))
    visitor_master.addVisitor(&visitor_printfrequency);
  visitor_master.addVisitor(&visitor_runinfo);
  visitor_master.addVisitor(&visitor_nsb);
  if(!no_l3)visitor_master.addVisitor(&visitor_l3pointing);

  // Dispatch packets to all visitors -----------------------------------------
  VSSimpleVBFDispatcher dispatcher(&visitor_master);
  if(no_l3)dispatcher.setDontUseL3Time(true);
  VSSimpleVBFDispatcher::catchSignalAndStopProcessingFile(SIGINT);
  dispatcher.openFile(vbf_filename.c_str());
  uint32_t flags = VSSimpleVBFDispatcher::DISPATCH_NONE;
  //flags |= VSSimpleVBFDispatcher::DISPATCH_DISCARD_INVALID;
  //flags |= VSSimpleVBFDispatcher::DISPATCH_IN_ORDER;
  //flags |= VSSimpleVBFDispatcher::DISPATCH_VERBOSE;
#ifndef NOTHREADS
  if(!no_threads)flags|=VSSimpleVBFDispatcher::DISPATCH_THREADED;
#endif
  if(do_use_overflow)
    flags |= VSSimpleVBFDispatcher::DISPATCH_INTEGRATE_OVERFLOW;
  unsigned nevent = 0;
  nevent = dispatcher.dispatchAllPackets(m_npackets,flags);
  dispatcher.resetVisitor(0);

  // 1. Run information -------------------------------------------------------

  visitor_runinfo.getData(data.run_info);
  data.sim_info = visitor_runinfo.getSimData();

  if(data.run_info.run_number == 0)
    data.run_info.run_number = 
      VSFileUtility::extractNumberFromFilename(vbf_filename);

  // 2. NSB summary information -----------------------------------------------

  visitor_nsb.getData(data.nsb);
  for(unsigned iscope=0;iscope<data.nsb.scope.size();iscope++)
    {
      if((m_nsb_suppress_dev.size()<=iscope)||(m_nsb_suppress_dev[iscope]==0))
	data.nsb.scope[iscope].suppress_dev 
	  = data.nsb.scope[iscope].suppress_dev_suggested;
      else 
	data.nsb.scope[iscope].suppress_dev 
	  = m_nsb_suppress_dev[iscope];
      
      std::cout << "Suppress dev " << iscope << ": " 
		<< data.nsb.scope[iscope].suppress_dev << std::endl;
    }

#if 0
  for(VSSimpleHist<double>::iterator ibin = data.nsb.all_dev.begin();
      ibin != data.nsb.all_dev.end(); ibin++)
    std::cout << ibin->val() << ' ' << ibin->count() << std::endl;
#endif

  // 3. Suppress channels -----------------------------------------------------

  unsigned nslice = visitor_nsb.nslice();
  data.suppress.resize(m_events_per_slice, nslice, data.run_info.nchan);
  data.suppress_event_slice =
    visitor_runinfo.getEventSlices(m_events_per_slice);

#if 0
  std::cerr << nslice << std::endl;
  for(unsigned islice=0;islice<data.suppress_event_slice.size();islice++)
    std::cerr << data.suppress_event_slice[islice].event_num_lo << ' '
	      << data.suppress_event_slice[islice].event_num_hi << ' '
	      << data.suppress_event_slice[islice].event_time_lo << ' '
	      << data.suppress_event_slice[islice].event_time_hi << std::endl;
#endif

  if(!m_no_nsb_suppress)
    {
      const unsigned nscope = data.nsb.scope.size();
      std::vector<double> suppress_dev(nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	suppress_dev[iscope] = data.nsb.scope[iscope].suppress_dev;
      visitor_nsb.suppressFromPedDev(data.suppress,suppress_dev);
    }

  if(!no_db)
    data.suppress.suppressFromIMon(data.suppress_event_slice);

#if 0
  for(unsigned iscope=0;iscope<data.suppress.nscope();iscope++)
    {
      unsigned nchan = data.suppress.nchan(iscope);
      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  unsigned count = data.suppress.getSuppressedCount(iscope,ichan);
	  if(count)std::cout << iscope << ' ' << ichan << ' ' << count
			     << std::endl;
	}
    }
#endif

  // 4. Pointing information from L3 ------------------------------------------

  if(!no_l3)
    {
      data.l3_pointing = new VSL3PointingData;
      visitor_l3pointing.getData(*data.l3_pointing);

      VSL3Pointing pointing(*data.l3_pointing);
      data.run_info.calculateMeanPointing(pointing, earth_position); 
    }
  
  // 5. Pointing information from database ------------------------------------

  if(!no_db)
    {
      VSTime lo_time = data.run_info.lo_event_time;
      VSTime hi_time = data.run_info.hi_event_time;
      lo_time -= INT64_C(60000000000);
      hi_time += INT64_C(60000000000);
      data.db_pointing = new VSDBPointingData;
      for(unsigned iscope=0;iscope<data.run_info.nchan.size(); iscope++)
	if(data.run_info.nchan[iscope])
	  VSDBPointingDataSource::loadPointings(iscope, lo_time, hi_time, 
						*data.db_pointing);

      VSDBPointing pointing(*data.db_pointing);
      data.run_info.calculateMeanPointing(pointing, earth_position); 
    }

  // 6. Dispatch pedestal events to time dependent pedestal calculator --------

  double slice_time = 
    visitor_runinfo.getSliceWidthForApproxWidth(m_ped_slice_time);
  data.pedestal_time_slice = visitor_runinfo.getTimeSlices(slice_time);

#if 0
  std::cerr << slice_time << std::endl;
  std::cerr << data.pedestal_time_slice.size() << std::endl;
  for(unsigned islice=0;islice<data.pedestal_time_slice.size();islice++)
    std::cerr << data.pedestal_time_slice[islice].event_num_lo << ' '
	      << data.pedestal_time_slice[islice].event_num_hi << ' '
	      << data.pedestal_time_slice[islice].event_time_lo << ' '
	      << data.pedestal_time_slice[islice].event_time_hi << std::endl;
#endif

  const VSTimeDepSupData* suppress = 0;
  if(!m_ped_no_suppress)suppress = &data.suppress;
  VBFTimeDepPeds visitor_peds(data.pedestal_time_slice, suppress, 
			      m_ped_min_window, m_ped_sample_separation);

  dispatcher.resetVisitor(&visitor_peds);

  //VBFDumper dumper(std::cout, true);
  //dispatcher.resetVisitor(&dumper);

  visitor_runinfo.dispatchPedestals(dispatcher);
  dispatcher.resetVisitor(0);

  visitor_peds.getData(data.pedestals);

  dispatcher.closeFile();

  // 7. Get the miscellaneous DB data -----------------------------------------

  if(!no_db)
    {
      std::vector<unsigned> telescopes;
      for(unsigned iscope=0;iscope<data.run_info.nchan.size(); iscope++)
	if(data.run_info.nchan[iscope])telescopes.push_back(iscope);
      data.misc_db = new VSMiscellaneousDBData;
      data.misc_db->fillFromDB(data.run_info.lo_event_time, 
			       data.run_info.hi_event_time, telescopes);

      VSTargetTable tt(data.misc_db->target_table);
      data.target_table = new VSTargetTable::Data(tt.data());
    }

  // --------------------------------------------------------------------------
  // Stop the time counter and print out the execution time
  // --------------------------------------------------------------------------

  analysis_data.stop();
  data.analyze = analysis_data;
  analysis_data.printTiming(vstream,vbf_filename,nevent);
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSAnalysisStage1::
configure(VSOptions& options, const std::string& profile,
	  const std::string& opt_prefix)
{
  options.addCatagory("s1_misc", "Stage 1 miscellaneous options.");

  options.addCatagory("s1_ped", 
		      "Stage 1 options conrolling pedestal and high "
		      "resolution NSB calculation, timeslice duration, etc.");

  options.findWithValue(OPTNAME(opt_prefix,"print_frequency"),
			s_default_print_frequency,
			"Set the frequency of printing event numbers.",
			"s1_misc");

  options.findWithValue(OPTNAME(opt_prefix,"npackets"),
			s_default_npackets,
			"Set the maximum number of packets to process.",
			"s1_misc");

  options.findBoolValue(OPTNAME(opt_prefix,"no_pedestals_in_core"),
			s_default_no_pedestals_in_core, true,
			"Do not keep VBF packets of pedestals in memory, "
			"instead make a second pass through the file to "
			"process pedestal events.",
			"s1_ped");
  
  options.findWithValue(OPTNAME(opt_prefix,"events_per_slice"),
			s_default_events_per_slice,
			"Set the number of events per time slice used for "
			"suppressing channels and calculating the high "
			"resolution NSB level.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"nsb_window_start"),
			s_default_nsb_window_start,
			"Set the starting sample number for integrating "
			"to find the NSB level. The NSB level is usually "
			"calculated from the first few samples of every "
			"event.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"nsb_window_width"),
			s_default_nsb_window_width,
			"Set the sampling window width for integrating "
			"to find the NSB level. The NSB level is usually "
			"calculated from the first few samples of every "
			"event.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"nsb_suppress_dev"),
			s_default_nsb_suppress_dev,
			"Set the RMS NSB level below which a channel is "
			"considered to be switched off. The channel will "
			"be set as suppressed for the time slice. A value "
			"of zero requests that the suppression threshold be "
			"chosen automatically based on the RMS NSB "
			"distribution.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"nsb_min_events_per_slice"),
			s_default_nsb_min_events_per_slice,
			"Set the minum number of events per slice required "
			"for the NSB estimate to be valid.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"no_nsb_suppress"),
			s_default_no_nsb_suppress,
			"Do not use the high resolution NSB estimate to "
			"suppress channels which appear to be switched off.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"ped_slice_time"),
			s_default_ped_slice_time,
			"Set the width of the pedestal time slices in "
			"seconds.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"ped_min_window"),
			s_default_ped_min_window,
			"Set the minimum window size over which to calculate "
			"pedestal variances.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"ped_sample_separation"),
			s_default_ped_sample_separation,
			"Set the number of samples to leave between seperate "
			"evaluations of the pedestal from a single trace. "
			"This should be long enough to break the correlation "
			"that would otherwise exist if multiple evaluations "
			"were made from on trace.",
			"s1_ped");

  options.findWithValue(OPTNAME(opt_prefix,"ped_no_suppress"),
			s_default_ped_no_suppress,
			"Do not suppress accumulation of pedestal information "
			"for a channel during slices in which they are "
			"suppressed.",
			"s1_ped");
}

VSAnalysisStage1::Data::~Data()
{
  delete l3_pointing;
  delete db_pointing;
  delete laser;
  delete target_table;
  delete misc_db;
  delete sim_info;
  delete hilo;      
}

void VSAnalysisStage1::Data::clear()
{
  analyze.clear();
  run_info.clear();
  nsb.clear();
  suppress.clear();
  pedestals.clear();
  suppress_event_slice.clear();
  pedestal_time_slice.clear();
  delete l3_pointing;  l3_pointing  = 0;
  delete db_pointing;  db_pointing  = 0;
  delete laser;        laser        = 0;
  delete target_table; target_table = 0;
  delete misc_db;      misc_db      = 0;
  delete sim_info;     sim_info     = 0;
  delete hilo;         hilo         = 0;
}

void VSAnalysisStage1::Data::load(VSOctaveH5ReaderStruct* reader)
{
  analyze.load(reader->readStruct("analyze"));
  run_info.load(reader->readStruct("run_info"));
  nsb.load(reader->readStruct("nsb"));
  suppress.load(reader->readStruct("suppress"));
  pedestals.load(reader->readStruct("pedestals"));
  VBFRunInfo::Slice::load(reader->readStruct("suppress_event_slice"),
			  suppress_event_slice);
  VBFRunInfo::Slice::load(reader->readStruct("pedestal_time_slice"),
			  pedestal_time_slice);

  delete(l3_pointing); l3_pointing=0;
  if(reader->isStruct("l3_pointing"))
    {
      l3_pointing=new VSL3PointingData;
      l3_pointing->load(reader->readStruct("l3_pointing"));
    }

  delete(db_pointing); db_pointing=0;
  if(reader->isStruct("db_pointing"))
    {
      db_pointing=new VSDBPointingData;
      db_pointing->load(reader->readStruct("db_pointing"));
    }

  delete(laser); laser=0;
  if(reader->isStruct("laser"))
    {
      laser=new VSLaserData;
      laser->load(reader->readStruct("laser"));
    }

  delete(target_table); target_table=0;
  if(reader->isStruct("target_table"))
    {
      target_table=new VSTargetTableData;
      target_table->load(reader->readStruct("target_table"));
    }

  delete(misc_db); misc_db=0;
  if(reader->isStruct("misc_db"))
    {
      misc_db=new VSMiscellaneousDBData;
      misc_db->load(reader->readStruct("misc_db"));
    }

  delete(sim_info); sim_info=0;
  if(reader->isStruct("sim_info"))
    {
      sim_info=new VSSimInfoData;
      sim_info->load(reader->readStruct("sim_info"));
    }

  delete(hilo); hilo=0;
  if(reader->isStruct("hilo"))
    {
      hilo=new VSHiLoData;
      hilo->load(reader->readStruct("hilo"));
    }
}

void VSAnalysisStage1::Data::save(VSOctaveH5WriterStruct* writer) const
{
  analyze.save(writer->writeStruct("analyze"));
  run_info.save(writer->writeStruct("run_info"));
  nsb.save(writer->writeStruct("nsb"));
  suppress.save(writer->writeStruct("suppress"));
  pedestals.save(writer->writeStruct("pedestals"));
  VBFRunInfo::Slice::save(writer->writeStruct("suppress_event_slice"),
			  suppress_event_slice);
  VBFRunInfo::Slice::save(writer->writeStruct("pedestal_time_slice"),
			  pedestal_time_slice);
  if(l3_pointing)l3_pointing->save(writer->writeStruct("l3_pointing"));
  if(db_pointing)db_pointing->save(writer->writeStruct("db_pointing"));
  if(laser)laser->save(writer->writeStruct("laser"));
  if(target_table)target_table->save(writer->writeStruct("target_table"));
  if(misc_db)misc_db->save(writer->writeStruct("misc_db"));
  if(sim_info)sim_info->save(writer->writeStruct("sim_info"));
  if(hilo)hilo->save(writer->writeStruct("hilo"));
}
