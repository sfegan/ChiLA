//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFAnalysisStage2.cpp

  Stage 2 analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 3.66 $
  \date       05/10/2006

  $Id: VBFAnalysisStage2.cpp,v 3.66 2010/03/23 00:49:12 matthew Exp $

*/

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// TODO ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// *) Add gain variance to spot bi-modal gain channels
// *) Proper logging class whose o/p goes to cout and to file
// *) Plot in diagnostics detailing motion of sources in FoV
// *) Add accumulation of ped histogram for 1 sample

// *) Add default parameter profiles
// *) Count events with lost scope data to deadtime (L3 threshold setting)
// *) Add size spectrum calculation directly to analyze
// *) Need to use tracking correction - nsbmovie or camera measurememnts?
// *) Fix help text for sample_start and sample_width
// *) Implement listTargets
// *) L3 bits histogram / Scope triggered / Scope triggered co-variance
//    Scope trigger types (full,extra)
// *) Fix method 3
// *) Pointing from RADEC target for when the DB is not near by
// *) No L3 flag to suppress requirement that L3 events be present

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// DONE ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// *) Test HIT bit in LeaveEvent
// *) Fix orientation problem when in event frame
// *) Always calculate theta with respect to mean tracking position
// *) Wobble, target location (for theta -- or list of them?)
// *) Simple theta cut to make file size managable
// *) X and Y in FoV (relative to mean array location perhaps?)
// *) Command line recording
// *) Diagonstics, times, livetime, RA, Dec, mean array position
// *) Livetimes, GPS delta t
// *) Only save scope parameters for scopes that are in the run
// *) Test integration window -- variable window - caution: high gain remnant
// *) Finish "tij" diagnostics
// *) Honor suppression of telescope in visitArrayTriggerScope
// *) Option to suppress writing of diagnostics
// *) Venn diagram of how many channels were suppressed by each cut
// *) Threaded dispatcher -- VASequentialPacketRetriever
// *) Impact parameter on telescope by telescope basis
// *) Photon arrival direction WRT mean array position
// *) Cleaning factory
// *) Test regional cleaner
// *) Have VSOptions remember all options as strings - write them to H5
// *) New telescope positions
// *) Write targets list as RA, DEC (number and stringified) and target name
// *) Pointing from L3
// *) Move events variables to "events" tree
// *) Add L1 calibration as matter of routine - gain 7.81 - looks good!!
// *) Packet linearizer - lookahead in VBF file to find next event
// *) Added elapsed_ticks and live_ticks to the events tree
// *) Pass more variables in ArrayData, removing potential threading error
// *) Fixed a shameful mistake in calculation of signal to noise which was 
//    reason why optimized cleaning thresholds were so low
// *) New "cleaning_mult" in merged info so cleaning uses multiplication -
//    this allows strict threshold in PEs if needed
// *) Added more variables to merged info for future time variable peds etc
// *) visitHitChannel now uses non-VSTimingCalc non-virtual API call 
// *) New "analyze" tree in H5 file with version, run execution information
// *) Added user data pass-through to VSSimpleVBF
// *) visitHitChannel now takes c-style arrays for speed
// *) Cleaner and pointing moved to seperate files
// *) Pointing data now logically separated from code
// *) MAJOR UPDATE -- CODE MOVED TO SEPERATE FILES AND THREADING MODEL CHANGED
// *) Pedestal variance should be a function of time
// *) Channel suppression from database and from high resolution NSB estimate
// *) Fix vector normalization problem when telescope coordinate is zero
// *) Move most channel diagnostics to visitHitChannel
// *) Have new PartialDiagnostics class with free-list to store diagnostics
//    data. Instantiate one per thread and store reference in IED
// *) Move all access to channel data to visitHitChannel and leaveScopeEvent.
//    Now stage2 analysis is faster than analyze - why are diagnostics so slow
// *) Histogram of nimage-ntrig_in_largest_region
// *) Move to version number 2.0
// *) Separate automatic tracking determination from production of sources
// *) Fix to threading code to handle writing of events and fetching of IED
//    buffers properly.
// *) Test database before processing file in Stage 1
// *) Add auto targeting based on mean position, commanded DB target etc
// *) Read commanded telescope target from DB in stage 1
// *) Switched over some histograms to VSLimitedHist for safety
// *) Better clock estimation
// *) Download currents, L1 rates, FIR
// *) Plots of median pedvar, L1 rates, currents
// *) Muon ring detection
// *) Auto determine suppression levels based on 95% x 1.2 formula
// *) Weighting of telescopes should be independent of reconstruction method
// *) Diagnostics and event data classes - reloadable
// *) Add Dave Hanna's gain calculation to LASER
// *) Implement supression of telescope
// *) Implement software L3 trigger
// *) Calculate mean scaled parameters
// *) Split functionality over multiple CPP files
// *) Rename class to VBFAnalysisStage2

/* 
   The threading model in this class has become very complicated. 
   Multiple threads run through most of the class. Event processing is
   linearized at the writing stage, and partially linearized at the
   beginning when fetching IED buffers from the free list. The
   seq_packet_number plays a large role at both points. 

   In the writer the seq_number is used to determine what event is
   next to write. In my first attempt at the writer the IEDs were
   registered on the "processing" list and written when the first one
   or more events on the list were complete. However the next event in
   sequence has not necessarily been placed on the processing list (if
   its thread was swapped out), so this method fails badly, especially
   under load. The correct way to do it to additionally ensure that
   the first IED on the processing list is the next one in sequence,
   using the seq_number.

   At the upper end of the code, IED buffers are passed out to the
   threads as they become available on the free list. Initially, it
   was "first come first served", and with the mistake in the writer
   there was no problem evident here. However with the new writer this
   could mean that the next-in-sequence event is delayed arbitrarily
   if they do not wake up from their sleep quickly enough to grab the
   buffer. Under load the system would lock up, with all IED buffers
   used by subsequent events. Now we keep track of all seq_numbers
   using the largest requested one to date. When free IED buffers are
   plentiful they are handed out on a first come first served basis,
   but as they become scarce, they are reserved for the
   next-in-sequence events, even though they may not yet have even
   been requested. This is done with a priority queue. As an IED is
   requested, the seq_number of all events between the highest one
   seen so far and the currently requested one are entered on the
   priority queue. When there are less free buffers than there are
   threads a buffer is only returned to threads handling events close
   to the start of the queue. So say there are 10 seq_numbers on the
   queue and 4 IED buffers remaining on the free list, then only
   threads with seq_numbers corresponding to the first 4 on the queue
   will get an IED buffer. All others will just wait on the condition
   variable. Entries are deleted from the queue when they receive a
   buffer.

   The model does not necessarily have to be so complicated. For
   example we could continually allocate and assign IEDs if the free
   list is empty. That might be better. It would certainly be easier!
   Or a hybrid approach could be used, where they are dynamically
   allocated up to a certain maximum number, after which the wait
   queue is used. That would be a little bit more complex :-)
*/

#include <sstream>
#include <iomanip>
#include <cmath>

#include <Astro.h>
#include <VSFileUtility.hpp>
#include <fast_alloc.hpp>
#include <WhippleCams.h>
#include <VSAAlgebra.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSALinearLeastSquares.hpp>

#include <VBFAnalysisStage2.hpp>

using namespace VERITAS;
using namespace SEphem;

template<typename T>
static inline T norm(const T& a, const T&b) { return std::sqrt(a*a+b*b); }

static double NEGINF = -std::numeric_limits<double>::infinity();
static double POSINF = std::numeric_limits<double>::infinity();

// ----------------------------------------------------------------------------
// __   _____ ___ _             _            _    ___ _                 ___
// \ \ / / _ ) __/_\  _ _  __ _| |____  _ __(_)__/ __| |_ __ _ __ _ ___|_  )
//  \ V /| _ \ _/ _ \| ' \/ _` | (_-< || (_-< (_-<__ \  _/ _` / _` / -_)/ /
//   \_/ |___/_/_/ \_\_||_\__,_|_/__/\_, /__/_/__/___/\__\__,_\__, \___/___|
//                                   |__/                     |___/
// ----------------------------------------------------------------------------

VBFAnalysisStage2::
VBFAnalysisStage2(VSCleaner* clean,
		  const SecondaryCleaning& secondary_clean,
		  VSPointing* pointing,
		  VSAReconstruction* reconstruction, 
		  std::vector<VSMuonAnalysis*> muon_analysis,
		  VSScaledParameterCalc* sp_calc,
		  VSEnergyCalc* energy_calc,
		  VSCutsCalc* cuts,
		  const VSAnalysisStage1Data& stage1,
		  const VSChannelMap& channel_map,
		  const Settings& settings,
		  VSOctaveH5Writer* io,
		  const VSAnalysisStage1Data* pad_stage1):
  m_clean(clean),
  m_secondary_clean(secondary_clean),
  m_pointing(pointing),
  m_reconstruction(reconstruction),
  m_muon_analysis(muon_analysis),
  m_sp_calc(sp_calc),
  m_energy_calc(energy_calc),
  m_cuts_calc(cuts),
  m_io(io),
  m_stage1(&stage1),
  m_pad_stage1(pad_stage1),
  m_channel_map(&channel_map),
  m_settings(settings),
  m_slow_diagnostics(),
  m_qcalc(),
  m_diagnostics(),
  m_edw(),
  m_acw(),
  m_asw(),
  m_s2ew(),
  m_e2sw(),
  m_cal(),
  m_muon_data(),
  m_sim_transform(),
  m_j2000_phi(),
  m_j2000_theta(),
  m_j2000_psi(),
  m_sim_ntable(),
  m_sim_ebinner(),
  m_sim_package(VBFRunInfo::SimData::SP_UNKNOWN),
  m_ied_maximum(),
  m_ied_created(),
  m_ied_free(),
  m_ied_processing(),
  m_ied_writer_active(),
  m_ied_writer_next_seq_packet_number(),
  m_tsd_list(),
#ifndef NOTHREADS
  m_threaded(),
  m_ied_list_mutex(),
  m_ied_list_wait_cond(),
  m_pointing_mutex(),
  m_tsd_key(),
  m_nthreads(),
  m_ied_wait_priority_next(),
  m_ied_wait_priority_list(),
  m_ied_free_stat(),
  m_ied_wait_count(),
#endif
  m_last_event_time(),
  m_last_event_time_sec(),
  m_last_event_time_hist(),
  m_last_event_time_ten_mhz_elap(),
  m_last_event_time_ten_mhz_live(),
  m_last_event_num(),
  m_last_l3_event_num(),
  m_last_ten_mhz_elapsed(),
  m_last_ten_mhz_veto_both(),
  m_last_ten_mhz_veto_vdaq(),
  m_last_ten_mhz_veto_l3(),
  m_last_counts_l2(),
  m_mean_have_ra_zero(),
  m_mean_ra_zero(),
  m_mean_ra(),
  m_mean_dec(),
  m_have_one_array_event(),
  m_have_one_array_event_with_l3(),
  m_first_time(),
  m_first_scope_has_event(),
  m_sim_header(),
  m_nevent_written(),
  m_nsim_written(),
  m_sim_ntel()
{
  const Settings& S(settings);

  const unsigned nscope = stage1.run_info.nchan.size();

  m_last_counts_l2.resize(nscope);
  m_first_scope_has_event.resize(nscope);

  vsassert(S.integration_lo_zero_sample.size() >= nscope);
  vsassert(S.integration_hi_zero_sample.size() >= nscope);

  m_qcalc = new VSConstantFractionCalc(S.integration_threshold_frac,
				       S.integration_window_start,
				       S.integration_window_width,
				       0, 0, 
				       S.integration_threshold_charge,
				       S.integration_window_start_increase,
				       S.integration_window_width_increase);

  m_cal = new ArrayMergedCal(stage1.run_info, 
			     &stage1.suppress,
			     &stage1.pedestals,
			     S.integration_hi_zero_sample,
			     S.integration_lo_zero_sample,
			     S.cleaning_scale,
			     S.integration_window_width,
			     m_channel_map,
			     reconstruction->array(),
			     stage1.laser,
			     S.permissive_laser,
			     S.unity_gain,
			     stage1.hilo,
			     S.scope_hi_lo_gain_ratio,
			     S.ped_event_count_min,
			     S.scope_suppress,
			     S.scope_gain,
			     S.auto_suppress_l2_chan,
			     S.chan_suppress,
			     S.chan_has_no_pmt,
			     pad_stage1?&pad_stage1->suppress:0,
			     pad_stage1?&pad_stage1->pedestals:0,
			     stage1.misc_db,
			     S.scope_dev_scaling);

  if(!m_muon_analysis.empty())m_muon_data = new VSMuonAnalysisData(nscope);

  for(unsigned iscope=0;iscope<m_cal->scope.size();iscope++)
    {
      unsigned nchan = m_cal->scope[iscope].nchan;
      for(unsigned ichan=0;ichan<nchan;ichan++)
	if(m_cal->scope[iscope].channel[ichan].l2chan == nchan)
	  m_settings.enable_l2_corrections = false;
    }

  if(!S.do_not_write_diagnostics)
    m_diagnostics = 
      new ArrayDiagnostics(stage1.run_info.config_mask,
			   stage1.run_info.nchan, 
			   stage1.run_info.nsample,
			   m_settings.limited_dt_range);

  m_slow_diagnostics = 
    (!S.do_not_write_diagnostics)&&(!S.no_slow_diagnostics);

  m_edw = new VSEventDataWriter(io->writeStruct("events"), 
				stage1.run_info.nchan,
				S.theta_targets.size()+1);

  if(m_cuts_calc)
    m_acw = 
      new VSArrayCutsWriter(io->writeStruct("cuts"), stage1.run_info.nchan);

  Astro::precessionalAngles(Astro::julianEpochToMJD(2000.0), 
			    stage1.run_info.lo_event_time.getMJDDbl(),
			    m_j2000_phi, m_j2000_theta, m_j2000_psi);

  if(stage1.sim_info)
    {
      m_sim_transform = new VSSimCoordTransform;
      SphericalCoords obs_radec = 
	SEphem::SphericalCoords::makeLatLongDeg(stage1.sim_info->obs_dec_deg, 
						stage1.sim_info->obs_ra_deg);
      m_sim_transform->setObsRADec(obs_radec);

      m_sim_package = stage1.sim_info->package;
      if(m_sim_package != VBFRunInfo::SimData::SP_CHILA)
	{
	  for(std::vector<VBFRunInfo::SimData::ParticleSpectrum>::iterator
		ispect = stage1.sim_info->particle_spectrum.begin();
	      ispect != stage1.sim_info->particle_spectrum.end(); ispect++)
	    if(ispect->dlog_energy>0)
	      {
		double e0 = 
		  ispect->discrete_binning_start-ispect->dlog_energy*0.5;
		double en = log10(ispect->max_energy_tev);
		EBinner eb(e0, en, ispect->dlog_energy, m_sim_ntable);
		m_sim_ebinner[ispect->corsika_particle_id] = eb;
		m_sim_ntable += eb.ntable();
	      }
	    else
	      {
		double e0 = log10(ispect->min_energy_tev);
		double en = log10(ispect->max_energy_tev);
		EBinner eb(e0, en, 0.0625, m_sim_ntable);
		m_sim_ebinner[ispect->corsika_particle_id] = eb;
		m_sim_ntable += eb.ntable();
	      }
	}
    }

  m_ied_maximum = 3;

  // --------------------------------------------------------------------------
  // Load the mean scaled parameter tables
  // --------------------------------------------------------------------------
  sp_calc->load(stage1.run_info.nchan,
		stage1.run_info.mean_zn_1dim_deg,
		stage1.run_info.az_mean_deg,
		m_cal->mean_scaled_dev);
  
  // --------------------------------------------------------------------------
  // Load the energy lookup tables
  // --------------------------------------------------------------------------
  energy_calc->load(stage1.run_info.nchan,
		    stage1.run_info.mean_zn_1dim_deg,
		    stage1.run_info.az_mean_deg,
		    m_cal->mean_scaled_dev);

#ifndef NOTHREADS
  vsassert(pthread_mutex_init(&m_ied_list_mutex, NULL) == 0);
  vsassert(pthread_cond_init(&m_ied_list_wait_cond, NULL) == 0);
  vsassert(pthread_mutex_init(&m_pointing_mutex, NULL) == 0);
  vsassert(pthread_key_create(&m_tsd_key, NULL) == 0);
#endif
}

VBFAnalysisStage2::~VBFAnalysisStage2()
{
  // --------------------------------------------------------------------------
  // Finish writing events
  // --------------------------------------------------------------------------

  processConstructedIEDs();

  // --------------------------------------------------------------------------
  // Write simulation header
  // --------------------------------------------------------------------------

  if(m_sim_header)
    {
      if(m_sim_package == VBFRunInfo::SimData::SP_KASCADE)
	{
	  for(std::vector<VSTableSimulationDatum>::iterator itable =
		m_sim_header->tables.begin(); 
	      itable != m_sim_header->tables.end();
	      itable++)
	    {
	      itable->sampling_radius_m = sqrt(itable->sampling_radius_m/M_PI);
	      itable->event_count = itable->num_events_read;
	    }
	}

      // Write header
      VSOctaveH5WriterStruct* s = m_io->writeStruct("sim_header");
      m_sim_header->save(s);
      delete s;
      delete m_sim_header;
    }

  // --------------------------------------------------------------------------
  // Calculate mean RA and DEC across all telescopes
  // --------------------------------------------------------------------------

  double ra = Angle::toDeg(m_mean_ra.mean() + m_mean_ra_zero);
  double dec = Angle::toDeg(m_mean_dec.mean());
  if(ra<0)ra+=360;
  SphericalCoords radec;
  radec.setLatLongDeg(dec,ra);

  SphericalCoords radec_j2000(radec);
  radec_j2000.rotate(m_j2000_phi, m_j2000_theta, m_j2000_psi);
  
  // --------------------------------------------------------------------------
  // Write targets
  // --------------------------------------------------------------------------

  VSOctaveH5WriterCellVector* tcell = 
    m_io->writeCellVector("targets",m_settings.theta_targets.size()+1);
    
  VSOctaveH5WriterStruct* tstruct = tcell->writeStruct(0);
  tstruct->writeString("name","Mean pointing");
  tstruct->writeScalar("ra",ra);
  tstruct->writeScalar("dec",dec);
  tstruct->writeScalar("ra_J2000",radec_j2000.longitudeDeg());
  tstruct->writeScalar("dec_J2000",radec_j2000.latitudeDeg());
  tstruct->writeString("ra_string",radec.longitude().hmsString(1));
  tstruct->writeString("dec_string",radec.latitude().dmsString(1));
  tstruct->writeString("ra_J2000_string",radec_j2000.longitude().hmsString(1));
  tstruct->writeString("dec_J2000_string",radec_j2000.latitude().dmsString(1));
  tstruct->writeScalar("offset",0.0);
  tstruct->writeScalar("direction",0.0);
  delete tstruct;

  for(unsigned itarget=0;itarget<m_settings.theta_targets.size();itarget++)
    {
      SphericalCoords radec_j2000(m_settings.theta_targets[itarget].coord);
      radec_j2000.rotate(m_j2000_phi, m_j2000_theta, m_j2000_psi);

      Angle s;
      Angle d;
      radec.separationAndDirectionTo(m_settings.theta_targets[itarget].coord,
				     s,d);
      d = M_PI/2 - d; // for convention with 0=N, 90=E
      VSOctaveH5WriterStruct* tstruct = tcell->writeStruct(itarget+1);
      tstruct->writeString("name",m_settings.theta_targets[itarget].name);
      tstruct->writeScalar("ra",
		      m_settings.theta_targets[itarget].coord.longitudeDeg());
      tstruct->writeScalar("dec",
		      m_settings.theta_targets[itarget].coord.latitudeDeg());
      tstruct->writeScalar("ra_J2000",radec_j2000.longitudeDeg());
      tstruct->writeScalar("dec_J2000",radec_j2000.latitudeDeg());
      tstruct->writeString("ra_string",
		      m_settings.theta_targets[itarget].coord.longitude().
			   hmsString(1));
      tstruct->writeString("dec_string",
		      m_settings.theta_targets[itarget].coord.latitude().
			   dmsString(1));
      tstruct->writeString("ra_J2000_string",
			   radec_j2000.longitude().hmsString(1));
      tstruct->writeString("dec_J2000_string",
			   radec_j2000.latitude().dmsString(1));
      tstruct->writeScalar("offset",s.deg());
      tstruct->writeScalar("direction",d.deg());
      delete tstruct;
    }

  delete tcell;

  // --------------------------------------------------------------------------
  // Write muons, diagnostics and channel info
  // --------------------------------------------------------------------------

  if(m_diagnostics)
    {
      std::list<PartialArrayDiagnostics*> diags;
      for(std::list<ThreadSpecificData*>::iterator itsd = m_tsd_list.begin();
	  itsd != m_tsd_list.end(); itsd++)
	diags.push_back((*itsd)->partial_diagnostics);
      m_diagnostics->integratePartial(diags);

      m_diagnostics->finalize(m_cal, m_stage1->run_info, m_stage1->misc_db,
			      &m_stage1->nsb, m_settings.earth_position);

      VSOctaveH5WriterStruct* s = m_io->writeStruct("diagnostics");
      m_diagnostics->save(s, m_slow_diagnostics);
      delete s;
    }

  if(m_muon_data)
    {
      VSOctaveH5WriterStruct* s = m_io->writeStruct("muon_analysis");
      m_muon_data->save(s);
      delete s;
    }

  // MergedCal's save() function is different than all others, we do not need
  // to create a separate structure for it to use. It creates a CellVector of
  // the correct size

  m_cal->save(m_io);  
  
  m_io->writeScalar("run_number",m_stage1->run_info.run_number);
  m_io->writeScalar("run_start_time",
		    m_stage1->run_info.first_event_time.getMJDDbl());
  m_io->writeString("run_start_time_string",
		    m_stage1->run_info.first_event_time.getString());

  // --------------------------------------------------------------------------
  // Write scaled parameter and energy lookup tables
  // --------------------------------------------------------------------------
  if(m_sp_calc)
    {
      VSOctaveH5WriterStruct* s = m_io->writeStruct("sp_calc");
      m_sp_calc->save(s);
      delete s;
    }

  if(m_energy_calc)
    {
      VSOctaveH5WriterStruct* s = m_io->writeStruct("egy_calc");
      m_energy_calc->save(s);
      delete s;
    }

  // --------------------------------------------------------------------------
  // Clean up
  // --------------------------------------------------------------------------

#ifndef NOTHREADS
  if((m_settings.verbose)&&(m_nthreads))
    vstream << "IED buffer statistics:" << std::endl
	    << "  Created:                      " 
	    << m_ied_created << '/' << m_ied_maximum << std::endl
	    << "  Average free list occupancy:  " 
	    << (double(m_ied_free_stat.sum())/double(m_ied_free_stat.count()))
	    << std::endl
	    << "  Number of thread waits:       "
	    << m_ied_wait_count << std::endl << std::endl;

  vsassert(pthread_key_delete(m_tsd_key) == 0);
  vsassert(pthread_cond_destroy(&m_ied_list_wait_cond) == 0);
  vsassert(pthread_mutex_destroy(&m_ied_list_mutex) == 0);
  vsassert(pthread_mutex_destroy(&m_pointing_mutex) == 0);
#endif

  for(std::list<ThreadSpecificData*>::iterator itsd = m_tsd_list.begin();
      itsd != m_tsd_list.end(); itsd++)
    delete (*itsd);
  
  vsassert(m_ied_processing.empty());
  vsassert(m_ied_free.size() == m_ied_created);

  for(std::list<IED*>::iterator iied = m_ied_processing.begin();
      iied != m_ied_processing.end(); iied++)delete *iied;
  m_ied_processing.clear();
  for(std::list<IED*>::iterator iied = m_ied_free.begin();
      iied != m_ied_free.end(); iied++)delete *iied;
  m_ied_free.clear();
  
  delete m_qcalc;
  delete m_diagnostics;
  delete m_edw;
  delete m_acw;
  delete m_asw;
  delete m_s2ew;
  delete m_e2sw;
  delete m_cal;
  delete m_muon_data;
  delete m_sim_transform;
}

#ifndef NOTHREADS
void VBFAnalysisStage2::
usingThreads(unsigned nthreads)
{
  m_nthreads = nthreads;

  if(nthreads>1)
    {
      m_ied_maximum = nthreads*50;
      m_threaded = true;
    }
  else
    {
      m_ied_maximum = 1;
      m_threaded = false;
    }
}
#endif

void VBFAnalysisStage2::
visitPacket(bool& veto_packet, void*& user_data,
	    uint32_t                   seq_packet_number,
	    uint32_t                   vbf_packet_number,
	    bool                       has_array_event,
	    bool                       has_sim_header,
	    bool                       has_sim_event,
	    uint32_t                   num_overflow_datum,
	    const VBFPacket*           packet)
{
  IED* ied = getEmptyIED(seq_packet_number, vbf_packet_number);
  user_data = static_cast<void*>(ied);
}

void VBFAnalysisStage2::
leavePacket(bool veto_packet, void* user_data)
{
  IED* ied = static_cast<IED*>(user_data);

  unsigned vbf_packet_number = ied->vbf_packet_number;
  if((m_settings.verbose)&&(vbf_packet_number)&&(m_settings.print_frequency)
     &&(vbf_packet_number%m_settings.print_frequency==0))
    {
      std::ostringstream stream;
      stream << std::left
	     << std::setw(5) << m_stage1->run_info.run_number << "  "
	     << std::setw(7) << vbf_packet_number;
      if(ied->has_array_event)stream << ' ' << ied->event_time;
      std::cout << stream.str() << std::endl;
    }
  lockIEDList();
  ied->ied_processing=false;
  unlockIEDList();
  processConstructedIEDs();
}

void VBFAnalysisStage2::
visitArrayEvent(bool& veto_array_event, void* user_data,
		uint32_t               num_scope_events,
		bool                   has_array_trigger,
		bool                   has_event_number,
		uint32_t               event_num,
		bool                   has_good_event_time,
		const VSTime&          best_event_time,
		EventType              event_type,
		uint32_t               l2_trigger_mask,
		const VArrayEvent*     array_event)
{
  IED* ied = static_cast<IED*>(user_data);

  ied->has_array_event            = true;
  ied->event_num                  = event_num;
  ied->event_type                 = event_type;
  ied->event_l2_trigger_mask      = l2_trigger_mask;

  if(!has_good_event_time)
    {
      ied->has_event_time         = false;
    }
  else
    {
      ied->has_event_time         = true;
      ied->event_time             = best_event_time;
      ied->mjd                    = ied->event_time.getMJDDbl();
      ied->lmst                   = 
	Astro::mjdToLMST(ied->mjd,m_settings.earth_position.longitudeRad());
    }

  if((event_type==ET_L2)||(m_settings.process_all_events))
    {
      ied->processed_array_event = true;
      if(!m_settings.software_trigger_masks.empty())
	{
	  ied->event_failed_software_trigger = true;
	  for(std::vector<unsigned>::const_iterator imask =
		m_settings.software_trigger_masks.begin(); 
	      imask != m_settings.software_trigger_masks.end(); imask++)
	    if((l2_trigger_mask & *imask) == *imask)
	      {
		ied->event_failed_software_trigger = false;
		break;
	      }
	  if(ied->event_failed_software_trigger)
	    ied->processed_array_event = false;
	}
    }

  unsigned nscope = ied->scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if((ied->scope[iscope])&&(l2_trigger_mask & 1<<iscope))
      ied->scope[iscope]->l3_l2_received = true;

  nscope = ied->images.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    ied->images[iscope].pixels.clear(), 
      ied->images[iscope].use_in_reconstruction=false;

  if(m_secondary_clean.first)
    for(unsigned iscope=0;iscope<nscope;iscope++)
      ied->secondary_images[iscope].pixels.clear(), 
	ied->secondary_images[iscope].use_in_reconstruction=false;

#if 0
  std::cout << event_type << ' ' << l2_trigger_mask;
  for(unsigned iscope=0;iscope<nscope;iscope++)
    std::cout << ' ' << ied->scope[iscope]->l3_l2_received;
  std::cout << ' ' << ied->processed_array_event << '\n';
#endif

  ied->tsd->cal.setEventNumber(event_num);
}

void VBFAnalysisStage2::
visitArrayTrigger(bool& veto_array_event, void* user_data,
		  uint32_t             event_num,
		  const VEventType&    event_type,
		  uint32_t             trigger_mask,
		  uint32_t             flags,
		  const VSTime&        raw_time,
		  uint32_t             at_flags,
		  uint32_t             config_mask,
		  uint32_t             num_telescopes,
		  uint32_t             num_trigger_telescopes,
		  uint32_t             run_number,
		  const uint32_t*      ten_mhz_clocks,
		  const uint32_t*      cal_count,
		  const uint32_t*      ped_count,
		  const VArrayTrigger* trigger)
{
  IED* ied = static_cast<IED*>(user_data);

  ied->has_l3              = true;
  ied->ten_mhz_elapsed     = ten_mhz_clocks[0];
  ied->ten_mhz_veto_both   = ten_mhz_clocks[1];
  ied->ten_mhz_veto_vdaq   = ten_mhz_clocks[2];
  ied->ten_mhz_veto_l3     = ten_mhz_clocks[3];

  unsigned nscope = ied->scope.size();
  if(nscope>8)nscope=8;

  ied->event_sent_trigger_mask = trigger_mask;
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(trigger_mask & 1<<iscope)
      {
	if(ied->scope[iscope])ied->scope[iscope]->l3_sent=true;
	ied->nscope_sent_l3++;
      }

  if(ied->has_event_time)
    {
      ied->gps_dt = raw_time-ied->event_time;
      //if(ied->gps_dt<INT64_C(-20000))ied->gps_dt=INT64_C(-20000);
      //else if(ied->gps_dt>INT64_C(20000))ied->gps_dt=INT64_C(20000);
    }
}
    
void VBFAnalysisStage2::
visitArrayTriggerScope(bool& veto_array_event, 
		       void* user_data,
		       uint32_t        telescope_num,
		       bool            triggered,
		       uint32_t        event_type,
		       float           altitude,
		       float           azimuth,
		       uint32_t        tdc,
		       uint32_t        shower_delay,
		       uint32_t        comp_delay,
		       const uint32_t* l2_counts,
		       const uint32_t* cal_counts)
{
  IED* ied = static_cast<IED*>(user_data);

  ScopeIED* sied(ied->scope[telescope_num]);
  sied->has_l3              = true;
  if(triggered)ied->nscope_trigger++;
  sied->l3_az               = Angle::frDeg(azimuth);
  sied->l3_el               = Angle::frDeg(altitude);
  sied->l3_counts_l2        = l2_counts[0];
  sied->l3_tdc              = tdc;
}

void VBFAnalysisStage2::
visitScopeEvent(bool& veto_scope_event, void* user_data,
		uint32_t               event_num,
		uint32_t               telescope_num, 
		const VEventType&      event_type,
		uint32_t               trigger_mask,
		uint32_t               flags,
		const VSTime&          raw_time,
		uint32_t               num_samples,
		uint32_t               num_channels_saved,
		uint32_t               num_channels_total,
		uint32_t               num_clock_trigger,
		const VEvent*          event)
{
  IED* ied = static_cast<IED*>(user_data);

  ScopeIED* sied(ied->scope[telescope_num]);

  if((ied->scope_ied==0)&&(!ied->has_l3))
    {
      unsigned nscope = ied->scope.size();
      if(nscope>8)nscope=8;
      ied->event_sent_trigger_mask = trigger_mask;
      for(unsigned iscope=0;iscope<nscope;iscope++)
	if(trigger_mask & 1<<iscope)
	  {
	    if(ied->scope[iscope])ied->scope[iscope]->l3_sent=true;
	    ied->nscope_sent_l3++;
	  }
    }

  ied->scope_num              = telescope_num;
  ied->scope_ied              = sied;
  ied->scope_cal              = &ied->tsd->cal.scope[telescope_num];
  sied->has_vdaq              = true;
  ied->nscope_has_event++;

  ied->event_has_vdaq_mask   |= 1<<telescope_num;

  ied->resetForScopeImage(telescope_num, num_channels_total);

  sied->processed_scope_event = 
    ied->processed_array_event && (!ied->scope_cal->suppress);

  if(ied->has_event_time)
    {
      sied->gps_dt = raw_time-ied->event_time;
      //if(sied->gps_dt<INT64_C(-20000))sied->gps_dt=INT64_C(-20000);
      //else if(sied->gps_dt>INT64_C(20000))sied->gps_dt=INT64_C(20000);
    }

  if((sied->processed_scope_event==false)&&(ied->event_type!=ET_PED))
    veto_scope_event=true;
}

void VBFAnalysisStage2::
visitChannel(bool& veto_channel, void* user_data,
	     uint32_t                  channel_num, 
	     bool                      hit, 
	     bool                      trigger)
{
  IED* ied = static_cast<IED*>(user_data);
  ScopeIED* sied(ied->scope_ied);

  if(!sied->processed_scope_event)return;

  ied->chan_has_signal[channel_num] = hit;
  ied->chan_hit[channel_num]        = hit;
  ied->chan_trigger[channel_num]    = trigger;
  ied->this_chan_trigger            = trigger;
  if((m_diagnostics)&&(trigger))
    {
      ied->scope_diag->recordTriggeredChannel(channel_num,
					      ied->scope_ied->l3_l2_received);
      ied->scope_nchan_trigger++;
    }
}

void VBFAnalysisStage2::
visitHitChannel(void* user_data,
		uint32_t               channel_num,
		uint32_t               charge, 
		uint32_t               pedestal,
		bool                   lo_gain,
		unsigned               nsample,
		const uint32_t*        samples,
		const uint32_t*        integrated)
{
  IED* ied = static_cast<IED*>(user_data);

  if((ied->event_type==ET_PED)&&(m_slow_diagnostics))
    {
      ied->scope_diag->recordPedestal(channel_num, nsample, samples);
      return;
    }

  ScopeIED* sied(ied->scope_ied);
  if(!sied->processed_scope_event)return;

  ied->chan_logain[channel_num]    = lo_gain;

  const ScopeMergedCal& scal(*ied->scope_cal);

  const unsigned sample_zero = lo_gain?scal.lo_sample_zero:scal.hi_sample_zero;
  double ped;
  if(lo_gain)
    {  
      ied->scope_nchan_logain++;      
      ped = scal.ch_ped_lo[channel_num];
      if(scal.ch_suppress[channel_num]||scal.ch_suppress_lo_gain[channel_num])
	ied->chan_has_signal[channel_num] = false;
    }
  else
    {
      ped = scal.ch_ped_hi[channel_num];
      if(scal.ch_suppress[channel_num])
	ied->chan_has_signal[channel_num] = false;
    }

  const double lo_gain_mult = scal.ch_hi_lo_gain_ratio[channel_num];

  double sij;
  double tij;
  unsigned window_width;
  double window_start;
  double total_signal;
  const uint32_t* peak_ptr;

  m_qcalc->calc(sij, tij, window_width, window_start, total_signal, peak_ptr,
		lo_gain, lo_gain_mult, nsample, samples, integrated, 
		ped, sample_zero);

  ied->chan_sij[channel_num]       = sij;
  ied->chan_sij_total[channel_num] = total_signal;
  ied->chan_tij[channel_num]       = tij;

  ied->scope_nchan_hit++;

  PartialScopeDiagnostics* sdiag = ied->scope_diag;
  if(sdiag)
    sdiag->recordHitChannel(sied->l3_l2_received,
			    channel_num, lo_gain, ied->this_chan_trigger,
			    *peak_ptr, sij, tij, nsample, samples,
			    m_slow_diagnostics);

  if(scal.ch_has_pmt[channel_num])
    {
      ied->scope_total_signal+=total_signal;
      
      if(sij > ied->chan_sij_max3_raw)
	if(sij > ied->chan_sij_max2_raw)
	  if(sij > ied->chan_sij_max1_raw)
	    ied->chan_sij_max3_raw     = ied->chan_sij_max2_raw, 
	      ied->chan_sij_max3_raw_j = ied->chan_sij_max2_raw_j,
	      ied->chan_sij_max2_raw   = ied->chan_sij_max1_raw, 
	      ied->chan_sij_max2_raw_j = ied->chan_sij_max1_raw_j,
	      ied->chan_sij_max1_raw   = sij, 
	      ied->chan_sij_max1_raw_j = channel_num;
	  else
	    ied->chan_sij_max3_raw     = ied->chan_sij_max2_raw, 
	      ied->chan_sij_max3_raw_j = ied->chan_sij_max2_raw_j,
	      ied->chan_sij_max2_raw   = sij,
	      ied->chan_sij_max2_raw_j = channel_num;
	else
	  ied->chan_sij_max3_raw       = sij,
	    ied->chan_sij_max3_raw_j   = channel_num;
    }
}

VBFAnalysisStage2::IEDRawDataFetcher::~IEDRawDataFetcher()
{
  // nothing to see here
}

bool VBFAnalysisStage2::IEDRawDataFetcher::
fetchRawData(unsigned npixel, double* raw_data, double* raw_total_data)
{
  static const double missing_val = std::numeric_limits<double>::quiet_NaN();
  const unsigned ncal = m_scal->channel.size();

  unsigned nloop = npixel;
  if(nloop>ncal)nloop=ncal;
  if(nloop>m_ied->scope_nchan)nloop=m_ied->scope_nchan;

  unsigned ipixel;

  if(raw_data)
    {
      for(ipixel=0;ipixel<nloop;ipixel++)
	{
	  if(m_ied->chan_has_signal[ipixel])
	    raw_data[ipixel] = 
	      m_ied->chan_sij[ipixel] * m_scal->channel[ipixel].gain;
	  else
	    raw_data[ipixel] = missing_val;
	}
      while(ipixel<npixel)raw_data[ipixel] = missing_val;
    }

  if(raw_total_data)
    {
      for(ipixel=0;ipixel<nloop;ipixel++)
	{
	  if(m_ied->chan_has_signal[ipixel])
	    raw_total_data[ipixel] = 
	      m_ied->chan_sij_total[ipixel] * m_scal->channel[ipixel].gain;
	  else
	    raw_total_data[ipixel] = missing_val;
	}
      while(ipixel<npixel)raw_total_data[ipixel] = missing_val;
    }

  return true;
}

void VBFAnalysisStage2::leaveScopeEvent(bool veto_scope_event, void* user_data)
{
  IED* ied = static_cast<IED*>(user_data);
  ScopeIED* sied(ied->scope_ied);

  if((veto_scope_event)||(!sied->processed_scope_event))
    {
      const unsigned iscope = ied->scope_num;
#warning TEMPORARY ASSERT
      vsassert(ied->images[iscope].pixels.empty());
      if(ied->processed_array_event)
	{
	  ied->images[iscope].pixels.clear();
	  if(m_secondary_clean.first)
	    ied->secondary_images[iscope].pixels.clear();
	}
      return;
    }

  const unsigned iscope = ied->scope_num;
  const unsigned nchan = ied->scope_nchan;

  if(nchan==0)
    {
#warning TEMPORARY ASSERT
      vsassert(0);
      ied->images[iscope].pixels.clear();
      if(m_secondary_clean.first)
	ied->secondary_images[iscope].pixels.clear();
      return;
    }

#warning TEMPORARY ASSERT
  vsassert(iscope < ied->tsd->cal.scope.size());
  ScopeMergedCal& scal(ied->tsd->cal.scope[iscope]);

  // --------------------------------------------------------------------------
  // Set up for cleaning
  // --------------------------------------------------------------------------
  
  double* level_ij = FASTCALLOC(double,nchan);
  unsigned* masked_ij = FASTCALLOC(unsigned,nchan);

  bool l2_correction_enabled = m_settings.enable_l2_corrections;

  for(unsigned ichan=0;ichan<nchan;ichan++)
    {
      if(!ied->chan_has_signal[ichan] && m_settings.pad_zero_suppressed_chan)
	{
	  ied->chan_sij[ichan] = 
	    ied->tsd->rng->Normal()*scal.channel[ichan].dev;
	}
      else if(!ied->chan_has_signal[ichan])
	{
	  masked_ij[ichan] = true;
	  continue;
	}

      if(l2_correction_enabled)
	{
	  unsigned l2chan = scal.ch_l2chan[ichan];
	  if((l2chan==nchan)||(!ied->chan_hit[l2chan]))
	    l2_correction_enabled=false;
#warning some diagnostic counter here
	}
	  
      if(scal.ch_dopad[ichan])
	ied->chan_sij[ichan] += ied->tsd->rng->Normal()*scal.ch_pad[ichan];
      level_ij[ichan] = ied->chan_sij[ichan] * scal.ch_clean_multiplier[ichan];
      masked_ij[ichan] = false;
    }

  PartialScopeDiagnostics* sdiag = ied->scope_diag;
  VSAReconstruction::ScopeImage& simage(ied->images[iscope]);
  
  // --------------------------------------------------------------------------
  // Clean image
  // --------------------------------------------------------------------------

  const unsigned nimage = m_clean->clean(nchan, masked_ij, level_ij);

  // --------------------------------------------------------------------------
  // Fill scope level reconstruction structure
  // --------------------------------------------------------------------------
  
  double sij_max1_sig = 0;
  double sij_max2_sig = 0;
  double sij_max3_sig = 0;
  unsigned sij_max1_sig_j = nchan;
  unsigned sij_max2_sig_j = nchan;
  unsigned sij_max3_sig_j = nchan;

  if(nimage >= m_settings.nimage_cut)
    {
      sied->has_image = true;
      ied->event_has_image_mask |= 1<<iscope;

      // Transfer to image & do some diagnostics ------------------------------

      simage.pixels.resize(nimage);

      unsigned ipix=0;
      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  if(!masked_ij[ichan])
	    {
	      double sij = ied->chan_sij[ichan] * scal.ch_gain[ichan];
	      double tij = ied->chan_tij[ichan];
	      if(l2_correction_enabled)
		{
		  tij -= ied->chan_tij[scal.ch_l2chan[ichan]];
		  tij -= scal.ch_chantime[ichan];
		}
	      tij *= 2.0; // from samples to NS

	      if(sij > sij_max3_sig)
		if(sij > sij_max2_sig)
		  if(sij > sij_max1_sig)
		    sij_max3_sig     = sij_max2_sig, 
		      sij_max3_sig_j = sij_max2_sig_j,
		      sij_max2_sig   = sij_max1_sig,
		      sij_max2_sig_j = sij_max1_sig_j,
		      sij_max1_sig   = sij, 
		      sij_max1_sig_j = ichan;
		  else 
		    sij_max3_sig     = sij_max2_sig,
		      sij_max3_sig_j = sij_max2_sig_j,
		      sij_max2_sig   = sij, 
		      sij_max2_sig_j = ichan;
		else 
		  sij_max3_sig       = sij, 
		    sij_max3_sig_j   = ichan;
	      
	      if(sdiag)sdiag->recordImageChannel(ichan, 
						 ied->chan_logain[ichan], 
						 sij, tij,
						 m_slow_diagnostics);

	      simage.pixels[ipix].j = ichan;
	      simage.pixels[ipix].nij = sij*scal.gain;
	      simage.pixels[ipix].tij = tij;

	      ipix++;
	    }
	}

      simage.use_in_reconstruction = true;
      ied->nscope_image++;
    }

  // --------------------------------------------------------------------------
  // Secondary cleaning
  // --------------------------------------------------------------------------

  if(m_secondary_clean.first)
    {
      // Reset channel mask ---------------------------------------------------
      for(unsigned ichan=0;ichan<nchan;ichan++)
	masked_ij[ichan] =
	  ((!ied->chan_hit[ichan])||(scal.ch_suppress[ichan]));
      
      // Clean image ----------------------------------------------------------
      const unsigned snimage = 
	m_secondary_clean.first->clean(nchan, masked_ij, level_ij);

      // Transfer to secondary image structure --------------------------------
      VSAReconstruction::ScopeImage& ssimage(ied->secondary_images[iscope]);

      if(snimage >= m_settings.nimage_cut)
	{
	  ssimage.pixels.resize(snimage);

	  unsigned ipix=0;
	  for(unsigned ichan=0;ichan<nchan;ichan++)
	    {
	      if(!masked_ij[ichan])
		{
		  double sij = ied->chan_sij[ichan] * scal.ch_gain[ichan];
		  double tij = ied->chan_tij[ichan];
		  if(l2_correction_enabled)
		    {
		      tij -= ied->chan_tij[scal.ch_l2chan[ichan]];
		      tij -= scal.ch_chantime[ichan];
		    }
		  tij *= 2.0; // from samples to NS
		  ssimage.pixels[ipix].j = ichan;
		  ssimage.pixels[ipix].nij = sij*scal.gain;
		  ssimage.pixels[ipix].tij = tij;
		  ipix++;
		}
	    }

	  ssimage.use_in_reconstruction = true;
	}
    }

  // --------------------------------------------------------------------------
  // Delete temporary arrays
  // --------------------------------------------------------------------------

  FASTFREE(level_ij);
  FASTFREE(masked_ij);

  // --------------------------------------------------------------------------
  // Pointing and astronomy
  // --------------------------------------------------------------------------

  lockPointing();

  // Pass the L3 Az and El to the L3Pointing (HACK)
  VSDirectL3Pointing::getInstance()->
    setAzElRad(iscope, sied->l3_az, sied->l3_el);

  // Get Az and Zn
  bool got_azel =
    m_pointing->getAzZn(sied->az_rad, sied->zn_rad, iscope, ied->event_time);

  unlockPointing();

  if(got_azel)
    {
      // Astronomy and mean pointing --------------------------------------

      SphericalCoords radec(sied->zn_rad,sied->az_rad);
      //Astro::azElToRaDec(lmst, m_earth_position, radec);
      Astro::azElToMeanRaDec(ied->lmst, ied->mjd, 
			     m_settings.earth_position, radec);
      SphericalCoords gal(radec);
      Astro::raDecToGal(ied->mjd, gal);

      sied->has_azzn         = true;
      sied->ra_rad           = radec.longitudeRad();
      sied->dec_rad          = radec.latitudeRad();
      sied->l_rad            = gal.longitudeRad();
      sied->b_rad            = gal.latitudeRad();

      simage.azimuth         = sied->az_rad;
      simage.zenith          = sied->zn_rad;
      simage.camera_rotation = 0;

      // Accumulate scope positions to get mean pointing
      if(!ied->got_one_position)
	{
	  ied->zn_zero = simage.zenith;
	  ied->az_zero = simage.azimuth;
	  ied->mean_position_x.accumulate(0);
	  ied->mean_position_y.accumulate(0);
	  ied->got_one_position = true;
	}
      else
	{
	  SphericalCoords c(simage.zenith,simage.azimuth);
	  c.rotate(0,-ied->zn_zero,-ied->az_zero);
	  ied->mean_position_x.accumulate(c.theta()*cos(c.phi()));
	  ied->mean_position_y.accumulate(c.theta()*sin(c.phi()));
	}
    }
  else
    {
      sied->has_azzn = false;
      std::cerr << ied->event_num << ": Could not get T" << iscope+1
		<< " positional information" << std::endl;
      simage.use_in_reconstruction = false;
    }

  // --------------------------------------------------------------------------
  // Do some diagnostics 
  // --------------------------------------------------------------------------

  sied->nchan_hit          = ied->scope_nchan_hit;
  sied->nchan_logain       = ied->scope_nchan_logain;
  sied->nchan_trigger      = ied->scope_nchan_trigger;
  sied->nchan_image        = nimage;
  sied->total_signal       = ied->scope_total_signal;
 
  sied->trigger_j          = nchan;

  if((m_slow_diagnostics)&&(sdiag))
    doTriggerDiagnostics(ied->chan_trigger, ied->chan_sij, 
			 sied->l3_l2_received,
			 sdiag, sied->largest_region_nchan, sied->trigger_j);

  sied->sij_max1_raw_j = ied->chan_sij_max1_raw_j;
  sied->sij_max2_raw_j = ied->chan_sij_max2_raw_j;
  sied->sij_max3_raw_j = ied->chan_sij_max3_raw_j;  
  sied->sij_max1_sig_j = sij_max1_sig_j;
  sied->sij_max2_sig_j = sij_max2_sig_j;
  sied->sij_max3_sig_j = sij_max3_sig_j;
  sied->sij_max1_sig   = sij_max1_sig;
  sied->sij_max2_sig   = sij_max2_sig;
  sied->sij_max3_sig   = sij_max3_sig;

  // --------------------------------------------------------------------------
  // Muon analysis
  // --------------------------------------------------------------------------

  if((iscope<m_muon_analysis.size())&&(m_muon_analysis[iscope]))
    {
      IEDRawDataFetcher fetcher(ied, &scal);

      sied->has_muon =
	m_muon_analysis[iscope]->
	analyze(sied->muon_data, ied->event_num, ied->images[iscope],
		&fetcher);
    }
}

void VBFAnalysisStage2::
visitSimulationHeader(void* user_data,
		      uint32_t           run_number,
		      const VSTime&      date_of_sims,
		      uint32_t           simulation_package,
		      const std::string& simulator_name,
		      const VSTime&      date_of_array,
		      uint32_t           corsika_atm_model,
		      float              obs_altitude_m,
		      const std::vector<VArrayConfiguration>& array,
		      const std::string& stringified_config)
{
  IED* ied = static_cast<IED*>(user_data);
  vsassert(ied->sim_header == 0);
  ied->sim_header = new VSHeaderSimulationDatum;
  VSHeaderSimulationDatum* header = ied->sim_header;

  header->run_number          = run_number;
  header->date_of_sims        = date_of_sims;
  header->simulation_package  = simulation_package;
  header->simulator_name      = simulator_name;
  header->date_of_array       = date_of_array;
  header->corsika_atm_model   = corsika_atm_model;
  header->obs_altitude_m      = obs_altitude_m;
  header->sim_config          = stringified_config;

  m_sim_ntel                  = array.size();

  if(!m_sim_ebinner.empty())
    {
      header->database_name = "not applicable";
      header->tables.resize(m_sim_ntable);
      for(std::map<unsigned,EBinner>::const_iterator ibinner = 
	    m_sim_ebinner.begin(); ibinner!=m_sim_ebinner.end(); ibinner++)
	{
	  const unsigned ntable = ibinner->second.ntable();
	  const unsigned ztable = ibinner->second.ztable();
	  for(unsigned itable=0;itable<ntable;itable++)
	    {
	      double etev = std::pow(10,ibinner->second.energy(itable));
	      unsigned zitable = ztable+itable;

	      VSTableSimulationDatum& vs(header->tables[itable]);
	      vs.table_index          = zitable;
	      vs.table_db_id          = zitable;
	      vs.primary_id           = ibinner->first;
	      vs.energy_tev           = etev;
	      vs.zenith_min_deg       = 0.0;
	      vs.zenith_max_deg       = 90.0;
	      vs.azimuth_min_deg      = 0.0;
	      vs.azimuth_max_deg      = 360.0;
	      vs.optics_id            = 0;
	      vs.sampling_radius_m    = 0;
	      vs.event_count          = 0;
	      vs.table_name           = 
		std::string("pseudo/")+VSDataConverter::toString(zitable);
	    }
	}
    }
}

void VBFAnalysisStage2::
visitChiLAHeader(void* user_data,
		 const std::string&    database_name,
		 const std::vector<VChiLASimParamTable>&
		                       sim_param_tables,
		 const std::vector<VChiLAOpticsConfiguration>&
		                       optics_configurations,
		 const std::vector<VChiLAElectronicsConfiguration>& 
		                       electronics_configurations)
{
  IED* ied = static_cast<IED*>(user_data);
  vsassert(ied->sim_header);
  VSHeaderSimulationDatum* header = ied->sim_header;
  vsassert(header);
  
  header->database_name       = database_name;
  unsigned ntable = sim_param_tables.size();
  header->tables.resize(ntable); 
  for(unsigned itable=0; itable<ntable; itable++)
    {
      VSTableSimulationDatum& vs(header->tables[itable]);
      const VChiLASimParamTable& vbf(sim_param_tables[itable]);
      vs.table_index          = itable;
      vs.table_db_id          = vbf.fTableID;
      vs.primary_id           = vbf.fPrimaryID;
      vs.energy_tev           = vbf.fEnergyTeV;
      vs.zenith_min_deg       = Angle::toDeg(vbf.fZenithMinRad);
      vs.zenith_max_deg       = Angle::toDeg(vbf.fZenithMaxRad);
      vs.azimuth_min_deg      = Angle::toDeg(vbf.fAzimuthMinRad);
      vs.azimuth_max_deg      = Angle::toDeg(vbf.fAzimuthMaxRad);
      vs.optics_id            = vbf.fOpticsID;
      vs.sampling_radius_m    = vbf.fSamplingRadiusM;
      vs.event_count          = vbf.fEventCount;
      vs.table_name           = vbf.fTableName;
    }
}

void VBFAnalysisStage2::
visitKascadeHeader(void* user_data,
		   uint32_t            corsika_particle_id,
		   float               energy_gev,
		   uint32_t            shower_id,
		   float               x_area_width_m,
		   float               y_area_width_m,
		   uint32_t            north_south_grid)
{
  // nothing to see here
}
 
void VBFAnalysisStage2::
visitSimulationEvent(bool& veto_packet, void* user_data,
		     uint32_t          run_number,
		     uint32_t          event_num,
		     uint32_t          corsika_particle_id,
		     float             energy_gev,
		     float             obs_zenith_deg,
		     float             obs_azimuth_deg,
		     float             primary_zenith_deg,
		     float             primary_azimuth_deg,
		     float             ref_zenith_deg,
		     float             ref_azimuth_deg,
		     float             ref_position_angle_deg,
		     float             core_east_m,
		     float             core_south_m,
		     float             core_elevation_asl_m)
{
  IED* ied = static_cast<IED*>(user_data);
  ied->has_sim = true;
  if(ied->sim == 0)ied->sim = new VSArraySimulationDatum;
  VSArraySimulationDatum* sim = ied->sim;

  sim->event_num               = event_num;
  sim->corsika_particle_id     = corsika_particle_id;
  sim->energy_tev              = energy_gev*0.001;
  sim->obs_zenith_deg          = obs_zenith_deg;
  sim->obs_azimuth_deg         = obs_azimuth_deg;
  sim->primary_zenith_deg      = primary_zenith_deg;
  sim->primary_azimuth_deg     = primary_azimuth_deg;
  sim->ref_zenith_deg          = ref_zenith_deg;
  sim->ref_azimuth_deg         = ref_azimuth_deg;
  sim->ref_position_angle_deg  = ref_position_angle_deg;
  sim->core_east_m             = core_east_m;
  sim->core_north_m            = -core_south_m;
  sim->core_elevation_asl_m    = core_elevation_asl_m;

  if(corsika_particle_id > 0)
    {
      VSAAlgebra::Vec3D e; 
      double szn = sin(sim->primary_zenith_deg*M_PI/180.0);
      double saz = sin(sim->primary_azimuth_deg*M_PI/180.0);
      double caz = cos(sim->primary_azimuth_deg*M_PI/180.0);
      e.set(szn*saz,szn*caz,cos(sim->primary_zenith_deg*M_PI/180.0));
      
      VSAAlgebra::Vec3D R(sim->core_east_m,sim->core_north_m,0);
      R = R - (R*e)*e;
      
      VSAAlgebra::Vec3D e2 = e^VSAAlgebra::Vec3D(0,0,1);      
      if(e2.norm() <= 0) e2 = VSAAlgebra::Vec3D(0,0,1);
      e2.normalize();
      R.rotate(e2*e.theta());
      
      sim->core_x_m                = R.x();
      sim->core_y_m                = R.y();
      sim->core_R_m                = R.norm();

      SphericalCoords primary_azzn =
	SEphem::SphericalCoords::makeDeg(sim->primary_zenith_deg,
					 sim->primary_azimuth_deg);
      SphericalCoords obs_azzn =
	SEphem::SphericalCoords::makeDeg(sim->obs_zenith_deg, 
					 sim->obs_azimuth_deg);
	  
      SphericalCoords primary_radec = primary_azzn;

      m_sim_transform->transform(primary_radec,sim->corsika_particle_id,
				 primary_azzn,obs_azzn);

      sim->primary_dec_deg = primary_radec.latitudeDeg();
      sim->primary_ra_deg = primary_radec.phiDeg();
    }
  else
    {
      sim->primary_dec_deg         = 0;
      sim->primary_ra_deg          = 0;
      sim->core_x_m                = 0;
      sim->core_y_m                = 0;
      sim->core_R_m                = 0;
    }

  sim->table_index             = 0xFFFFFFFFU;
  sim->electronics_id          = 0;
  sim->table_event_index       = 0;

  if(m_sim_ntable)
    {
      std::map<unsigned,EBinner>::const_iterator ibin = 
	m_sim_ebinner.find(corsika_particle_id);
      if(ibin != m_sim_ebinner.end())
	sim->table_index = ibin->second.table(log10(sim->energy_tev));
    }

  ied->sim_aomega              = 0;
}

void VBFAnalysisStage2::
visitChiLAEvent(bool& veto_packet, void* user_data,
		uint32_t               table_index,
		uint32_t               electronics_id,
		uint32_t               table_event_index,
		float                  a_omega,
		float                  a_omega_var)
{
  IED* ied = static_cast<IED*>(user_data);
  vsassert(ied->has_sim);
  VSArraySimulationDatum* sim = ied->sim;
  vsassert(sim);

  sim->table_index           = table_index;
  sim->electronics_id        = electronics_id;
  sim->table_event_index     = table_event_index;
}

void VBFAnalysisStage2::
visitKascadeEvent(bool& veto_packet, void* user_data,
		  int32_t              nx,
		  int32_t              ny,
		  int32_t              direction_index,
		  float                emission_altitude_m,
		  float                emission_altitude_sigma,
		  float                muon_ratio,
		  float                a_omega,
		  float                rel_tel_trigger_time_ns,
		  float                differential_rate_per_event_hz,
		  float                integral_rate_per_event_hz)
{
  IED* ied = static_cast<IED*>(user_data);
  vsassert(ied->has_sim);
  VSArraySimulationDatum* sim = ied->sim;
  vsassert(sim);

  ied->sim_aomega            = a_omega;
}

void VBFAnalysisStage2::
leaveArrayEvent(bool veto_array_event, void* user_data)
{
  IED* ied = static_cast<IED*>(user_data);
  if((veto_array_event)||(ied->processed_array_event==false))return;

  const unsigned nscope = ied->scope.size();

  const VSAReconstruction::ArrayMoments* am = &ied->recon.moments;

  // --------------------------------------------------------------------------
  // Reconstruction
  // --------------------------------------------------------------------------

  if(ied->nscope_image >= m_settings.nscope_cut)
    {
      ied->reconstruction_attempted = true;

      ied->reconstruction_successful = 
	m_reconstruction->
	reconstruct(ied->recon, m_settings.method, m_settings.weighting,
		    ied->images);

      if(ied->reconstruction_successful)
	{
	  for(unsigned iscope=0;iscope<nscope;iscope++)
	    if(am->scopes[iscope].use_in_reconstruction)
	      ied->event_used_in_reconstruction_mask |= 1<<iscope;
	  
	  if(!m_settings.qc_masks.empty())
	    {
	      ied->reconstruction_successful = false;

	      for(std::vector<unsigned>::const_iterator imask =
		    m_settings.qc_masks.begin(); 
		  imask != m_settings.qc_masks.end(); imask++)
		if((ied->event_used_in_reconstruction_mask & *imask) == *imask)
		  {
		    ied->reconstruction_successful = true;
		    break;
		  }
	    }
	}

      ied->nscope_quality = ied->recon.moments.nuse_in_reconstruction;

      // ----------------------------------------------------------------------
      // Secondary moments
      // ----------------------------------------------------------------------

      if(m_secondary_clean.first)
	{

	  m_reconstruction->
	    calculateArrayMoments(ied->secondary_moments,
				  m_settings.weighting,
				  ied->secondary_images,
				  m_secondary_clean.second);
	  
	  am = &ied->secondary_moments;
	}

      // ----------------------------------------------------------------------
      // Calculate per telescope Hillas parameters (length/width/psi)
      // ----------------------------------------------------------------------

      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  const VSAReconstruction::ScopeMoments& sm(am->scopes[iscope]);

	  ScopeIED* sied(ied->scope[iscope]);

	  sied->used_in_reconstruction = 
	    ied->reconstruction_successful &&
	    ied->recon.moments.scopes[iscope].use_in_reconstruction;

	  if(sm.Ni > 0)
	    {
	      sied->has_moments   = true;
	      sied->fp_Ni         = sm.Ni;
	      sied->fp_xc         = Angle::toDeg(sm.Vi_cam.x());
	      sied->fp_yc         = Angle::toDeg(sm.Vi_cam.y());

	      const double fp_dist2 = 
		sied->fp_xc*sied->fp_xc+sied->fp_yc*sied->fp_yc;
	      const double fp_dist = std::sqrt(fp_dist2);

	      sied->fp_dist       = fp_dist;

	      const VSAAlgebra::Symmetric2D 
		m2(sm.Ti_cam - VSAAlgebra::Symmetric2D(sm.Vi_cam));

	      VSAAlgebra::Eigen2D e;
	      m2.eigen(e);

	      // This is tricky - the order of the coefficients in the
	      // std::max is important when one of them is a NAN,
	      // because max(NAN,0.0)=NAN while max(0.0,NAN)=0.0

	      sied->fp_length     = Angle::toDeg(sqrt(std::max(0.0,e.val[1])));
	      sied->fp_width      = Angle::toDeg(sqrt(std::max(0.0,e.val[0])));

	      if(sm.Vi_cam*e.vec[1] >= 0)
		sied->fp_psi      = 
		  Angle::toDeg(atan2(e.vec[1].y(), e.vec[1].x()));
	      else 
		sied->fp_psi      = 
		  Angle::toDeg(atan2(-e.vec[1].y(), -e.vec[1].x()));

	      if(!m_settings.psf_poly_radial.empty())
		{
		  double psf_r2 = 
		    VSAMath::PolyFit::val(m_settings.psf_poly_radial,
					  fp_dist2);
		  psf_r2 *= psf_r2;

		  double psf_l2 = psf_r2;
		  double psf_w2 = psf_r2;
		  
		  if(!m_settings.psf_poly_tangential.empty())
		    {
		      double psf_t2 = 
			VSAMath::PolyFit::val(m_settings.psf_poly_tangential, 
					      fp_dist2);
		      psf_t2 *= psf_t2;

		      const double exr = sied->fp_xc/fp_dist;
		      const double eyr = sied->fp_yc/fp_dist;
		      const double exl = e.vec[1].x();
		      const double eyl = e.vec[1].y();

		      double erl2 = exl*exr + eyl*eyr;
		      erl2 *= erl2;

		      double etl2 = eyl*exr - exl*eyr;
		      etl2 *= etl2;

		      psf_l2 = psf_r2*erl2 + psf_t2*etl2;
		      psf_w2 = psf_r2*etl2 + psf_t2*erl2;
		    }

#if 0
		  std::cout << psf_r2 << ' ' << psf_l2 << ' ' << psf_r2 << ' '
			    << sied->fp_length*sied->fp_length << ' ' 
			    << sied->fp_width*sied->fp_width
			    << '\n';
#endif

		  sied->intrinsic_length = 
		    sied->fp_length*sied->fp_length - psf_l2;
		  sied->intrinsic_width  = 
		    sied->fp_width*sied->fp_width - psf_w2;
		  
		  if(sied->intrinsic_length>0)
		    sied->intrinsic_length = sqrt(sied->intrinsic_length);
		  else
		    sied->intrinsic_length = 0;

		  if(sied->intrinsic_width>0)
		    sied->intrinsic_width = sqrt(sied->intrinsic_width);
		  else
		    sied->intrinsic_width = 0;
		}
	      else
		{
		  sied->intrinsic_length = sied->fp_length;
		  sied->intrinsic_width  = sied->fp_width;
		}
	    }
	  else
	    {
	      sied->has_moments   = false;
	      sied->fp_Ni         = 0;
	      sied->fp_xc         = 0;
	      sied->fp_yc         = 0;
	      sied->fp_dist       = 0;
	      sied->fp_length     = 0;
	      sied->fp_width      = 0;
	      sied->fp_psi        = 0;
	      sied->intrinsic_length = 0;
	      sied->intrinsic_width  = 0;
	    }
	}
    }
  else
    {
      ied->nscope_quality = 0;
    }

  // --------------------------------------------------------------------------
  // Calculate array parameters if reconstruction succeeded
  // --------------------------------------------------------------------------

  if(ied->reconstruction_successful)
    {
      m_reconstruction->
	calculateArrayParameters(ied->param, ied->images, ied->recon);

      // ----------------------------------------------------------------------
      // Calculate reconstructed direction in all interesting geometries
      // ----------------------------------------------------------------------

      const double mpx = ied->mean_position_x.mean();
      const double mpy = ied->mean_position_y.mean();
      SphericalCoords mean_position(norm(mpy,mpx),atan2(mpy,mpx));
      mean_position.rotate(ied->az_zero,ied->zn_zero,0);

      const double mean_zn = mean_position.theta();
      const double mean_az = mean_position.phi();

      ied->mean_zn = mean_zn;
      ied->mean_az = mean_az;

      const double smean_zn = sin(mean_zn);
      const VSAAlgebra::Vec3D mean_e(-smean_zn*sin(mean_az),
				     -smean_zn*cos(mean_az),
				     -cos(mean_zn));

      const double theta0 = acos(ied->recon.e*mean_e)/M_PI*180.0;

      const double zn = atan2(norm(ied->recon.e.x(), ied->recon.e.y()), 
			      -ied->recon.e.z());
      
      const double az = fmod(atan2(-ied->recon.e.x(),-ied->recon.e.y())
			     +2.0*M_PI,2.0*M_PI);

      ied->recon_zn = zn;
      ied->recon_az = az;
      
      SphericalCoords fov;
      fov.setThetaPhi(zn,M_PI-az);
      fov.rotate(0,-mean_zn,mean_az-M_PI);

      ied->recon_fov_x = fov.thetaRad()*sin(fov.phiRad());
      ied->recon_fov_y = fov.thetaRad()*cos(fov.phiRad());

      SphericalCoords radec(zn,az);
      SphericalCoords mean_radec(mean_zn, mean_az);
 
      if(ied->has_sim)
	{
	  VSArraySimulationDatum* sim = ied->sim;

	  SphericalCoords primary_azzn =
 	    SEphem::SphericalCoords::makeDeg(sim->primary_zenith_deg,
 					     sim->primary_azimuth_deg);
	  SphericalCoords obs_azzn =
	    SEphem::SphericalCoords::makeDeg(sim->obs_zenith_deg, 
					     sim->obs_azimuth_deg);
	  
	  m_sim_transform->transform(radec,sim->corsika_particle_id,
				     primary_azzn,obs_azzn);
	  m_sim_transform->transform(mean_radec,sim->corsika_particle_id,
				     primary_azzn,obs_azzn);
	  
	}
      else
	{
	  Astro::azElToMeanRaDec(ied->lmst, ied->mjd, 
				 m_settings.earth_position, radec);
	  Astro::azElToMeanRaDec(ied->lmst, ied->mjd, 
				 m_settings.earth_position, mean_radec);  
	}	  

      ied->recon_ra = radec.longitudeRad();
      ied->recon_dec = radec.latitudeRad();

      const double mean_codec = mean_radec.theta();
      const double mean_ra = mean_radec.phi();

      SphericalCoords derotated_fov;
      derotated_fov.setLatLongRad(ied->recon_dec,ied->recon_ra);
      derotated_fov.rotate(0,-mean_codec,-mean_ra);

      ied->recon_derotated_fov_x 
	= -derotated_fov.thetaRad()*sin(derotated_fov.phiRad());
      ied->recon_derotated_fov_y 
	= -derotated_fov.thetaRad()*cos(derotated_fov.phiRad());

      VSAAlgebra::Vec3D R0;
      if(ied->recon.e.z() >= -DBL_EPSILON)
	R0 = VSAAlgebra::Vec3D(1e8, 1e8, 0);
      else
	R0 = ied->recon.R-ied->recon.e*(ied->recon.R.z()/ied->recon.e.z());

      ied->recon_r0_x = R0.x();
      ied->recon_r0_y = R0.y();

      // ----------------------------------------------------------------------
      // Calculate theta^2 for all targets
      // ----------------------------------------------------------------------

      ied->theta[0] = theta0;

      for(unsigned itheta=0;itheta<m_settings.theta_targets.size(); itheta++)
	{
	  ied->theta[itheta+1] = 
	    radec.separation(m_settings.theta_targets[itheta].coord).deg();
	}

      if(m_settings.theta_cut>0)
	{
	  for(unsigned itheta=0;itheta<ied->theta.size();itheta++)
	    if(ied->theta[itheta] <= m_settings.theta_cut)
	      {
		ied->write_event = true;
		break;
	      }
	}
      else ied->write_event = true;

      // ----------------------------------------------------------------------
      // Calculate J2000 RA/Dec
      // ----------------------------------------------------------------------
      
      if(!ied->has_sim)
	radec.rotate(m_j2000_phi, m_j2000_theta, m_j2000_psi);

      ied->recon_ra_j2000 = radec.longitudeRad();
      ied->recon_dec_j2000 = radec.latitudeRad();
      
      // ----------------------------------------------------------------------
      // Calculate size2
      // ----------------------------------------------------------------------
      std::vector< double > N;
      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  const VSAReconstruction::ScopeMoments& sm(am->scopes[iscope]);
	  if(sm.use_in_reconstruction) N.push_back(sm.Ni);
	}
      std::sort( N.begin(), N.end(), std::greater<double>() );
      if(N.size() >= 2) ied->N2 = N[1];
      else ied->N2 = 0;

      // ----------------------------------------------------------------------
      // Calculate scaled parameters and lookup table energy
      // ----------------------------------------------------------------------

      VSScaledParameterCalc::State sp_state;
      VSEnergyCalc::State energy_state;

      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  const VSAReconstruction::ScopeMoments& sm(am->scopes[iscope]);
	  
	  ScopeIED* sied(ied->scope[iscope]);
	      
	  VSAAlgebra::Vec3D fp_e = 
	    ied->recon.moments.U_earth_to_event.transFwd(ied->recon.e);
	  fp_e = 
	    ied->recon.moments.scopes[iscope].U_cam_to_event.transBwd(fp_e);

	  sied->fp_ex         = Angle::toDeg(fp_e.x());
	  sied->fp_ey         = Angle::toDeg(fp_e.y());

	  if((sm.Ni > 0)&&(sm.use_in_reconstruction))
	    {
	      double dx = sied->fp_ex-sied->fp_xc;
	      double dy = sied->fp_ey-sied->fp_yc;

	      sied->fp_disp   = norm(dx, dy);

#if 0
	      std::cout << sied->fp_ex << ' ' << sied->fp_ey << ' '
			<< sied->fp_xc << ' ' << sied->fp_xc << ' '
			<< dx << ' ' << dy << ' '
			<< sied->fp_disp << std::endl;
#endif

	      const VSAReconstruction::ScopeParameters&
		sp(ied->param.scopes[iscope]);

	      if(m_sp_calc)
		m_sp_calc->getScopeSP(sp_state, iscope, sp.Ri, sm.Ni, 
				      sied->intrinsic_width,
				      sied->intrinsic_length,
				      sied->fp_disp,
				      sied->fp_dist,
				      sied->sc_width, 
				      sied->sc_length,
				      sied->sc_disp);

	      if(m_energy_calc)
		m_energy_calc->getScopeEnergy(energy_state, iscope,
					      sp.Ri, sm.Ni, 
					      sied->fp_disp,
					      sied->fp_dist,
					      sied->lt_log10_energy,
					      sied->lt_log10_energy_err);
	      
	    }
	  else
	    {
	      sied->fp_disp                 = 0;
	      sied->sc_width                = 0;
	      sied->sc_length               = 0;
	      sied->sc_disp                 = 0;
	      sied->lt_log10_energy         = 0;
	      sied->lt_log10_energy_err     = 0;
	    }
	}
      
      if(m_sp_calc)
	m_sp_calc->getArraySP(sp_state, 
			      ied->msc_width, ied->msc_length, ied->msc_disp);
      else
	{
	  ied->msc_width = POSINF;
	  ied->msc_length = POSINF;
	  ied->msc_disp = POSINF;
	}

      if(m_energy_calc)
	m_energy_calc->getArrayEnergy(energy_state, 
				      ied->mlt_log10_energy,
				      ied->mlt_log10_energy_chi2);
      else 
	{
	  ied->mlt_log10_energy = NEGINF;
      	  ied->mlt_log10_energy_chi2 = NEGINF;
	}

      // ----------------------------------------------------------------------
      // Transfer results to event data structure
      // ----------------------------------------------------------------------

      ied->transferDataToEAD();

      // ----------------------------------------------------------------------
      // Apply cuts
      // ----------------------------------------------------------------------

      if(m_cuts_calc)
	{
	  m_cuts_calc->getCutResults(*ied->acd, ied->ead);
	}
    }
}

void VBFAnalysisStage2::
doTriggerDiagnostics(//const bool* trigger, const double* signal,
		     const std::vector<bool>& trigger, const double* signal,
		     bool l2_trigger_received,
		     PartialScopeDiagnostics* sdiag,
		     unsigned& largest_region_nchan, unsigned& trigger_ichan)
{
  // This complex function does the following trigger diagnostics:
  // 1) Assigns each channel to a region of continuous triggered channels 
  //    using the same methodology as the regional cleaner
  // 2) Find the largest region of channels
  // 3) Find the "triggering channel" in the largest region, that channel
  //    that has the m_settings.l2_trigger_threshold'th largest signal
  // 4) Accumulate various counts in the partial scope diagnostics

  unsigned nnchan = m_settings.neighbor_nchan;

  enum TriggerState { TS_NO_TRIGGER, TS_TRIGGERED, TS_VISITED };
  TriggerState* state = FASTCALLOC(TriggerState,nnchan);

  for(unsigned int ichan=0;ichan<nnchan;ichan++)
    state[ichan] = trigger[ichan]?TS_TRIGGERED:TS_NO_TRIGGER;

  // (1) ----------------------------------------------------------------------

  // Go through all the tubes and assign them to a "region" of
  // continuous tubes that trigger. The "0" region is a pretend one
  // for channels that do not trigger. Hence the initialisation of
  // region_nchan with 1... which skips the "0" entry in the vector

  unsigned* region_nchan = FASTCALLOC(unsigned,nnchan);
  unsigned* channel_stack = FASTCALLOC(unsigned,nnchan);
  unsigned* channel_region = FASTCALLOC(unsigned,nnchan);

  unsigned nregion = 0;
  region_nchan[nregion] = 0; // For tubes that are not part of a region

  for(unsigned ichan=0;ichan<nnchan;ichan++)
    {
      if(state[ichan] == TS_VISITED)continue;
      unsigned channel_stack_occupancy = 0;
      if(state[ichan] == TS_TRIGGERED)
	{
	  region_nchan[++nregion] = 0;
	  channel_stack[channel_stack_occupancy++] = ichan;
	  state[ichan] = TS_VISITED;
	}
      else
	{
	  channel_region[ichan]=0;
	  state[ichan] = TS_VISITED;
	  continue;
	}

      while(channel_stack_occupancy)
	{
	  unsigned jchan = channel_stack[--channel_stack_occupancy];
	  region_nchan[nregion]++;
	  channel_region[jchan]=nregion;

          for(unsigned jneighbor=0;jneighbor<NUM_NEIGHBORS;jneighbor++)
            {
              int neighbor_ichan=m_settings.neighbors[jchan][jneighbor];
	      if(neighbor_ichan==-1)break;
	      if(state[neighbor_ichan]==TS_TRIGGERED)
		{
		  channel_stack[channel_stack_occupancy++] = neighbor_ichan;
		  state[neighbor_ichan] = TS_VISITED;
		}
	      else if(state[neighbor_ichan]==TS_NO_TRIGGER)
		{
		  channel_region[neighbor_ichan]=0;		  
		  state[neighbor_ichan] = TS_VISITED;
		}
	    }
	}
    }
  nregion++;
  
  // (2) ----------------------------------------------------------------------

  largest_region_nchan = 0;
  for(unsigned iregion=1;iregion<nregion;iregion++)
    if(region_nchan[iregion] > largest_region_nchan)
      largest_region_nchan = region_nchan[iregion];
  
  // (3) ----------------------------------------------------------------------

  if((l2_trigger_received)
     &&(largest_region_nchan>=m_settings.l2_trigger_threshold))
    {
      typedef std::pair<double,unsigned> PDU;

      // Bad C++ here!!! I really need to write an allocator
      PDU* signal_vec = FASTCALLOC(PDU,largest_region_nchan);
      unsigned nsignal_vec;
      
      bool found_one = false;
      PDU max_signal(0,0);

      for(unsigned iregion=1;iregion<nregion;iregion++)
	if(region_nchan[iregion] == largest_region_nchan)
	  {
	    nsignal_vec = 0;
	    for(unsigned ichan=0;
		(ichan<nnchan)&&(nsignal_vec<largest_region_nchan);ichan++)
	      if(channel_region[ichan] == iregion)
		{
		  signal_vec[nsignal_vec].first = signal[ichan];
		  signal_vec[nsignal_vec].second = ichan;
		  nsignal_vec++;
		}
	    std::sort(signal_vec, signal_vec+nsignal_vec);
	    PDU& this_max(signal_vec[largest_region_nchan -
				     m_settings.l2_trigger_threshold]);
	    if(!found_one)
	      max_signal = this_max, found_one = true;
	    else if(max_signal < this_max)
	      max_signal = this_max;
	  }

      trigger_ichan = max_signal.second;

      FASTFREE(signal_vec);
    }

  // (4) ----------------------------------------------------------------------

  if(largest_region_nchan)
    {
      for(unsigned ichan=0;ichan<nnchan;ichan++)
	{
	  unsigned iregion = channel_region[ichan];
	  unsigned iregion_nchan = region_nchan[iregion];
	  if(iregion_nchan == 1)
	    if(l2_trigger_received)
	      sdiag->recordChannelIsTriggeredIsolated(ichan);
	    else
	      sdiag->recordChannelIsTriggeredIsolatedNoL2(ichan);
	  else if(iregion_nchan == largest_region_nchan)
	    if(l2_trigger_received)
	      sdiag->recordChannelIsTriggeredInLargestRegion(ichan);
	    else
	      sdiag->recordChannelIsTriggeredInLargestRegionNoL2(ichan);
	}

      if(largest_region_nchan<m_settings.l2_trigger_threshold)
	{
	  if(l2_trigger_received)
	    {
	      for(unsigned ichan=0;ichan<nnchan;ichan++)
		if(trigger[ichan])
		  sdiag->recordChannelIsTriggeredInSubThresholdEvent(ichan);
	    }
	  else
	    {
	      for(unsigned ichan=0;ichan<nnchan;ichan++)
		if(trigger[ichan])
		  sdiag->
		    recordChannelIsTriggeredInSubThresholdEventNoL2(ichan);
	    }
	}
    }

  FASTFREE(region_nchan);
  FASTFREE(channel_stack);
  FASTFREE(channel_region);
  FASTFREE(state);
}

VBFAnalysisStage2::IED* VBFAnalysisStage2::
getEmptyIED(unsigned seq_packet_number, unsigned vbf_packet_number)
{
  lockIEDListAndWait(seq_packet_number);

#warning TEMPORARY ASSERT
  vsassert(!m_ied_free.empty());

  IED* ied = m_ied_free.front();
  m_ied_free.pop_front();

  ied->reset();
  ied->ied_processing = true;
  ied->seq_packet_number = seq_packet_number;
  ied->vbf_packet_number = vbf_packet_number;

#ifndef NOTHREADS
  ied->tsd = static_cast<ThreadSpecificData*>(pthread_getspecific(m_tsd_key));
#else
  ied->tsd = m_tsd_list.empty()?0:m_tsd_list.front();
#endif

  if(ied->tsd == 0)
    {
      bool construct_rng = false;
      if(m_pad_stage1 || m_settings.pad_zero_suppressed_chan)
	construct_rng = true;

      ied->tsd = 
	new ThreadSpecificData(m_stage1, m_settings, *m_cal, m_diagnostics,
			       construct_rng);

#ifndef NOTHREADS
      vsassert(pthread_setspecific(m_tsd_key, ied->tsd) == 0);
#endif
      m_tsd_list.push_front(ied->tsd);
    }

  std::list<IED*>::iterator iied=m_ied_processing.begin();
  while((iied!=m_ied_processing.end())
	&&(ied->seq_packet_number>(*iied)->seq_packet_number))iied++;
  m_ied_processing.insert(iied,ied);

  unlockIEDList();

  return ied;
}

void VBFAnalysisStage2::releaseIED(IED* ied)
{
  lockIEDList();
  m_ied_free.push_front(ied);
  unlockIEDListAndBroadcast();
}

void VBFAnalysisStage2::processConstructedIEDs()
{
  lockIEDList();
  if(!m_ied_writer_active)
    {
      m_ied_writer_active = true;

      IED* ied = m_ied_processing.empty()?0:m_ied_processing.front();
      while((ied)&&(!ied->ied_processing)&&
	    (ied->seq_packet_number==m_ied_writer_next_seq_packet_number))
	{
	  m_ied_writer_next_seq_packet_number++;
	  m_ied_processing.pop_front();
	  unlockIEDList();
	  processConstructedIED(ied);
	  lockIEDList();
	  ied = m_ied_processing.empty()?0:m_ied_processing.front();
	}
      
      m_ied_writer_active = false;
    }
  unlockIEDList();
}

template<typename T> static inline T nearestPowerOfTwo(const T x)
{
  const unsigned nbit = 8*sizeof(T);
  const T mask = ~T(0)>>2 | T(1)<<(nbit-1);
  if((x>mask)||(x==T(0)))return T(0);

  unsigned nshift = nbit>>2;
  unsigned ishift = nbit>>1;

  while(nshift)
    {
      if(x>(mask>>ishift))ishift -= nshift;      
      else ishift += nshift;
      nshift>>=1;
    }

  if(x>(mask>>ishift))ishift--;
  return T(1)<<(nbit-ishift-1);
}

template<typename T> static inline void 
checkBitPattern(const T c0, T& c1, T& c2, const T c3)
{
  // Try to correct bit patterns

  // Assuming a single bit error the sequence of counts in four
  // consecutive packets will be either of the following:

  // N0                             N0
  // N0 +   dN + ERROR              N0 +   dN
  // N0 + 2 dN          ***         N0 + 2 dN - ERROR ***
  // N0 + 3 dN                      N0 + 3 dN

  // N0
  // N0 +   dN
  // N0 + 2 dN - ROLLOVER ***
  // N0 + 3 dN - ROLLOVER

  // Correct the middle two using information from the first and last
  // packets.

  if(c2<c1)
    if((c3>c1)&&(c0>=c2))
      {
	const T dN1 = c3-c1;
	const T dN2 = c0-c2;
	const T correction = nearestPowerOfTwo(dN1+dN2);
	c2+=correction;
      }
    else if((c2>c0)&&(c1>=c3))
      {
	const T dN1 = c2-c0;
	const T dN2 = c1-c3;
	const T correction = nearestPowerOfTwo(dN1+dN2);
	c1-=correction;
      }
}

static inline int64_t ticks32(uint32_t _last, uint32_t _this)
{
  int64_t x = int64_t(_this) - int64_t(_last);
  if(x<INT64_C(-2147483648))x+=INT64_C(4294967296);
  return x;
}

void VBFAnalysisStage2::processConstructedIED(IED* ied)
{
  if(m_diagnostics)m_diagnostics->packetFound();

  if(ied->sim_header)
    {
      vsassert(ied->vbf_packet_number == 0);
      if(!m_sim_header)
	{
	  m_sim_header = ied->sim_header;

	  // Construct writers for simulation events
	  m_asw = new VSArraySimulationWriter(m_io->writeStruct("sim_event"),
					      m_sim_ntel);
	  m_e2sw = m_io->writeExpandableVector<unsigned>("event_to_sim");
	  m_s2ew = m_io->writeExpandableVector<unsigned>("sim_to_event");
	}
      else
	delete ied->sim_header;
      ied->sim_header = 0;
    }

  if(ied->has_sim)
    {
      vsassert(m_asw);

      bool is_pedestal = ied->has_array_event && ied->event_type==ET_PED;
      bool write_sim_event = !is_pedestal;

      if(write_sim_event)
	{
	  ied->sim->has_array_event = ied->has_array_event;
	  m_asw->append(*ied->sim);
      
	  bool valid_table_index = 
	    ied->sim->table_index<m_sim_header->tables.size();

	  if(valid_table_index)
	    {
	      m_sim_header->tables[ied->sim->table_index].num_events_read++;
	      if((ied->has_array_event)&&(ied->processed_array_event))
	       m_sim_header->tables[ied->sim->table_index]
		 .num_events_processed++;

	      if(m_sim_package == VBFRunInfo::SimData::SP_KASCADE)
		m_sim_header->tables[ied->sim->table_index]
		  .sampling_radius_m += ied->sim_aomega;
	    }
      
	  if((ied->has_array_event)&&(ied->write_event))
	    {
	      m_e2sw->append(m_nsim_written);
	      m_s2ew->append(m_nevent_written);

	      if(valid_table_index)
		m_sim_header->tables[ied->sim->table_index]
		  .num_events_written++;
	    }
	  else 
	    {
	      m_s2ew->append(0xFFFFFFFFU);
	    }

	  m_nsim_written++;      
	}
      else if((ied->has_array_event)&&(ied->write_event))
	{
	  m_e2sw->append(0xFFFFFFFFU);
	}
    }
  else if((m_asw)&&(ied->has_array_event)&&(ied->write_event))
    {
      m_e2sw->append(0xFFFFFFFFU);
    }

  if(ied->has_array_event)
    {
      if(ied->write_event)m_nevent_written++;
      if((m_e2sw)&&(!ied->has_sim))m_e2sw->append(0xFFFFFFFFU);

      const unsigned nscope = ied->scope.size();

      if(!m_have_one_array_event)
	{
	  m_first_time                        = 
	    m_stage1->run_info.lo_event_time;
	  //m_first_time                        = ied->event_time;
	  m_last_event_time_sec               = 0;
	  m_last_event_time_hist              = 0;
	  m_last_event_num                    = 0;
	  m_have_one_array_event              = true;
	}

      if(ied->has_event_time)
	{
	  ied->event_time_sec = double(ied->event_time-m_first_time)/1e9;

	  if((ied->event_time_sec>=0)&&(ied->event_time_sec<3600.0))
	    ied->event_time_hist = ied->event_time_sec/60.0;
	  else
	    ied->event_time_hist = -1;
	}
      else
	{
	  ied->event_time_sec = -3;
	  ied->event_time_hist = -3;
	}

      for(unsigned iscope=0;iscope<nscope;iscope++)
	if(ied->scope[iscope]->has_azzn)
	  {
	    if(!m_mean_have_ra_zero)
	      {
		m_mean_ra_zero = ied->scope[iscope]->ra_rad;
		m_mean_have_ra_zero = true;
	      }
	    
	    Angle _ra(ied->scope[iscope]->ra_rad-m_mean_ra_zero);
	    m_mean_ra.accumulate(_ra.radPM());
	    m_mean_dec.accumulate(ied->scope[iscope]->dec_rad);
	  }
      
      if(ied->has_l3)
	{
	  if(!m_have_one_array_event_with_l3)
	    {
	      m_last_l3_event_num             = 0;
	      m_last_event_time_ten_mhz_elap  = 0;
	      m_last_event_time_ten_mhz_live  = 0;
	      m_last_ten_mhz_elapsed          = 0;//ied->ten_mhz_elapsed;
	      m_last_ten_mhz_veto_both        = 0;//ied->ten_mhz_veto_both;
	      m_last_ten_mhz_veto_vdaq        = 0;//ied->ten_mhz_veto_vdaq;
	      m_last_ten_mhz_veto_l3          = 0;//ied->ten_mhz_veto_l3;
	      m_have_one_array_event_with_l3  = true;
	    }      

	  for(unsigned iscope=0;iscope<nscope;iscope++)
	    if(ied->scope[iscope]->has_l3)
	      if(!m_first_scope_has_event[iscope])
		m_last_counts_l2[iscope] = ied->scope[iscope]->l3_counts_l2,
		  m_first_scope_has_event[iscope] = true;
	}

      finalizeIED(ied);
    }
  else
    {
      releaseIED(ied);
    }
}

void VBFAnalysisStage2::finalizeIED(IED* ied)
{
  ied->last_event_time_sec            = m_last_event_time_sec;
  ied->last_event_time_hist           = m_last_event_time_hist;
  ied->last_event_num                 = m_last_event_num;
  ied->last_ten_mhz_elapsed           = m_last_ten_mhz_elapsed;
  ied->last_ten_mhz_veto_both         = m_last_ten_mhz_veto_both;
  ied->last_ten_mhz_veto_vdaq         = m_last_ten_mhz_veto_vdaq;
  ied->last_ten_mhz_veto_l3           = m_last_ten_mhz_veto_l3;

  if(ied->has_l3)
    {
      ied->ticks_elap = 
	ticks32(m_last_ten_mhz_elapsed, ied->ten_mhz_elapsed);
      ied->ticks_both                 =
	ticks32(m_last_ten_mhz_veto_both, ied->ten_mhz_veto_both);
      ied->ticks_vdaq                 =
	ticks32(m_last_ten_mhz_veto_vdaq, ied->ten_mhz_veto_vdaq);
      ied->ticks_lev3                 =
	ticks32(m_last_ten_mhz_veto_l3, ied->ten_mhz_veto_l3);

      // Missed L3 events - figure all time since last L3 event as DEAD
      if(ied->event_num - m_last_l3_event_num > 1)
	ied->ticks_both               = ied->ticks_elap; 
      m_last_l3_event_num             = ied->event_num;

      m_last_ten_mhz_elapsed          = ied->ten_mhz_elapsed;
      m_last_ten_mhz_veto_both        = ied->ten_mhz_veto_both;
      m_last_ten_mhz_veto_vdaq        = ied->ten_mhz_veto_vdaq;
      m_last_ten_mhz_veto_l3          = ied->ten_mhz_veto_l3;	  
    }
  else
    {
      ied->ticks_elap                 = 0;
      ied->ticks_both                 = 0;
      ied->ticks_vdaq                 = 0;
      ied->ticks_lev3                 = 0;
    }

  m_last_event_time_ten_mhz_elap     += ied->ticks_elap;
  m_last_event_time_ten_mhz_live     += ied->ticks_elap-ied->ticks_both;

  for(unsigned iscope=0;iscope<ied->scope.size();iscope++)
    {
      ScopeIED* sied = ied->scope[iscope];

      if(sied->has_l3)
	{
	  sied->counts_l2 =  
	    ticks32(m_last_counts_l2[iscope], sied->l3_counts_l2);
	  m_last_counts_l2[iscope]    = sied->l3_counts_l2;
	}
      else
	{
	  sied->counts_l2             = 0;
	}

      if(sied->has_muon)
	m_muon_data->addDatum(iscope,sied->muon_data);
    }

  if(ied->has_event_time)
    {
      m_last_event_time_sec           = ied->event_time_sec;
      m_last_event_time_hist          = ied->event_time_hist;
    }
  m_last_event_num                    = ied->event_num;

  // --------------------------------------------------------------------------
  // Final transfer results to event data structure
  // --------------------------------------------------------------------------
  
  ied->finalTransferDataToEAD(m_last_event_time_ten_mhz_elap,
			      m_last_event_time_ten_mhz_live);

  // --------------------------------------------------------------------------
  // Do diagnostics, write data and release IED
  // --------------------------------------------------------------------------

  if(m_diagnostics)m_diagnostics->integrateEventIED(ied);
  if(ied->write_event)writeIED(ied);
  releaseIED(ied);
}

bool VBFAnalysisStage2::writeIED(const IED* ied)
{
  bool good;
  if(m_cuts_calc)good = m_acw->append(*ied->acd);
  else good = true;
  good &= m_edw->append(ied->ead); 
  return  good;
}


#ifndef NOTHREADS
void VBFAnalysisStage2::doIEDListWaitForSeqNumber(unsigned seq_number)
{
  vsassert(pthread_mutex_lock(&m_ied_list_mutex) == 0);

  std::list<unsigned>::iterator iseq_number = m_ied_wait_priority_list.begin();

  if(seq_number >= m_ied_wait_priority_next)
    {
      for(unsigned iadd=m_ied_wait_priority_next;
	  iadd<=seq_number; iadd++)
	{
	  while((iseq_number!=m_ied_wait_priority_list.end())
		&&(iadd>(*iseq_number)))iseq_number++;
	  m_ied_wait_priority_list.insert(iseq_number,iadd);
	}
      m_ied_wait_priority_next = seq_number+1;
    }

#ifdef TESTING
  std::cout << "W0: " << seq_number << " T" << pthread_self()
	    << ' ' << m_ied_wait_priority_list.size();
  for(iseq_number = m_ied_wait_priority_list.begin();
      iseq_number!=m_ied_wait_priority_list.end();iseq_number++)
    std::cout << ' ' << *iseq_number;		
  std::cout << std::endl;
#endif

  m_ied_free_stat.accumulate(m_ied_free.size());

  bool waited = false;

  if((m_ied_free.size()<=m_nthreads)
     ||(!m_ied_wait_priority_list.empty()))
    {
      while(1)
	{
	  unsigned nfree = m_ied_free.size();
#ifdef TESTING
	  std::cout << "NFREE: " << nfree << std::endl;
#endif		
	  if((nfree)||(m_ied_created!=m_ied_maximum))
	    {
	      unsigned ipriority = 0;
	      for(iseq_number=m_ied_wait_priority_list.begin();
		  *iseq_number!=seq_number;iseq_number++)ipriority++;

#ifdef TESTING
	      std::cout << "W1: " << seq_number << ' ' 
			<< ipriority << " T" << pthread_self() 
			<< std::endl;
#endif

	      if(ipriority<nfree)break;
	      unsigned ialloc = allocateIED(ipriority-nfree+1);
	      if(ialloc)
		{
		  vsassert(pthread_cond_broadcast(&m_ied_list_wait_cond) == 0);
		  nfree += ialloc;
		  if(ipriority<nfree)break;
		}
	    }

	  waited = true;
	  vsassert(pthread_cond_wait(&m_ied_list_wait_cond,
				   &m_ied_list_mutex)==0);
	}
    }

  if(waited)m_ied_wait_count++;

  for(iseq_number = m_ied_wait_priority_list.begin();
      iseq_number!=m_ied_wait_priority_list.end();iseq_number++)
    if(*iseq_number == seq_number)
      {
	m_ied_wait_priority_list.erase(iseq_number);
	break;
      }

#ifdef TESTING
  std::cout << "W2: " << seq_number << " T" << pthread_self()
	    << ' ' << m_ied_wait_priority_list.size();
  for(iseq_number = m_ied_wait_priority_list.begin();
      iseq_number!=m_ied_wait_priority_list.end();iseq_number++)
    std::cout << ' ' << *iseq_number;		
  std::cout << std::endl;
#endif
}
#endif

