#include <VSPixelCleaningVisitor.hpp>
#include <VSFileUtility.hpp>
#include <fast_alloc.hpp>

using namespace VERITAS;

VSPixelCleaningVisitor::Options::Options():
  integration_lo_zero_sample(4,6),
  integration_hi_zero_sample(4,0),
  integration_threshold_frac(0.5),
  integration_window_start(0.5),
  integration_window_width(5),
  integration_apply_increase(true),
  integration_threshold_charge(100.0),
  integration_window_start_increase(0.5),
  integration_window_width_increase(3),
  cleaning_args(),
  cleaning_scale("pedrms"),
  camera("veritas/499"),
  include_trigger_channels(true),
  include_logain_channels(true),
  output_file("")
{
  cleaning_args.push_back("regional");
  cleaning_args.push_back("2");
  cleaning_args.push_back("3");
  cleaning_args.push_back("4");
}

VSPixelCleaningVisitor::Options VSPixelCleaningVisitor::s_default_options;

VSPixelCleaningVisitor::
VSPixelCleaningVisitor(const VSAnalysisStage1Data& stage1): 
  VSSimpleVBFVisitor(), m_stream(std::cout), 
  m_stage1(&stage1), m_sup(), m_ped(&stage1.pedestals), 
  m_qcalc(), m_clean(), m_channel_map(),
  m_iscope(), m_scope_nchan_total(), m_scope_nchan_hit(), m_nchan(),
  m_chan_has_signal(), m_chan_hit(), m_chan_trigger(), m_chan_logain(),
  m_chan_sij(), m_chan_sij_total(), m_chan_tij(),
  m_chan_clean_multiplier(),
  m_nsup_stat(),
  m_scope_nsup_stat(), m_scope_ninc_stat(),
  m_event(), m_array_event(), m_packet(), m_bank_writer(),
  m_options(s_default_options)
{ 
  // --------------------------------------------------------------------------
  // Set up the camera definition
  // --------------------------------------------------------------------------
  unsigned camera_nchan = 0;
  const float* camera_xcoord = 0;
  const float* camera_ycoord = 0;

  if(m_options.camera == "veritas/499")
    {
      camera_nchan = sizeof(VC499GroundYcoord)/sizeof(*VC499GroundYcoord);
      camera_xcoord = VC499GroundXcoord;
      camera_ycoord = VC499GroundYcoord;
      m_camera_neighbors = VC499Neighbors;
    }
  else
    {
      std::cerr << "Unknown camera \"" << m_options.camera
		<< "\" .. aborting" << std::endl;
      exit(EXIT_FAILURE);
    }

  unsigned max_nchan = 0;
  for(unsigned iscope=0;iscope<stage1.run_info.nchan.size();iscope++)
    if(stage1.run_info.nchan[iscope]>max_nchan)
      max_nchan=stage1.run_info.nchan[iscope];

  m_chan_has_signal.resize(max_nchan);
  m_chan_hit.resize(max_nchan);
  m_chan_trigger.resize(max_nchan);
  m_chan_logain.resize(max_nchan);

  m_chan_sij        = new double[max_nchan];
  m_chan_sij_total  = new double[max_nchan];
  m_chan_tij        = new double[max_nchan];

  m_chan_clean_multiplier.resize(max_nchan);

  m_nchan = stage1.run_info.nchan;
  m_scope_nsup_stat.resize(m_nchan.size());
  m_scope_ninc_stat.resize(m_nchan.size());

  m_clean = 
    VSCleanerFactory::getCleaner(m_options.cleaning_args,
				 camera_nchan, m_camera_neighbors);

  VSPixelCleaningVisitor::CleaningScale cleaning_scale;
  if((m_options.cleaning_scale == "pedrms")
     ||(m_options.cleaning_scale == "pedvar"))
    cleaning_scale=VSPixelCleaningVisitor::CS_PEDRMS;
  else if(m_options.cleaning_scale == "gain")
    cleaning_scale=VSPixelCleaningVisitor::CS_GAIN;
  else if(m_options.cleaning_scale == "unity")
    cleaning_scale=VSPixelCleaningVisitor::CS_UNITY;
  else
    {
      std::cerr << "Unknown cleaning scaling: " << m_options.cleaning_scale
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  m_qcalc = 
    new VSConstantFractionCalc(m_options.integration_threshold_frac,
			       m_options.integration_window_start,
			       m_options.integration_window_width,
			       0, 0, 
			       m_options.integration_threshold_charge,
			       m_options.integration_window_start_increase,
			       m_options.integration_window_width_increase);

  m_channel_map = new VSChannelMap(m_stage1->run_info.lo_event_time);

  const unsigned nscope = m_nchan.size();
  m_chan_l2.resize(nscope);
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      m_chan_l2[iscope].resize(max_nchan);
      for(unsigned ichan=0;ichan<max_nchan;ichan++)
	{
	  const unsigned icrate = 
	    m_channel_map->crateForChannel(iscope,ichan);
	  const unsigned il2 = 
	    m_channel_map->l2ChannelForCrate(iscope, icrate);

	  if(il2 == ichan) m_chan_l2[iscope][ichan] = true;
	}
    }

  // Create the bank writer ---------------------------------------------------
  std::string output_file = m_options.output_file;

  if(output_file.empty())output_file="?.cvbf";

  if(output_file.find('?') != output_file.npos)
    {
      std::string::size_type iquestion = output_file.find('?');
      output_file.
	replace(iquestion,1,
		VSDataConverter::toString(m_stage1->run_info.run_number));
    }

  if(VSFileUtility::exists(output_file))
    {
      std::cerr << "Error: File " << output_file << " already exists."
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  m_bank_writer = new VBankFileWriter(output_file, 
				      m_stage1->run_info.run_number, 
				      m_stage1->run_info.config_mask);
}

VSPixelCleaningVisitor::~VSPixelCleaningVisitor()
{
  delete[] m_chan_sij;
  delete[] m_chan_sij_total;
  delete[] m_chan_tij;

  delete m_clean;
  delete m_qcalc;
  delete m_channel_map;

  m_bank_writer->finish();

  delete m_bank_writer;


  std::cout 
    << "SUMMARY -------------------------------------------------------------" 
    << std::endl << std::endl;

//   std::cout << "TOTAL AVG SUP: " 
// 	    << m_nsup_stat.mean() << " "
// 	    << m_nsup_stat.dev() << std::endl;

  std::cout << std::setw(10) << "SCOPE"
	    << std::setw(20) << "AVG PIX INCLUDED"
	    << std::setw(20) << "AVG PIX SUPPRESSED"
	    << std::endl;

  for(unsigned iscope = 0; iscope < m_scope_nsup_stat.size(); iscope++)
    {
      std::cout << std::setw(10) << iscope 
		<< std::setw(20) << m_scope_ninc_stat[iscope].mean() 
		<< std::setw(20) << m_scope_nsup_stat[iscope].mean() 
		<< std::endl;
    }

}

void VSPixelCleaningVisitor::visitPacket(bool& veto_packet, 
					 void*& user_data,
					 uint32_t  seq_packet_number,
					 uint32_t  vbf_packet_number,
					 bool      has_array_event,
					 bool      has_sim_header,
					 bool      has_sim_event,
					 uint32_t  num_overflow_datum,
					 const VBFPacket* packet)
{
  m_packet = new VPacket;
  m_packet_number = vbf_packet_number;
}

void VSPixelCleaningVisitor::leavePacket(bool veto_packet, void* user_data)
{
  if(m_array_event)
    {
      m_packet->putArrayEvent(m_array_event);
      m_array_event = NULL;
    }

  m_bank_writer->writePacket(m_packet_number,m_packet);
  delete m_packet;
  m_packet = NULL;
}

void VSPixelCleaningVisitor::
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
  if(event_num % 10000 == 0) std::cout << event_num << std::endl;

  m_event_type = event_type;
  m_array_event = new VArrayEvent;
}

void VSPixelCleaningVisitor::
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
  m_array_event->setTrigger(trigger->copyAT());

//   std::cout << cal_count[0] << " " << cal_count[1] << std::endl;

//   std::cout << m_array_event->getTrigger()->getOptCalCountArray()[0] << " "
// 	    << m_array_event->getTrigger()->getOptCalCountArray()[1] 
// 	    << std::endl;
}

void VSPixelCleaningVisitor::
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
  m_iscope = telescope_num;
  m_scope_nchan_total = num_channels_total;
  m_scope_nchan_hit = num_channels_saved;

  std::fill(m_chan_has_signal.begin(),m_chan_has_signal.end(),false);
  std::fill(m_chan_hit.begin(),m_chan_hit.end(),false);
  std::fill(m_chan_trigger.begin(),m_chan_trigger.end(),false);
  std::fill(m_chan_logain.begin(),m_chan_logain.end(),false);

//   std::cout << std::setw(10) << event_num
// 	    << std::setw(10) << telescope_num
// 	    << std::setw(15) << num_channels_saved
// 	    << std::setw(15) << num_channels_total
// 	    << std::setw(15) << m_nchan[telescope_num]
// 	    << std::endl;

  m_event = event;
}

void VSPixelCleaningVisitor::leaveScopeEvent(bool veto_scope_event, 
					     void* user_data)
{
  VEvent *new_event = m_event->copyEvent();
  new_event->setCompressedBit(true);

  if(m_event_type == ET_PED)
    {
      m_array_event->addEvent(new_event);
      return;
    }

  const unsigned iscope = m_iscope;
  const unsigned nchan = m_scope_nchan_total;

  // --------------------------------------------------------------------------
  // Set up for cleaning
  // --------------------------------------------------------------------------
  
  double* level_ij = FASTCALLOC(double,nchan);
  unsigned* masked_ij = FASTCALLOC(unsigned,nchan);

  //  bool l2_correction_enabled = m_settings.enable_l2_corrections;

  for(unsigned ichan=0;ichan<nchan;ichan++)
    {
      if(!m_chan_has_signal[ichan] || m_chan_l2[m_iscope][ichan])
	{
	  masked_ij[ichan] = true;
	  continue;
	}

      double ped_dev = m_ped->dev(m_iscope, ichan, 
				  m_options.integration_window_width);

      if(ped_dev > 0)
	{
	  level_ij[ichan] = m_chan_sij[ichan]/ped_dev;
	  masked_ij[ichan] = false;
	}
      else
	{	  
	  level_ij[ichan] = 0;
	  masked_ij[ichan] = true;	  
	}
    }
  
  // --------------------------------------------------------------------------
  // Clean image
  // --------------------------------------------------------------------------
  std::vector<bool> suppress(nchan);
  unsigned nsuppressed = 0;
  unsigned nincluded = 0;

  const unsigned nimage = m_clean->clean(nchan, masked_ij, level_ij);

  for(unsigned ichan=0;ichan<nchan;ichan++)
    {



      if(m_chan_l2[m_iscope][ichan] ||
	 (m_chan_trigger[ichan] && m_options.include_trigger_channels) ||
	 (m_chan_logain[ichan] && m_options.include_logain_channels))
	{
	  nincluded++;
	  suppress[ichan] = false;	  
	}
      else if(masked_ij[ichan] == true)
	{
	  nsuppressed++;
	  suppress[ichan] = true;
	}
      else
	{
	  nincluded++;
	  suppress[ichan] = false;	  
	}

//       if(level_ij[ichan] < 0 && masked_ij[ichan] == 0)
// 	{
// 	  std::cout << ichan << " " << level_ij[ichan]
// 		    << " " << masked_ij[ichan] << " ----------- "
// 		    << std::endl;
// 	  for(unsigned jneighbor=0;jneighbor<NUM_NEIGHBORS;jneighbor++)
//             {
// 	      int jchan = m_camera_neighbors[ichan][jneighbor];
// 	      if(jchan == -1) continue;
// 	      std::cout 
// 		<< jneighbor << " " << jchan << " "
// 		<< masked_ij[jchan] << " " << level_ij[jchan] << " "
// 		<< std::endl;
// 	    }
// 	}

    }

  new_event->resizeChannelData(new_event->getNumSamples(), nchan-nsuppressed);
  for (unsigned k=0,l=0,m=0;k<new_event->getMaxNumChannels();++k) 
    {
      if(m_event->getHitBit(k)) 
	{
	  if(!suppress[k]) 
	    {
	      new_event->setCharge(m,m_event->getCharge(l));
	      new_event->setPedestalAndHiLo(m,m_event->getPedestalAndHiLo(l));
	      memcpy(new_event->getSamplePtr(m,0),
		     m_event->getSamplePtr(l,0),
		     m_event->getNumSamples());
	      ++m;
	    }
	  ++l;
	}
      new_event->setHitBit(k,!suppress[k]);
    }

  m_scope_nsup_stat[m_iscope].accumulate(nsuppressed);
  m_scope_ninc_stat[m_iscope].accumulate(nincluded);
  m_nsup_stat.accumulate(nsuppressed);

  FASTFREE(level_ij);
  FASTFREE(masked_ij);

  m_array_event->addEvent(new_event);
}

void VSPixelCleaningVisitor::visitChannel(bool& veto_channel, 
					  void* user_data,
					  uint32_t channel_num, 
					  bool hit, 
					  bool trigger)
{
  m_chan_has_signal[channel_num] = hit;
  m_chan_hit[channel_num] = hit;
  m_chan_trigger[channel_num] = trigger;
}

void VSPixelCleaningVisitor::visitHitChannel(void* user_data,
				      uint32_t               channel_num,
				      uint32_t               charge, 
				      uint32_t               pedestal,
				      bool                   lo_gain,
				      unsigned               nsample,
				      const uint32_t*        samples,
				      const uint32_t*        integrated)
{
  m_chan_logain[channel_num] = lo_gain;

  unsigned sample_zero = 0;
  if(lo_gain) sample_zero = m_options.integration_lo_zero_sample[m_iscope];
  else sample_zero = m_options.integration_hi_zero_sample[m_iscope];

  double ped;
  if(lo_gain)
    {
      ped = m_ped->loPed(m_iscope,channel_num);
// if(scal.ch_suppress[channel_num]||scal.ch_suppress_lo_gain[channel_num])
//  ied->chan_has_signal[channel_num] = false;
    }
  else
    {
      ped = m_ped->hiPed(m_iscope,channel_num);
//       if(scal.ch_suppress[channel_num])
// 	ied->chan_has_signal[channel_num] = false;
    }

  const double lo_gain_mult = 6;

  double sij;
  double tij;
  unsigned window_width;
  double window_start;
  double total_signal;
  const uint32_t* peak_ptr;

  m_qcalc->calc(sij, tij, window_width, window_start, total_signal, peak_ptr,
		lo_gain, lo_gain_mult, nsample, samples, integrated, 
		ped, sample_zero);

  m_chan_sij[channel_num]       = sij;
  m_chan_sij_total[channel_num] = total_signal;
  m_chan_tij[channel_num]       = tij;
}

void VSPixelCleaningVisitor::
visitSimulationHeader(void* user_data,
		      uint32_t         run_number,
		      const VSTime&    date_of_sims,
		      uint32_t         simulation_package,
		      const std::string& simulator_name,
		      const VSTime&    date_of_array,
		      uint32_t         corsika_atm_model,
		      float            obs_altitude_m,
		      const std::vector<VArrayConfiguration>&
		      array,
		      const std::string& stringified_config)
{
  VSimulationHeader* sh = new VSimulationHeader(date_of_sims,
						simulation_package,
						simulator_name,
						date_of_array,
						corsika_atm_model,
						obs_altitude_m,
						array,
						stringified_config);

  m_packet->putSimulationHeader(sh);
}

void VSPixelCleaningVisitor::
visitSimulationEvent(bool& veto_packet, 
		     void* user_data,
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
  VSimulationData* sd = new VSimulationData(corsika_particle_id,
					    energy_gev,
					    obs_zenith_deg,
					    obs_azimuth_deg,
					    primary_zenith_deg,
					    primary_azimuth_deg,
					    ref_zenith_deg,
					    ref_azimuth_deg,
					    ref_position_angle_deg,
					    core_east_m,
					    core_south_m,
					    core_elevation_asl_m);
  
  m_packet->putSimulationData(sd);
}

void VSPixelCleaningVisitor::
visitChiLAHeader(void* user_data,
		 const std::string& database_name,
		 const std::vector<VChiLASimParamTable>& sim_param_tables,
		 const std::vector<VChiLAOpticsConfiguration>&
		 optics_configurations,
		 const std::vector<VChiLAElectronicsConfiguration>& 
		 electronics_configurations)
{
  VChiLASimulationHeader* csh = 
    new VChiLASimulationHeader(database_name,
			       sim_param_tables,
			       optics_configurations,
			       electronics_configurations);

  m_packet->put(VGetChiLASimulationHeaderBankName(),csh);
}

void VSPixelCleaningVisitor::
visitChiLAEvent(bool& veto_packet, void* user_data,
		uint32_t               table_index,
		uint32_t               electronics_id,
		uint32_t               table_event_index,
		float                  a_omega,
		float                  a_omega_var)
{
  VChiLASimulationData* csd = 
    new VChiLASimulationData(table_index,
			     electronics_id,
			     table_event_index,
			     a_omega, a_omega_var);

  m_packet->put(VGetChiLASimulationDataBankName(),csd);
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSPixelCleaningVisitor::configure(VSOptions& options,
				       const std::string& profile,
				       const std::string& opt_prefix)
{

  options.findWithValue(OPTNAME(opt_prefix,"cleaning_scale"),
			s_default_options.cleaning_scale,
			"Set the scale of the cleaning algorithm. Valid "
			"settings are \"pedrms\", \"gain\", and \"unity\".");

  options.findWithValue(OPTNAME(opt_prefix,"cleaning"),
			s_default_options.cleaning_args,
			"Configure the primary cleaning algorithm. Valid "
			"forms are \"picbnd,piclevel,bndlevel\" or "
			"\"regional,rlevel,rsize,ilevel\".");

  
  options.findWithValue(OPTNAME(opt_prefix,"hi_gain_zero_sample"),
			s_default_options.integration_hi_zero_sample,
			"Set the minimum sample number on each telescope from "
			"which HIGH gain channels can be integrated. This can "
			"usually be set at ZERO.");

  options.findWithValue(OPTNAME(opt_prefix,"lo_gain_zero_sample"),
			s_default_options.integration_lo_zero_sample,
			"Set the minimum sample number on each telescope from "
			"from which LOW gain channels can be integrated. This "
			"should usually be set to a value which excludes the "
			"remnants of the HIGH gain signal which may remain at "
			"the edge of the integration window.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_start"),
			s_default_options.integration_window_start,
			"Sample number to integrate from. The window start "
			"is defined relative to T-zero, the point at which "
			"the trace reaches some fraction of the peak height. "
			"Positive values means the integration starts before "
			"the T-zero point.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_width"),
			s_default_options.integration_window_width,
			"Number of samples to integrate when calculating "
			"charge in each channel.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"peak_threshold_fraction"),
			s_default_options.integration_threshold_frac,
			"Fraction of peak of trace at which T0 is defined.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"apply_window_increase"),
			s_default_options.integration_apply_increase,
			"Apply increase in window size for large pulses.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_increase_threshold"),
			s_default_options.integration_threshold_charge,
			"Total charge in channel above which the integration "
			"window is increased. The window size is increased "
			"and shifted earlier in the trace logarithmically "
			"with the charge above this threshold.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_start_increase"),
			s_default_options.integration_window_start_increase,
			"Amount of advance to apply to the start of the "
			"integration window (in samples) per decade that the "
			"total charge exceeds the threshold.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_width_increase"),
			s_default_options.integration_window_width_increase,
			"Amount to increase the width of the integration "
			"window (in samples) per decade that the total charge "
			"exceeds the threshold.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"include_trigger_channels"),
			s_default_options.include_trigger_channels,
			"Always include triggered channels.");

  options.findWithValue(OPTNAME(opt_prefix,"include_logain_channels"),
			s_default_options.include_logain_channels,
			"Always include lo gain channels.");

  options.findWithValue(OPTNAME(opt_prefix,"o"),
			s_default_options.output_file,
			"Set the output file name.  If the output file is "
			"empty a name defined by the run number will be used, "
			"e.g. 12345.cvbf.");
}
