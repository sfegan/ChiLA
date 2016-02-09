#include <VSRVBFVisitor.hpp>
#include <fast_alloc.hpp>

#include <VSRCanvas.hpp>

using namespace VERITAS;

VSRVBFVisitor::Options::Options():
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

VSRVBFVisitor::Options VSRVBFVisitor::s_default_options;

VSRVBFVisitor::VSRVBFVisitor(const VSAnalysisStage1Data& stage1):
  VSSimpleVBFVisitor(), m_stream(std::cout),
  m_stage1(&stage1), m_sup(), m_ped(&stage1.pedestals), 
  m_qcalc(), m_clean(), m_channel_map(),
  m_iscope(), m_scope_nchan_total(), m_scope_nchan_hit(), m_nchan(),
  m_chan_sij(), m_chan_sij_total(), m_chan_tij(),
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

  m_clean = 
    VSCleanerFactory::getCleaner(m_options.cleaning_args,
				 camera_nchan, m_camera_neighbors);

  VSRVBFVisitor::CleaningScale cleaning_scale;
  if((m_options.cleaning_scale == "pedrms")
     ||(m_options.cleaning_scale == "pedvar"))
    cleaning_scale=VSRVBFVisitor::CS_PEDRMS;
  else if(m_options.cleaning_scale == "gain")
    cleaning_scale=VSRVBFVisitor::CS_GAIN;
  else if(m_options.cleaning_scale == "unity")
    cleaning_scale=VSRVBFVisitor::CS_UNITY;
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
}

VSRVBFVisitor::~VSRVBFVisitor()
{
  // nothing to see here
}

void VSRVBFVisitor::visitPacket(bool& veto_packet, 
				  void*& user_data,
				  uint32_t  seq_packet_number,
				  uint32_t  vbf_packet_number,
				  bool      has_array_event,
				  bool      has_sim_header,
				  bool      has_sim_event,
				  uint32_t  num_overflow_datum,
				  const VBFPacket* packet)
{
  std::cout << "PACKET " << vbf_packet_number << std::endl;
}

void VSRVBFVisitor::
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
  m_event_num = event_num;
  m_iscope = telescope_num;
  m_scope_nchan_total = num_channels_total;
  m_scope_nchan_hit = num_channels_saved;
  //  m_scope_event_type = event_type.trigger;

  
}

void VSRVBFVisitor::leaveScopeEvent(bool veto_scope_event, 
					     void* user_data)
{
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

  //  const unsigned nimage = m_clean->clean(nchan, masked_ij, level_ij);
  
  VSRHistogramCamera* hist = new VSRHistogramCamera;

  for(unsigned i = 0; i < m_scope_nchan_total; i++)
    {
      if(m_chan_hit[i] && !m_chan_l2[m_iscope][i]) 
	hist->setPixel(i,m_chan_sij[i]);
    }

  FASTFREE(level_ij);
  FASTFREE(masked_ij);

  m_scope_hist.push_back(hist);
}

void VSRVBFVisitor::visitChannel(bool& veto_channel, 
					  void* user_data,
					  uint32_t channel_num, 
					  bool hit, 
					  bool trigger)
{
  m_chan_has_signal[channel_num] = hit;
  m_chan_hit[channel_num] = hit;
  m_chan_trigger[channel_num] = trigger;
}

void VSRVBFVisitor::visitHitChannel(void* user_data,
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


  //  std::cout << channel_num << " " << sij << std::endl;

}

void VSRVBFVisitor::
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
  std::cout << m_event_num << " " << energy_gev/1000. << std::endl;
}

void VSRVBFVisitor::draw()
{
  VSRCanvas* canvas = new VSRCanvas("event");

  canvas->setDimensions(1000,1000);

  for(unsigned iscope = 0; iscope < m_scope_hist.size(); iscope++)
    canvas->add(iscope,m_scope_hist[iscope]);

  canvas->draw();
}
