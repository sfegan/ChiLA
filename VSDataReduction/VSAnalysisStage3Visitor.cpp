//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage3Visitor.cpp
  

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       12/17/2009
*/

#include <SphericalCoords.h>

#include <VSAnalysisStage3Visitor.hpp>
#include <VSStage3SimCalc.hpp>

using namespace VERITAS;
using namespace SEphem;

// ============================================================================
// VSAnalysisStage3Visitor
// ============================================================================

VSAnalysisStage3Visitor::
VSAnalysisStage3Visitor(VSCutsEvaluator* cuts_evaluator,
			VSOctaveH5WriterStruct* writer,
			double sky_bin_width_deg,
			double offset_max,
			const std::string& coord_system,
			const std::string& src_name,
			const SEphem::SphericalCoords& origin_radec,
			const SEphem::SphericalCoords& src_radec,
			const std::vector< SEphem::SphericalCoords >& 
			ptg_radec,
			const VSExclusionRegion& exclusion_region,
			const std::map< unsigned, Run >& run_list,
			bool no_run_results,
			bool write_run_hists,
			unsigned rng_seed):
  VSEventDataVisitor(),
  m_integral_analysis(),
  m_cuts_evaluator(cuts_evaluator),
  m_writer(writer),
  m_sp_calc(), m_egy_calc(),
  m_egywt_calc(),
  m_source_injector(),
  m_spectrum_calc(),
  m_rng(),
  m_data(), m_run_data(),
  m_sim_data(),
  m_sim_run_data(),
  m_sim_table_data(),
  m_event_data(),
  m_sim(),
  m_nchan(),
  m_sky_bin_width_deg(sky_bin_width_deg),
  m_offset_max(offset_max),
  m_theta_cut(),
  m_spectrum_theta_cut(),
  m_ring_cut(),
  m_source_spectrum(),
  m_coord_system(coord_system),
  m_src_name(src_name),
  m_origin_radec(origin_radec),
  m_src_radec(src_radec),
  m_ptg_radec(ptg_radec),
  m_run_list(run_list),
  m_no_run_results(no_run_results),
  m_write_run_hists(write_run_hists),
  m_rng_seed(rng_seed),
  m_log10_egy_min(), m_log10_egy_max(), m_egy_bin_width(),
  m_has_sim(false), m_event_code()
{
  m_integral_analysis = 
    VSIntegralAnalysisFactory::getInstance()->create(offset_max,
						     sky_bin_width_deg);

  m_integral_analysis->setExclusionRegion(exclusion_region);
  m_theta_cut = m_integral_analysis->thetaCut();
  m_ring_cut = m_integral_analysis->ringCut();
  m_source_spectrum = m_integral_analysis->spmodel();

  m_integral_analysis->setSourcePosition(m_origin_radec,m_src_radec);

  m_rng = new RandomNumbers(rng_seed);
  m_egywt_calc = new VSSimEnergyWeightCalc;
  m_source_injector = new VSSourceInjector(m_rng,m_theta_cut,m_ring_cut);

  VSSpectrumCalcFactory* spcf = VSSpectrumCalcFactory::getInstance();

  m_spectrum_calc = spcf->create();
  m_egy_bin_width = spcf->egyBinWidth();
  m_log10_egy_min = spcf->egyMin();
  m_log10_egy_max = spcf->egyMax();
  m_spectrum_theta_cut = spcf->thetaCut();
}

VSAnalysisStage3Visitor::~VSAnalysisStage3Visitor()
{
  std::cout << std::string(79,'-') << std::endl;
  
  // --------------------------------------------------------------------------
  // Create the combined dataset
  // --------------------------------------------------------------------------
  VSAnalysisStage3Data data(m_src_name,m_src_radec,m_ptg_radec,
			    m_origin_radec);

  data.addData(m_run_data);

//   const unsigned nruns = m_run_data.size();
//   for(unsigned irun = 0; irun < nruns; irun++)
//     data.addData(*m_run_data[irun]);

  VSSpectrumFitData fit_data;

  VSIntegralAnalysisData integral_data;

  VSAcceptanceData* acceptance_data = NULL;

  if(!m_has_sim)
    {
      acceptance_data = m_integral_analysis->fitAcceptance(data);
      m_integral_analysis->analyze(m_intanl_run_data,
				   data,*acceptance_data,integral_data);

      // Save the results to a file -------------------------------------------
      integral_data.save(m_writer->writeStruct("results"));

      m_spectrum_calc->reconstruct(data,m_spectrum_data);
      m_spectrum_calc->fit(data,fit_data);
    }

  if(!m_has_sim && !m_no_run_results)
    {
      std::vector< VSIntegralAnalysisRunDatum > run_results;
      std::vector< VSIntegralAnalysisRunDatum > cumulative_run_results;

      VSAnalysisStage3Data cumulative_run_data(m_src_name,m_src_radec,
					       m_origin_radec);

      std::vector<VSIntegralAnalysis::Data> intanl_cum_data;

      // Calculate results for each run individually --------------------------
      const unsigned nrun = m_run_data.size();
      for(unsigned irun = 0; irun < nrun; irun++)
	{
	  VSIntegralAnalysisRunDatum run_datum(m_run_data[irun]->run_number);
	  VSIntegralAnalysisRunDatum 
	    cumulative_run_datum(m_run_data[irun]->run_number);

	  std::vector<VSIntegralAnalysis::Data> intanl_data;
	  intanl_cum_data.push_back(m_intanl_run_data[irun]);
	  intanl_data.push_back(m_intanl_run_data[irun]);

	  run_datum.run_start_time_string = 
	    m_run_data[irun]->lo_event_time.getString();
	  run_datum.run_stop_time_string = 
	    m_run_data[irun]->hi_event_time.getString();
	  run_datum.run_start_time_mjd = 
	    m_run_data[irun]->lo_event_time.getMJDDbl();
	  run_datum.run_stop_time_mjd = 
	    m_run_data[irun]->hi_event_time.getMJDDbl();
	  
	  VSAnalysisStage3Data run_data(m_src_name,m_src_radec,m_origin_radec);
	  run_data.addData(*m_run_data[irun]);
	  cumulative_run_data.addData(*m_run_data[irun]);

	  std::cout << "Analyzing Run " << m_run_data[irun]->run_number
		    << std::endl;
	  
	  m_integral_analysis->analyze(intanl_data,
				       run_data,*acceptance_data,run_datum);

	  std::cout << "Computing cumulative results... " << std::endl;

	  m_integral_analysis->analyze(intanl_cum_data,
				       cumulative_run_data,*acceptance_data,
				       cumulative_run_datum);


	  run_results.push_back(run_datum);
	  cumulative_run_results.push_back(cumulative_run_datum);
	}

      m_writer->writeStructCellVector("run_results",run_results);
      m_writer->writeStructCellVector("cumulative_run_results",
				      cumulative_run_results);
    }

  data.save(m_writer->writeStruct("stage3"),m_write_run_hists);
  
  m_spectrum_data.save(m_writer->writeStruct("spectrum"));
  fit_data.save(m_writer->writeStruct("spectrum_fit"));

  // Write the cuts -----------------------------------------------------------
  m_cuts_evaluator->save(m_writer);

  // Write fake source data ---------------------------------------------------
  m_source_injector->getData(m_egy_bin_width,
			     m_log10_egy_min,
			     m_log10_egy_max).
    save(m_writer->writeStruct("fake_src"));
  

  // Cleanup ------------------------------------------------------------------
  const unsigned nrun = m_run_data.size();
  for(unsigned irun = 0; irun < nrun; irun++)
    delete m_run_data[irun];
  for(unsigned irun = 0; irun < nrun; irun++)
    delete m_sim_run_data[irun];

  delete m_egywt_calc;
  delete m_source_injector;
  delete m_spectrum_calc;
  delete m_integral_analysis;
}

void VSAnalysisStage3Visitor::visitRun(const VSAnalysisStage1Data& stage1,
				       const VSTargetTable::Observation& obs,
				       const VSArrayMergedCalibrationData& cal)
{
  m_nchan = stage1.run_info.nchan;

  if(m_sim_data) m_sim_data->setNChan(m_nchan);

  // --------------------------------------------------------------------------
  // Create stage3 data structure
  // --------------------------------------------------------------------------
  m_data = 
    new VSAnalysisStage3Data::RunData(m_offset_max,m_sky_bin_width_deg,
				      m_egy_bin_width, m_log10_egy_min,
				      m_log10_egy_max,
				      m_src_radec,
				      m_origin_radec, obs);


  m_data->run_number = stage1.run_info.run_number;
  m_data->src_name = obs.name;
  m_data->mode = obs.mode;
  m_data->wobble_phi_deg = obs.wobble_phi_rad*180./M_PI;
  m_data->wobble_theta_deg = obs.wobble_theta_rad*180./M_PI;

  m_data->zn_mean_deg = stage1.run_info.zn_mean_deg;
  m_data->zn_rms_deg = stage1.run_info.zn_rms_deg;
  m_data->az_mean_deg = stage1.run_info.az_mean_deg;
  m_data->az_rms_deg = stage1.run_info.az_rms_deg;
  m_data->mean_scaled_dev = cal.mean_scaled_dev;

  std::map< unsigned, Run >::iterator itr = 
    m_run_list.find(m_data->run_number);
  vsassert(itr != m_run_list.end());

  m_data->setPointing(itr->second.iptg,itr->second.ptg_radec);

  m_intanl_data = 
    m_integral_analysis->create(m_data->ptg_xy(),m_data->obs_xy());


  // Calculate parallactic angle ----------------------------------------------
  m_data->lo_event_time = stage1.run_info.lo_event_time;
  m_data->lo_event_time_string = stage1.run_info.lo_event_time.getString();
  m_data->hi_event_time = stage1.run_info.hi_event_time;
  m_data->hi_event_time_string = stage1.run_info.hi_event_time.getString();

  calcParallacticAngle(m_data->lo_event_time,m_data->hi_event_time,
		       m_data->obs_radec(), m_data->pangle_mean_deg,
		       m_data->pangle_rms_deg);

  // --------------------------------------------------------------------------
  // Load Scaled Parameter and Energy lookup tables
  // --------------------------------------------------------------------------
  m_data->sp_calc()->load(m_data->lo_event_time,
			  stage1.run_info.nchan,stage1.run_info.zn_mean_deg,
			  stage1.run_info.az_mean_deg,cal.mean_scaled_dev);

  m_data->egy_calc()->load(m_data->lo_event_time,
			   stage1.run_info.nchan,stage1.run_info.zn_mean_deg,
			   stage1.run_info.az_mean_deg,cal.mean_scaled_dev);

  m_data->irf_calc().load(m_data->lo_event_time,
			  stage1.run_info.nchan,stage1.run_info.zn_mean_deg,
			  stage1.run_info.az_mean_deg,cal.mean_scaled_dev);

  m_cuts_evaluator->setDate(stage1.run_info.lo_event_time);

  // --------------------------------------------------------------------------
  // Load Effective Area, PSF Model, and Energy Kernel
  // --------------------------------------------------------------------------

  m_integral_analysis->accumulate(*m_data,m_intanl_data);

      
  m_data->irf_calc().
    getEffectiveAreaAperture(m_source_spectrum,       
			     m_theta_cut,
			     m_data->int_effarea_aperture_nspace);

  m_data->irf_calc().getEffectiveAreaOffset(m_source_spectrum,
					    m_data->int_effarea_nspace,
					    m_data->int_effarea_psf_nspace);

  double offset = m_data->src_xy().d(m_data->obs_xy());

  m_data->effarea_nspace = 
    m_data->irf_calc().getEffectiveArea(m_egy_bin_width,
					m_log10_egy_min,
					m_log10_egy_max,
					offset,m_spectrum_theta_cut);    

  // Get Off regions ----------------------------------------------------------
  VSOffRegion::getRegions(m_data->src_xy(), m_data->obs_xy(),
			  m_theta_cut, m_integral_analysis->maxOffRegions(),
			  m_data->off_xy());

  VSOffRegion::getRegions(m_data->src_xy(), m_data->obs_xy(),
			  m_spectrum_calc->thetaCut(),
			  m_spectrum_calc->maxOffRegions(),
			  m_data->sp_off_xy());

  if(!stage1.sim_info) return;

  // Do sim related tasks here ------------------------------------------------
  double egy_binsize = 0;
  double egy_min = 0;
  double egy_max = 0;

  for(std::vector<VSSimInfoData::ParticleSpectrum>::iterator itr =
	stage1.sim_info->particle_spectrum.begin(); itr != 
	stage1.sim_info->particle_spectrum.end(); ++itr)
    {
      if(itr->corsika_particle_id >= 1 && itr->dlog_energy != 0)
	{
	  egy_binsize = itr->dlog_energy;
	  egy_min = std::log10(itr->min_energy_tev) - egy_binsize/2.;
	  egy_max = std::log10(itr->max_energy_tev) + egy_binsize/2.;
	  break;
	}
    }

  vsassert(egy_binsize);

  m_sim_data->src_ra_deg = stage1.sim_info->src_ra_deg;
  m_sim_data->src_dec_deg = stage1.sim_info->src_dec_deg;
  m_sim_data->obs_ra_deg = stage1.sim_info->obs_ra_deg;
  m_sim_data->obs_dec_deg = stage1.sim_info->obs_dec_deg;
  m_sim_data->log10_egy_binsize = egy_binsize;
  m_sim_data->log10_egy_min = egy_min;
  m_sim_data->log10_egy_max = egy_max;

  typedef VSLimitedErrorsHist<double,double> Hist1D;
  typedef VSSimple2DHist<double,double> Hist2D;

  for(std::vector<VSStage3SimArrayTableDatum*>::iterator itr =
 	m_sim_table_data.begin(); itr != m_sim_table_data.end(); ++itr)
    (*itr)->setEnergyBinSize(egy_binsize, egy_min, egy_max);
  
  m_sim_data->setEnergyBinSize(egy_binsize, egy_min, egy_max);
}

void VSAnalysisStage3Visitor::leaveRun()
{
  m_data->livetime_min = m_data->livetime_ticks/(1E7*60.);
  m_data->elaptime_min = m_data->elaptime_ticks/(1E7*60.);

  m_intanl_data.livetime_min = m_data->livetime_min;
  m_intanl_data.elaptime_min = m_data->elaptime_min;

  for(VSSimple2DHist<double, double>::iterator itr = 
  	m_data->sky_livetime_hist.begin(); 
      itr != m_data->sky_livetime_hist.end(); ++itr)
    {
      VSAAlgebra::Vec2D xy(itr->x(),itr->y());

      if(m_data->ptg_xy().d(xy) < m_offset_max)
	{
	  m_data->sky_livetime_hist.setBin(itr->bin(),m_data->livetime_min);
	  //	  m_data->sky_exposure_hist.setBin(itr->bin(),1);
	}
    }

  // Injector the source ------------------------------------------------------
  m_source_injector->initialize(m_data);
  m_source_injector->fill(m_data);

  m_data->calcExcess();

  m_intanl_run_data.push_back(m_intanl_data);

  m_run_data.push_back(m_data);
  m_data = NULL;

  m_sim_run_data.push_back(m_sim_data);
  m_sim_data = NULL;
}

void VSAnalysisStage3Visitor::visitEvent(const VSEventArrayDatum& event)
{
  m_data->elaptime_ticks = event.elapsed_ticks;
  m_data->livetime_ticks = event.live_ticks;

  m_event_code |= VSAnalysisStage3Data::EC_RECONSTRUCTED;
  m_event_data = event;

  // Transform to XY coordinates ----------------------------------------------
  if(!std::isfinite(m_event_data.dec_J2000) || 
     !std::isfinite(m_event_data.ra_J2000))
    return;

  VSAAlgebra::Vec2D event_xy;
  SEphem::SphericalCoords event_radec = 
    SEphem::SphericalCoords::makeLatLongDeg(m_event_data.dec_J2000,
					    m_event_data.ra_J2000);

  if(m_coord_system == "radec")
    {
      event_radec.rotate(0,
			 -m_origin_radec.thetaRad(),
			 -m_origin_radec.phiRad());
      event_xy.set(-event_radec.thetaDeg()*sin(event_radec.phiRad()),
		   -event_radec.thetaDeg()*cos(event_radec.phiRad()));
    }
  else if(m_coord_system == "galactic")
    {
      SEphem::SphericalCoords origin_galactic = m_origin_radec;

      SEphem::Astro::raDecToGal(2000.0,origin_galactic);
      SEphem::Astro::raDecToGal(2000.0,event_radec);
      
      event_radec.rotate(0,-origin_galactic.thetaRad(),
			 -origin_galactic.phiRad());
      event_xy.set(-event_radec.thetaDeg()*sin(event_radec.phiRad()),
		   -event_radec.thetaDeg()*cos(event_radec.phiRad()));
    }
  else
    {
      std::cerr 
	<< std::string(__PRETTY_FUNCTION__) << ": "
	<< "Unrecognized coordinate system: " << m_coord_system
	<< std::endl;
      exit(EXIT_FAILURE);
    }

  // Recalculate energy estimator ---------------------------------------------
  m_data->egy_calc()->calcEnergy(m_event_data);
  
  // Recalculate scaled parameters --------------------------------------------
  m_data->sp_calc()->calcSP(m_event_data);
  
  if(m_cuts_evaluator->isSelected(m_event_data) && 
     m_data->ptg_xy().d(event_xy) < m_offset_max)
    m_event_code |= VSAnalysisStage3Data::EC_SELECTED;

  double wt = 1.0;
  if(m_has_sim) wt = m_egywt_calc->calcWeight(m_sim.energy_tev);

  // Do all accumulation here -------------------------------------------------
  if(m_event_code & VSAnalysisStage3Data::EC_SELECTED) 
    {
      m_spectrum_calc->accumulate(m_event_data,event_xy,*m_data);

      m_integral_analysis->accumulate(m_event_data,event_xy,
				      m_intanl_data);

      m_data->sky_counts_hist.accumulate(event_xy.x(),event_xy.y());
      m_data->fov_counts_hist.accumulate(m_event_data.mean_derotated_fov_x,
					 m_event_data.mean_derotated_fov_y);
      m_data->cam_counts_hist.accumulate(m_event_data.mean_fov_x,
					 m_event_data.mean_fov_y);
      m_data->fovr2_counts_hist.accumulate(std::pow(m_event_data.theta0,2));

      m_data->log10_N2_hist.accumulate(std::log10(m_event_data.N2));
      m_data->msc_width_hist.accumulate(m_event_data.msc_width);
      m_data->msc_length_hist.accumulate(m_event_data.msc_length);
      m_data->msc_disp_hist.accumulate(m_event_data.msc_disp);

      m_data->triggered_hist.accumulate(event.trigger_mask);
      m_data->has_image_hist.accumulate(event.has_image_mask);
      m_data->used_in_reconstruction_hist.
	accumulate(event.used_in_reconstruction_mask);

      m_data->accumulate(m_event_data,wt);

      // Accumulate on counts -------------------------------------------------
      double dr_src = m_data->src_xy().d(event_xy);

      if(dr_src < m_theta_cut)
	{
	  m_event_code |= VSAnalysisStage3Data::EC_ON;
	  m_data->on_counts++;
	  m_data->accumulateOn(m_event_data,wt);
	}

      // Accumulate ring counts -----------------------------------------------
      if(event_xy.d(m_data->src_xy()) > m_ring_cut.first &&
	 event_xy.d(m_data->src_xy()) < m_ring_cut.second)
	m_data->ring_counts++;
      
      m_data->th2_on_counts_hist.accumulate(std::pow(dr_src,2));

      // Accumulate off counts ------------------------------------------------
      for(std::vector< VSAAlgebra::Vec2D >::const_iterator itr =
	    m_data->off_xy().begin(); itr != m_data->off_xy().end(); ++itr)
	{
	  double dr_off = (*itr - event_xy).norm();
	  if(dr_off < m_theta_cut)
	    {
	      m_event_code |= VSAnalysisStage3Data::EC_OFF;
	      m_data->accumulateOff(m_event_data,wt);
	      m_data->th2_off_counts_hist.accumulate(std::pow(dr_off,2));
	      m_data->off_counts++;
	    }
	}
    }
}

void VSAnalysisStage3Visitor::leaveEvent()
{
  m_event_code = 0;
}

void VSAnalysisStage3Visitor::visitSimEvent(const VSArraySimulationDatum& sim)
{
  m_sim = sim;
  if(sim.has_array_event) m_event_code |= VSAnalysisStage3Data::EC_TRIGGERED;
  double wt = m_egywt_calc->calcWeight(sim.energy_tev);

  m_sim_data->accumulate(sim,m_event_data,m_event_code,wt);

  for(std::vector<VSStage3SimArrayTableDatum*>::iterator itr =
	m_sim_table_data.begin(); itr != m_sim_table_data.end(); ++itr)
    (*itr)->accumulate(sim,m_event_data,m_event_code,wt);
}

void VSAnalysisStage3Visitor::leaveSimEvent()
{

}

void VSAnalysisStage3Visitor::
visitSimHeader(const VSHeaderSimulationDatum& sim_header)
{
  m_sim_data = new VSStage3SimDatum;
  m_has_sim = true;

  m_egywt_calc->calcWeighting(sim_header);

  // Fill sampling area histogram ---------------------------------------------
  for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
	sim_header.tables.begin(); itr != sim_header.tables.end(); ++itr)
    {
      double loge = log10(itr->energy_tev);

      if(m_sim_data->egymc_sampling_area_hist.countForVal(loge) == 0)
	m_sim_data->egymc_sampling_area_hist.
	  accumulate(loge,M_PI*std::pow(itr->sampling_radius_m,2),0);

      m_sim_data->ntotal += itr->event_count;
      m_sim_data->egymc_count_hist.accumulate(loge,itr->event_count,
					    itr->event_count);

      double wt = m_egywt_calc->calcWeight(itr->energy_tev);

      m_sim_data->egymc_total_hist.
	accumulate(loge,wt*itr->event_count,wt*wt*itr->event_count);
    }
  
  for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
 	sim_header.tables.begin(); itr != sim_header.tables.end(); ++itr)
    {
      if(itr->table_index >= m_sim_table_data.size())
	m_sim_table_data.resize(itr->table_index + 1,NULL);

      if(m_sim_table_data[itr->table_index] == NULL)
	{
	  VSStage3SimArrayTableDatum* table_datum = 
	    new VSStage3SimArrayTableDatum(m_nchan);
	  table_datum->sampling_area = M_PI*std::pow(itr->sampling_radius_m,2);
	  table_datum->table_index = itr->table_index;
	  table_datum->egy_tev = itr->energy_tev;
	  m_sim_table_data[itr->table_index] = table_datum;
	}
      else
	vsassert(itr->energy_tev == 
		 m_sim_table_data[itr->table_index]->egy_tev);

      m_sim_table_data[itr->table_index]->ntotal += itr->event_count;
    }
}

void VSAnalysisStage3Visitor::leaveSimHeader()
{

}

void VSAnalysisStage3Visitor::finalizeSim()
{  
  m_sim_data = m_sim_run_data.front();

  for(std::vector< VSStage3SimDatum* >::iterator itr = 
	m_sim_run_data.begin()+1; itr != m_sim_run_data.end(); ++itr)
    m_sim_data->merge(*itr);

//   for(std::vector<VSStage3SimArrayTableDatum*>::iterator itr =
// 	m_sim_table_data.begin(); itr != m_sim_table_data.end(); ++itr)
//     (*itr)->finalize();

  m_sim_data->finalize();
  
  VSStage3SimCalc sim_calc;

  sim_calc.analyze(*m_sim_data);
  sim_calc.analyze(m_sim_table_data);


  // Compute containment radius for thetasq -----------------------------------
  for(std::vector<VSStage3SimArrayTableDatum*>::iterator itr =
	m_sim_table_data.begin(); itr != m_sim_table_data.end(); ++itr)
    {
      double loge = log10((*itr)->egy_tev);

      m_sim_data->egymc_thsq68_hist.
	accumulate(loge, (*itr)->thsq68, std::pow((*itr)->thsq68_err,2));
      m_sim_data->egymc_thsq90_hist.
	accumulate(loge,(*itr)->thsq90, std::pow((*itr)->thsq90_err,2));
      m_sim_data->egymc_thsq95_hist.
	accumulate(loge, (*itr)->thsq95, std::pow((*itr)->thsq95_err,2));
    }
}



void VSAnalysisStage3Visitor::
calcParallacticAngle(const VSTime& lo_time, const VSTime& hi_time,
		     const SEphem::SphericalCoords& radec,
		     double& pangle_mean, double& pangle_rms)
{
  SEphem::Angle latitude;
  SEphem::Angle longitude;
  latitude.setFromDMSString("31d40m29.04s");
  longitude.setFromDMSString("-110d57m10.08s");
  SEphem::SphericalCoords earth_position;
  earth_position.setLatLong(latitude,longitude);

  VSSimpleStat2<double> pangle_stat;

  for(VSTime t = lo_time; t < hi_time; t+= (int64_t)1E9)
    {
      double mjd  = t.getMJDDbl();
      SEphem::Angle lmst = 
	SEphem::Astro::mjdToLMST(mjd,earth_position.longitudeRad());
      
      SEphem::SphericalCoords zn_radec(0,0);      
      SEphem::Astro::azElToMeanRaDec(lmst, mjd, earth_position, zn_radec);
      double pangle = radec.directionTo(zn_radec).deg();
      pangle_stat.accumulate(pangle);
    }

  pangle_mean = pangle_stat.mean();
  pangle_rms = pangle_stat.dev();
}
