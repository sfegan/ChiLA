//-*-mode:c++; mode:font-lock;-*-

/*! \file VSStage3SimData.cpp

  Class to store results of the stage 3 analysis

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.11 $
  \date       08/16/2008

  $Id: VSResultsSimData.cpp,v 3.11 2009/11/12 23:58:40 matthew Exp $

*/

#include <VSResultsSimData.hpp>
#include <VSSimulationData.hpp>
#include <VSRunInfoData.hpp>
#include <VSNSpaceOctaveH5IO.hpp>

using namespace VERITAS;

// ============================================================================
// VSStage3SimArrayDatumBase
// ============================================================================
VSStage3SimArrayDatumBase::
VSStage3SimArrayDatumBase(const std::vector<unsigned>& nchan):
  ntotal(), ntriggered(), nreconstructed(), 
  ncuts_selected(), nselected(),
  th68(), th68_err(), th90(), th90_err(), th95(), th95_err(),
  thsq68(), thsq68_err(), thsq90(), thsq90_err(), thsq95(), thsq95_err(),
  thetasq_hist(0.001,0.,0.1),
  thetasq_core_hist(0.001,0.,0.1,10,0,500),
  thetasq_coremc_hist(0.001,0.,0.1,10,0,500),
  core_log10_egyerr_hist(10,0,500,0.02,-1,1),
  coremc_log10_egyerr_hist(10,0,500,0.02,-1,1),
  log10_N2_log10_egyerr_hist(0.1,0,5,0.02,-1,1),
  thetasq_log10_egyerr_hist(0.001,0.,0.1,0.02,-1,1),
  fov_total_hist(0.05,-2,2),
  fov_triggered_hist(0.05,-2,2),
  fov_reconstructed_hist(0.05,-2,2),
  fov_cuts_selected_hist(0.05,-2,2),
  fov_selected_hist(0.05,-2,2),
  core_triggered_hist(10.,-400.,400.),
  core_reconstructed_hist(10.,-400.,400.),
  core_reconstructed_err_x_hist(5.,-400.,400.),
  core_reconstructed_err_y_hist(5.,-400.,400.),
  core_reconstructed_err_r_hist(5.,0.,400.),
  core_cuts_selected_hist(10.,-400.,400.),
  core_cuts_selected_err_x_hist(5.,-400.,400.),
  core_cuts_selected_err_y_hist(5.,-400.,400.),
  core_cuts_selected_err_r_hist(5.,0.,400.),
  core_selected_hist(10.,-400.,400.),
  core_selected_err_x_hist(5.,-400.,400.),
  core_selected_err_y_hist(5.,-400.,400.),
  core_selected_err_r_hist(5.,0.,400.)
{
  const unsigned nscope = nchan.size();
  scope.resize(nscope);
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    if(nchan[iscope]) scope[iscope]=new VSStage3SimScopeDatum;
}

VSStage3SimArrayDatumBase::~VSStage3SimArrayDatumBase()
{
  const unsigned nscope = scope.size();
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    delete scope[iscope];
}

void VSStage3SimArrayDatumBase::setNChan(const std::vector<unsigned>& nchan)
{
  const unsigned nscope = nchan.size();
  scope.resize(nscope);
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    if(nchan[iscope] && !scope[iscope]) 
      scope[iscope]=new VSStage3SimScopeDatum;
}

void VSStage3SimArrayDatumBase::merge(VSStage3SimArrayDatumBase* sim)
{
  // Event Counts -------------------------------------------------------------
  ntotal += sim->ntotal;
  ntriggered += sim->ntriggered;
  nreconstructed += sim->nreconstructed;
  ncuts_selected += sim->ncuts_selected;
  nselected += sim->nselected;

  // Theta Histograms ---------------------------------------------------------
  thetasq_hist += sim->thetasq_hist;

  // Bias histograms ----------------------------------------------------------
  core_log10_egyerr_hist += sim->core_log10_egyerr_hist;
  coremc_log10_egyerr_hist += sim->coremc_log10_egyerr_hist;
  log10_N2_log10_egyerr_hist += sim->log10_N2_log10_egyerr_hist;
  thetasq_log10_egyerr_hist += sim->thetasq_log10_egyerr_hist;

  // Direction Histograms -----------------------------------------------------
  fov_total_hist += sim->fov_total_hist;
  fov_triggered_hist += sim->fov_triggered_hist;
  fov_reconstructed_hist += fov_reconstructed_hist;
  fov_cuts_selected_hist += fov_cuts_selected_hist;
  fov_selected_hist += sim->fov_selected_hist;

  // Core Position Histograms -------------------------------------------------
  core_triggered_hist += sim->core_triggered_hist;
  core_reconstructed_hist += sim->core_reconstructed_hist;
  core_reconstructed_err_x_hist += sim->core_reconstructed_err_x_hist;
  core_reconstructed_err_y_hist += sim->core_reconstructed_err_y_hist;
  core_reconstructed_err_r_hist += sim->core_reconstructed_err_r_hist;
  core_cuts_selected_hist += sim->core_cuts_selected_hist;
  core_cuts_selected_err_x_hist += sim->core_cuts_selected_err_x_hist;
  core_cuts_selected_err_y_hist += sim->core_cuts_selected_err_y_hist;
  core_cuts_selected_err_r_hist += sim->core_cuts_selected_err_r_hist;
  core_selected_hist += sim->core_selected_hist;
  core_selected_err_x_hist += sim->core_selected_err_x_hist;
  core_selected_err_y_hist += sim->core_selected_err_y_hist;
  core_selected_err_r_hist += sim->core_selected_err_r_hist;

  // Temporary data ---------------------------------------------------------
  m_thetasq.insert(m_thetasq.end(),
		   sim->m_thetasq.begin(),sim->m_thetasq.end());

  const unsigned nscope = scope.size();
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(scope[iscope] != NULL && iscope < sim->scope.size()) 
	scope[iscope]->merge(sim->scope[iscope]);
    }
}

void VSStage3SimArrayDatumBase::finalize()
{

}

void VSStage3SimArrayDatumBase::
accumulate(const VSArraySimulationDatum& sim,
	   const VSEventArrayDatum& event,
	   unsigned event_code, double wt)
{
  double thetasq = std::pow(event.theta1,2);
  double err_x = sim.core_east_m - event.Rx;
  double err_y = sim.core_north_m - event.Ry;
  double err_r = sqrt(std::pow(err_x,2) + std::pow(err_y,2));

  double log10_egy = log10(sim.energy_tev);
  double log10_egyerr = event.mlt_log10_energy-log10_egy;

  SEphem::SphericalCoords azzn = 
    SEphem::SphericalCoords::makeDeg(sim.primary_zenith_deg,
				    sim.primary_azimuth_deg);
  SEphem::SphericalCoords obs_azzn = 
    SEphem::SphericalCoords::makeDeg(sim.obs_zenith_deg,
				     sim.obs_azimuth_deg);

  azzn.rotate(0,-obs_azzn.thetaRad(),-obs_azzn.phiRad());
  double sim_x = -azzn.thetaDeg()*sin(azzn.phiRad());
  double sim_y = -azzn.thetaDeg()*cos(azzn.phiRad());

  fov_total_hist.accumulate(sim_x,sim_y,wt);
    
  // --------------------------------------------------------------------------
  // Triggered
  // --------------------------------------------------------------------------
  if(sim.has_array_event)
    {
      ntriggered+=wt;
      fov_triggered_hist.accumulate(sim_x,sim_y,wt);
      core_triggered_hist.accumulate(sim.core_x_m,sim.core_y_m,wt);
    }

  // --------------------------------------------------------------------------
  // Reconstructed
  // --------------------------------------------------------------------------
  if(event_code & VSAnalysisStage3Data::EC_RECONSTRUCTED)
    {
      nreconstructed+=wt;
      fov_reconstructed_hist.accumulate(sim_x,sim_y,wt);
      core_reconstructed_hist.accumulate(sim.core_x_m,sim.core_y_m,wt);
      core_reconstructed_err_r_hist.accumulate(err_r,wt,wt*wt);
      core_reconstructed_err_x_hist.accumulate(err_x,wt,wt*wt);
      core_reconstructed_err_y_hist.accumulate(err_y,wt,wt*wt);
    }

  // --------------------------------------------------------------------------
  // Cuts Selected
  // --------------------------------------------------------------------------
  if(event_code & VSAnalysisStage3Data::EC_SELECTED)
    {
      ncuts_selected+=wt;
      fov_cuts_selected_hist.accumulate(sim_x,sim_y,wt);
      core_cuts_selected_hist.accumulate(sim.core_x_m,sim.core_y_m,wt);
      core_cuts_selected_err_r_hist.accumulate(err_r,wt,wt*wt);
      core_cuts_selected_err_x_hist.accumulate(err_x,wt,wt*wt);
      core_cuts_selected_err_y_hist.accumulate(err_y,wt,wt*wt);
      thetasq_hist.accumulate(thetasq,wt,wt*wt);

      coremc_log10_egyerr_hist.accumulate(sim.core_R_m,log10_egyerr,wt);
      thetasq_log10_egyerr_hist.accumulate(thetasq,log10_egyerr,wt);
      log10_N2_log10_egyerr_hist.accumulate(log10(event.N2),log10_egyerr,wt);

      if(thetasq < 0.1) m_thetasq.push_back(thetasq);
    }

  // --------------------------------------------------------------------------
  // Selected
  // --------------------------------------------------------------------------
  if(event_code & VSAnalysisStage3Data::EC_ON)
    {
      nselected+=wt;
    }

  // Scope Datums -------------------------------------------------------------
  const unsigned nscope = scope.size();
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(!event.scope[iscope] || 
	 !(event_code & VSAnalysisStage3Data::EC_RECONSTRUCTED)) 
	continue;
      
      scope[iscope]->accumulate(sim,*event.scope[iscope],
				event_code,wt);
    }

}

void VSStage3SimArrayDatumBase::load(VSOctaveH5ReaderStruct* reader)
{
  reader->readCompositeHere(*this);

  thetasq_hist.load(reader->readStruct("thetasq_hist"));

  core_log10_egyerr_hist.load(reader->readStruct("core_log10_egyerr_hist"));
  coremc_log10_egyerr_hist.load(reader->readStruct("core_log10_egyerr_hist"));
  log10_N2_log10_egyerr_hist.
    load(reader->readStruct("log10_N2_log10_egyerr_hist"));
  thetasq_log10_egyerr_hist.
    load(reader->readStruct("thetasq_log10_egyerr_hist"));

  fov_total_hist.load(reader->readStruct("fov_total_hist"));
  fov_triggered_hist.load(reader->readStruct("fov_triggered_hist"));
  fov_reconstructed_hist.
    load(reader->readStruct("fov_reconstructed_hist"));
  fov_cuts_selected_hist.
    load(reader->readStruct("fov_cuts_selected_hist"));
  fov_selected_hist.load(reader->readStruct("fov_selected_hist"));

  core_triggered_hist.load(reader->readStruct("core_triggered_hist"));
  core_reconstructed_hist.load(reader->readStruct("core_reconstructed_hist"));
  core_reconstructed_err_x_hist.
    load(reader->readStruct("core_reconstructed_err_x_hist"));
  core_reconstructed_err_y_hist.
    load(reader->readStruct("core_reconstructed_err_y_hist"));
  core_reconstructed_err_r_hist.
    load(reader->readStruct("core_reconstructed_err_r_hist"));
  core_cuts_selected_hist.
    load(reader->readStruct("core_cuts_selected_hist"));
  core_cuts_selected_err_x_hist.
    load(reader->readStruct("core_cuts_selected_err_x_hist"));
  core_cuts_selected_err_y_hist.
    load(reader->readStruct("core_cuts_selected_err_y_hist"));
  core_cuts_selected_err_r_hist.
    load(reader->readStruct("core_cuts_selected_err_r_hist"));
  core_selected_hist.
    load(reader->readStruct("core_selected_hist"));
  core_selected_err_x_hist.
    load(reader->readStruct("core_selected_err_x_hist"));
  core_selected_err_y_hist.
    load(reader->readStruct("core_selected_err_y_hist"));
  core_selected_err_r_hist.
    load(reader->readStruct("core_selected_err_r_hist"));

  // Scope Datums -------------------------------------------------------------
  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("t");
  vsassert(wc);
  const unsigned nscope = wc->dimensions();
  scope.resize(nscope);
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(!wc->isStruct(iscope)) continue;
      VSOctaveH5ReaderStruct* ws = wc->readStruct(iscope);	
      vsassert(ws);
      scope[iscope] = new VSStage3SimScopeDatum;
      scope[iscope]->load(ws);
      delete ws;
    }
  delete wc;
}
    
void VSStage3SimArrayDatumBase::save(VSOctaveH5WriterStruct* writer,
				      bool write_hists) const
{  
  writer->writeCompositeHere(*this);

  thetasq_hist.save(writer->writeStruct("thetasq_hist"));

  if(write_hists)
    {
      core_log10_egyerr_hist.
	save(writer->writeStruct("core_log10_egyerr_hist"));
      coremc_log10_egyerr_hist.
	save(writer->writeStruct("coremc_log10_egyerr_hist"));
      log10_N2_log10_egyerr_hist.
	save(writer->writeStruct("log10_N2_log10_egyerr_hist"));
      thetasq_log10_egyerr_hist.
	save(writer->writeStruct("thetasq_log10_egyerr_hist"));

      fov_total_hist.save(writer->writeStruct("fov_total_hist"));
      fov_triggered_hist.save(writer->writeStruct("fov_triggered_hist"));
      fov_reconstructed_hist.
	save(writer->writeStruct("fov_reconstructed_hist"));
      fov_cuts_selected_hist.
	save(writer->writeStruct("fov_cuts_selected_hist"));
      fov_selected_hist.save(writer->writeStruct("fov_selected_hist"));

      core_triggered_hist.save(writer->writeStruct("core_triggered_hist"));
      core_reconstructed_hist.
	save(writer->writeStruct("core_reconstructed_hist"));
      core_reconstructed_err_x_hist.
	save(writer->writeStruct("core_reconstructed_err_x_hist"));
      core_reconstructed_err_y_hist.
	save(writer->writeStruct("core_reconstructed_err_y_hist"));
      core_reconstructed_err_r_hist.
	save(writer->writeStruct("core_reconstructed_err_r_hist"));
      core_cuts_selected_hist.
	save(writer->writeStruct("core_cuts_selected_hist"));
      core_cuts_selected_err_x_hist.
	save(writer->writeStruct("core_cuts_selected_err_x_hist"));
      core_cuts_selected_err_y_hist.
	save(writer->writeStruct("core_cuts_selected_err_y_hist"));
      core_cuts_selected_err_r_hist.
	save(writer->writeStruct("core_cuts_selected_err_r_hist"));
      core_selected_hist.
	save(writer->writeStruct("core_selected_hist"));
      core_selected_err_x_hist.
	save(writer->writeStruct("core_selected_err_x_hist"));
      core_selected_err_y_hist.
	save(writer->writeStruct("core_selected_err_y_hist"));
      core_selected_err_r_hist.
	save(writer->writeStruct("core_selected_err_r_hist"));
    }

  // Scope Datums -------------------------------------------------------------
  const unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("t",nscope);  
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(!scope[iscope]) continue;
      VSOctaveH5WriterStruct* s = wc->writeStruct(iscope);
      scope[iscope]->save(s,write_hists);
      delete s;
    }
  delete wc;
}

// ============================================================================
// VSStage3SimArrayTableDatum
// ============================================================================
VSStage3SimArrayTableDatum::
VSStage3SimArrayTableDatum(const std::vector<unsigned>& nchan):
  VSStage3SimArrayDatumBase(nchan), egy_tev(), 
  sampling_area(),
  effarea_triggered(), effarea_triggered_err(),
  effarea_reconstructed(), effarea_reconstructed_err(),
  effarea_cuts_selected(), effarea_cuts_selected_err(),
  effarea_selected(), effarea_selected_err(),
  egy_counts_hist()
{

}

VSStage3SimArrayTableDatum::~VSStage3SimArrayTableDatum()
{

}

void VSStage3SimArrayTableDatum::merge(VSStage3SimArrayTableDatum* sim)
{
  VSStage3SimArrayDatumBase::merge(sim);
}

void VSStage3SimArrayTableDatum::finalize()
{


}

void VSStage3SimArrayTableDatum::
accumulate(const VSArraySimulationDatum& sim,
	   const VSEventArrayDatum& event,
	   unsigned event_code, double wt)
{
  if(table_index != sim.table_index) return;

  VSStage3SimArrayDatumBase::accumulate(sim,event,event_code,wt);

  if(event_code & VSAnalysisStage3Data::EC_SELECTED)
    {
      egy_counts_hist.accumulate(event.mlt_log10_energy,wt,wt*wt);
    }
}

void VSStage3SimArrayTableDatum::setEnergyBinSize(double egy_binsize,
						  double egy_lo,
						  double egy_hi)
{
  egy_counts_hist = 
    VSLimitedErrorsHist<double, double>(egy_binsize*0.5,egy_lo,egy_hi);
}

void VSStage3SimArrayTableDatum::load(VSOctaveH5ReaderStruct* reader)
{
  reader->readCompositeHere(*this);

  egy_counts_hist.load(reader->readStruct("egy_counts_hist"));

  VSStage3SimArrayDatumBase::load(reader);
}
    
void VSStage3SimArrayTableDatum::save(VSOctaveH5WriterStruct* writer,
				       bool write_hists) const
{
  writer->writeCompositeHere(*this);
  
  egy_counts_hist.save(writer->writeStruct("egy_counts_hist"));

  VSStage3SimArrayDatumBase::save(writer,write_hists);
}

// ============================================================================
// VSStage3SimScopeDatum
// ============================================================================
VSStage3SimScopeDatum::VSStage3SimScopeDatum():
  log10_N_hist(0.1,0,5),
  fp_width_hist(0.01,0.0,0.5),
  fp_length_hist(0.01,0.0,0.5),
  fp_dist_hist(0.05,0.0,1.7),
  fp_centroid_hist(0.1,-1.7,1.7),
  core_log10_egyerr_hist(10,0,500,0.02,-1,1),
  log10_N_log10_egyerr_hist(0.1,0,5,0.02,-1,1),
  fp_dist_log10_egyerr_hist(0.05,0,2,0.02,-1,1),
  fp_disp_log10_egyerr_hist(0.05,0,2,0.02,-1,1),
  nimage_log10_egyerr_hist(1,0,40,0.02,-1,1)  
{

}

VSStage3SimScopeDatum::~VSStage3SimScopeDatum()
{

}

void VSStage3SimScopeDatum::merge(VSStage3SimScopeDatum* sim)
{
  // Histograms of Focal Plane Parameters -------------------------------------
  log10_N_hist += sim->log10_N_hist;
  fp_width_hist += sim->fp_width_hist;
  fp_length_hist += sim->fp_length_hist;
  fp_dist_hist += sim->fp_dist_hist;
  fp_centroid_hist += sim->fp_centroid_hist;

  // Energy Bias Histograms ---------------------------------------------------
  core_log10_egyerr_hist += sim->core_log10_egyerr_hist;
  log10_N_log10_egyerr_hist += sim->log10_N_log10_egyerr_hist;
  fp_dist_log10_egyerr_hist += sim->fp_dist_log10_egyerr_hist;
  fp_disp_log10_egyerr_hist += sim->fp_disp_log10_egyerr_hist;
  nimage_log10_egyerr_hist += sim->nimage_log10_egyerr_hist;
}

void VSStage3SimScopeDatum::
accumulate(const VSArraySimulationDatum& sim,
	   const VSEventScopeDatum& scope_data,
	   unsigned event_code, double wt)
{
  if(!scope_data.used_in_reconstruction) return;

  if(event_code & VSAnalysisStage3Data::EC_ON)
    {
      log10_N_hist.accumulate(log10(scope_data.N),wt,wt*wt);
      fp_width_hist.accumulate(scope_data.fp_width,wt,wt*wt);
      fp_length_hist.accumulate(scope_data.fp_length,wt,wt*wt);
      fp_dist_hist.accumulate(scope_data.fp_dist,wt,wt*wt);
      fp_centroid_hist.accumulate(scope_data.fp_xc,scope_data.fp_yc,wt);

      double log10_egy = log10(sim.energy_tev);
      double log10_egyerr = scope_data.lt_log10_energy - log10_egy;
	  
      core_log10_egyerr_hist.accumulate(scope_data.R,log10_egyerr,wt);
      log10_N_log10_egyerr_hist.accumulate(std::log10(scope_data.N),
					 log10_egyerr,wt);
      fp_dist_log10_egyerr_hist.accumulate(scope_data.fp_dist,log10_egyerr,wt);
      fp_disp_log10_egyerr_hist.accumulate(scope_data.fp_disp,log10_egyerr,wt);
      nimage_log10_egyerr_hist.accumulate(scope_data.nimage,log10_egyerr,wt);
    }      
}

void VSStage3SimScopeDatum::load(VSOctaveH5ReaderStruct* reader)
{
  log10_N_hist.load(reader->readStruct("log10_N_hist"));
  fp_width_hist.load(reader->readStruct("fp_width_hist"));
  fp_length_hist.load(reader->readStruct("fp_length_hist"));
  fp_dist_hist.load(reader->readStruct("fp_dist_hist"));
  fp_centroid_hist.load(reader->readStruct("fp_centroid_hist"));
  core_log10_egyerr_hist.load(reader->readStruct("core_log10_egyerr_hist"));
  log10_N_log10_egyerr_hist.
    load(reader->readStruct("log10_N_log10_egyerr_hist"));
  fp_dist_log10_egyerr_hist.
    load(reader->readStruct("fp_dist_log10_egyerr_hist"));
  fp_disp_log10_egyerr_hist.
    load(reader->readStruct("fp_disp_log10_egyerr_hist"));
  nimage_log10_egyerr_hist.
    load(reader->readStruct("nimage_log10_egyerr_hist"));
}
    
void VSStage3SimScopeDatum::save(VSOctaveH5WriterStruct* writer,
				  bool write_hists) const
{
  if(write_hists)
    {
      log10_N_hist.save(writer->writeStruct("log10_N_hist"));
      fp_width_hist.save(writer->writeStruct("fp_width_hist"));
      fp_length_hist.save(writer->writeStruct("fp_length_hist"));
      fp_dist_hist.save(writer->writeStruct("fp_dist_hist"));
      fp_centroid_hist.save(writer->writeStruct("fp_centroid_hist"));
      core_log10_egyerr_hist.
	save(writer->writeStruct("core_log10_egyerr_hist"));
      log10_N_log10_egyerr_hist.
	save(writer->writeStruct("log10_N_log10_egyerr_hist"));
      fp_dist_log10_egyerr_hist.
	save(writer->writeStruct("fp_dist_log10_egyerr_hist"));
      fp_disp_log10_egyerr_hist.
	save(writer->writeStruct("fp_disp_log10_egyerr_hist"));
      nimage_log10_egyerr_hist.
	save(writer->writeStruct("nimage_log10_egyerr_hist"));
    }
}

// ============================================================================
// VSStage3SimDatum
// ============================================================================
VSStage3SimDatum::VSStage3SimDatum(const std::vector<unsigned>& nchan):
  VSStage3SimArrayDatumBase(nchan),
  src_ra_deg(), src_dec_deg(), obs_ra_deg(), obs_dec_deg(),
  gamma_rate_triggered(), gamma_rate_triggered_err(), 
  gamma_rate_reconstructed(), gamma_rate_reconstructed_err(),
  gamma_rate_cuts_selected(), gamma_rate_cuts_selected_err(),
  gamma_rate_selected(), gamma_rate_selected_err(),
  crab_rate_triggered(), crab_rate_triggered_err(), 
  crab_rate_reconstructed(), crab_rate_reconstructed_err(),
  crab_rate_cuts_selected(), crab_rate_cuts_selected_err(),
  crab_rate_selected(), crab_rate_selected_err(),
  proton_rate_triggered(), proton_rate_triggered_err(), 
  proton_rate_reconstructed(), proton_rate_reconstructed_err(), 
  proton_rate_cuts_selected(), proton_rate_cuts_selected_err(), 
  proton_rate_selected(), proton_rate_selected_err(),
  egy_threshold_triggered(),
  egy_threshold_reconstructed(),
  egy_threshold_cuts_selected(),
  egy_threshold_selected(),
  log10_egy_binsize(), log10_egy_min(), log10_egy_max(),
  table_id_hist(1), 
  egymc_count_hist(), 
  egymc_total_hist(), 
  egymc_triggered_hist(), 
  egymc_reconstructed_hist(), 
  egymc_cuts_selected_hist(),
  egymc_selected_hist(),
  egymc_fluence_hist(),
  egymc_sampling_area_hist(),
  effarea_triggered_hist(),
  effarea_reconstructed_hist(),
  effarea_cuts_selected_hist(),
  effarea_selected_hist(),
  diffrate_triggered_hist(),
  diffrate_reconstructed_hist(),
  diffrate_cuts_selected_hist(),
  diffrate_selected_hist(),
  egymc_core_total_hist(),
  egymc_core_triggered_hist(),
  egymc_core_triggered_norm_hist(),
  egymc_core_triggered_eff_hist(),
  egymc_core_reconstructed_hist(),
  egymc_core_reconstructed_norm_hist(),
  egymc_core_reconstructed_eff_hist(),
  egymc_core_selected_hist(),
  egymc_core_selected_norm_hist(),
  egymc_core_selected_eff_hist(),
  egymc_thetasq_hist(), egymc_thetasq_norm_hist(), egymc_theta_hist(),
  egymc_thsq68_hist(), egymc_thsq90_hist(), egymc_thsq95_hist(),
  egy_kernel(), effarea(),
  egymc_bias_hist(), egymc_rms_hist(),
  egymc_log10_bias_hist(), egymc_log10_rms_hist()
{

}

VSStage3SimDatum::~VSStage3SimDatum()
{

}

void VSStage3SimDatum::merge(VSStage3SimDatum* sim)
{
  // Event Counts -------------------------------------------------------------
  table_id_hist += sim->table_id_hist;
  egymc_count_hist += sim->egymc_count_hist;
  egymc_total_hist += sim->egymc_total_hist;
  egymc_triggered_hist += sim->egymc_triggered_hist;
  egymc_reconstructed_hist += sim->egymc_reconstructed_hist;
  egymc_cuts_selected_hist += sim->egymc_cuts_selected_hist;
  egymc_selected_hist += sim->egymc_selected_hist;

  // Trigger and Selection efficiencies ---------------------------------------
  egymc_core_total_hist += sim->egymc_core_total_hist;
  egymc_core_triggered_hist += sim->egymc_core_triggered_hist;
  egymc_core_reconstructed_hist += sim->egymc_core_reconstructed_hist;
  egymc_core_selected_hist += sim->egymc_core_selected_hist;

  // TeV resolution -----------------------------------------------------------
  egymc_thetasq_hist += sim->egymc_thetasq_hist;
  egymc_thetasq_norm_hist += sim->egymc_thetasq_norm_hist;
  egymc_theta_hist += sim->egymc_theta_hist;
  
  // Energy Reconstruction ----------------------------------------------------
  egy_kernel += sim->egy_kernel;

  VSStage3SimArrayDatumBase::merge(sim);
}

void VSStage3SimDatum::finalize()
{

  // Fill the theta opt hist --------------------------------------------------
//   const unsigned nenergy = egy_theta_hist.nXBins();
//   for(unsigned ienergy = 0; ienergy < nenergy; ienergy++)
//     {
//       double log10_egy = egy_theta_hist.xBinToCenter(ienergy);

//       const unsigned ntheta = egy_theta_hist.nYBins();

//       double sum = 0;
//       double theta = 0;
//       for(unsigned itheta = 0; itheta < ntheta; itheta++)
// 	{
// 	  theta = egy_theta_hist.yBinToVal(itheta) + egy_theta_hist.yBinSize();
// 	  sum += egy_theta_hist.count(ienergy,itheta);
// 	  if(theta >  thetaopt) break;
// 	}

//       egy_thetaopt_hist.accumulate(log10_egy,sum,sum);
//     }

//   VSLimitedErrorsHist<double, double> theta_cum_hist =
//     theta_hist.getCumulativeHist();

//   double tot = theta_cum_hist.count(theta_cum_hist.onePastLastBin()-1);
//   for(VSLimitedErrorsHist<double, double>::iterator itr = 
// 	theta_cum_hist.begin(); itr != theta_cum_hist.end(); ++itr)
//     {
//       if(!th68 && itr->count()/tot > 0.68) 
// 	th68 = itr->val() + theta_cum_hist.binSize();
//       if(!th90 && itr->count()/tot > 0.90) 
// 	th90 = itr->val() + theta_cum_hist.binSize();
//       if(!th95 && itr->count()/tot > 0.95) 
// 	th95 = itr->val() + theta_cum_hist.binSize();
//     }

//   thsq68 = std::pow(th68,2);
//   thsq90 = std::pow(th90,2);
//   thsq95 = std::pow(th95,2);

  // --------------------------------------------------------------------------
  // Calculate rates for a Crab-like source
  // --------------------------------------------------------------------------
//   const double crab_flux_constant = 3.2E-7; // m^-2 s^-1 TeV^-1
//   const double crab_index = 2.5;

//   for(VSLimitedErrorsHist<double, double>::iterator itr = 
// 	egy_total_hist.begin(); itr != egy_total_hist.end(); ++itr)
//     {
//       double log10_egy = itr->center();      
//       double egy_tev = std::pow(10,log10_egy);
//       VSNSpace::Point p(1);
//       p.x[0] = log10_egy;

//       double crab25_diff_flux = crab_flux_constant*pow(egy_tev,-crab_index);

//       if(egy_total_hist.countForVal(log10_egy))
// 	effarea[p] = selected_effarea_hist.countForVal(log10_egy);
//     }
  
  // Calculate efficiency as a function of projected distance -----------------
  egymc_core_triggered_eff_hist.clear();
  egymc_core_reconstructed_eff_hist.clear();
  egymc_core_selected_eff_hist.clear();
  egymc_core_triggered_norm_hist.clear();
  egymc_core_reconstructed_norm_hist.clear();
  egymc_core_selected_norm_hist.clear();

  for(VSSimple2DHist<double,double>::iterator itr = 
	egymc_core_total_hist.begin();
      itr != egymc_core_total_hist.end(); ++itr)
    {
      int xbin = egymc_core_total_hist.indexToXBin(itr->bin());
      double x =  egymc_core_total_hist.xBinToCenter(xbin);

      double ntot = egymc_core_total_hist.countForIndex(itr->bin());
      double ntrig = egymc_core_triggered_hist.countForIndex(itr->bin());
      double nrec = egymc_core_reconstructed_hist.countForIndex(itr->bin());
      double nsel = egymc_core_selected_hist.countForIndex(itr->bin());

      double ntrig_tot = egymc_triggered_hist.countForVal(x);
      double nrec_tot = egymc_reconstructed_hist.countForVal(x);
      double nsel_tot = egymc_selected_hist.countForVal(x);

      if(ntot != 0)
	{
	  egymc_core_triggered_eff_hist.accumulateBin(itr->bin(),ntrig/ntot);
	  egymc_core_reconstructed_eff_hist.
	    accumulateBin(itr->bin(),nrec/ntot);
	  egymc_core_selected_eff_hist.
	    accumulateBin(itr->bin(),nsel/ntot);	 
	}

      egymc_core_triggered_norm_hist.accumulateBin(itr->bin(),ntrig/ntrig_tot);
      egymc_core_reconstructed_norm_hist.
	accumulateBin(itr->bin(),nrec/nrec_tot);
      egymc_core_selected_norm_hist.accumulateBin(itr->bin(),nsel/nsel_tot);
    }
}

void VSStage3SimDatum::
accumulate(const VSArraySimulationDatum& sim,
	   const VSEventArrayDatum& event,
	   unsigned event_code, double wt)
{
  double thetasq = std::pow(event.theta1,2);
  double log10_egy = log10(sim.energy_tev);

  egymc_core_total_hist.accumulate(log10_egy,sim.core_R_m,wt);
  table_id_hist.accumulate(sim.table_index);

  // --------------------------------------------------------------------------
  // Triggered
  // --------------------------------------------------------------------------
  if(event_code & VSAnalysisStage3Data::EC_TRIGGERED)
    {
      egymc_core_triggered_hist.accumulate(log10_egy,sim.core_R_m,wt);
      egymc_triggered_hist.accumulate(log10_egy,wt,wt*wt);
    }

  // --------------------------------------------------------------------------
  // Reconstructed
  // --------------------------------------------------------------------------
  if(event_code & VSAnalysisStage3Data::EC_RECONSTRUCTED)
    {
      egymc_core_reconstructed_hist.accumulate(log10_egy,sim.core_R_m,wt);
      egymc_reconstructed_hist.accumulate(log10_egy,wt,wt*wt);
    }
  
  // --------------------------------------------------------------------------
  // Cuts Selected
  // --------------------------------------------------------------------------
  if(event_code & VSAnalysisStage3Data::EC_SELECTED)
    {
      egymc_cuts_selected_hist.accumulate(log10_egy,wt,wt*wt);      
      egymc_thetasq_hist.accumulate(log10_egy,thetasq,wt);
      egymc_thetasq_norm_hist.accumulate(log10_egy,thetasq,wt);
      egymc_theta_hist.accumulate(log10_egy,event.theta1,wt);

      egymc_egy_hist.accumulate(log10_egy,event.mlt_log10_energy);
    }

  // --------------------------------------------------------------------------
  // Selected
  // --------------------------------------------------------------------------
  if(event_code & VSAnalysisStage3Data::EC_ON)
    {
      double log10_egyerr = event.mlt_log10_energy-log10_egy;

      VSNSpace::Point p_kernel(2);
      p_kernel.x[0] = log10_egy;
      p_kernel.x[1] = event.mlt_log10_energy;

      egymc_selected_hist.accumulate(log10_egy,wt,wt*wt);
      egy_kernel.accumulate(p_kernel);      
      egymc_core_selected_hist.accumulate(log10_egy,sim.core_R_m,wt);
      egymc_log10_egyerr_hist.accumulate(log10_egy,log10_egyerr,wt);       
    }


  // Table Datums -------------------------------------------------------------
//   table[sim.table_index]->accumulate(sim,event,
// 				     is_reconstructed,
// 				     is_cuts_selected,
// 				     is_selected,wt);

  // Scope Datums -------------------------------------------------------------
//   const unsigned nscope = scope.size();
//   for(unsigned iscope = 0; iscope < nscope; iscope++)
//     {
//       if(!event.scope[iscope] || !is_reconstructed) continue;
      
//       scope[iscope]->accumulate(sim,*event.scope[iscope],
// 				is_reconstructed,
// 				is_cuts_selected,
// 				is_selected,wt);
//     }

  VSStage3SimArrayDatumBase::accumulate(sim,event,event_code,wt);
}

void VSStage3SimDatum::setEnergyBinSize(double ebin, double elo, double ehi)
{
  vsassert(ehi > elo && ebin > 0);

  // 1D Histograms ------------------------------------------------------------
  typedef VSLimitedErrorsHist<double,double> Hist1D;

  egymc_count_hist = Hist1D(ebin, elo, ehi);
  egymc_total_hist = Hist1D(ebin, elo, ehi);
  egymc_triggered_hist = Hist1D(ebin, elo, ehi);  
  egymc_reconstructed_hist = Hist1D(ebin, elo, ehi);
  egymc_cuts_selected_hist = Hist1D(ebin, elo, ehi);  
  egymc_selected_hist = Hist1D(ebin, elo, ehi);    
  egymc_fluence_hist = Hist1D(ebin, elo, ehi);    
  egymc_sampling_area_hist = Hist1D(ebin, elo, ehi);   

  egymc_bias_hist = Hist1D(ebin, elo, ehi);
  egymc_rms_hist = Hist1D(ebin, elo, ehi);
  egymc_log10_bias_hist = Hist1D(ebin, elo, ehi);
  egymc_log10_rms_hist = Hist1D(ebin, elo, ehi);

  // 2D Histograms ------------------------------------------------------------
  typedef VSSimple2DHist<double,double> Hist2D;

  const double core_dist_max = 600;

  egymc_core_total_hist = Hist2D(ebin, elo,ehi, 10.,0.,core_dist_max);
  egymc_core_triggered_hist = Hist2D(ebin, elo, ehi, 10.,0.,core_dist_max);
  egymc_core_triggered_norm_hist = Hist2D(ebin, elo,ehi, 10.,0.,core_dist_max);
  egymc_core_triggered_eff_hist = Hist2D(ebin, elo, ehi, 10.,0.,core_dist_max);
  egymc_core_reconstructed_hist = Hist2D(ebin, elo, ehi, 10.,0.,core_dist_max);
  egymc_core_reconstructed_norm_hist = 
    Hist2D(ebin, elo,ehi, 10.,0.,core_dist_max);
  egymc_core_reconstructed_eff_hist = 
    Hist2D(ebin, elo, ehi, 10.,0.,core_dist_max);
  egymc_core_selected_hist = Hist2D(ebin, elo, ehi, 10.,0.,core_dist_max);
  egymc_core_selected_norm_hist = Hist2D(ebin, elo, ehi, 10.,0.,core_dist_max);
  egymc_core_selected_eff_hist = Hist2D(ebin, elo, ehi, 10.,0.,core_dist_max);
  egymc_thetasq_hist = Hist2D(ebin, elo, ehi, 0.001,0.,0.06);  
  egymc_thetasq_norm_hist = Hist2D(ebin, elo, ehi, 0.001,0.,0.06);
  egymc_theta_hist = Hist2D(ebin, elo, ehi, 0.0005,0.0,0.32);  

  unsigned nbin = lround((ehi-elo)/ebin);

  effarea = VSNSpace(1,elo, ehi, nbin);
  egy_kernel = VSNSpace(2,elo, ehi, nbin);

  egymc_egy_hist = Hist2D(ebin,elo,ehi, 0.5*ebin,elo, ehi);
}

void VSStage3SimDatum::load(VSOctaveH5ReaderStruct* reader)
{
  reader->readCompositeHere(*this);

  // Event counts -------------------------------------------------------------
  table_id_hist.load(reader->readStruct("table_id_hist"));
  egymc_total_hist.load(reader->readStruct("egymc_total_hist"));
  egymc_triggered_hist.load(reader->readStruct("egymc_triggered_hist"));
  egymc_reconstructed_hist.
    load(reader->readStruct("egymc_reconstructed_hist"));
  egymc_cuts_selected_hist.
    load(reader->readStruct("egymc_cuts_selected_hist"));
  egymc_selected_hist.load(reader->readStruct("egymc_selected_hist"));
  egymc_fluence_hist.load(reader->readStruct("egymc_fluence_hist"));
  egymc_sampling_area_hist.
    load(reader->readStruct("egymc_sampling_area_hist"));

  // Effective Area Hists ----------------------------------------------------
  effarea_triggered_hist.load(reader->readStruct("effarea_triggered_hist"));
  effarea_reconstructed_hist.
    load(reader->readStruct("effarea_reconstructed_hist"));
  effarea_cuts_selected_hist.
    load(reader->readStruct("effarea_cuts_selected_hist"));
  effarea_selected_hist.load(reader->readStruct("effarea_selected_hist"));

  // Differential Rate Hists -------------------------------------------------
  diffrate_triggered_hist.
    load(reader->readStruct("diffrate_triggered_hist"));
  diffrate_reconstructed_hist.
    load(reader->readStruct("diffrate_reconstructed_hist"));
  diffrate_cuts_selected_hist.
    load(reader->readStruct("diffrate_cuts_selected_hist"));
  diffrate_selected_hist.load(reader->readStruct("diffrate_selected_hist"));

  // Trigger and Selection efficiencies ---------------------------------------
  egymc_core_total_hist.load(reader->readStruct("egymc_core_total_hist"));
  egymc_core_triggered_hist.
    load(reader->readStruct("egymc_core_triggered_hist"));
  egymc_core_triggered_norm_hist.
    load(reader->readStruct("egymc_core_triggered_eff_hist"));
  egymc_core_triggered_eff_hist.
    load(reader->readStruct("egymc_core_triggered_norm_hist"));
  egymc_core_reconstructed_hist.
    load(reader->readStruct("egymc_core_reconstructed_hist"));
  egymc_core_reconstructed_norm_hist.
    load(reader->readStruct("egymc_core_reconstructed_eff_hist"));
  egymc_core_reconstructed_eff_hist.
    load(reader->readStruct("egymc_core_reconstructed_norm_hist"));
  egymc_core_selected_hist.
    load(reader->readStruct("egymc_core_selected_hist"));
  egymc_core_selected_norm_hist.
    load(reader->readStruct("egymc_core_selected_norm_hist"));
  egymc_core_selected_eff_hist.
    load(reader->readStruct("egymc_core_selected_eff_hist"));

  // TeV resolution -----------------------------------------------------------
  egymc_thetasq_hist.load(reader->readStruct("egymc_thetasq_hist"));
  egymc_thetasq_norm_hist.load(reader->readStruct("egymc_thetasq_norm_hist"));
  egymc_theta_hist.load(reader->readStruct("egymc_theta_hist"));
  egymc_thsq68_hist.load(reader->readStruct("egymc_thsq68_hist"));
  egymc_thsq90_hist.load(reader->readStruct("egymc_thsq90_hist"));
  egymc_thsq95_hist.load(reader->readStruct("egymc_thsq95_hist"));

  // Energy Reconstruction ----------------------------------------------------
  egymc_egy_hist.load(reader->readStruct("egymc_egy_hist"));

  VSNSpaceOctaveH5IO io;
  io.readHistogram(reader->readStruct("egy_kernel"),egy_kernel);
  io.readHistogram(reader->readStruct("effarea"),effarea);

  egymc_log10_egyerr_hist.
    load(reader->readStruct("egymc_log10_egyerr_hist"));
  
  egymc_bias_hist.load(reader->readStruct("egymc_bias_hist"));
  egymc_rms_hist.load(reader->readStruct("egymc_rms_hist"));
  egymc_log10_bias_hist.load(reader->readStruct("egymc_log10_bias_hist"));
  egymc_log10_rms_hist.load(reader->readStruct("egymc_log10_rms_hist"));


  VSStage3SimArrayDatumBase::load(reader);  
}

void VSStage3SimDatum::save(VSOctaveH5WriterStruct* writer, 
			     bool write_hists) const
{
  writer->writeCompositeHere(*this);
  
  if(write_hists)
    {
      // Event counts ---------------------------------------------------------
      table_id_hist.save(writer->writeStruct("table_id_hist"));
      egymc_total_hist.save(writer->writeStruct("egymc_total_hist"));
      egymc_triggered_hist.save(writer->writeStruct("egymc_triggered_hist"));
      egymc_reconstructed_hist.
	save(writer->writeStruct("egymc_reconstructed_hist"));
      egymc_cuts_selected_hist.
	save(writer->writeStruct("egymc_cuts_selected_hist"));
      egymc_selected_hist.save(writer->writeStruct("egymc_selected_hist"));
      egymc_fluence_hist.save(writer->writeStruct("egymc_fluence_hist"));
      egymc_sampling_area_hist.
	save(writer->writeStruct("egymc_sampling_area_hist"));

      // Effective Area Hists ------------------------------------------------
      effarea_triggered_hist.
	save(writer->writeStruct("effarea_triggered_hist"));
      effarea_reconstructed_hist.
	save(writer->writeStruct("effarea_reconstructed_hist"));
      effarea_cuts_selected_hist.
	save(writer->writeStruct("effarea_cuts_selected_hist"));
      effarea_selected_hist.
	save(writer->writeStruct("effarea_selected_hist"));

      // Differential Rate Hists ---------------------------------------------
      diffrate_triggered_hist.
	save(writer->writeStruct("diffrate_triggered_hist"));
      diffrate_reconstructed_hist.
	save(writer->writeStruct("diffrate_reconstructed_hist"));
      diffrate_cuts_selected_hist.
	save(writer->writeStruct("diffrate_cuts_selected_hist"));
      diffrate_selected_hist.
	save(writer->writeStruct("diffrate_selected_hist"));

      // Trigger and Selection efficiencies -----------------------------------
      egymc_core_total_hist.save(writer->writeStruct("egymc_core_total_hist"));
      egymc_core_triggered_hist.
	save(writer->writeStruct("egymc_core_triggered_hist"));
      egymc_core_triggered_norm_hist.
	save(writer->writeStruct("egymc_core_triggered_norm_hist"));
      egymc_core_triggered_eff_hist.
	save(writer->writeStruct("egymc_core_triggered_eff_hist"));
      egymc_core_reconstructed_hist.
	save(writer->writeStruct("egymc_core_reconstructed_hist"));
      egymc_core_reconstructed_norm_hist.
	save(writer->writeStruct("egymc_core_reconstructed_norm_hist"));
      egymc_core_reconstructed_eff_hist.
	save(writer->writeStruct("egymc_core_reconstructed_eff_hist"));
      egymc_core_selected_hist.
	save(writer->writeStruct("egymc_core_selected_hist"));
      egymc_core_selected_norm_hist.
	save(writer->writeStruct("egymc_core_selected_norm_hist"));
      egymc_core_selected_eff_hist.
	save(writer->writeStruct("egymc_core_selected_eff_hist"));

      // TeV resolution -------------------------------------------------------
      egymc_thetasq_hist.save(writer->writeStruct("egymc_thetasq_hist"));
      egymc_thetasq_norm_hist.
	save(writer->writeStruct("egymc_thetasq_norm_hist"));
      egymc_theta_hist.save(writer->writeStruct("egymc_theta_hist"));
      egymc_thsq68_hist.save(writer->writeStruct("egymc_thsq68_hist"));
      egymc_thsq90_hist.save(writer->writeStruct("egymc_thsq90_hist"));
      egymc_thsq95_hist.save(writer->writeStruct("egymc_thsq95_hist"));

      // Energy Reconstruction ------------------------------------------------
      egymc_egy_hist.save(writer->writeStruct("egymc_egymc_hist"));

      VSNSpaceOctaveH5IO io;
      io.writeHistogram(writer->writeStruct("egy_kernel"),egy_kernel);
      io.writeHistogram(writer->writeStruct("effarea"),effarea);

      egymc_log10_egyerr_hist.
	save(writer->writeStruct("egymc_log10_egyerr_hist"));

      egymc_log10_bias_hist.save(writer->writeStruct("egymc_log10_bias_hist"));
      egymc_log10_rms_hist.save(writer->writeStruct("egymc_log10_rms_hist"));
      egymc_bias_hist.save(writer->writeStruct("egymc_bias_hist"));
      egymc_rms_hist.save(writer->writeStruct("egymc_rms_hist"));
    }

  VSStage3SimArrayDatumBase::save(writer,write_hists);  
}

