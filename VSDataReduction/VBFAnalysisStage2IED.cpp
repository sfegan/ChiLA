//-*-mode:c++; mode:font-lock;-*-

/*! \file VBFAnalysisStage2IED.cpp

  Stage 2 analysis (interim event data)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 3.19 $
  \date       05/10/2006

  $Id: VBFAnalysisStage2IED.cpp,v 3.19 2009/11/05 22:52:57 matthew Exp $

*/

#include <VSAnalysisStage2.hpp>

using namespace VERITAS;
using namespace SEphem;

// ----------------------------------------------------------------------------
//  ___     _           _       ___             _   ___       _
// |_ _|_ _| |_ ___ _ _(_)_ __ | __|_ _____ _ _| |_|   \ __ _| |_ __ _
//  | || ' \  _/ -_) '_| | '  \| _|\ V / -_) ' \  _| |) / _` |  _/ _` |
// |___|_||_\__\___|_| |_|_|_|_|___|\_/\___|_||_\__|___/\__,_|\__\__,_|
// 
// ----------------------------------------------------------------------------

VBFAnalysisStage2::ThreadSpecificData::
ThreadSpecificData(const VSAnalysisStage1Data* stage1,
		   const Settings& settings, const ArrayMergedCal& _cal,
		   const ArrayDiagnostics* diagnostics, bool construct_rng):
  partial_diagnostics(), trace_nsample(stage1->run_info.nsample), 
  cal(_cal), rng(), m_diagnostics(diagnostics)
{
  if(diagnostics)
    partial_diagnostics = 
      new PartialArrayDiagnostics(stage1->run_info.config_mask, 
				  stage1->run_info.nchan, 
				  stage1->run_info.nsample);

  if(construct_rng)
    rng = new RandomNumbers(RandomNumbers::defaultFilename());
}

VBFAnalysisStage2::ThreadSpecificData::~ThreadSpecificData()
{
  delete partial_diagnostics;
  delete rng;
}

void VBFAnalysisStage2::ThreadSpecificData::destructor(void* ptr)
{
  delete static_cast<ThreadSpecificData*>(ptr);
}

VBFAnalysisStage2::ScopeIED::
ScopeIED(const IED* _ied, unsigned _iscope, unsigned _nchan):
  has_l3(),
  l3_sent(),
  l3_l2_received(),
  l3_az(),
  l3_el(),
  l3_counts_l2(),
  l3_tdc(),
  counts_l2(),
  has_vdaq(),
  processed_scope_event(),
  gps_dt(),
  nchan_hit(),
  sij_max1_raw_j(),
  sij_max2_raw_j(),
  sij_max3_raw_j(),
  trigger_j(),
  total_signal(),
  nchan_logain(),
  nchan_trigger(),
  nchan_image(),
  has_image(),
  sij_max1_sig_j(),
  sij_max2_sig_j(),
  sij_max3_sig_j(),
  sij_max1_sig(),
  sij_max2_sig(),
  sij_max3_sig(),
  largest_region_nchan(),
  az_rad(),
  zn_rad(),
  ra_rad(),
  dec_rad(),
  l_rad(),
  b_rad(),
  has_moments(),
  fp_Ni(),
  fp_xc(),
  fp_yc(),
  fp_dist(),
  fp_length(),
  fp_width(),
  fp_psi(),
  fp_ex(),
  fp_ey(),
  fp_disp(),
  intrinsic_length(),
  intrinsic_width(),
  used_in_reconstruction(),
  sc_width(),
  sc_length(),
  sc_disp(),
  lt_log10_energy(),
  lt_log10_energy_err(),
  has_muon(false),
  muon_data(),
  ied(_ied),
  iscope(_iscope),
  nchan(_nchan)
{
  // nothing to see here
}

VBFAnalysisStage2::ScopeIED::~ScopeIED()
{
  // nothing to see here
}

void VBFAnalysisStage2::ScopeIED::
transferDataToESD(VSEventScopeDatum* esd,
		  const VSAReconstruction::Reconstruction& recon,
		  const VSAReconstruction::ScopeParameters& sp,
		  const VSAReconstruction::ScopeMoments& sm)
{
  if(processed_scope_event)
    {
      esd->has_image              = has_image;
      esd->used_in_reconstruction = used_in_reconstruction;
      esd->nimage                 = nchan_image;
      esd->ntrig                  = nchan_trigger;
    }
  else
    {
      esd->has_image              = false;
      esd->used_in_reconstruction = false;
      esd->nimage                 = 0;
      esd->ntrig                  = 0;
    }

#define SQR(x) ((x)*(x))

  if((processed_scope_event)&&(has_moments)&&(ied->reconstruction_attempted))
    {
      esd->N                      = sp.Ni;
      esd->fp_trigger_ichan       = trigger_j;
      esd->fp_N                   = fp_Ni;
      esd->fp_xc                  = fp_xc;
      esd->fp_yc                  = fp_yc;
      esd->fp_dist                = fp_dist;
      esd->fp_length              = fp_length;
      esd->fp_width               = fp_width;
      esd->fp_psi                 = fp_psi;
      esd->intrinsic_length       = intrinsic_length;
      esd->intrinsic_width        = intrinsic_width;
    }
  else
    {
      esd->N                      = 0;
      esd->fp_trigger_ichan       = trigger_j;
      esd->fp_N                   = 0;
      esd->fp_xc                  = 0;
      esd->fp_yc                  = 0;
      esd->fp_dist                = 0;
      esd->fp_length              = 0;
      esd->fp_width               = 0;
      esd->fp_psi                 = 0;
      esd->intrinsic_length       = 0;
      esd->intrinsic_width        = 0;
    }

  if((processed_scope_event)&&(used_in_reconstruction)&&(sp.Ni > 0)
     &&(ied->reconstruction_successful)&&(recon.e.z() < 0))
    {
      VSAAlgebra::Eigen3D delta2i_eigen;
      sp.delta2i.eigen(delta2i_eigen);

      esd->R                      = sp.Ri;
      esd->d1                     = sp.d1i;
      esd->d2                     = sp.d2i;
      esd->theta1                 = Angle::toDeg(sp.theta1i);
      esd->theta2                 = Angle::toDeg(sp.theta2i);
      esd->delta1                 = sp.delta1i.norm();
      esd->delta2l                = sqrt(delta2i_eigen.val[2]);
      esd->delta2m                = sqrt(delta2i_eigen.val[1]);
      esd->G                      = sp.G1i;
      esd->t01                    = sp.t01i;
      esd->t02                    = sp.t02i;
      esd->lambdad                = sp.lambdadi;
      esd->lambdac                = sp.lambdaci;

      esd->fp_ex                  = fp_ex;
      esd->fp_ey                  = fp_ey;
      esd->fp_disp                = fp_disp;

      esd->sc_width               = sc_width;
      esd->sc_length              = sc_length;
      esd->sc_disp                = sc_disp;
      esd->lt_log10_energy        = lt_log10_energy;
      esd->lt_log10_energy_err    = lt_log10_energy_err;
    }
  else
    {
      esd->R                      = 0;
      esd->d1                     = 0;
      esd->d2                     = 0;
      esd->theta1                 = 0;
      esd->theta2                 = 0;
      esd->delta1                 = 0;
      esd->delta2l                = 0;
      esd->delta2m                = 0;
      esd->G                      = 0;
      esd->t01                    = 0;
      esd->t02                    = 0;
      esd->lambdad                = 0;
      esd->lambdac                = 0;

      esd->fp_ex                  = 0;
      esd->fp_ey                  = 0;
      esd->fp_disp                = 0;

      esd->sc_width               = 0;
      esd->sc_length              = 0;
      esd->sc_disp                = 0;
      esd->lt_log10_energy        = -std::numeric_limits<double>::infinity();
      esd->lt_log10_energy_err    = 0;
    }
}

VBFAnalysisStage2::IED::
IED(const VSAnalysisStage1Data* stage1, unsigned ntheta, 
    VSCutsCalc *cuts_calc):
  event_type(),
  event_failed_software_trigger(),
  event_l2_trigger_mask(),
  event_sent_trigger_mask(),
  event_has_vdaq_mask(),
  scope_num(),
  scope_ied(),
  scope_cal(),
  scope_diag(),
  scope_nchan(),
  scope_nchan_hit(),
  scope_nchan_logain(),
  scope_nchan_trigger(),
  scope_total_signal(),
  this_chan_trigger(),
  chan_has_signal(),
  chan_hit(),
  chan_sij(),
  chan_sij_total(),
  chan_tij(),
  chan_trigger(),
  chan_logain(),
  chan_sij_max1_raw(),
  chan_sij_max2_raw(),
  chan_sij_max3_raw(),
  chan_sij_max1_raw_j(),
  chan_sij_max2_raw_j(),
  chan_sij_max3_raw_j(),
  got_one_position(),
  zn_zero(),
  az_zero(),
  mean_position_x(),
  mean_position_y(),
  scope(stage1->run_info.nchan.size()),
  ied_processing(),
  tsd(),
  vbf_packet_number(),
  seq_packet_number(),
  has_array_event(),
  processed_array_event(),
  event_num(),
  has_event_time(),
  event_time(),
  event_time_sec(),
  event_time_hist(),
  ticks_elap(),
  ticks_both(),
  ticks_vdaq(),
  ticks_lev3(),
  has_l3(),
  ten_mhz_elapsed(),
  ten_mhz_veto_both(),
  ten_mhz_veto_vdaq(),
  ten_mhz_veto_l3(),
  nscope_trigger(),
  nscope_sent_l3(),
  gps_dt(),
  nscope_has_event(),
  event_has_image_mask(),
  nscope_image(),
  nscope_quality(),
  event_used_in_reconstruction_mask(),
  reconstruction_attempted(),
  reconstruction_successful(),
  images(stage1->run_info.nchan.size()),
  recon(),
  param(),
  secondary_images(stage1->run_info.nchan.size()),
  secondary_moments(),
  theta(ntheta),
  mjd(),
  lmst(),
  mean_az(),
  mean_zn(),
  recon_az(),
  recon_zn(),
  recon_ra(),
  recon_dec(),
  recon_ra_j2000(),
  recon_dec_j2000(),
  recon_fov_x(),
  recon_fov_y(),
  recon_derotated_fov_x(),
  recon_derotated_fov_y(),
  recon_r0_x(),
  recon_r0_y(),
  N2(),
  msc_width(),
  msc_length(),
  msc_disp(),
  mlt_log10_energy(),
  mlt_log10_energy_chi2(),
  write_event(),
  ead(stage1->run_info.nchan, ntheta),
  acd(),
  has_sim(false),
  sim(),
  sim_header(),
  sim_num(),
  sim_aomega(),
  last_event_time_sec(),
  last_event_time_hist(),
  last_event_num(),
  last_ten_mhz_elapsed(),
  last_ten_mhz_veto_both(),
  last_ten_mhz_veto_vdaq(),
  last_ten_mhz_veto_l3()
{
  unsigned max_nchan = 0;
  for(unsigned iscope=0;iscope<stage1->run_info.nchan.size();iscope++)
    if(stage1->run_info.nchan[iscope]>max_nchan)
      max_nchan=stage1->run_info.nchan[iscope];

  for(unsigned iscope=0;iscope<stage1->run_info.nchan.size();iscope++)
    scope[iscope] = new ScopeIED(this, iscope, stage1->run_info.nchan[iscope]);

  chan_has_signal .resize(max_nchan);
  chan_hit        .resize(max_nchan);
  chan_sij        = new double[max_nchan];
  chan_sij_total  = new double[max_nchan];
  chan_tij        = new double[max_nchan];
  chan_trigger    .resize(max_nchan);
  chan_logain     .resize(max_nchan);

  if(cuts_calc)acd = new VSArrayCutsDatum;
}

VBFAnalysisStage2::IED::~IED() 
{
  for(unsigned iscope=0;iscope<scope.size();iscope++)delete scope[iscope];
  delete[] chan_sij;
  delete[] chan_sij_total;
  delete[] chan_tij;
  delete acd;
  delete sim;
  delete sim_header;
}

void VBFAnalysisStage2::IED::
finalTransferDataToEAD(uint64_t event_time_ten_mhz_elap,
		       uint64_t event_time_ten_mhz_live)
{
  if(has_event_time)
    {
      ead.abs_event_time    = event_time;
      ead.mjd               = mjd;
      ead.event_time        = event_time_sec;
    }
  else
    {
      ead.abs_event_time    = VSTime();
      ead.mjd               = -1;
      ead.event_time        = -1;
    }

  ead.elapsed_ticks         = event_time_ten_mhz_elap;
  ead.live_ticks            = event_time_ten_mhz_live;
}

void VBFAnalysisStage2::IED::transferDataToEAD()
{
  ead.event_num             = event_num;

  ead.trigger_mask          = event_l2_trigger_mask;
  ead.l3_sent_mask          = event_sent_trigger_mask;
  ead.has_datum_mask        = event_has_vdaq_mask;
  ead.has_image_mask        = event_has_image_mask;
  ead.used_in_reconstruction_mask = event_used_in_reconstruction_mask;
  ead.mean_array_zn         = mean_zn;
  ead.mean_array_az         = mean_az;
  ead.nscope_image          = nscope_image;
  ead.chi2e                 = recon.chi2e;
  ead.chi2R                 = recon.chi2R;
  ead.mean_fov_x            = Angle::toDeg(recon_fov_x);
  ead.mean_fov_y            = Angle::toDeg(recon_fov_y);
  ead.mean_derotated_fov_x  = Angle::toDeg(recon_derotated_fov_x);
  ead.mean_derotated_fov_y  = Angle::toDeg(recon_derotated_fov_y);
  ead.zn                    = Angle::toDeg(recon_zn);
  ead.az                    = Angle::toDeg(recon_az);
  ead.ra                    = Angle::toDeg(recon_ra);
  ead.dec                   = Angle::toDeg(recon_dec);
  ead.ra_J2000              = Angle::toDeg(recon_ra_j2000);
  ead.dec_J2000             = Angle::toDeg(recon_dec_j2000);
  ead.Rx                    = recon_r0_x;
  ead.Ry                    = recon_r0_y;
  ead.R                     = sqrt(ead.Rx*ead.Rx + ead.Ry*ead.Ry);
  ead.deltael               = Angle::toDeg(recon.deltael);
  ead.deltaew               = Angle::toDeg(recon.deltaew);
  ead.deltaRl               = recon.deltaRl;
  ead.deltaRw               = recon.deltaRw;
  ead.N2                    = N2;
  ead.msc_width             = msc_width;
  ead.msc_length            = msc_length;
  ead.msc_disp              = msc_disp;
  ead.mlt_log10_energy      = mlt_log10_energy;
  ead.mlt_log10_energy_chi2 = mlt_log10_energy_chi2;

  unsigned ntheta = theta.size();
  if(ead.theta.size() != ntheta)
    ead.theta.resize(ntheta);

  for(unsigned itheta=0;itheta<ntheta;itheta++)
    ead.theta[itheta]       = theta[itheta]; // already in degrees

  if(ntheta > 1)
    {
      ead.theta0 = theta[0];
      ead.theta1 = theta[1];
    }
  else if(ntheta == 1)
    {
      ead.theta0 = theta[0];
      ead.theta1 = 0;
    }   
  else
    {
      ead.theta0 = 0;
      ead.theta1 = 0;
    }
  
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(ead.scope[iscope] && scope[iscope])
      scope[iscope]->transferDataToESD(ead.scope[iscope],
				       recon,
				       param.scopes[iscope],
				       recon.moments.scopes[iscope]);
}

