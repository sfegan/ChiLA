//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage3Data.cpp

  Stage 3 data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/12/2006

  $Id: VSAnalysisStage3Data.cpp,v 3.25 2010/10/20 03:59:31 matthew Exp $

*/


#include <VSNSpaceOctaveH5IO.hpp>
#include <VSAnalysisStage3Data.hpp>

using namespace VERITAS;
using namespace SEphem;

// ============================================================================
// VSAnalysisStage3Datum::Scope
// ============================================================================
// VSAnalysisStage3Datum::Scope::Scope():
//   fp_dist(0.05,0.,1.7),
//   fp_disp(0.05,0.,1.7),
//   fp_width(0.01,0.,0.3),
//   fp_length(0.01,0.,0.5),
//   log10_N(0.1,1,6),
//   G(10,0,1000),
//   delta2l(0.5,0,50),
//   delta2m(0.5,0,25),
//   log10_lambdac(0.1,0,10),
//   log10_lambdad(0.1,0,10),
//   sc_disp(0.1,-4,10),
//   sc_width(0.1,-4,10),
//   sc_length(0.1,-4,10)
// {
  
// }

// void 
// VSAnalysisStage3Datum::Scope::accumulate(const VSEventScopeDatum& scope_data,
// 					 bool is_on_event,
// 					 bool is_off_event)
// {
//   if(!scope_data.used_in_reconstruction) return;

//   fp_dist.accumulate(scope_data.fp_dist,is_on_event,is_off_event);
//   fp_disp.accumulate(scope_data.fp_disp,is_on_event,is_off_event);
//   fp_width.accumulate(scope_data.fp_width,is_on_event,is_off_event);
//   fp_length.accumulate(scope_data.fp_length,is_on_event,is_off_event);
//   log10_N.accumulate(log10(scope_data.N),is_on_event,is_off_event);
//   G.accumulate(scope_data.G,is_on_event,is_off_event);
//   delta2l.accumulate(scope_data.delta2l,is_on_event,is_off_event);
//   delta2m.accumulate(scope_data.delta2m,is_on_event,is_off_event);
//   log10_lambdac.accumulate(log10(scope_data.lambdac),is_on_event,is_off_event);
//   log10_lambdad.accumulate(log10(scope_data.lambdad),is_on_event,is_off_event);
//   sc_disp.accumulate(scope_data.sc_disp,is_on_event,is_off_event);
//   sc_width.accumulate(scope_data.sc_width,is_on_event,is_off_event);
//   sc_length.accumulate(scope_data.sc_length,is_on_event,is_off_event);
// }

// void VSAnalysisStage3Datum::Scope::calcExcess(double alpha)
// {
// //   const unsigned nhist1d = hist1d.size();
// //   for(unsigned ihist = 0; ihist < nhist1d; ihist++)
// //     hist1d[ihist].calcExcess(alpha);

// //   const unsigned nhist2d = hist2d.size();
// //   for(unsigned ihist = 0; ihist < nhist2d; ihist++)
// //     hist2d[ihist].calcExcess(alpha);

//   fp_dist.calcExcess(alpha);
//   fp_disp.calcExcess(alpha);
//   fp_width.calcExcess(alpha);
//   fp_length.calcExcess(alpha);
//   log10_N.calcExcess(alpha);
//   G.calcExcess(alpha);
//   delta2l.calcExcess(alpha);
//   delta2m.calcExcess(alpha);
//   log10_lambdac.calcExcess(alpha);
//   log10_lambdad.calcExcess(alpha);
//   sc_disp.calcExcess(alpha);
//   sc_width.calcExcess(alpha);
//   sc_length.calcExcess(alpha);
// }

// void VSAnalysisStage3Datum::Scope::merge(VSAnalysisStage3Datum::Scope* results)
// {
// //   const unsigned nhist1d = results->hist1d.size();
// //   for(unsigned ihist = 0; ihist < nhist1d; ihist++)
// //     hist1d[ihist] += results->hist1d[ihist];

// //   const unsigned nhist2d = results->hist2d.size();
// //   for(unsigned ihist = 0; ihist < nhist2d; ihist++)
// //     hist2d[ihist] += results->hist2d[ihist];

//   fp_dist += results->fp_dist;
//   fp_disp += results->fp_disp;
//   fp_width += results->fp_width;
//   fp_length += results->fp_length;
//   log10_N += results->log10_N;
//   G += results->G;
//   delta2l += results->delta2l;
//   delta2m += results->delta2m;
//   log10_lambdac += results->log10_lambdac;
//   log10_lambdad += results->log10_lambdad;
//   sc_disp += results->sc_disp;
//   sc_width += results->sc_width;
//   sc_length += results->sc_length;
// }

// void VSAnalysisStage3Datum::Scope::load(VSOctaveH5ReaderStruct* reader)
// {
//   fp_dist.load(reader->readStruct("fp_dist"));
//   fp_disp.load(reader->readStruct("fp_disp"));
//   fp_width.load(reader->readStruct("fp_width"));
//   fp_length.load(reader->readStruct("fp_length"));
//   log10_N.load(reader->readStruct("log10_N"));
//   G.load(reader->readStruct("G"));
//   delta2l.load(reader->readStruct("delta2l"));
//   delta2m.load(reader->readStruct("delta2m"));
//   log10_lambdac.load(reader->readStruct("log10_lambdac"));
//   log10_lambdad.load(reader->readStruct("log10_lambdad"));
//   sc_disp.load(reader->readStruct("sc_disp"));
//   sc_width.load(reader->readStruct("sc_width"));
//   sc_length.load(reader->readStruct("sc_length"));
// }

// void VSAnalysisStage3Datum::Scope::save(VSOctaveH5WriterStruct* writer, 
// 					bool write_hists) const
// {
//   fp_dist.save(writer->writeStruct("fp_dist"));
//   fp_disp.save(writer->writeStruct("fp_disp"));
//   fp_width.save(writer->writeStruct("fp_width"));
//   fp_length.save(writer->writeStruct("fp_length"));
//   log10_N.save(writer->writeStruct("log10_N"));
//   G.save(writer->writeStruct("G"));
//   delta2l.save(writer->writeStruct("delta2l"));
//   delta2m.save(writer->writeStruct("delta2m"));
//   log10_lambdac.save(writer->writeStruct("log10_lambdac"));
//   log10_lambdad.save(writer->writeStruct("log10_lambdad"));
//   sc_disp.save(writer->writeStruct("sc_disp"));
//   sc_width.save(writer->writeStruct("sc_width"));
//   sc_length.save(writer->writeStruct("sc_length"));
// }

VSAnalysisStage3OnOffHistDatum::VSAnalysisStage3OnOffHistDatum():
  msc_width("msc_width",0.1,-3,3),
  msc_length("msc_length",0.1,-3,3),
  msc_disp("msc_disp",0.1,-3,3),
  log10_N2("log10_N2",0.1,1,6)
{

}

void VSAnalysisStage3OnOffHistDatum::calcExcess(double alpha)
{
  msc_width.calcExcess(alpha);
  msc_length.calcExcess(alpha);
  msc_disp.calcExcess(alpha);
  log10_N2.calcExcess(alpha);
}

void VSAnalysisStage3OnOffHistDatum::accumulate(const VSEventArrayDatum& event,
						double wt)
{
  msc_width.accumulate(event.msc_width,wt);
  msc_length.accumulate(event.msc_length,wt);
  msc_disp.accumulate(event.msc_disp,wt);
  log10_N2.accumulate(log10(event.N2),wt);
}

void VSAnalysisStage3OnOffHistDatum::
accumulateOn(const VSEventArrayDatum& event, double wt)
{
  msc_width.accumulateOn(event.msc_width,wt);
  msc_length.accumulateOn(event.msc_length,wt);
  msc_disp.accumulateOn(event.msc_disp,wt);
  log10_N2.accumulateOn(log10(event.N2),wt);
}

void VSAnalysisStage3OnOffHistDatum::
accumulateOff(const VSEventArrayDatum& event, double wt)
{
  msc_width.accumulateOff(event.msc_width,wt);
  msc_length.accumulateOff(event.msc_length,wt);
  msc_disp.accumulateOff(event.msc_disp,wt);
  log10_N2.accumulateOff(log10(event.N2),wt);
}

VSAnalysisStage3OnOffHistDatum& VSAnalysisStage3OnOffHistDatum::operator+= 
(const VSAnalysisStage3OnOffHistDatum& o)
{
  msc_width += o.msc_width;
  msc_length += o.msc_length;
  msc_disp += o.msc_disp;
  log10_N2 += o.log10_N2;

  return *this;
}

bool VSAnalysisStage3OnOffHistDatum::load(VSOctaveH5ReaderStruct* reader)
{
  msc_width.load(reader);
  msc_length.load(reader);
  msc_disp.load(reader);
  log10_N2.load(reader);

  return true;
}

void VSAnalysisStage3OnOffHistDatum::save(VSOctaveH5WriterStruct* writer) const
{
  msc_width.save(writer);
  msc_length.save(writer);
  msc_disp.save(writer);
  log10_N2.save(writer);
}

// ============================================================================
// VSAcceptanceData
// ============================================================================
VSAcceptanceData::VSAcceptanceData(double offset_max,
				   const std::string& model,
				   const VSAAlgebra::VecND& param,
				   const VSAAlgebra::MatrixND& cov): 
  model(model),
  chi2(),
  param(param),
  param_err(),
  param_cov(cov),
  fov_acceptance_hist(0.02,-(offset_max+0.3),offset_max+0.3),
  fov_bkgnd_hist(0.02,-(offset_max+0.3),offset_max+0.3),
  fovr2_acceptance_hist(),
  fovr2_bkgnd_hist(),
  sky_bkgnd_hist(),
  th2_bkgnd_hist(0.002,0.,0.15),
  lnl_hist(), lnl_diff_hist()
{ 
  param_err.resize(param.ndim());
  for(unsigned ip = 0; ip < param.ndim(); ip++)
    param_err(ip) = sqrt(param_cov(ip,ip));
}


VSAcceptanceData* VSAcceptanceData::clone()
{
  return new VSAcceptanceData(*this);
}

bool VSAcceptanceData::load(VSOctaveH5ReaderStruct* reader)
{
  return true;
}

void VSAcceptanceData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeString("model",model);
  writer->writeScalar("chi2",chi2);
  param.save(writer,"param");
  param_err.save(writer,"param_err");
  param_cov.save(writer,"param_cov");
  fov_acceptance_hist.save(writer->writeStruct("fov_acceptance_hist"));
  fov_bkgnd_hist.save(writer->writeStruct("fov_bkgnd_hist"));
  fovr2_acceptance_hist.save(writer->writeStruct("fovr2_acceptance_hist"));
  fovr2_bkgnd_hist.save(writer->writeStruct("fovr2_bkgnd_hist"));
  sky_bkgnd_hist.save(writer->writeStruct("sky_bkgnd_hist"));
  th2_bkgnd_hist.save(writer->writeStruct("th2_bkgnd_hist"));

  lnl_hist.save(writer->writeStruct("lnl_hist"));
  lnl_diff_hist.save(writer->writeStruct("lnl_diff_hist"));
}

// ============================================================================
// VSAnalysisStage3DataBase
// ============================================================================
VSAnalysisStage3DataBase::
VSAnalysisStage3DataBase(const SEphem::SphericalCoords& origin_radec,
			 double offset_max, double bin_size,
			 double log10_egy_binsize, double log10_egy_min,
			 double log10_egy_max):
  origin_ra_rad(origin_radec.longitudeRad()),
  origin_dec_rad(origin_radec.latitudeRad()),
  origin_ra_hms(origin_radec.longitude().hmsString(1)),
  origin_dec_dms(origin_radec.latitude().dmsString(1)), 
  sky_counts_hist(), 
  sky_livetime_hist(),
  cam_counts_hist(bin_size,-(offset_max+0.3),offset_max+0.3), 
  fov_counts_hist(bin_size,-(offset_max+0.3),offset_max+0.3),
  fovr2_counts_hist(0.05,0,std::pow(offset_max,2)),
  th2_on_counts_hist(0.002,0,0.15),
  th2_off_counts_hist(0.002,0,0.15),
  egy_on_hist(0.0625, -1.34375, 2.15625),
  egy_off_hist(0.0625, -1.34375, 2.15625),
  egy_ring_hist(0.0625, -1.34375, 2.15625),
  egy_th2_on_hist(0.0625, -1.34375, 2.15625,0.002,0,0.15),
  egy_th2_off_hist(0.0625, -1.34375, 2.15625,0.002,0,0.15),
  egy_on_np_hist(log10_egy_binsize, log10_egy_min, log10_egy_max),
  egy_off_np_hist(log10_egy_binsize, log10_egy_min, log10_egy_max),
  egy_ring_np_hist(log10_egy_binsize, log10_egy_min, log10_egy_max),
  egy_th2_on_np_hist(log10_egy_binsize, log10_egy_min, log10_egy_max,
  		     0.002,0,0.15),
  egy_th2_off_np_hist(log10_egy_binsize, log10_egy_min, log10_egy_max,
  		      0.002,0,0.15),
  log10_N2_hist(0.1,1,6),
  msc_width_hist(0.1,-3,3),
  msc_length_hist(0.1,-3,3),
  msc_disp_hist(0.1,-3,3),
  triggered_hist(1,0,16),
  has_image_hist(1,0,16),
  used_in_reconstruction_hist(1,0,16),
  hist_datum(),
  int_effarea_nspace(), 
  int_effarea_aperture_nspace(),
  int_effarea_psf_nspace(), 
  effarea_nspace(),
  m_origin_radec(origin_radec)
{
  int_effarea_nspace.clear(1);
}

VSAnalysisStage3DataBase::
VSAnalysisStage3DataBase(const SEphem::SphericalCoords& origin_radec):
  origin_ra_rad(origin_radec.longitudeRad()),
  origin_dec_rad(origin_radec.latitudeRad()),
  origin_ra_hms(origin_radec.longitude().hmsString(1)),
  origin_dec_dms(origin_radec.latitude().dmsString(1)), 
  sky_counts_hist(), 
  sky_livetime_hist(),
  cam_counts_hist(),
  fov_counts_hist(),
  fovr2_counts_hist(),
  th2_on_counts_hist(0.002,0,0.15),
  th2_off_counts_hist(0.002,0,0.15),
  egy_on_hist(),
  egy_off_hist(),
  egy_ring_hist(),
  egy_th2_on_hist(),
  egy_th2_off_hist(),
  egy_on_np_hist(),
  egy_off_np_hist(),
  egy_ring_np_hist(),
  egy_th2_on_np_hist(),
  egy_th2_off_np_hist(),
  log10_N2_hist(0.1,1,6),
  msc_width_hist(0.1,-3,3),
  msc_length_hist(0.1,-3,3),
  msc_disp_hist(0.1,-3,3),
  triggered_hist(1,0,16),
  has_image_hist(1,0,16),
  used_in_reconstruction_hist(1,0,16),
  hist_datum(),
  int_effarea_nspace(), 
  int_effarea_aperture_nspace(),
  int_effarea_psf_nspace(), 
  effarea_nspace(),
  m_origin_radec(origin_radec)
{

}

VSAnalysisStage3DataBase::VSAnalysisStage3DataBase():
  origin_ra_rad(),
  origin_dec_rad(),
  origin_ra_hms(),
  origin_dec_dms(),
  sky_counts_hist(), 
  sky_livetime_hist(),
  cam_counts_hist(),
  fov_counts_hist(),
  fovr2_counts_hist(),
  th2_on_counts_hist(0.002,0,0.15),
  th2_off_counts_hist(0.002,0,0.15),
  egy_on_hist(),
  egy_off_hist(),
  egy_ring_hist(),
  egy_th2_on_hist(),
  egy_th2_off_hist(),
  egy_on_np_hist(),
  egy_off_np_hist(),
  egy_ring_np_hist(),
  egy_th2_on_np_hist(),
  egy_th2_off_np_hist(),
  log10_N2_hist(0.1,1,6),
  msc_width_hist(0.1,-3,3),
  msc_length_hist(0.1,-3,3),
  msc_disp_hist(0.1,-3,3),
  triggered_hist(1,0,16),
  has_image_hist(1,0,16),
  used_in_reconstruction_hist(1,0,16),
  hist_datum(),
  int_effarea_nspace(), 
  int_effarea_aperture_nspace(),
  int_effarea_psf_nspace(), 
  effarea_nspace(),
  m_origin_radec()
{

}

void VSAnalysisStage3DataBase::accumulate(const VSEventArrayDatum& event,
					  double wt)
{
  hist_datum.accumulate(event,wt);
}

void VSAnalysisStage3DataBase::accumulateOn(const VSEventArrayDatum& event,
					    double wt)
{
  hist_datum.accumulateOn(event,wt);
}

void VSAnalysisStage3DataBase::accumulateOff(const VSEventArrayDatum& event,
					     double wt)
{
  hist_datum.accumulateOff(event,wt);
}

void VSAnalysisStage3DataBase::load(VSOctaveH5ReaderStruct* reader)
{
  reader->readCompositeHere(*this);

  m_origin_radec = 
    SEphem::SphericalCoords::makeLatLong(origin_dec_rad,origin_ra_rad);

  // Distributions in Sky/FoV/Camera Coordinate Systems -----------------------
  sky_counts_hist.load(reader->readStruct("sky_counts_hist"));
  sky_livetime_hist.load(reader->readStruct("sky_livetime_hist"));
  cam_counts_hist.load(reader->readStruct("cam_counts_hist"));
  fov_counts_hist.load(reader->readStruct("fov_counts_hist"));
  
  fovr2_counts_hist.load(reader->readStruct("fovr2_counts_hist"));
  th2_on_counts_hist.load(reader->readStruct("th2_on_counts_hist"));
  th2_off_counts_hist.load(reader->readStruct("th2_off_counts_hist"));
  
  VSNSpaceOctaveH5IO io;
  io.readHistogram(reader->readStruct("int_effarea_nspace"),
		   int_effarea_nspace);
  io.readHistogram(reader->readStruct("int_effarea_aperture_nspace"),
		   int_effarea_aperture_nspace);
  io.readHistogram(reader->readStruct("int_effarea_psf_nspace"),
		   int_effarea_psf_nspace);

  // Energy Histograms --------------------------------------------------------
  egy_on_hist.load(reader->readStruct("egy_on_hist"));
  egy_off_hist.load(reader->readStruct("egy_off_hist"));
  egy_ring_hist.load(reader->readStruct("egy_ring_hist"));

  egy_th2_on_hist.load(reader->readStruct("egy_th2_on_hist"));
  egy_th2_off_hist.load(reader->readStruct("egy_th2_off_hist"));

  egy_on_np_hist.load(reader->readStruct("egy_on_np_hist"));
  egy_off_np_hist.load(reader->readStruct("egy_off_np_hist"));
  egy_ring_np_hist.load(reader->readStruct("egy_ring_np_hist"));

  egy_th2_on_np_hist.load(reader->readStruct("egy_th2_on_np_hist"));
  egy_th2_off_np_hist.load(reader->readStruct("egy_th2_off_np_hist"));

  // Parameter Histograms -----------------------------------------------------
  log10_N2_hist.load(reader->readStruct("log10_N2_hist"));
  msc_width_hist.load(reader->readStruct("msc_width_hist"));
  msc_length_hist.load(reader->readStruct("msc_length_hist"));
  msc_disp_hist.load(reader->readStruct("msc_disp_hist"));

  triggered_hist.load(reader->readStruct("triggered_hist"));
  has_image_hist.load(reader->readStruct("has_image_hist"));
  used_in_reconstruction_hist.
    load(reader->readStruct("used_in_reconstruction_hist"));

  hist_datum.load(reader->readStruct("hist"));
}

void VSAnalysisStage3DataBase::save(VSOctaveH5WriterStruct* writer, 
				    bool write_hists) const
{
  writer->writeCompositeHere(*this);

  sky_counts_hist.save(writer->writeStruct("sky_counts_hist"));

  int_effarea_nspace.save(writer->writeStruct("int_effarea_nspace"));
  int_effarea_aperture_nspace.
    save(writer->writeStruct("int_effarea_aperture_nspace"));
  int_effarea_psf_nspace.save(writer->writeStruct("int_effarea_psf_nspace"));

  if(write_hists)
    {
      // Distributions in Sky/FoV/Camera Coordinate Systems -------------------
      sky_livetime_hist.save(writer->writeStruct("sky_livetime_hist"));
      cam_counts_hist.save(writer->writeStruct("cam_counts_hist"));
      fov_counts_hist.save(writer->writeStruct("fov_counts_hist"));
      
      fovr2_counts_hist.save(writer->writeStruct("fovr2_counts_hist"));
      th2_on_counts_hist.save(writer->writeStruct("th2_on_counts_hist"));
      th2_off_counts_hist.save(writer->writeStruct("th2_off_counts_hist"));
    }

  // Energy Histograms --------------------------------------------------------
  egy_on_hist.save(writer->writeStruct("egy_on_hist"));
  egy_off_hist.save(writer->writeStruct("egy_off_hist"));
  egy_ring_hist.save(writer->writeStruct("egy_ring_hist"));

  egy_th2_on_hist.save(writer->writeStruct("egy_th2_on_hist"));
  egy_th2_off_hist.save(writer->writeStruct("egy_th2_off_hist"));

  egy_on_np_hist.save(writer->writeStruct("egy_on_np_hist"));
  egy_off_np_hist.save(writer->writeStruct("egy_off_np_hist"));
  egy_ring_np_hist.save(writer->writeStruct("egy_ring_np_hist"));

  egy_th2_on_np_hist.save(writer->writeStruct("egy_th2_on_np_hist"));
  egy_th2_off_np_hist.save(writer->writeStruct("egy_th2_off_np_hist"));

  // Parameter Histograms -----------------------------------------------------
  log10_N2_hist.save(writer->writeStruct("log10_N2_hist"));
  msc_width_hist.save(writer->writeStruct("msc_width_hist"));
  msc_length_hist.save(writer->writeStruct("msc_length_hist"));
  msc_disp_hist.save(writer->writeStruct("msc_disp_hist"));

  triggered_hist.save(writer->writeStruct("triggered_hist"));
  has_image_hist.save(writer->writeStruct("has_image_hist"));
  used_in_reconstruction_hist.
    save(writer->writeStruct("used_in_reconstruction_hist"));

  hist_datum.save(writer->writeStruct("hist"));
}

// ============================================================================
// VSAnalysisStage3Data
// ============================================================================
VSAnalysisStage3Data::VSAnalysisStage3Data():
  VSAnalysisStage3DataBase(),
  src_name(),
  src_ra_rad(),
  src_dec_rad(),
  src_ra_hms(),
  src_dec_dms(),  
  pangle_mean_deg(),
  pangle_rms_deg(),
  livetime_ticks(),
  elaptime_ticks(),
  livetime_min(),
  elaptime_min(),
  src_radec(),
  ptg_radec(),
  m_ptg_xy(),
  m_acceptance_data(),
  m_run_data()
{
  m_acceptance_data = new VSAcceptanceData;
}

VSAnalysisStage3Data::
VSAnalysisStage3Data(const std::string& src_name,
		     const SEphem::SphericalCoords& src_radec,
		     const SEphem::SphericalCoords& origin_radec):
  VSAnalysisStage3DataBase(origin_radec),
  src_name(src_name),
  src_ra_rad(src_radec.longitudeRad()),
  src_dec_rad(src_radec.latitudeRad()),
  src_ra_hms(src_radec.longitude().hmsString(1)),
  src_dec_dms(src_radec.latitude().dmsString(1)),  
  pangle_mean_deg(),
  pangle_rms_deg(),
  livetime_ticks(),
  elaptime_ticks(),
  livetime_min(),
  elaptime_min(),
  src_radec(src_radec),
  ptg_radec(),
  m_ptg_xy(),
  m_acceptance_data(),
  m_run_data()
{
  m_acceptance_data = new VSAcceptanceData;
}

VSAnalysisStage3Data::
VSAnalysisStage3Data(const std::string& src_name,
		     const SEphem::SphericalCoords& src_radec,
		     const std::vector<SEphem::SphericalCoords>& ptg_radec,
		     const SEphem::SphericalCoords& origin_radec):
  VSAnalysisStage3DataBase(origin_radec),
  src_name(src_name),
  src_ra_rad(src_radec.longitudeRad()),
  src_dec_rad(src_radec.latitudeRad()),
  src_ra_hms(src_radec.longitude().hmsString(1)),
  src_dec_dms(src_radec.latitude().dmsString(1)),  
  pangle_mean_deg(),
  pangle_rms_deg(),
  livetime_ticks(),
  elaptime_ticks(),
  livetime_min(),
  elaptime_min(),
  src_radec(src_radec),
  ptg_radec(ptg_radec),
  m_ptg_xy(),
  m_acceptance_data(),
  m_run_data()
{
  const unsigned nptg = ptg_radec.size();
  m_ptg_xy.resize(nptg);
  for(unsigned iptg = 0; iptg < nptg; iptg++)
    {
      std::pair< Angle, Angle > xy;
      Astro::raDecToXY(ptg_radec[iptg],m_origin_radec,xy);  
      m_ptg_xy[iptg].set(xy.first.degPM(),xy.second.degPM());
    }
  m_acceptance_data = new VSAcceptanceData;
}

VSAnalysisStage3Data::~VSAnalysisStage3Data()
{
  delete m_acceptance_data;
}

void VSAnalysisStage3Data::addData(const std::vector< RunData* >& run_data)
{
  const unsigned nruns = run_data.size();
  for(unsigned irun = 0; irun < nruns; irun++) addData(*run_data[irun]);
}

void VSAnalysisStage3Data::addData(const RunData& data)
{
  m_run_data.push_back(data);

  unsigned iptg = 0;
  for(iptg = 0; iptg < m_ptg_xy.size(); iptg++)
    if(data.ptg_xy() == m_ptg_xy[iptg]) break;	

  if(iptg == m_ptg_xy.size()) m_ptg_xy.push_back(data.ptg_xy());

  hist_datum += data.hist_datum;

  sky_counts_hist.merge(data.sky_counts_hist);
  sky_livetime_hist.merge(data.sky_livetime_hist);
  cam_counts_hist.merge(data.cam_counts_hist);
  fov_counts_hist.merge(data.fov_counts_hist);

  fovr2_counts_hist.merge(data.fovr2_counts_hist);
  th2_on_counts_hist.merge(data.th2_on_counts_hist);
  th2_off_counts_hist.merge(data.th2_off_counts_hist);

  egy_on_hist.merge(data.egy_on_hist);
  egy_off_hist.merge(data.egy_off_hist);
  egy_ring_hist.merge(data.egy_ring_hist);

  egy_th2_on_hist.merge(data.egy_th2_on_hist);
  egy_th2_off_hist.merge(data.egy_th2_off_hist);

  egy_on_np_hist.merge(data.egy_on_np_hist);
  egy_off_np_hist.merge(data.egy_off_np_hist);
  egy_ring_np_hist.merge(data.egy_ring_np_hist);

  egy_th2_on_np_hist.merge(data.egy_th2_on_np_hist);
  egy_th2_off_np_hist.merge(data.egy_th2_off_np_hist);

  log10_N2_hist.merge(data.log10_N2_hist);
  msc_width_hist.merge(data.msc_width_hist);
  msc_length_hist.merge(data.msc_length_hist);
  msc_disp_hist.merge(data.msc_disp_hist);

  triggered_hist.merge(data.triggered_hist);
  has_image_hist.merge(data.has_image_hist);
  used_in_reconstruction_hist.merge(data.used_in_reconstruction_hist);

  livetime_ticks += data.livetime_ticks;
  elaptime_ticks += data.elaptime_ticks;
  livetime_min += data.livetime_min;
  elaptime_min += data.elaptime_min;

  // Average Parallactic Angle ------------------------------------------------
  VSSimpleStat2<double> pangle_stat;

  int_effarea_nspace = 
    VSNSpace(m_run_data.front().int_effarea_nspace.space());
  int_effarea_aperture_nspace = 
    VSNSpace(m_run_data.front().int_effarea_aperture_nspace.space());
  int_effarea_psf_nspace = 
    VSNSpace(m_run_data.front().int_effarea_psf_nspace.space());

  for(unsigned irun = 0; irun < m_run_data.size(); irun++)
    {
      double livetime = m_run_data[irun].livetime_min;
      int_effarea_nspace += 
	m_run_data[irun].int_effarea_nspace*livetime;
      int_effarea_aperture_nspace += 
	m_run_data[irun].int_effarea_aperture_nspace*livetime;
      int_effarea_psf_nspace += 
	m_run_data[irun].int_effarea_psf_nspace*livetime;

      pangle_stat.accumulate(m_run_data[irun].pangle_mean_deg);
    }

  pangle_mean_deg = pangle_stat.mean();
  pangle_rms_deg = pangle_stat.dev();

  int_effarea_nspace *= (1./livetime_min); 
  int_effarea_aperture_nspace *= (1./livetime_min); 
  int_effarea_psf_nspace *= (1./livetime_min); 
}

void VSAnalysisStage3Data::load(VSOctaveH5ReaderStruct* reader)
{
  VSAnalysisStage3DataBase::load(reader);
  reader->readCompositeHere(*this);

  src_radec = SEphem::SphericalCoords::makeLatLong(src_dec_rad,src_ra_rad);


  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("run");
  vsassert(wc);
  const unsigned nrun = wc->dimensions();
  m_run_data.resize(nrun);
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      VSOctaveH5ReaderStruct* ws = wc->readStruct(irun);
      vsassert(ws);
      m_run_data[irun].load(ws);
      delete ws;
    }

  delete wc;

  wc = reader->readCellVector("ptg_xy");
  vsassert(wc);
  const unsigned nptg = wc->dimensions();
  m_ptg_xy.resize(nptg);
  for(unsigned iptg = 0; iptg < nptg; iptg++) m_ptg_xy[iptg].load(wc,iptg);
  delete wc;

  delete m_acceptance_data;
  m_acceptance_data = new VSAcceptanceData;
}

void VSAnalysisStage3Data::save(VSOctaveH5WriterStruct* writer, 
				bool write_hists) const
{
  VSAnalysisStage3DataBase::save(writer);
  writer->writeCompositeHere(*this);

  const unsigned nrun = m_run_data.size();
  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("run", nrun);
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(irun);
      m_run_data[irun].save(ws,write_hists);
      delete ws;
    }
  delete wc;

  const unsigned nptg = m_ptg_xy.size();
  wc = writer->writeCellVector("ptg_xy", nptg);
  vsassert(wc);
  for(unsigned iptg = 0; iptg < nptg; iptg++) m_ptg_xy[iptg].save(wc,iptg);    
  delete wc;

  m_acceptance_data->save(writer->writeStruct("acceptance"));
}

void VSAnalysisStage3Data::clear()
{
  m_run_data.clear();
  sky_counts_hist.clear();
  sky_livetime_hist.clear();
  cam_counts_hist.clear();
  fov_counts_hist.clear();
}

VSAnalysisStage3Data::VSAnalysisStage3Data(const VSAnalysisStage3Data& o):
  VSAnalysisStage3DataBase(o)
{
  m_acceptance_data = o.m_acceptance_data->clone();
  copy(o);
}

VSAnalysisStage3Data& 
VSAnalysisStage3Data::operator=(const VSAnalysisStage3Data& o)
{
  delete m_acceptance_data;
  m_acceptance_data = o.m_acceptance_data->clone();
  copy(o);

  // Call base class assignment operator --------------------------------------
  VSAnalysisStage3DataBase::operator= (o);
  return *this;
}

void VSAnalysisStage3Data::copy(const VSAnalysisStage3Data& o)
{
  src_name = o.src_name;
  src_ra_rad = o.src_ra_rad;
  src_dec_rad = o.src_dec_rad;
  src_ra_hms = o.src_ra_hms;
  src_dec_dms = o.src_dec_dms;

  pangle_mean_deg = o.pangle_mean_deg;
  pangle_rms_deg = o.pangle_rms_deg;

  livetime_ticks = o.livetime_ticks;
  elaptime_ticks = o.elaptime_ticks;
  livetime_min = o.livetime_min;
  elaptime_min = o.elaptime_min;
  
  src_radec = o.src_radec;
  ptg_radec = o.ptg_radec;
  
  m_ptg_xy = o.m_ptg_xy;

  m_run_data = o.m_run_data;
}

// ============================================================================
// VSAnalysisStage3Data::RunData
// ============================================================================
VSAnalysisStage3Data::RunData::RunData(): 
  VSAnalysisStage3DataBase(),
  src_name(), 
  src_ra_rad(), src_dec_rad(), src_ra_hms(), src_dec_dms(),
  obs_ra_rad(), obs_dec_rad(), obs_ra_hms(), obs_dec_dms(),
  ptg_ra_rad(), ptg_dec_rad(), ptg_ra_hms(), ptg_dec_dms(),
  run_number(), mode(), wobble_theta_deg(),
  wobble_phi_deg(), zn_mean_deg(), zn_rms_deg(),
  az_mean_deg(), az_rms_deg(), pangle_mean_deg(), pangle_rms_deg(),
  mean_scaled_dev(),
  livetime_ticks(), elaptime_ticks(), 
  livetime_min(), elaptime_min(), 
  ring_counts(), on_counts(), off_counts(),
  m_src_xy(), m_obs_xy(), m_ptg_xy(), m_off_xy(), m_sp_off_xy(), 
  m_src_radec(), m_obs_radec(), m_ptg_radec(),
  m_irf_calc(), m_sp_calc(), m_egy_calc()
{
  m_sp_calc = new VSScaledParameterCalc;
  m_egy_calc = new VSEnergyCalcLT;
}

VSAnalysisStage3Data::RunData::
RunData(double offset_max, double bin_size,
	double log10_egy_binsize, double log10_egy_min,
	double log10_egy_max,
	const SEphem::SphericalCoords& src_radec,
	const SEphem::SphericalCoords& origin_radec,
	const VSTargetTable::Observation& obs): 
  VSAnalysisStage3DataBase(origin_radec,
			   offset_max, bin_size, log10_egy_binsize,
			   log10_egy_min, log10_egy_max),
  src_name(), 
  src_ra_rad(), src_dec_rad(), src_ra_hms(), src_dec_dms(),
  obs_ra_rad(), obs_dec_rad(), obs_ra_hms(), obs_dec_dms(),
  ptg_ra_rad(), ptg_dec_rad(), ptg_ra_hms(), ptg_dec_dms(),
  run_number(), mode(), wobble_theta_deg(),
  wobble_phi_deg(), zn_mean_deg(), zn_rms_deg(),
  az_mean_deg(), az_rms_deg(), pangle_mean_deg(), pangle_rms_deg(),
  mean_scaled_dev(),
  livetime_ticks(), elaptime_ticks(), 
  livetime_min(), elaptime_min(), 
  ring_counts(), on_counts(), off_counts(),
  m_src_xy(), m_obs_xy(), m_ptg_xy(), m_off_xy(),  m_sp_off_xy(), 
  m_src_radec(src_radec), m_obs_radec(obs.obs_radec_J2000), m_ptg_radec(),
  m_irf_calc(), m_sp_calc(), m_egy_calc()
{
  m_sp_calc = new VSScaledParameterCalc;
  m_egy_calc = new VSEnergyCalcLT;

  std::pair< Angle, Angle > xy;
  Astro::raDecToXY(m_obs_radec,m_origin_radec,xy);  
  m_obs_xy.set(xy.first.degPM(),xy.second.degPM());
  Astro::raDecToXY(m_src_radec,m_origin_radec,xy);  
  m_src_xy.set(xy.first.degPM(),xy.second.degPM());

  src_ra_rad = m_src_radec.longitudeRad();
  src_dec_rad = m_src_radec.latitudeRad();
  src_ra_hms = m_src_radec.longitude().hmsString(1);
  src_dec_dms = m_src_radec.latitude().dmsString(1);

  obs_ra_rad = m_obs_radec.longitudeRad();
  obs_dec_rad = m_obs_radec.latitudeRad();
  obs_ra_hms = m_obs_radec.longitude().hmsString(1);
  obs_dec_dms = m_obs_radec.latitude().dmsString(1);

  double x1 = (lround((m_obs_xy.x()-offset_max-0.2)/bin_size)+0.5)*bin_size;
  double x2 = (lround((m_obs_xy.x()+offset_max+0.2)/bin_size)+0.5)*bin_size;
  double y1 = (lround((m_obs_xy.y()-offset_max-0.2)/bin_size)+0.5)*bin_size;
  double y2 = (lround((m_obs_xy.y()+offset_max+0.2)/bin_size)+0.5)*bin_size;

  sky_counts_hist = 
    VSSimple2DHist<double, double>(bin_size,x1,x2,bin_size,y1,y2); 
  sky_livetime_hist =
    VSSimple2DHist<double, double>(bin_size,x1,x2,bin_size,y1,y2);
}

VSAnalysisStage3Data::RunData::~RunData()
{
  delete m_sp_calc;
  delete m_egy_calc;
}

void VSAnalysisStage3Data::RunData::
setPointing(unsigned iptg, const SEphem::SphericalCoords& radec)
{
  m_ptg_radec = radec;
  std::pair< Angle, Angle > xy;
  Astro::raDecToXY(m_ptg_radec,m_origin_radec,xy);  
  m_ptg_xy.set(xy.first.degPM(),xy.second.degPM());
}

void VSAnalysisStage3Data::RunData::load(VSOctaveH5ReaderStruct* reader)
{
  VSAnalysisStage3DataBase::load(reader);
  reader->readCompositeHere(*this);

  m_obs_xy.load(reader,"obs_xy");
  m_ptg_xy.load(reader,"ptg_xy");
  m_src_xy.load(reader,"src_xy");

  m_src_radec = SEphem::SphericalCoords::makeLatLong(src_dec_rad,src_ra_rad);
  m_obs_radec = SEphem::SphericalCoords::makeLatLong(obs_dec_rad,obs_ra_rad);
  m_ptg_radec = SEphem::SphericalCoords::makeLatLong(ptg_dec_rad,ptg_ra_rad);

  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("off_xy");
  vsassert(wc);
  const unsigned noff = wc->dimensions();
  m_off_xy.resize(noff);
  for(unsigned ioff = 0; ioff < noff; ioff++) m_off_xy[ioff].load(wc,ioff);
  delete wc;

  wc = reader->readCellVector("sp_off_xy");
  vsassert(wc);
  const unsigned sp_noff = wc->dimensions();
  m_sp_off_xy.resize(sp_noff);
  for(unsigned ioff = 0; ioff < sp_noff; ioff++) 
    m_sp_off_xy[ioff].load(wc,ioff);
  delete wc;

  m_irf_calc.load(reader->readStruct("irf"));
}

void VSAnalysisStage3Data::RunData::save(VSOctaveH5WriterStruct* writer, 
					 bool write_hists) const
{
  VSAnalysisStage3DataBase::save(writer,write_hists);
  writer->writeCompositeHere(*this);

  m_obs_xy.save(writer,"obs_xy");
  m_ptg_xy.save(writer,"ptg_xy");
  m_src_xy.save(writer,"src_xy");

  const unsigned noff = m_off_xy.size();
  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("off_xy", noff);
  vsassert(wc);
  for(unsigned ioff = 0; ioff < noff; ioff++) m_off_xy[ioff].save(wc,ioff);    
  delete wc;

  const unsigned sp_noff = m_sp_off_xy.size();
  wc = writer->writeCellVector("sp_off_xy", sp_noff);
  vsassert(wc);
  for(unsigned ioff = 0; ioff < sp_noff; ioff++) 
    m_sp_off_xy[ioff].save(wc,ioff);    
  delete wc;

  m_irf_calc.save(writer->writeStruct("irf"));
}

VSAnalysisStage3Data::RunData::RunData(const RunData& o):
  VSAnalysisStage3DataBase(o)
{
  m_sp_calc = new VSScaledParameterCalc(*o.m_sp_calc);
  m_egy_calc = new VSEnergyCalcLT(*o.m_egy_calc);

  copy(o);
}

VSAnalysisStage3Data::RunData& 
VSAnalysisStage3Data::RunData::operator=(const RunData& o)
{
  delete m_sp_calc;
  delete m_egy_calc;
  m_sp_calc = new VSScaledParameterCalc(*o.m_sp_calc);
  m_egy_calc = new VSEnergyCalcLT(*o.m_egy_calc);

  copy(o);

  // Call base class assignment operator --------------------------------------
  VSAnalysisStage3DataBase::operator= (o);
  return *this;
}

void VSAnalysisStage3Data::RunData::copy(const RunData& o)
{
  src_name = o.src_name;

  src_ra_rad = o.src_ra_rad;
  src_dec_rad = o.src_dec_rad;
  src_ra_hms = o.src_ra_hms;
  src_dec_dms = o.src_dec_dms;

  obs_ra_rad = o.obs_ra_rad;
  obs_dec_rad = o.obs_dec_rad;
  obs_ra_hms = o.obs_ra_hms;
  obs_dec_dms = o.obs_dec_dms;

  ptg_ra_rad = o.ptg_ra_rad;
  ptg_dec_rad = o.ptg_dec_rad;
  ptg_ra_hms = o.ptg_ra_hms;
  ptg_dec_dms = o.ptg_dec_dms;

  run_number = o.run_number;
  mode = o.mode;
  wobble_theta_deg = o.wobble_theta_deg;
  wobble_phi_deg = o.wobble_phi_deg;

  zn_mean_deg = o.zn_mean_deg;
  zn_rms_deg = o.zn_rms_deg;
  az_mean_deg = o.az_mean_deg;
  az_rms_deg = o.az_rms_deg;
  pangle_mean_deg = o.pangle_mean_deg;
  pangle_rms_deg = o.pangle_rms_deg;
  mean_scaled_dev = o.mean_scaled_dev;
 
  lo_event_time = o.lo_event_time;
  lo_event_time_string = o.lo_event_time_string;
  hi_event_time = o.hi_event_time;
  hi_event_time_string = o.hi_event_time_string;

  livetime_ticks = o.livetime_ticks;
  elaptime_ticks = o.elaptime_ticks;
  livetime_min = o.livetime_min;
  elaptime_min = o.elaptime_min;
  ring_counts = o.ring_counts;
  on_counts = o.on_counts;
  off_counts = o.off_counts;
  
  m_src_xy = o.m_src_xy;
  m_obs_xy = o.m_obs_xy;
  m_ptg_xy = o.m_ptg_xy;
  m_off_xy = o.m_off_xy;
  m_sp_off_xy = o.m_sp_off_xy;
  m_src_radec = o.m_src_radec;
  m_obs_radec = o.m_obs_radec;
  m_ptg_radec = o.m_ptg_radec;

  m_irf_calc = o.m_irf_calc;
}
