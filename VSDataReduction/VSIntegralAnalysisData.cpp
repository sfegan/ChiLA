//-*-mode:c++; mode:font-lock;-*-

/*! \file VSIntegralAnalysisData.cpp

  Class to store results of the stage 3 analysis

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.6 $
  \date       07/29/2007

  $Id: VSIntegralAnalysisData.cpp,v 3.6 2010/10/20 03:46:18 matthew Exp $

*/

#include <VSIntegralAnalysisData.hpp>
#include <VSRunInfoData.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSNSpace.hpp>

using namespace VERITAS;

// ============================================================================
// VSIntegralAnalysisDatum
// ============================================================================
VSIntegralAnalysisDatum::VSIntegralAnalysisDatum():
  method(), src_dec_J2000(), src_ra_J2000(), src_xy(),
  elaptime_min(), livetime_min(),
  on_counts(), on_counts_err(), on_rate(), on_rate_err(),
  off_counts(), off_counts_err(), off_rate(), off_rate_err(),
  bkgnd(), bkgnd_err(), bkgnd_rate(), bkgnd_rate_err(),
  bkgnd_density(), bkgnd_density_err(), 
  bkgnd_density_rate(), bkgnd_density_rate_err(),
  excess(), excess_err(), excess_rate(), excess_rate_err(),
  flux(), flux_err(),
  flux100(), flux100_err(),
  flux316(), flux316_err(),
  flux1000(), flux1000_err(),
  flux_eth(), flux_eth_err(),
  dfde(), dfde_err(),
  dfde100(), dfde100_err(),
  dfde316(), dfde316_err(),
  dfde1000(), dfde1000_err(),
  dfde_eth(), dfde_eth_err(),
  alpha(), significance(), sigma_sqrthr(),
  excess_rate_ul95(), excess_rate_ul99(), 
  dfde_ul95(), dfde_ul99(),
  dfde100_ul95(), dfde100_ul99(),
  dfde316_ul95(), dfde316_ul99(),
  dfde1000_ul95(), dfde1000_ul99(),
  dfde_eth_ul95(), dfde_eth_ul99(),
  e2dfde_ul95(), e2dfde_ul99(),
  e2dfde100_ul95(), e2dfde100_ul99(),
  e2dfde316_ul95(), e2dfde316_ul99(),
  e2dfde1000_ul95(), e2dfde1000_ul99(),
  e2dfde_eth_ul95(), e2dfde_eth_ul99(),
  th2_bkgnd_counts_model_hist(0.002,0.,0.15),
  th2_source_counts_model_hist(0.002,0.,0.15),
  th2_source_model_hist(0.002,0.,0.15),
  th2_counts_model_hist(0.002,0.,0.15),
  effarea(), psf(), drde(), edrde()
{

}

void VSIntegralAnalysisDatum::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);
  src_xy.save(writer,"src_xy");
  th2_bkgnd_counts_model_hist.
    save(writer->writeStruct("th2_bkgnd_counts_model_hist"));
  th2_source_counts_model_hist.
    save(writer->writeStruct("th2_source_counts_model_hist"));
  th2_source_model_hist.
    save(writer->writeStruct("th2_source_model_hist"));
  th2_counts_model_hist.save(writer->writeStruct("th2_counts_model_hist"));

  VSNSpaceOctaveH5IO io;
  io.writeHistogram(writer->writeStruct("effarea"),effarea);
  io.writeHistogram(writer->writeStruct("psf"),psf);
  io.writeHistogram(writer->writeStruct("drde"),drde);
  io.writeHistogram(writer->writeStruct("edrde"),edrde);
}

void VSIntegralAnalysisDatum::load(VSOctaveH5ReaderStruct* reader) 
{
  reader->readCompositeHere(*this);
  src_xy.load(reader,"src_xy");
  th2_bkgnd_counts_model_hist.
    load(reader->readStruct("th2_bkgnd_counts_model_hist"));
  th2_source_counts_model_hist.
    load(reader->readStruct("th2_source_counts_model_hist"));
  th2_source_model_hist.load(reader->readStruct("th2_source_model_hist"));
  th2_counts_model_hist.load(reader->readStruct("th2_counts_model_hist"));

  VSNSpaceOctaveH5IO io;
  io.readHistogram(reader->readStruct("effarea"),effarea);
  io.readHistogram(reader->readStruct("psf"),psf);
  io.readHistogram(reader->readStruct("drde"),drde);
  io.readHistogram(reader->readStruct("edrde"),edrde);
}

// ============================================================================
// VSIntegralAnalysisData
// ============================================================================
VSIntegralAnalysisData::VSIntegralAnalysisData(): 
  sky_on_counts_hist(), 
  sky_off_counts_hist(), 
  sky_bkgnd_hist(),
  sky_bkgnd_rate_hist(),
  sky_bkgnd_density_rate_hist(),
  sky_excess_hist(), 
  sky_excess_rate_hist(), 
  sky_source_hist(),
  sky_flux_hist(),
  sky_flux_ul95_hist(),
  sky_alpha_hist(), 
  sky_domega_hist(), 
  sky_significance_hist(), 
  sky_livetime_hist(),
  sky_acceptance_hist(),
  significance_hist(0.2,-10.,10.), 
  significance_excluded_hist(0.2,-10.,10.),
  src_data(),
  maxsig_data()
{

}

VSIntegralAnalysisData::~VSIntegralAnalysisData() 
{

}

void VSIntegralAnalysisData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);

  // Sky maps -------------------------------------------------------------
  sky_on_counts_hist.save(writer->writeStruct("sky_on_counts_hist"));
  sky_off_counts_hist.save(writer->writeStruct("sky_off_counts_hist"));
  sky_bkgnd_hist.save(writer->writeStruct("sky_bkgnd_hist"));
  sky_bkgnd_rate_hist.save(writer->writeStruct("sky_bkgnd_rate_hist"));
  sky_bkgnd_density_rate_hist.
    save(writer->writeStruct("sky_bkgnd_density_rate_hist"));
  sky_excess_hist.save(writer->writeStruct("sky_excess_hist"));
  sky_excess_rate_hist.save(writer->writeStruct("sky_excess_rate_hist"));
  sky_source_hist.save(writer->writeStruct("sky_source_hist"));
  sky_flux_hist.save(writer->writeStruct("sky_flux_hist"));
  sky_flux_ul95_hist.save(writer->writeStruct("sky_flux_ul95_hist"));
  sky_alpha_hist.save(writer->writeStruct("sky_alpha_hist"));
  sky_domega_hist.save(writer->writeStruct("sky_domega_hist"));
  sky_significance_hist.save(writer->writeStruct("sky_significance_hist"));
  sky_livetime_hist.save(writer->writeStruct("sky_livetime_hist"));
  sky_acceptance_hist.save(writer->writeStruct("sky_acceptance_hist"));
  
  significance_hist.save(writer->writeStruct("significance_hist"));
  significance_excluded_hist.
    save(writer->writeStruct("significance_excluded_hist"));

  src_data.save(writer);
}  

void VSIntegralAnalysisData::load(VSOctaveH5ReaderStruct* reader) 
{
  reader->readCompositeHere(*this);

  // Sky maps -----------------------------------------------------------------
  sky_on_counts_hist.load(reader->readStruct("sky_on_counts_hist"));
  sky_off_counts_hist.load(reader->readStruct("sky_off_counts_hist"));
  sky_bkgnd_hist.load(reader->readStruct("sky_bkgnd_hist"));
  sky_bkgnd_rate_hist.load(reader->readStruct("sky_bkgnd_rate_hist"));
  sky_bkgnd_density_rate_hist.
    load(reader->readStruct("sky_bkgnd_density_rate_hist"));
  sky_excess_hist.load(reader->readStruct("sky_excess_hist"));
  sky_excess_rate_hist.load(reader->readStruct("sky_excess_rate_hist"));
  sky_flux_hist.load(reader->readStruct("sky_flux_hist"));
  sky_flux_ul95_hist.load(reader->readStruct("sky_flux_ul95_hist"));
  sky_alpha_hist.load(reader->readStruct("sky_alpha_hist"));
  sky_domega_hist.load(reader->readStruct("sky_domega_hist"));
  sky_significance_hist.load(reader->readStruct("sky_significance_hist"));
  sky_livetime_hist.load(reader->readStruct("sky_livetime_hist"));
  sky_acceptance_hist.load(reader->readStruct("sky_acceptance_hist"));

  significance_hist.load(reader->readStruct("significance_hist"));
  significance_excluded_hist.
    load(reader->readStruct("significance_excluded_hist"));

  src_data.load(reader);
}

// ----------------------------------------------------------------------------
// VSIntegralAnalysisRunData
// ----------------------------------------------------------------------------
VSIntegralAnalysisRunDatum::VSIntegralAnalysisRunDatum(unsigned _run_number):
  VSIntegralAnalysisDatum(),
  run_number(_run_number), run_start_time_string(), run_stop_time_string(), 
  run_start_time_mjd(), run_stop_time_mjd(), obs_dec_J2000(), obs_ra_J2000()
{

}

VSIntegralAnalysisRunDatum::~VSIntegralAnalysisRunDatum()
{

}

void VSIntegralAnalysisRunDatum::load(VSOctaveH5ReaderStruct* reader) 
{
  reader->readCompositeHere(*this);
  VSIntegralAnalysisDatum::load(reader);
}

void VSIntegralAnalysisRunDatum::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);
  VSIntegralAnalysisDatum::save(writer);
}

