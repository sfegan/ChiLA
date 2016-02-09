#include <VSStage3SimCalc.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSSimpleStat.hpp>

using namespace VERITAS;




void VSStage3SimCalc::analyze(VSStage3SimDatum& data)
{
  // Effarea ------------------------------------------------------------------
  calc(data.egymc_triggered_hist,data.egymc_total_hist,
       data.egymc_sampling_area_hist, data.effarea_triggered_hist, 
       data.diffrate_triggered_hist);
  calc(data.egymc_reconstructed_hist,data.egymc_total_hist,
       data.egymc_sampling_area_hist, data.effarea_reconstructed_hist, 
       data.diffrate_reconstructed_hist);
  calc(data.egymc_cuts_selected_hist,data.egymc_total_hist,
       data.egymc_sampling_area_hist, data.effarea_cuts_selected_hist, 
       data.diffrate_cuts_selected_hist);

  data.egymc_fluence_hist = 
    data.egymc_total_hist/data.egymc_sampling_area_hist;

  // Find Energy Threshold ----------------------------------------------------
  data.egy_threshold_triggered = fitThreshold(data.diffrate_triggered_hist);
  data.egy_threshold_reconstructed = 
    fitThreshold(data.diffrate_reconstructed_hist);
  data.egy_threshold_cuts_selected = 
    fitThreshold(data.diffrate_cuts_selected_hist);
  data.egy_threshold_selected = fitThreshold(data.diffrate_selected_hist);

  // Calculate Various Rates --------------------------------------------------
  const double crab_flux_constant = 3.2E-7; // m^-2 s^-1 TeV^-1
  const double crab_index = -2.5;
  //  const double proton_flux_constant = 0.143; // m^-2 s^-1 TeV^-1 sr^-1
  const double proton_flux_constant = 1.2314E-3; // m^-2 s^-1 TeV^-1
  const double proton_index = -2.7;

  VSSpectrumFn* crab_spectrum = 
    new VSSpectrumFnPowerLaw(crab_index,crab_flux_constant);
  VSSpectrumFn* proton_spectrum = 
    new VSSpectrumFnPowerLaw(proton_index,proton_flux_constant);
  //  VSSpectrumFn* gamma_spectrum = m_egywt_calc->getSpectrum();
  
//   calcRate(data.gamma_rate_triggered, data.gamma_rate_triggered_err,
// 	   data.effarea_triggered_hist,gamma_spectrum);
//   calcRate(data.gamma_rate_reconstructed,data.gamma_rate_reconstructed_err,
// 	   data.effarea_reconstructed_hist,gamma_spectrum);
//   calcRate(data.gamma_rate_cuts_selected,data.gamma_rate_cuts_selected_err,
// 	   data.effarea_cuts_selected_hist,gamma_spectrum);
//   calcRate(data.gamma_rate_selected,data.gamma_rate_selected_err,
// 	   data.effarea_selected_hist,gamma_spectrum);

  calcRate(data.crab_rate_triggered,data.crab_rate_triggered_err,
	   data.effarea_triggered_hist,crab_spectrum);
  calcRate(data.crab_rate_reconstructed,data.crab_rate_reconstructed_err,
	   data.effarea_reconstructed_hist,crab_spectrum);
  calcRate(data.crab_rate_cuts_selected,data.crab_rate_cuts_selected_err,
	   data.effarea_cuts_selected_hist,crab_spectrum);
  calcRate(data.crab_rate_selected,data.crab_rate_selected_err,
	   data.effarea_selected_hist,crab_spectrum);
  
  calcRate(data.proton_rate_triggered,data.proton_rate_triggered_err,
	   data.effarea_triggered_hist,proton_spectrum);
  calcRate(data.proton_rate_reconstructed,data.proton_rate_reconstructed_err,
	   data.effarea_reconstructed_hist,proton_spectrum);
  calcRate(data.proton_rate_cuts_selected,data.proton_rate_cuts_selected_err,
	   data.effarea_cuts_selected_hist,proton_spectrum);
  calcRate(data.proton_rate_selected,data.proton_rate_selected_err,
	   data.effarea_selected_hist,proton_spectrum);

  delete crab_spectrum;
  delete proton_spectrum;

  // Energy Stuff -------------------------------------------------------------
  calcKernel(data.egymc_egy_hist,data.egy_kernel);
  calcEnergyBias(data.egymc_egy_hist,data.egymc_bias_hist,
		 data.egymc_rms_hist, data.egymc_log10_bias_hist,
		 data.egymc_log10_rms_hist);

  // PSF Fitting --------------------------------------------------------------
  
}

void VSStage3SimCalc::analyze(std::vector<VSStage3SimArrayTableDatum*>& tables)
{
  for(std::vector<VSStage3SimArrayTableDatum*>::iterator itr =
	tables.begin(); itr != tables.end(); ++itr)
    {
      calcEffarea(**itr);
    }
}

void VSStage3SimCalc::calcEffarea(VSStage3SimArrayTableDatum& d)
{
  // Compute containment radius for thetasq -----------------------------------
  containmentIntervalUpper(d.m_thetasq,0.68,d.thsq68,d.thsq68_err);
  containmentIntervalUpper(d.m_thetasq,0.90,d.thsq90,d.thsq90_err);
  containmentIntervalUpper(d.m_thetasq,0.95,d.thsq95,d.thsq95_err);

  d.th68 = sqrt(d.thsq68);
  d.th68_err = 0.5*(d.thsq68_err/d.th68);
  d.th90 = sqrt(d.thsq90);
  d.th90_err = 0.5*(d.thsq90_err/d.th90);
  d.th95 = sqrt(d.thsq95);
  d.th95_err = 0.5*(d.thsq95_err/d.th95);

  double eff_triggered         = d.ntriggered/d.ntotal;
  d.effarea_triggered     = eff_triggered*d.sampling_area;
  d.effarea_triggered_err = 
    sqrt(eff_triggered*(1-eff_triggered)/d.ntotal)*d.sampling_area;
  
  double eff_reconstructed  = d.nreconstructed/d.ntotal;
  d.effarea_reconstructed     = eff_reconstructed*d.sampling_area;
  d.effarea_reconstructed_err = 
    sqrt(eff_reconstructed*(1-eff_reconstructed)/d.ntotal)*d.sampling_area;
  
  double eff_cuts_selected  = d.ncuts_selected/d.ntotal;
  d.effarea_cuts_selected     = eff_cuts_selected*d.sampling_area;
  d.effarea_cuts_selected_err = 
    sqrt(eff_cuts_selected*(1-eff_cuts_selected)/d.ntotal)*d.sampling_area;
  
  double eff_selected  = d.nselected/d.ntotal;
  d.effarea_selected     = eff_selected*d.sampling_area;
  d.effarea_selected_err = 
    sqrt(eff_selected*(1-eff_selected)/d.ntotal)*d.sampling_area;
}

void
VSStage3SimCalc::calc(const VSLimitedErrorsHist<double, double>& h,
		      const VSLimitedErrorsHist<double, double>& htot,
		      const VSLimitedErrorsHist<double, double>& harea,
		      VSLimitedErrorsHist<double,double>& effarea_hist,
		      VSLimitedErrorsHist<double,double>& diffrate_hist)
{
  effarea_hist = htot; effarea_hist.clear();
  diffrate_hist = htot; diffrate_hist.clear();

  const double crab_flux_constant = 3.2E-7; // m^-2 s^-1 TeV^-1
  const double crab_index = 2.5;

  for(VSLimitedErrorsHist<double, double>::iterator itr = htot.begin(); 
      itr != htot.end(); ++itr)
    {
      double log10_egy = itr->center();      
      double egy_tev = std::pow(10,log10_egy);

      double flux = crab_flux_constant*pow(egy_tev,-crab_index);

      double ntot = itr->count();
      double n = h.countForVal(itr->center());

      double eff = n/ntot;
      double eff_var = eff*(1-eff)/ntot;
      double area = harea.countForVal(itr->center());
      double effarea = eff*area;
      double effarea_err = sqrt(eff_var)*area;

      effarea_hist.accumulate(log10_egy,effarea,std::pow(effarea_err,2));
      diffrate_hist.
	accumulate(log10_egy,effarea*flux,std::pow(effarea_err*flux,2));
    }
}

double VSStage3SimCalc::
fitThreshold(const VSLimitedErrorsHist<double, double>& h)
{

  VSLimitedErrorsHist<double, double>::iterator max_itr = h.begin();
  double max = max_itr->count();

  for(VSLimitedErrorsHist<double, double>::iterator itr = h.begin(); 
      itr != h.end(); ++itr)
    {
      if(itr->count() > max)
	{
	  max = itr->count();
	  max_itr = itr;
	}
    }

  VSLimitedErrorsHist<double, double>::iterator itr1 = 
    std::max(h.begin(),max_itr-2);

  VSLimitedErrorsHist<double, double>::iterator itr2 = 
    std::min(h.end()-1,max_itr+2);

  VSAMath::Data<double> data;

  while(itr1 <= itr2)
    {
      if(itr1->count() == 0 || itr1->err() <= 0) return 0;
      data.insert(VSAMath::DataPoint<double>
		  (itr1->center(),itr1->count(),itr1->err()));
      itr1++;
    }

  VSAAlgebra::MatrixND cov;
  VSAAlgebra::VecND param;

  try
    {
      VSAMath::PolyFit::fit(2,data,param,&cov);
    }
  catch(const std::string& s)
    {
      std::cerr << s << std::endl;
    }

  return std::pow(10,-0.5*param[1]/param[2]);
}

void VSStage3SimCalc::
calcRate(double& rate, double& rate_err,
	 const VSLimitedErrorsHist<double, double>& effarea_hist,
	 VSSpectrumFn* spectrum)
{
  rate = 0;
  rate_err = 0;
  double dx = (effarea_hist.begin()+1)->val() - effarea_hist.begin()->val();

  for(VSLimitedErrorsHist<double, double>::iterator itr = 
	effarea_hist.begin(); itr != effarea_hist.end(); ++itr)
    {
      double log10_egy_tev = itr->center();
      double egy_tev = std::pow(10,log10_egy_tev);

      rate += 
	egy_tev*spectrum->diffFlux(log10_egy_tev)*
	dx*log(10.)*60*itr->count();
      rate_err += 
	egy_tev*spectrum->diffFlux(log10_egy_tev)*dx*log(10.)*60*itr->err();
    }

  return;
}

void VSStage3SimCalc::calcKernel(const VSSimple2DHist<double, double>& h,
				 VSNSpace& egy_kernel)
{
  // Normalize kernel ---------------------------------------------------------
  VSNSpace tmp = egy_kernel;
  tmp.project(0);

  VSNSpace::Cell c(2);
  for(c.i[0] = 0; c.i[0] < egy_kernel.space().axes[0].nbin; c.i[0]++)
    {
      VSNSpace::Weight w = 0;
      tmp.getWeight(c.i[0],w);

      if(w==0) continue;

      for(c.i[1] = 0; c.i[1] < egy_kernel.space().axes[1].nbin; c.i[1]++)
	egy_kernel.setWeight(c,egy_kernel[c]/w);
    }
}

void VSStage3SimCalc::
calcEnergyBias(const VSSimple2DHist<double, double>& h,
	       VSLimitedErrorsHist<double, double>& egy_bias_hist,
	       VSLimitedErrorsHist<double, double>& egy_rms_hist,
	       VSLimitedErrorsHist<double, double>& egy_log10_bias_hist,
	       VSLimitedErrorsHist<double, double>& egy_log10_rms_hist)
{
  // Compute bias and rms of energy estimator ---------------------------------
  const unsigned nbinx = h.nXBins();
  const unsigned nbiny = h.nYBins();

  for(unsigned ix = 0; ix < nbinx; ix++)
    {
      double log10_emc = h.xBinToCenter(ix);
      double emc = std::pow(10,log10_emc);

      VSSimpleStat2<double,double> stat_log;
      VSSimpleStat2<double,double> stat_lin;
      for(unsigned iy = 0; iy < nbiny; iy++)
	{
	  double log10_erec = h.yBinToCenter(iy);
	  double n = h.count(ix,iy);
	  stat_log.accumulate(log10_erec,n);
	  stat_lin.accumulate(std::pow(10,log10_erec),n);
	}

      if(stat_log.count() == 0)
	continue;

      double count = stat_log.count();

      double bias_log = stat_log.mean()-log10_emc;
      double bias_log_err = stat_log.dev()/sqrt(count);
      double rms_log = stat_log.dev();
      double rms_log_err = stat_log.dev()/sqrt(2*count);

      double bias = (stat_lin.mean()-emc)/emc;
      double bias_err = stat_lin.dev()/(emc*sqrt(count));
      double rms = stat_lin.dev()/emc;
      double rms_err = sqrt(2.)*stat_lin.dev()/(emc*sqrt(count));

//       std::cout << std::setw(15) << emc 
// 		<< std::setw(15) << stat_lin.mean()
// 		<< std::setw(15) << stat_lin.dev()
// 		<< std::setw(15) << stat_lin.dev()/emc
// 		<< std::endl;


      egy_bias_hist.accumulate(log10_emc,bias,std::pow(bias_err,2));
      egy_rms_hist.accumulate(log10_emc,rms,std::pow(rms_err,2));

      egy_log10_bias_hist.accumulate(log10_emc,bias_log,
				     std::pow(bias_log_err,2));
      egy_log10_rms_hist.accumulate(log10_emc,rms_log,std::pow(rms_log_err,2));
    }

}
