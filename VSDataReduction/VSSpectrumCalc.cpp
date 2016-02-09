//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSpectrumCalc.cpp

  Class which reconstructs the true energy spectrum.

  \author     Tim Arlen                   \n
              UCLA                        \n
              arlen@astro.ucla.edu        \n

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       10/03/2008

  $Id: VSSpectrumCalc.cpp,v 3.25 2010/10/20 04:02:22 matthew Exp $

*/

#include <VSNSpaceOctaveH5IO.hpp>
#include <VSSpectrumCalc.hpp>
#include <VSAStatistics.hpp>
#include <VSANonlinearFitting.hpp>
#include <VSALinearLeastSquares.hpp>
#include <VSABracketedMonotonicRootFinder.hpp>
#include <VSAConfidenceInterval.hpp>
#include <VSDataModel.hpp>

using namespace VERITAS;

// ============================================================================
// VSSpectrumFitData
// ============================================================================
VSSpectrumFitData::VSSpectrumFitData():
  model(),
  log10_enorm(),
  param(),
  param_err(),
  param_cov(),
  dfde(),
  dfde_err(),
  e2dfde(),
  e2dfde_err(),
  dfde100(),
  dfde100_err(),
  dfde200(),
  dfde200_err(),
  dfde316(),
  dfde316_err(),
  dfde1000(),
  dfde1000_err(),
  flux100(),
  flux200(),
  flux316(),
  flux1000(),
  chi2_pearson(),
  chi2_pearson_pval(),
  chi2_ml(),
  chi2_ml_pval(),
  ndf(),
  kernel_nspace(),
  on_hist(),
  off_hist(),
  mu_on_hist(),
  mu_off_hist(),
  mu_bkgnd_hist()
{ 

}

void VSSpectrumFitData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeString("model",model);
  writer->writeScalar("log10_enorm",log10_enorm);

  param.save(writer,"param");
  param_err.save(writer,"param_err");
  param_cov.save(writer,"param_cov");

  writer->writeScalar("flux100",flux100);
  writer->writeScalar("flux200",flux200);
  writer->writeScalar("flux316",flux316);
  writer->writeScalar("flux1000",flux1000);

  writer->writeScalar("dfde",dfde);
  writer->writeScalar("dfde_err",dfde_err);
  writer->writeScalar("e2dfde",e2dfde);
  writer->writeScalar("e2dfde_err",e2dfde_err);

  writer->writeScalar("dfde100",dfde100);
  writer->writeScalar("dfde100_err",dfde100_err);
  writer->writeScalar("dfde200",dfde200);
  writer->writeScalar("dfde200_err",dfde200_err);
  writer->writeScalar("dfde316",dfde316);
  writer->writeScalar("dfde316_err",dfde316_err);
  writer->writeScalar("dfde1000",dfde1000);
  writer->writeScalar("dfde1000_err",dfde1000_err);

  writer->writeScalar("chi2_pearson",chi2_pearson);
  writer->writeScalar("chi2_pearson_pval",chi2_pearson_pval);
  writer->writeScalar("chi2_ml",chi2_ml);
  writer->writeScalar("chi2_ml_pval",chi2_ml_pval);
  writer->writeScalar("ndf",ndf);

  VSNSpaceOctaveH5IO io;
  io.writeHistogram(writer->writeStruct("kernel_nspace"),kernel_nspace);

  lhist.save(writer->writeStruct("lhist"));

  on_hist.save(writer->writeStruct("on_hist"));
  off_hist.save(writer->writeStruct("off_hist"));
  mu_on_hist.save(writer->writeStruct("mu_on_hist"));
  mu_off_hist.save(writer->writeStruct("mu_off_hist"));
  mu_bkgnd_hist.save(writer->writeStruct("mu_bkgnd_hist"));

  dfde_butt68.save(writer->writeStruct("dfde_butt68"));
  e2dfde_butt68.save(writer->writeStruct("e2dfde_butt68"));

  writer->writeStructCellVector("lcont68",lcont68);
  writer->writeStructCellVector("lcont90",lcont90);
}

// ============================================================================
// VSSpectrumData
// ============================================================================
VSSpectrumData::VSSpectrumData():
  chi2(), lambda()
{ 

}

void VSSpectrumData::load(VSOctaveH5ReaderStruct* reader)
{

}

void VSSpectrumData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeScalar("chi2",chi2);
  writer->writeScalar("lambda",lambda);

  gcc_hist.save(writer->writeStruct("gcc_hist"));

  VSNSpaceOctaveH5IO io;
  io.writeHistogram(writer->writeStruct("krn_nspace"),krn_nspace);
  
  egy_on_hist.save(writer->writeStruct("egy_on_hist"));
  egy_off_hist.save(writer->writeStruct("egy_off_hist"));
  egy_excess_hist.save(writer->writeStruct("egy_excess_hist"));
  egy_nu_hist.save(writer->writeStruct("egy_nu_hist"));
  egy_resid_hist.save(writer->writeStruct("egy_resid_hist"));
  resid_hist.save(writer->writeStruct("resid_hist"));
  dfde_hist.save(writer->writeStruct("dfde_hist"));
  edfde_hist.save(writer->writeStruct("edfde_hist"));
  e2dfde_hist.save(writer->writeStruct("e2dfde_hist"));
  flux_hist.save(writer->writeStruct("flux_hist"));
  flux_cov_hist.save(writer->writeStruct("flux_cov_hist"));
  mu_cov_hist.save(writer->writeStruct("mu_cov_hist"));
  nu_cov_hist.save(writer->writeStruct("nu_cov_hist"));
  k_hist.save(writer->writeStruct("k_hist"));
  ksq_hist.save(writer->writeStruct("ksq_hist"));
  ksq_u_hist.save(writer->writeStruct("ksq_u_hist"));
  ksq_w_hist.save(writer->writeStruct("ksq_w_hist"));
  ksq_sqw_hist.save(writer->writeStruct("ksq_sqw_hist"));
  ksq_c_hist.save(writer->writeStruct("ksq_c_hist"));
  ksq_nc_hist.save(writer->writeStruct("ksq_nc_hist"));
  ksqr_hist.save(writer->writeStruct("ksqr_hist"));
  ksqr_c_hist.save(writer->writeStruct("ksqr_c_hist"));
  ksqr_u_hist.save(writer->writeStruct("ksqr_u_hist"));
  ksqr_f_hist.save(writer->writeStruct("ksqr_f_hist"));
  flux_bias_hist.save(writer->writeStruct("flux_bias_hist"));
  flux_rbias_hist.save(writer->writeStruct("flux_rbias_hist"));

  lambda_gcc_hist.save(writer->writeStruct("lambda_gcc_hist"));
  lambda_chi2_hist.save(writer->writeStruct("lambda_chi2_hist"));
  lambda_chi2b_hist.save(writer->writeStruct("lambda_chi2b_hist"));
  lambda_chi2_ml_hist.save(writer->writeStruct("lambda_chi2_ml_hist"));

  lambda_msbias_hist.save(writer->writeStruct("lambda_msbias_hist"));
  lambda_mvar_hist.save(writer->writeStruct("lambda_mvar_hist"));
  lambda_mse_hist.save(writer->writeStruct("lambda_mse_hist"));
  chi2_mse_graph.save(writer->writeStruct("chi2_mse_graph"));
  chi2_ml_mse_graph.save(writer->writeStruct("chi2_ml_mse_graph"));
  lambda_wmse_hist.save(writer->writeStruct("lambda_wmse_hist"));
  chi2_wmse_graph.save(writer->writeStruct("chi2_wmse_graph"));
  chi2_ml_wmse_graph.save(writer->writeStruct("chi2_ml_wmse_graph"));
}

// ============================================================================
// VSSpectrumFit
// ============================================================================
VSSpectrumFit::VSSpectrumFit(double theta_cut,
			     const std::string& sp_model,
			     double log10_enorm):
  m_theta_cut(theta_cut),
  m_sp_fn()
{
  m_sp_fn = VSSpectrumFn::create(sp_model);
  m_sp_fn->setNormEnergy(log10_enorm);
}

VSSpectrumFit::~VSSpectrumFit()
{
  delete m_sp_fn;
}

void VSSpectrumFit::fit(const VSAnalysisStage3Data& stage3_data,
			VSSpectrumFitData& data)
{
  std::cout << "VSSpectrumFit::fit(): "
	    << "Calculating forward-folded spectrum" << std::endl;

  data.model = m_sp_fn->name();
  data.log10_enorm = m_sp_fn->normEnergy();

  std::vector<Data> fit_data;

  double on_counts = 0;
  double off_counts = 0;

  const unsigned nrun = stage3_data.nrun();
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      const VSAnalysisStage3Data::RunData& rd = stage3_data.run_data(irun);
      double alpha = 1./rd.sp_noff();
      Data* fd = NULL;

      for(std::vector<Data>::iterator itr = fit_data.begin(); itr !=
	    fit_data.end(); ++itr)
	{
	  if(fabs(itr->alpha-alpha) < 1E-3)
	    {
	      fd = &(*itr);
	      break;
	    }	  
	}

      if(fd == NULL)
	{
	  fit_data.resize(fit_data.size()+1);
	  fd = &fit_data.back();
	}

      //      if(fit_data.size() <= iptg) fit_data.resize(iptg+1);

      double offset = rd.src_xy().d(rd.obs_xy());
      double livetime = rd.livetime_min;
      double ebin_size = rd.egy_on_hist.binSize();
      double elo = rd.egy_on_hist.loLimit();
      double ehi = rd.egy_on_hist.hiLimit();
      fd->alpha = alpha;

      VSNSpace k = 
	rd.irf_calc().getKernel(ebin_size, elo, ehi, offset, m_theta_cut);

      k *= livetime*60.;

      on_counts += rd.egy_on_hist.sum();
      off_counts += rd.egy_off_hist.sum();

      if(fd->on.ndim() == 0)
	{
	  fd->on = VSSpectrumCalc::toVecND(rd.egy_on_hist.contentsHist());
	  fd->off = VSSpectrumCalc::toVecND(rd.egy_off_hist.contentsHist());
	  fd->krn = k;
	  data.on_hist = rd.egy_on_hist;
	  data.off_hist = rd.egy_off_hist;
	  data.kernel_nspace = k;
	}
      else
	{
	  fd->on += VSSpectrumCalc::toVecND(rd.egy_on_hist.contentsHist());
	  fd->off += 
	    VSSpectrumCalc::toVecND(rd.egy_off_hist.contentsHist());
	  fd->krn += k;
	  data.on_hist += rd.egy_on_hist;
	  data.off_hist += rd.egy_off_hist;
	  data.kernel_nspace += k;
	}
    }

  if(on_counts == 0 && off_counts == 0) return;

  m_sp_fn->setParam(0,3.2E-8);

  VSSpectrumFit::LnLFn fn(fit_data,m_sp_fn);

  VSAMath::NLFitterLM<VSSpectrumFit::LnLFn>* fitter = 
    new VSAMath::NLFitterLM<VSSpectrumFit::LnLFn>(fn);
  fitter->setScaling(0,1E9);
  
  fitter->setLoBound(0,0);
  fitter->setLoBound(1,0);
  fitter->setHiBound(1,5);

  try
    {
      fitter->setTolerance(0.001);
      fitter->initialize();  
      std::cout << "VSSpectrumFit::fit(): Scanning" << std::endl;
      fitter->scan(0,1E-8,1E-7,1E-8,1,1.4,5.,0.2);
      std::cout << "VSSpectrumFit::fit(): Fitting" << std::endl;
      fitter->fit();

    }
  catch(const std::exception& e)
    {
      std::cout << e.what() << std::endl;
    }

  VSAAlgebra::VecND param = fitter->param();
  VSAAlgebra::VecND err = fitter->err();
  VSAAlgebra::MatrixND cov = fitter->cov();

  for(unsigned ip = 0; ip < fitter->nparam(); ip++)
    std::cout << std::setw(20) << param(ip) << " +/- "
	      << std::setw(20) << err(ip)
	      << std::endl;

  fn.gof(param,data.chi2_pearson,data.chi2_ml,data.ndf);

  data.chi2_ml_pval = VSAMath::chi2p(data.chi2_ml,data.ndf);
  data.chi2_pearson_pval = VSAMath::chi2p(data.chi2_pearson,data.ndf);

  std::cout << "VSSpectrumFit::fit(): CHI2/NDF " 
	    << std::setw(15) << data.chi2_ml 
	    << "/" << std::setw(5) << data.ndf 
	    << " (p=" << std::setw(15) << data.chi2_ml_pval << ")"
	    << std::endl;

  VSAAlgebra::VecND mu_on;
  VSAAlgebra::VecND mu_off;
  VSAAlgebra::VecND mu_bkgnd;

  fn.mu(param,mu_on,mu_off,mu_bkgnd);

  data.mu_on_hist = data.on_hist;
  data.mu_on_hist.clear();
  data.mu_on_hist.fill(0);

  data.mu_off_hist = data.on_hist;
  data.mu_off_hist.clear();
  data.mu_off_hist.fill(0);

  data.mu_bkgnd_hist = data.on_hist;
  data.mu_bkgnd_hist.clear();
  data.mu_bkgnd_hist.fill(0);

  for(unsigned iegy = 0; iegy < mu_on.ndim(); iegy++)
    {
      data.mu_on_hist.setBin(iegy,mu_on(iegy),0);
      data.mu_off_hist.setBin(iegy,mu_off(iegy),0);
      data.mu_bkgnd_hist.setBin(iegy,mu_bkgnd(iegy),0);
    }

  data.param = param;
  data.param_err = err;
  data.param_cov = cov;

  data.dfde = data.param(0);
  data.dfde_err = data.param_err(0);

  data.e2dfde = std::pow(10,2*m_sp_fn->normEnergy())*data.param(0);
  data.e2dfde_err = std::pow(10,2*m_sp_fn->normEnergy())*data.param_err(0);

  data.flux100  = m_sp_fn->integralFlux(-1,3,data.param);
  data.flux200  = m_sp_fn->integralFlux(-0.69897,3,data.param);
  data.flux316  = m_sp_fn->integralFlux(-0.5,3,data.param);
  data.flux1000 = m_sp_fn->integralFlux(0.,3,data.param);


  VSAAlgebra::VecND dfde100_da;
  VSAAlgebra::VecND dfde200_da;
  VSAAlgebra::VecND dfde316_da;
  VSAAlgebra::VecND dfde1000_da;

  m_sp_fn->dyda(-1,param,dfde100_da);
  m_sp_fn->dyda(-0.69897,param,dfde200_da);
  m_sp_fn->dyda(-0.5,param,dfde316_da);
  m_sp_fn->dyda(0.0,param,dfde1000_da);

  double dfde100_var = dfde100_da*(cov*dfde100_da);
  double dfde200_var = dfde200_da*(cov*dfde200_da);
  double dfde316_var = dfde316_da*(cov*dfde316_da);
  double dfde1000_var = dfde1000_da*(cov*dfde1000_da);

  data.dfde100  = m_sp_fn->val(-1,data.param);
  data.dfde100_err  = sqrt(dfde100_var);
  data.dfde200  = m_sp_fn->val(-0.69897,data.param);
  data.dfde200_err  = sqrt(dfde200_var);
  data.dfde316  = m_sp_fn->val(-0.5,data.param);
  data.dfde316_err  = sqrt(dfde316_var);
  data.dfde1000 = m_sp_fn->val(0.,data.param);
  data.dfde1000_err  = sqrt(dfde1000_var);


  std::vector< std::pair<double,double> > dfde_hilo;
  std::vector< std::pair<double,double> > dfde_hi;
  std::vector< std::pair<double,double> > dfde_lo;

  std::vector< std::pair<double,double> > e2dfde_hilo;
  std::vector< std::pair<double,double> > e2dfde_hi;
  std::vector< std::pair<double,double> > e2dfde_lo;

  for(double e = -1; e <= 1.01; e += 0.05)
    {
      VSAAlgebra::VecND dfde_da;
      m_sp_fn->dyda(e,param,dfde_da);
      double err = sqrt(dfde_da*cov*dfde_da);

      double ylo = m_sp_fn->val(e,param)-err;
      double yhi = m_sp_fn->val(e,param)+err;

      dfde_hi.push_back(std::make_pair(e,yhi));
      dfde_lo.push_back(std::make_pair(e,ylo));

      e2dfde_hi.push_back(std::make_pair(e,yhi*std::pow(10,2*e)));
      e2dfde_lo.push_back(std::make_pair(e,ylo*std::pow(10,2*e)));
    }

  dfde_hilo.insert(dfde_hilo.end(),dfde_lo.begin(),dfde_lo.end());
  std::reverse(dfde_hi.begin(),dfde_hi.end());
  dfde_hilo.insert(dfde_hilo.end(),dfde_hi.begin(),dfde_hi.end());
  dfde_hilo.insert(dfde_hilo.end(),dfde_lo.begin(),dfde_lo.begin()+1);

  e2dfde_hilo.insert(e2dfde_hilo.end(),e2dfde_lo.begin(),e2dfde_lo.end());
  std::reverse(e2dfde_hi.begin(),e2dfde_hi.end());
  e2dfde_hilo.insert(e2dfde_hilo.end(),e2dfde_hi.begin(),e2dfde_hi.end());
  e2dfde_hilo.insert(e2dfde_hilo.end(),e2dfde_lo.begin(),e2dfde_lo.begin()+1);

  data.dfde_butt68 = VSSimpleGraph<double,double>(dfde_hilo);
  data.e2dfde_butt68 = VSSimpleGraph<double,double>(e2dfde_hilo);

  double p0_step = 0.01;
  if(err(0) > 0) p0_step = 0.1*err(0);

  data.lhist = 
     VSSimple2DHist<double,double>(p0_step,
				   0, param(0)+100*p0_step,
				   0.1,0.,5.);

  for(VSSimple2DHist<double,double>::iterator itr = data.lhist.begin();
      itr != data.lhist.end(); ++itr)
    {
      VSAAlgebra::VecND a(2);

      a(0) = itr->x();
      a(1) = itr->y();

      data.lhist.setBin(itr->bin(),fn.val(a));
    }

  const unsigned np = fn.nparm();

  VSAAlgebra::VecND lo(np);
  VSAAlgebra::VecND hi(np);


  // VSAMath::ConfidenceInterval<VSSpectrumFit::LnLFn> ci(fn);

  // std::vector< std::pair<double,double> > xy_test;

  // ci.limits(0.68,param,0,1E-8,1,0.1,xy_test);

  // data.lcont68.push_back(VSSimpleGraph<double,double>(xy_test));

  for(unsigned ip1 = 1; ip1 < np; ip1++)
    {
      std::vector< std::pair<double,double> > xy68;
      std::vector< std::pair<double,double> > xy90;

      const unsigned ip2 = 0;

      try
	{
	  fitter->generateContour(2.3,ip1,0.,5.,ip2,0.,0.,xy68);
	  fitter->generateContour(4.61,ip1,0.,5.,ip2,0.,0.,xy90);
	}
      catch(const VSAMath::FunctionNotMonotonic& x)
	{
	  std::cout << "Function not monotonic" << std::endl;
	}
      catch(const std::exception& e)
	{
	  std::cout << e.what() << std::endl;
	}


      data.lcont68.push_back(VSSimpleGraph<double,double>(xy68));
      data.lcont90.push_back(VSSimpleGraph<double,double>(xy90));

    }

  delete fitter;
}

// ============================================================================
// VSSpectrumFit::LnLFn
// ============================================================================
VSSpectrumFit::LnLFn::LnLFn(): m_data(), m_lny(), m_negy(), m_sp_fn()
{

}

VSSpectrumFit::LnLFn::LnLFn(const std::vector<Data>& data,
			    VSSpectrumFn* sp_fn):
  m_data(data),
  m_lny(),
  m_negy(),
  m_sp_fn(sp_fn->clone())
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      vsassert(m_negy == 0 || m_negy == m_data[idata].on.ndim());
      m_negy = m_data[idata].on.ndim();
     
      vsassert(m_data[idata].on.ndim() == m_data[idata].off.ndim());
      for(unsigned iegy = 0; iegy < m_negy; iegy++)
	m_lny += VSAMath::lgamma(m_data[idata].on(iegy)+1) + 
	  VSAMath::lgamma(m_data[idata].off(iegy)+1);
    }
}

VSSpectrumFit::LnLFn::~LnLFn()
{
  delete m_sp_fn;
}

double VSSpectrumFit::LnLFn::val(const VSAAlgebra::VecND& a) const
{
  m_sp_fn->setParam(a);

  double lnl = 0;
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++) 
    lnl += val(a,m_data[idata],0,m_data[idata].on.ndim());

  return -2*(lnl - m_lny);
}

void VSSpectrumFit::LnLFn::gof(const VSAAlgebra::VecND& a,
			       double& chi2_pearson, double& chi2_ml,
			       unsigned& ndf) const
{
  m_sp_fn->setParam(a);

  ndf = 0;
  chi2_pearson = 0;
  chi2_ml = 0;
  std::vector< std::pair< unsigned,unsigned> > iegy_bins;

  const unsigned min_counts = 1;

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      const unsigned negy = m_data[idata].on.ndim();
      vsassert(m_data[idata].on.ndim() == m_data[idata].off.ndim());
      //      vsassert(iegylo < negy && iegyhi <= negy);

      const VSAAlgebra::VecND& on = m_data[idata].on;
      const VSAAlgebra::VecND& off = m_data[idata].off;
      double alpha = m_data[idata].alpha;
  
      VSAAlgebra::VecND mus;
      VSAAlgebra::VecND mub;

      VSSpectrumFit::LnLFn::mu_signal(m_data[idata].krn,mus);
      VSSpectrumFit::LnLFn::mu_bkgnd(on,off,alpha,mus,mub);

      // ----------------------------------------------------------------------
      // Find the energy range encompassing on source bins that have
      // have at least N counts.
      // ----------------------------------------------------------------------
      unsigned iegylo = 0;
      unsigned iegyhi = negy;

      for(unsigned iegy = 0; iegy < negy; iegy++)
	{
	  if(m_data[idata].on(iegy) >= min_counts)
	    {
	      iegylo = iegy;
	      break;
	    }
	}

      for(unsigned iegy = negy-1; iegy > iegylo; iegy--)
	{
	  if(m_data[idata].on(iegy) >= min_counts)
	    {
	      iegyhi = iegy;
	      break;
	    }
	}

      //      std::cout << "IEGYLO " << iegylo << " " << iegyhi << std::endl;
      ndf += iegyhi-iegylo;

      for(unsigned iegy = iegylo; iegy < iegyhi; iegy++)
	{
	  double mu_on = mus(iegy)+alpha*mub(iegy);
	  //	  double mu_off = mub(iegy);

	  chi2_pearson += std::pow(on(iegy)-mu_on,2)/mu_on;
	  //	  chi2_pearson += std::pow(off(iegy)-mu_off,2)/mu_off;
	  
	  if(on(iegy) > 0)
	    chi2_ml += 2*on(iegy)*log(on(iegy)/mu_on);

	  chi2_ml += 2*(mu_on - on(iegy));

// 	  if(off(iegy) > 0)
// 	    chi2_ml += 2*off(iegy)*log(off(iegy)/mu_off);

// 	  chi2_ml += 2*(mu_off - off(iegy));
	}

      iegy_bins.push_back(std::make_pair(iegylo,iegyhi));
    }

  ndf -= a.ndim();
}

void VSSpectrumFit::LnLFn::mu(const VSAAlgebra::VecND& a,
			      VSAAlgebra::VecND& mu_on,
			      VSAAlgebra::VecND& mu_off,
			      VSAAlgebra::VecND& mu_bkgnd)
{
  m_sp_fn->setParam(a);
  mu_on.resize(m_negy);
  mu_off.resize(m_negy);
  mu_bkgnd.resize(m_negy);

  mu_on.set(0.0);
  mu_off.set(0.0);
  mu_bkgnd.set(0.0);

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      const VSAAlgebra::VecND& on = m_data[idata].on;
      const VSAAlgebra::VecND& off = m_data[idata].off;
      double alpha = m_data[idata].alpha;
  
      VSAAlgebra::VecND mus;
      VSAAlgebra::VecND mub;

      VSSpectrumFit::LnLFn::mu_signal(m_data[idata].krn,mus);
      VSSpectrumFit::LnLFn::mu_bkgnd(on,off,alpha,mus,mub);

      mu_on += mus + alpha*mub;
      mu_off += mub;
      mu_bkgnd += alpha*mub;
    }
}

void VSSpectrumFit::LnLFn::mu_signal(const VSNSpace& krn_nspace,
				     VSAAlgebra::VecND& mus) const
{  
  VSAAlgebra::MatrixND krn = VSSpectrumCalc::toMatrixND(krn_nspace);
  //  krn.transpose();

  VSAAlgebra::VecND flux(krn.ncol());
  for(unsigned iegy = 0; iegy < flux.ndim(); iegy++)
    {
      double elo = krn_nspace.axis(1).minCoordUnchecked(iegy);
      double ehi = krn_nspace.axis(1).maxCoordUnchecked(iegy);

      flux(iegy) = m_sp_fn->integralFlux(elo,ehi);
    } 
  mus = krn*flux;
}

void VSSpectrumFit::LnLFn::mu_bkgnd(const VSAAlgebra::VecND& on,
				    const VSAAlgebra::VecND& off,
				    double alpha,
				    const VSAAlgebra::VecND& mus,
				    VSAAlgebra::VecND& mub) 
{
  const unsigned negy = mus.ndim();
  mub.resize(negy);
  mub.set(0.0);
  for(unsigned iegy = 0; iegy < negy; iegy++)
    {
      if(off(iegy) == 0) mub(iegy) = 0;
      else if(on(iegy) == 0) mub(iegy) = off(iegy)/(1+alpha);
      else
	{
	  const double a1 = alpha*(off(iegy)+on(iegy)) - mus(iegy)*(1+alpha);
	  vsassert(4*alpha*(1+alpha)*off(iegy)*mus(iegy) + a1*a1 > 0);
	  mub(iegy) = 0.5*std::pow(alpha+alpha*alpha,-1)*
	    (a1+sqrt(4*alpha*(1+alpha)*off(iegy)*mus(iegy) + a1*a1));
	}
    }
}

double VSSpectrumFit::LnLFn::val(const VSAAlgebra::VecND& a,
				 const Data& data, 
				 unsigned iegylo, unsigned iegyhi) const
{
  m_sp_fn->setParam(a);

  const VSAAlgebra::VecND& on = data.on;
  const VSAAlgebra::VecND& off = data.off;
  double alpha = data.alpha;
  
  VSAAlgebra::VecND mus;
  VSAAlgebra::VecND mub;
  VSSpectrumFit::LnLFn::mu_signal(data.krn,mus);
  VSSpectrumFit::LnLFn::mu_bkgnd(on,off,alpha,mus,mub);

  double lnl = 0;

  for(unsigned iegy = iegylo; iegy < iegyhi; iegy++)
    {
      lnl += - mus(iegy) - alpha*mub(iegy) - mub(iegy);

      if(on(iegy) > 0) lnl += on(iegy)*log(mus(iegy)+alpha*mub(iegy));
      if(off(iegy) > 0) lnl += off(iegy)*log(mub(iegy));
    }

  return lnl;
}

void VSSpectrumFit::LnLFn::val(const VSAAlgebra::VecND& a, 
			       double& lnl,
			       VSAAlgebra::MatrixND& beta,
			       VSAAlgebra::MatrixND& alpha) const
{
  vsassert(a.ndim() == nparm());
  m_sp_fn->setParam(a);
  beta.resize(nparm(),1);
  alpha.resize(nparm(),nparm());
  
  alpha.set(0.0);
  beta.set(0.0);
  lnl = 0;

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    val(a,lnl,beta,alpha,m_data[idata]);

  lnl -= m_lny;
  lnl *= -2;
}

void VSSpectrumFit::LnLFn::val(const VSAAlgebra::VecND& a, 
			       double& lnl,
			       VSAAlgebra::MatrixND& beta,
			       VSAAlgebra::MatrixND& cm,
			       const Data& data) const
{
  vsassert(a.ndim() == nparm());
  m_sp_fn->setParam(a);

//   std::cout << std::setw(15) << a(0) 
// 	    << std::setw(15) << a(1) 
// 	    << std::endl;

  const VSAAlgebra::VecND& on = data.on;
  const VSAAlgebra::VecND& off = data.off;
  double alpha = data.alpha;

  VSAAlgebra::MatrixND krn = VSSpectrumCalc::toMatrixND(data.krn);
  //  krn.transpose();

  VSAAlgebra::VecND phi(krn.ncol());
  VSAAlgebra::MatrixND dphida(krn.ncol(),nparm());
  for(unsigned iegy = 0; iegy < phi.ndim(); iegy++)
    {
      double elo = data.krn.axis(0).minCoordUnchecked(iegy);
      double ehi = data.krn.axis(0).maxCoordUnchecked(iegy);

      phi(iegy) = m_sp_fn->integralFlux(elo,ehi);

      VSAAlgebra::VecND dyda(nparm());
      m_sp_fn->integralFluxDerivative(elo,ehi,dyda);

      for(unsigned ip = 0; ip < nparm(); ip++)
	dphida(iegy,ip) = dyda(ip);

    }

  VSAAlgebra::VecND mus = krn*phi;
  VSAAlgebra::MatrixND dmusda = krn*dphida;

  VSAAlgebra::VecND mub;
  VSSpectrumFit::LnLFn::mu_bkgnd(on,off,alpha,mus,mub);

  const unsigned negy = on.ndim();
  for(unsigned iegy = 0; iegy < negy; iegy++)
    {
      const double a1 = alpha*(off(iegy)+on(iegy)) - mus(iegy)*(1+alpha);
      const double c1 = alpha+alpha*alpha;

      lnl += - mus(iegy) - alpha*mub(iegy) - mub(iegy);

      if(on(iegy) > 0) lnl += on(iegy)*log(mus(iegy)+alpha*mub(iegy));
      if(off(iegy) > 0) lnl += off(iegy)*log(mub(iegy));

      VSAAlgebra::VecND dmubda(nparm());
      for(unsigned ip1 = 0; ip1 < nparm(); ip1++)
	{
	  const double da1 = -dmusda(iegy,ip1)*(1+alpha);

	  if(on(iegy) == 0 || off(iegy) == 0)
	    dmubda(ip1) = 0;
	  else
	    dmubda(ip1) = 
	      0.5/c1*(da1+0.5*(4*c1*off(iegy)*dmusda(iegy,ip1) + 2*da1*a1)*
		      std::pow(4*c1*off(iegy)*mus(iegy) + a1*a1,-0.5));
	}

      for(unsigned ip1 = 0; ip1 < nparm(); ip1++)
	{
	  double ds1 = dmusda(iegy,ip1) + alpha*dmubda(ip1);
	  double db1 = dmubda(ip1);

	  beta(ip1,0) += - ds1 - db1;
	  
	  if(on(iegy) > 0) 
	    beta(ip1,0) += on(iegy)*ds1/(mus(iegy)+alpha*mub(iegy));
	  if(off(iegy) > 0) 
	    beta(ip1,0) += off(iegy)*db1/mub(iegy);

	  for(unsigned ip2 = 0; ip2 < nparm(); ip2++)
	    {
	      double ds2 = dmusda(iegy,ip2) + alpha*dmubda(ip2);
	      double db2 = dmubda(ip2);

	      if(on(iegy) > 0) 
		cm(ip1,ip2) += 
		  on(iegy)/std::pow(mus(iegy) + alpha*mub(iegy),2)*ds1*ds2;

	      if(off(iegy) > 0)
		cm(ip1,ip2) += off(iegy)/std::pow(mub(iegy),2)*db1*db2;
	    }	     

// 	  std::cout << std::setw(5) << iegy
//  		    << std::setw(5) << ip1
//  		    << std::setw(15) << dmubda(ip1)
//  		    << std::setw(15) << mub(iegy)
//  		    << std::setw(15) << mus(iegy)
//  		    << std::setw(15) << on(iegy)
//  		    << std::setw(15) << off(iegy)
//  		    << std::setw(15) << beta(ip1,0)
//  		    << std::setw(15) << cm(ip1,ip1)
// 		    << std::endl;
	}     
    }
}


// ============================================================================
// VSSpectrumCalc
// ============================================================================
VSSpectrumCalc::VSSpectrumCalc(double theta_cut, const std::string& sp_model,
			       double log10_enorm):
  m_theta_cut(theta_cut), m_max_off_regions(10)
{
  m_sp_fit = new VSSpectrumFit(theta_cut,sp_model,log10_enorm);
  m_sp_unfold = 
    VSSpectrumCalcFactory::getInstance()->createSpectrumUnfolding();
}

VSSpectrumCalc::~VSSpectrumCalc()
{
  delete m_sp_fit;
  delete m_sp_unfold;
}

void VSSpectrumCalc::accumulate(const VSEventArrayDatum& evt,
				const VSAAlgebra::Vec2D& evt_xy,
				VSAnalysisStage3Data::RunData& data)
{
  // Accumulate on counts -----------------------------------------------------
  double dr_src = data.src_xy().d(evt_xy);
  
  if(dr_src < m_theta_cut && std::isfinite(evt.mlt_log10_energy))
    {
      data.egy_on_hist.accumulate(evt.mlt_log10_energy);
      data.egy_on_np_hist.accumulate(evt.mlt_log10_energy);
    }

  data.egy_th2_on_hist.accumulate(evt.mlt_log10_energy,std::pow(dr_src,2));
  data.egy_th2_on_np_hist.accumulate(evt.mlt_log10_energy,std::pow(dr_src,2));

  // Accumulate off counts ----------------------------------------------------
  for(std::vector< VSAAlgebra::Vec2D >::const_iterator itr =
	data.sp_off_xy().begin(); itr != data.sp_off_xy().end(); ++itr)
    {
      double dr_off = itr->d(evt_xy);

      if(dr_off < m_theta_cut && std::isfinite(evt.mlt_log10_energy))
	{
	  data.egy_off_hist.accumulate(evt.mlt_log10_energy);
	  data.egy_off_np_hist.accumulate(evt.mlt_log10_energy);
	}

      data.egy_th2_off_hist.
	accumulate(evt.mlt_log10_energy,std::pow(dr_off,2));
      data.egy_th2_off_np_hist.
       	accumulate(evt.mlt_log10_energy,std::pow(dr_off,2));
    }
}

void VSSpectrumCalc::fit(const VSAnalysisStage3Data& stage3_data,
			 VSSpectrumFitData& data)
{
  m_sp_fit->fit(stage3_data,data);
}

void VSSpectrumCalc::reconstruct(const VSAnalysisStage3Data& stage3_data,
				 VSSpectrumData& data)
{
  m_sp_unfold->reconstruct(stage3_data,data);
}

// void VSSpectrumCalc::setRun(const VSAnalysisStage3Data::RunData& data)
// {
//   m_off_xy.clear();

//   // Get Off regions ----------------------------------------------------------
//   VSOffRegion::getRegions(data.src_xy(), data.obs_xy(),
// 			  m_theta_cut, m_off_regions, m_off_xy);

// //   for(unsigned i = 0; i < m_off_xy.size(); i++)
// //     {
// //       std::cout << m_off_xy[i].x() << " " << m_off_xy[i].y() << " "
// // 		<< m_off_xy[i].d(data.src_xy()) << " "
// // 		<< m_off_xy[i].d(data.obs_xy()) << std::endl;
// //     }
// }

VSAAlgebra::VecND 
VSSpectrumCalc::toVecND(const VSLimitedErrorsHist<double, double>& h)
{
  VSAAlgebra::VecND v(h.nBins());

  unsigned n = 0;
  for(VSLimitedErrorsHist<double, double>::iterator itr = h.begin();
      itr != h.end(); ++itr)
    v[n++] = itr->count();

  return v;
}

VSAAlgebra::VecND 
VSSpectrumCalc::toVecND(const VSLimitedHist<double, double>& h)
{
  VSAAlgebra::VecND v(h.nBins());

  unsigned n = 0;
  for(VSLimitedHist<double, double>::iterator itr = h.begin();
      itr != h.end(); ++itr)
    v[n++] = itr->count();

  return v;
}

VSAAlgebra::VecND VSSpectrumCalc::toVecND(const VSNSpace& nspace)
{
  vsassert(nspace.ndim() == 1);
  VSAAlgebra::VecND v(nspace.axis(0).nbin);

  for(unsigned i = 0; i < v.ndim(); i++)
    {
      VSNSpace::Cell c(1);
      c.i[0] = i;
      v[i] = nspace[c];
    }

  return v;
}

VSAAlgebra::MatrixND VSSpectrumCalc::toMatrixND(const VSNSpace& nspace)
{
  vsassert(nspace.ndim() == 2);
  VSAAlgebra::MatrixND m(nspace.axis(0).nbin,nspace.axis(1).nbin);

  for(unsigned irow = 0; irow < m.nrow(); irow++)
    {
      for(unsigned icol = 0; icol < m.ncol(); icol++)
	{
	  VSNSpace::Cell c(2);
	  c.i[0] = irow;
	  c.i[1] = icol;
	  m(irow,icol) = nspace[c];
	}
    }

  return m;
}

// ============================================================================
// VSSpectrumUnfolding
// ============================================================================
VSSpectrumUnfolding::VSSpectrumUnfolding(double theta_cut, 
					 const std::string& sp_model,
					 double ebin, double emin, double emax):
  m_ebin(ebin), m_emin(emin), m_emax(emax),
  m_theta_cut(theta_cut), m_max_off_regions(10)
{

}

VSSpectrumUnfolding::~VSSpectrumUnfolding()
{

}

void VSSpectrumUnfolding::finalize(VSSpectrumData& data)
{
  data.dfde_hist = data.flux_hist;
  data.edfde_hist = data.flux_hist;
  data.e2dfde_hist = data.flux_hist;

  for(VSLimitedErrorsHist<double,double>::iterator itr = 
	data.dfde_hist.begin(); itr != data.dfde_hist.end(); ++itr)
    {
      double ebin = std::pow(10,itr->center());
      double elo = std::pow(10,itr->val());
      double ehi = std::pow(10,itr->val()+data.dfde_hist.binSize());
      double flux = itr->count();
      double flux_var = itr->var();
      double de = ehi-elo;

      data.dfde_hist.setBin(itr->bin(),flux/de,flux_var*std::pow(1/de,2));
      data.edfde_hist.setBin(itr->bin(),ebin*flux/de,
			     flux_var*std::pow(ebin/de,2));
      data.e2dfde_hist.setBin(itr->bin(),ebin*ebin*flux/de,
			      flux_var*std::pow(ebin*ebin/de,2));
    }
}

// ============================================================================
// VSSpectrumUnfoldingCFM
// ============================================================================
VSSpectrumUnfoldingCFM::
VSSpectrumUnfoldingCFM(double theta_cut,
		       const std::string& sp_model,
		       double ebin, double emin, double emax,
		       double log10_enorm,
		       double ul_threshold,
		       unsigned niter):
  VSSpectrumUnfolding(theta_cut,sp_model,ebin,emin,emax),
  m_theta_cut(theta_cut), m_ul_threshold(ul_threshold), m_niter(niter)
{
  
}

void VSSpectrumUnfoldingCFM::
reconstruct(const VSAnalysisStage3Data& stage3_data,
	    VSSpectrumData& data) 
{
  VSAAlgebra::VecND excess;
  VSAAlgebra::VecND var;
  VSAAlgebra::VecND on;
  VSAAlgebra::VecND off;
  

  VSLimitedErrorsHist<double,double> excess_hist;

  VSNSpace kernel;

  const unsigned nrun = stage3_data.nrun();
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      const VSAnalysisStage3Data::RunData& rd = stage3_data.run_data(irun);
      double alpha = 1./rd.sp_noff();
      double offset = rd.src_xy().d(rd.obs_xy());

      VSNSpace k = rd.irf_calc().getKernel(m_ebin, m_emin, m_emax,
					   m_ebin, m_emin, m_emax,
					   offset, m_theta_cut);
      k *= rd.livetime_min*60;

      if(on.ndim() == 0)
	{
	  on = VSSpectrumCalc::toVecND(rd.egy_on_np_hist.contentsHist());
	  off = VSSpectrumCalc::toVecND(rd.egy_off_np_hist.contentsHist());
	  kernel = k;
	}
      else
	{
	  on += VSSpectrumCalc::toVecND(rd.egy_on_np_hist.contentsHist());
	  off += VSSpectrumCalc::toVecND(rd.egy_off_np_hist.contentsHist());
	  kernel += k;
	}

      if(irun == 0) 
	excess_hist = rd.egy_on_np_hist;
      else
	excess_hist += rd.egy_on_np_hist;

      excess_hist -= rd.egy_off_np_hist*alpha;
    }

  excess = VSSpectrumCalc::toVecND(excess_hist.contentsHist());
  var = VSSpectrumCalc::toVecND(excess_hist.varianceHist());

  data.egy_excess_hist = excess_hist;
  data.egy_nu_hist = excess_hist;
  data.egy_resid_hist = excess_hist;
  data.egy_resid_hist.clear();

  data.resid_hist = VSLimitedErrorsHist<double,double>(0.2,-10.,10.);

  data.flux_hist = VSLimitedErrorsHist<double,double>(m_ebin,m_emin,m_emax);
  data.flux_bias_hist = data.flux_hist;
  data.flux_rbias_hist = data.flux_hist;
  data.krn_nspace = kernel;

  data.flux_cov_hist = VSSimple2DHist<double,double>(m_ebin,m_emin,m_emax);

  VSSpectrumFn* sp_fn = new VSSpectrumFnPowerLaw(2.5);

  VSAAlgebra::MatrixND krn = VSSpectrumCalc::toMatrixND(kernel);

  const unsigned negy = krn.ncol();
  // const double erec_lo = kernel.axis(0).lo_bound;
  // const double erec_hi = kernel.axis(0).hi_bound;
  // const double erec_bin = kernel.axis(0).bin_size;

  VSAAlgebra::VecND phi(negy);
  VSAAlgebra::VecND nu;
  for(unsigned iegy = 0; iegy < negy; iegy++)
    {
      double elo = kernel.axis(1).minCoordUnchecked(iegy);
      double ehi = kernel.axis(1).maxCoordUnchecked(iegy);
      phi(iegy) = sp_fn->integralFlux(elo,ehi);
    }

  VSAAlgebra::VecND flux(negy);
  VSAAlgebra::VecND flux_var(negy);
  for(unsigned i = 0; i < m_niter; i++)
    {
      //std::cout << "-------------------------------------------" << std::endl;
      double chi2 = 0;
      nu = krn*phi;

      for(unsigned irow = 0; irow < negy; irow++)
	{
	  if(var(irow))
	    chi2 += std::pow(nu(irow)-excess(irow),2)/var(irow);

	  if(nu(irow) > 0)
	    {
	      flux(irow) = excess(irow)*phi(irow)/nu(irow);
	      flux_var(irow) = var(irow)*std::pow(phi(irow)/nu(irow),2);
	    }
	  else
	    {
	      flux(irow) = 0;
	      flux_var(irow) = 0;
	    }
	}

      VSAMath::Data<double> fdata;
      for(unsigned icol = 0; icol < negy; icol++)
	{
	  if(flux_var(icol) && flux(icol) > 0)
	    {	  
	      double y = std::log(flux(icol));
	      double yerr = sqrt(flux_var(icol))/flux(icol);
	      double egy = kernel.axis(1).midCoordUnchecked(icol);

	      fdata.insert(VSAMath::DataPoint<double>(egy,y,yerr));
	    }
	}

      VSAAlgebra::VecND p;
      VSAMath::PolyFit::fit(2,fdata,p);
      for(unsigned icol = 0; icol < negy; icol++)
	{
	  double egy = kernel.axis(1).midCoordUnchecked(icol);
	  phi(icol) = std::exp(VSAMath::PolyFit::val(p,egy));
	}      

      std::cout << "VSSpectrumUnfoldingCFM::reconstruct(): "
		<< "CHI2 " << chi2 << std::endl;
    }

  data.egy_nu_hist.set(nu.data());

  double chi2a = 0;
  double chi2b = 0;
  for(unsigned irow = 0; irow < negy; irow++)
  {
    double egy = kernel.axis(0).midCoordUnchecked(irow);
    if(var(irow))
      {
	data.egy_resid_hist.accumulate(egy,nu(irow)-excess(irow),var(irow));
	data.resid_hist.accumulate((nu(irow)-excess(irow))/sqrt(var(irow)));
	chi2a += std::pow(nu(irow)-excess(irow),2)/var(irow);
	chi2b += std::pow(nu(irow)-excess(irow),2)/var(irow);
      }
    else
      chi2b += std::pow(nu(irow)-excess(irow),2);
  }

  // std::cout << "CHI2A " << chi2a << std::endl;
  // std::cout << "CHI2B " << chi2b << std::endl;

  for(unsigned icol = 0; icol < negy; icol++)
    {
      double sigma = flux(icol)/sqrt(flux_var(icol));
      double egy = kernel.axis(1).midCoordUnchecked(icol);
	  
      if(sigma > m_ul_threshold)
	data.flux_hist.accumulate(egy, flux(icol), flux_var(icol));
      else
	{
	  double ul = VSAStatistics::heleneUL(flux(icol), 
					      sqrt(flux_var(icol)), 0.95);
	  data.flux_hist.accumulate(egy,ul,0);
	}      
    }

  data.chi2 = chi2a;

  finalize(data);
}

// ============================================================================
// VSSpectrumUnfoldingChi2
// ============================================================================
VSSpectrumUnfoldingChi2::
VSSpectrumUnfoldingChi2(double theta_cut, 
			const std::string& sp_model,
			double ebin, 
			double emin, double emax,
			double log10_enorm,
			double ul_threshold,
			const std::string& reg_method,
			double lambda, double index): 
  VSSpectrumUnfolding(theta_cut,sp_model,ebin,emin,emax),
  m_theta_cut(theta_cut), m_ul_threshold(ul_threshold), 
  m_reg_method(reg_method),  m_lambda(lambda), m_scaling_index(index)
{

}

void VSSpectrumUnfoldingChi2::
reconstruct(const VSAnalysisStage3Data& stage3_data,
	    VSSpectrumData& data) 
{
  VSAAlgebra::VecND excess;
  VSAAlgebra::VecND var;
  VSAAlgebra::VecND on;
  VSAAlgebra::VecND off;
  
  VSLimitedErrorsHist<double,double> excess_hist;

  VSNSpace krn;

  double alpha = 1;

  const unsigned nrun = stage3_data.nrun();
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      const VSAnalysisStage3Data::RunData& rd = stage3_data.run_data(irun);
      double alpha = 1./rd.sp_noff();
      double offset = (rd.obs_xy()-rd.src_xy()).norm();

      double erec_bin = stage3_data.run_data(irun).egy_on_hist.binSize();
      double erec_lo = stage3_data.run_data(irun).egy_on_hist.loLimit();
      double erec_hi = stage3_data.run_data(irun).egy_on_hist.hiLimit();

      VSNSpace k = rd.irf_calc().getKernel(m_ebin, m_emin, m_emax,
					   erec_bin, erec_lo, erec_hi,
					   offset, m_theta_cut);

      k *= rd.livetime_min*60;

      if(on.ndim() == 0)
	{
	  data.egy_on_hist = stage3_data.run_data(irun).egy_on_hist;
	  data.egy_off_hist = stage3_data.run_data(irun).egy_off_hist;
	  on = VSSpectrumCalc::toVecND(rd.egy_on_hist.contentsHist());
	  off = VSSpectrumCalc::toVecND(rd.egy_off_hist.contentsHist());
	  krn = k;
	}
      else
	{
	  data.egy_on_hist += stage3_data.run_data(irun).egy_on_hist;
	  data.egy_off_hist += stage3_data.run_data(irun).egy_off_hist;

	  on += VSSpectrumCalc::toVecND(stage3_data.run_data(irun).
					egy_on_hist.contentsHist());
	  off += VSSpectrumCalc::toVecND(stage3_data.run_data(irun).
					 egy_off_hist.contentsHist());
	  krn += k;
	}

      if(irun == 0) 
	excess_hist = stage3_data.run_data(irun).egy_on_hist;
      else
	excess_hist += stage3_data.run_data(irun).egy_on_hist;

      excess_hist -= stage3_data.run_data(irun).egy_off_hist*alpha;
    }

  // Create scaling matrix ----------------------------------------------------
  // VSSpectrumFit* spfit = new VSSpectrumFit(m_theta_cut,"powerlaw");
  // VSSpectrumFitData spfit_data;
  // spfit->fit(stage3_data,spfit_data);
  // delete spfit;

  excess = VSSpectrumCalc::toVecND(excess_hist.contentsHist());
  var = VSSpectrumCalc::toVecND(excess_hist.varianceHist());

  data.egy_excess_hist = excess_hist;
  data.egy_nu_hist = excess_hist;
  data.egy_resid_hist = excess_hist;
  data.resid_hist = VSLimitedErrorsHist<double,double>(0.2,-10.,10.);

  data.flux_hist = VSLimitedErrorsHist<double,double>(m_ebin,m_emin,m_emax);
  data.flux_bias_hist = data.flux_hist;
  data.flux_rbias_hist = data.flux_hist;
  data.flux_cov_hist = VSSimple2DHist<double,double>(m_ebin,m_emin,m_emax);

  data.krn_nspace = krn;
  reconstruct(krn,excess,var,on,off,alpha,data);

  finalize(data);
}

void VSSpectrumUnfoldingChi2::
reconstruct(const VSNSpace& kernel_nspace,
	    const VSAAlgebra::VecND& excess,
	    const VSAAlgebra::VecND& var,
	    const VSAAlgebra::VecND& on,
	    const VSAAlgebra::VecND& off,
	    double alpha,
	    VSSpectrumData& data) 
{ 
  VSAAlgebra::MatrixND kernel = VSSpectrumCalc::toMatrixND(kernel_nspace);
  data.krn = kernel;

  const unsigned nrow = kernel.nrow();
  const unsigned ncol = kernel.ncol();

  VSAAlgebra::MatrixND s = calcSmoothingMatrix(m_reg_method,kernel_nspace);
  VSAAlgebra::MatrixND k = kernel;
  VSAAlgebra::MatrixND v(nrow,nrow);
  VSAAlgebra::MatrixND vi(nrow,nrow);

  double ndf = 0;

  for(unsigned irow = 0; irow < nrow; irow++)
    {
      if(var[irow] > 0)
	{
	  vi(irow,irow) = 1/var[irow];
	  v(irow,irow) = var[irow];
	  ndf++;
	}
      else
	{
	  vi(irow,irow) = 1;
	  v(irow,irow) = 1;
	}
    }
  
  VSAAlgebra::MatrixND kt = k.getTranspose();
  VSAAlgebra::MatrixND ksq = kt*vi*k;

  VSAAlgebra::SVD svd(ksq);

 
  VSAAlgebra::VecND param(2);
  param(0) = 1.0;
  param(1) = m_scaling_index;

  VSSpectrumFn* spfn = VSSpectrumFn::create("powerlaw",param);

  const unsigned nbin = kernel_nspace.axis(1).nbin;
  const double emin = kernel_nspace.axis(1).lo_bound;
  const double ebin = kernel_nspace.axis(1).bin_size;

  VSAAlgebra::MatrixND sm(nbin,nbin);
  for(unsigned irow = 0; irow < nbin; irow++)
    {
      double elo = kernel_nspace.axis(1).minCoordUnchecked(irow);
      double ehi = kernel_nspace.axis(1).maxCoordUnchecked(irow);
      sm(irow,irow) = 
	spfn->integralFlux(emin,emin+ebin)/spfn->integralFlux(elo,ehi);
    }

  std::cout << sm << std::endl;

  s = sm*s*sm;

  double lambda = svd.w()(0);


  // double nbin = 0;
  // for(unsigned ibin = 0; ibin < excess.ndim(); ibin++)
  //   if(var(ibin) > 0) nbin++;
  
  Data spd;
  data.lambda_gcc_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);
  data.lambda_chi2_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);
  data.lambda_chi2b_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);
  data.lambda_chi2_ml_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);
  data.lambda_msbias_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);
  data.lambda_mvar_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);
  data.lambda_mse_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);
  data.lambda_wmse_hist = VSLimitedErrorsHist<double,double>(0.1,-14,0);

  //#if 0
  for(VSLimitedErrorsHist<double,double>::iterator itr = 
	data.lambda_mse_hist.begin(); itr != data.lambda_mse_hist.end(); ++itr)
    {

      double lam = lambda*std::pow(10,itr->center());

      Data d;
      calcSpectrum(k,s,vi,v,excess,lam,d);

      VSAAlgebra::VecND mub;
      VSSpectrumFit::LnLFn::mu_bkgnd(on,off,alpha,d.nu,mub);

      unsigned ml_ndf = 0;

      for(unsigned irow = 0; irow < nrow; irow++)
	{
	  if(on(irow)<5) continue;

	  ml_ndf++;

	  double mu = mub(irow);
	  double mu_on = d.nu(irow)+alpha*mu;

	  if(on(irow)) d.chi2_ml += 2*on(irow)*log(on(irow)/mu_on);
	  d.chi2_ml += 2*(mu_on-on(irow));
	}

      std::cout << "LAMBDA " 
		<< std::setw(15) << lam/lambda
		<< std::setw(15) << d.ksq_ndf
		<< std::setw(15) << d.ksqr_ndf
		<< std::setw(15) << d.chi2/ndf
		<< std::setw(15) << d.chi2_ml/ml_ndf
		<< std::setw(15) << d.mean_gcc
		<< std::endl;

      data.lambda_gcc_hist.accumulate(itr->center(),d.mean_gcc,0);
      data.lambda_chi2_hist.accumulate(itr->center(),d.chi2,0);
      data.lambda_chi2b_hist.accumulate(itr->center(),d.chi2b,0);
      data.lambda_chi2_ml_hist.accumulate(itr->center(),d.chi2_ml/ndf,0);
      data.lambda_msbias_hist.accumulate(itr->center(),d.msbias,0);
      data.lambda_mvar_hist.accumulate(itr->center(),d.mvar,0);
      data.lambda_mse_hist.accumulate(itr->center(),d.mse,0);
      data.lambda_wmse_hist.accumulate(itr->center(),d.wmse,0);
				      
      data.chi2_mse_graph.addVertex(d.chi2,d.mse);
      data.chi2_ml_mse_graph.addVertex(d.chi2_ml,d.mse);
      data.chi2_wmse_graph.addVertex(d.chi2,d.wmse);
      data.chi2_ml_wmse_graph.addVertex(d.chi2_ml,d.wmse);


      if(itr == data.lambda_mse_hist.begin()) spd = d;
      // else if(fabs(d.ksqr_ndf-d.ksq_ndf) < fabs(spd.ksqr_ndf-spd.ksq_ndf))
      else if(fabs(1-d.chi2/ndf) < fabs(1-spd.chi2/ndf))
	{
	  spd = d;
	  spd.lambda = std::pow(10,itr->center());
	}

      // else if(d.wmse < spd.wmse && d.chi2_ml/ndf < 1)
      // 	{
      // 	  spd = d;
      // 	  spd.lambda = std::pow(10,itr->center());
      // 	}
    }
  //#endif

  if(m_lambda > 0)
    calcSpectrum(k,s,vi,v,excess,m_lambda*lambda,spd);

  std::cout << "LAMBDA " << std::setw(25) << spd.lambda << std::endl;
  std::cout << "CHI2   " 
	    << std::setw(25) << spd.chi2 
	    << std::setw(25) << spd.chi2/ndf 
	    << std::setw(25) << ndf
	    << std::endl;

  data.chi2 = spd.chi2;
  data.lambda = spd.lambda;
  data.egy_nu_hist.set(spd.nu.data());

  data.egy_resid_hist.clear();
  for(unsigned irow = 0; irow < nrow; irow++)
  {
    double egy = kernel_nspace.axis(0).midCoordUnchecked(irow);
    if(var(irow))
      {
	data.egy_resid_hist.accumulate(egy,spd.nu(irow)-excess(irow),var(irow));
	data.resid_hist.accumulate((spd.nu(irow)-excess(irow))/sqrt(var(irow)));
      }
  }

  data.flux_cov_hist.set(spd.flux_cov.data());
  data.mu_cov_hist = 
    VSSimple2DHist<double,double>(1.0,0,spd.c.ncol(),1.0,0,spd.c.nrow());
  data.mu_cov_hist.set(spd.c.data());
  data.nu_cov_hist = 
    VSSimple2DHist<double,double>(1.0,0,spd.kc.ncol(),1.0,0,spd.kc.nrow());
  data.nu_cov_hist.set(spd.kc.data());

  data.k_hist = 
    VSSimple2DHist<double,double>(1.0,0,spd.k.ncol(),1.0,0,spd.k.nrow());
  data.k_hist.set(spd.k.data());

  data.ksq_hist = 
    VSSimple2DHist<double,double>(1.0,0,spd.ksq.ncol(),1.0,0,spd.ksq.nrow());
  data.ksq_hist.set(spd.ksq.data());

  data.ksqr_hist = 
    VSSimple2DHist<double,double>(1.0,0,spd.ksqr.ncol(),1.0,0,spd.ksqr.nrow());
  data.ksqr_hist.set(spd.ksqr.data());

  data.ksq_u_hist = 
    VSSimple2DHist<double,double>(1.0,0,spd.ksq_u.ncol(),
				  1.0,0,spd.ksq_u.nrow());
  data.ksq_u_hist.set(spd.ksq_u.data());

  data.ksq_w_hist = 
    VSLimitedHist<double,double>(1.0,0,spd.ksq_w.ndim(),spd.ksq_w.data());
  data.ksq_sqw_hist = 
    VSLimitedHist<double,double>(1.0,0,spd.ksq_sqw.ndim(),spd.ksq_sqw.data());

  data.ksq_c_hist = 
    VSLimitedErrorsHist<double,double>(1.0,0,spd.ksq_c.ndim());
  data.ksq_c_hist.fill(0);
  for(unsigned ibin = 0; ibin < data.ksq_c_hist.nBins(); ibin++)
    data.ksq_c_hist.setBin(ibin,spd.ksq_c(ibin),1);

  data.ksq_nc_hist = 
    VSLimitedErrorsHist<double,double>(1.0,0,spd.ksq_nc.ndim());
  data.ksq_nc_hist.fill(0);
  for(unsigned ibin = 0; ibin < data.ksq_nc_hist.nBins(); ibin++)
    data.ksq_nc_hist.setBin(ibin,spd.ksq_nc(ibin),
			    std::pow(spd.ksq_nc_err(ibin),2));

  data.ksqr_c_hist = 
    VSLimitedErrorsHist<double,double>(1.0,0,spd.ksqr_c.ndim());
  data.ksqr_c_hist.fill(0);
  for(unsigned ibin = 0; ibin < data.ksqr_c_hist.nBins(); ibin++)
    data.ksqr_c_hist.setBin(ibin,spd.ksqr_c(ibin),1);

  data.ksqr_f_hist = 
    VSLimitedErrorsHist<double,double>(1.0,0,spd.ksqr_f.ndim());
  data.ksqr_f_hist.set(spd.ksqr_f.data());

  data.ksqr_u_hist = 
    VSSimple2DHist<double,double>(1.0,0,spd.ksqr_u.ncol(),
				  1.0,0,spd.ksqr_u.nrow());
  data.ksqr_u_hist.set(spd.ksqr_u.data());

  for(unsigned icol = 0; icol < ncol; icol++)
    {
      double sigma = spd.flux[icol]/sqrt(spd.flux_cov(icol,icol));
      double egy = kernel_nspace.axis(1).midCoordUnchecked(icol);

      if(sigma > m_ul_threshold)
	data.flux_hist.accumulate(egy, spd.flux[icol], 
				  spd.flux_cov(icol,icol));
      else
	{
	  double ul = 
	    VSAStatistics::heleneUL(spd.flux[icol], 
				    sqrt(spd.flux_cov(icol,icol)), 0.95);
	  data.flux_hist.accumulate(egy,ul,0);
	}

    }


  for(unsigned icol = 0; icol < ncol; icol++)
    {
      double egy = kernel_nspace.axis(1).midCoordUnchecked(icol);
      data.flux_bias_hist.accumulate(egy, spd.bias[icol], 
				     spd.bias_cov(icol,icol));

      data.flux_rbias_hist.
	accumulate(egy, spd.bias[icol]/spd.flux[icol], 
		   spd.bias_cov(icol,icol)/std::pow(spd.flux[icol],2));
    }

}

void VSSpectrumUnfoldingChi2::
calcSpectrum(const VSAAlgebra::MatrixND& k,
	     const VSAAlgebra::MatrixND& s,
	     const VSAAlgebra::MatrixND& vi,
	     const VSAAlgebra::MatrixND& v,
	     const VSAAlgebra::VecND& excess,
	     double lambda,
	     VSSpectrumUnfolding::Data& spd)
{
  // Compute Spectrum ---------------------------------------------------------
  VSAAlgebra::MatrixND kt = k.getTranspose();
  VSAAlgebra::MatrixND ksq = kt*vi*k;
  
  VSAAlgebra::MatrixND h = ksq + lambda*s;
  VSAAlgebra::MatrixND b = h.inverse()*kt*vi;
  VSAAlgebra::MatrixND bsq = b*b.getTranspose();

  spd.flux = b*excess;
  spd.flux_cov = b*v*b.getTranspose();
  spd.nu = k*spd.flux;


  //  VSAAlgebra::VecND bp = kt*ep;
  VSAAlgebra::VecND nu = k*spd.flux;

  spd.resid = nu-excess;
  spd.chi2 = (nu-excess)*vi*(nu-excess);
  spd.rnorm = (nu-excess)*(nu-excess);
  spd.lambda = lambda;
  spd.k = k;
  spd.ksq = ksq;
  spd.ksqr = h;

  // Compute bias -------------------------------------------------------------
  VSAAlgebra::MatrixND ab = ksq + lambda*s;
  VSAAlgebra::MatrixND bb = -kt*vi;
  VSAAlgebra::MatrixND c = ab.inverse()*bb;
  VSAAlgebra::MatrixND kc = k*c;

  spd.c = c;
  spd.kc = kc;
  spd.bias = c*(nu-excess);
  spd.bias_cov = (c*k*c-c)*(c*k*c-c).getTranspose();

  // Fourier Decomposition of Kernel Matrix -----------------------------------
  VSAAlgebra::SVD svd(ksq);
  VSAAlgebra::MatrixND lam(ksq.nrow(),ksq.nrow());

  spd.ksq_w = svd.w();
  spd.ksq_sqw = svd.w();
  spd.ksq_u = svd.u();

  for(unsigned irow = 0; irow < ksq.nrow(); irow++)
    {
      spd.ksq_sqw(irow) = sqrt(spd.ksq_sqw(irow));
      lam(irow,irow) = 1/sqrt(svd.w()(irow));
    }

  spd.ksq_c = lam*svd.u().getTranspose()*kt*vi*excess;
  spd.ksq_nc = spd.ksq_c*lam;
  spd.ksq_nc_err = spd.ksq_nc;

  

  // Find NDF for unregularized solution --------------------------------------
  spd.ksq_ndf = 0;
  for(unsigned i = 0; i < spd.ksq_c.ndim(); i++)
    { 
      if(fabs(spd.ksq_c(i)) > 2) spd.ksq_ndf++;
      spd.ksq_nc_err(i) = lam(i,i);
    }

  VSAAlgebra::MatrixND u1 = spd.ksq_u;
  VSAAlgebra::MatrixND u1t = u1.getTranspose();

  VSAAlgebra::MatrixND s2 = lam*u1t*s*u1*lam;
  VSAAlgebra::SVD svds2(s2);
  
  VSAAlgebra::MatrixND u2 = svds2.u();
  VSAAlgebra::MatrixND u2t = u2.getTranspose();

  VSAAlgebra::MatrixND d = u2t*lam*u1t*s*u1*lam*u2;

  spd.ksqr_u = u1*lam*u2;

  // std::cout << "S" << std::endl;
  // std::cout << s << std::endl;
  // std::cout << "U" << std::endl;
  // std::cout << spd.ksq_u << std::endl;
  // std::cout << "W" << std::endl;
  // std::cout << spd.ksq_w << std::endl;
  // std::cout << "S2" << std::endl;
  // std::cout << s2 << std::endl;
  // std::cout << "D" << std::endl;
  // std::cout << d << std::endl;


  spd.ksqr_c = 
    svds2.u().getTranspose()*lam*svd.u().getTranspose()*kt*vi*excess;

  spd.ksqr_f = spd.ksq_c;
  spd.ksqr_ndf = 0;
  for(unsigned i = 0; i < spd.ksqr_f.ndim(); i++)
    {

      double w = svds2.w()(i);
      // double w2 = d(i,i);
      // std::cout << std::setw(10) << i 
      // 		<< std::setw(20) << w
      // 		<< std::setw(20) << w2
      // 		<< std::setw(20) << 1/w2
      // 		<< std::setw(20) << w/(w+lambda)
      // 		<< std::setw(20) << 1/(1+lambda*w2)
      // 		<< std::endl;



      spd.ksqr_f(i) = 1/(1+lambda*w);
      spd.ksqr_ndf += spd.ksqr_f(i);
    }

  VSAAlgebra::VecND f2(ksq.nrow());

  for(unsigned irow = 0; irow < ksq.nrow(); irow++)
    {
      VSAAlgebra::VecND u = spd.ksq_u.columnVector(irow);

      f2 += spd.ksqr_f(irow)*lam(irow,irow)*spd.ksq_c(irow)*u;
    }

  // for(unsigned irow = 0; irow < ksq.nrow(); irow++)
  //   {
  //     std::cout << std::setw(20) << f2(irow) 
  // 		<< std::setw(20) << spd.flux(irow) 
  // 		<< std::endl;
  //   }


  //  std::cout << "NDF " << spd.ksq_ndf << " " << spd.ksqr_ndf << std::endl;

  spd.msbias = 0;
  spd.mvar = 0;
  spd.mse = 0;
  spd.wmse = 0;
  spd.mean_gcc = 0;
  spd.gcc = VSAAlgebra::VecND(ksq.nrow());

  VSAAlgebra::MatrixND flux_covi = VSAAlgebra::SVD::inverse(spd.flux_cov);
  VSAAlgebra::VecND rho(ksq.nrow());

  const unsigned nrow = ksq.nrow();
  for(unsigned irow = 0; irow < nrow; irow++)
    {
      // Mean Variance
      spd.mvar += (1./(double)nrow)*spd.flux_cov(irow,irow);

      // Mean Squared Bias
      spd.msbias += (1./(double)nrow)*std::pow(spd.bias(irow),2);

      // Mean Squared Error
      spd.mse += (1./(double)nrow)*(spd.flux_cov(irow,irow) + 
				    std::pow(spd.bias(irow),2));

      // Chi-Squared of Bias
      spd.chi2b += 
	std::pow(spd.bias(irow),2)/spd.bias_cov(irow,irow);


      double vx = 
	std::min(1.,std::pow(spd.flux_cov(irow,irow)*flux_covi(irow,irow),-1));

      rho(irow) = sqrt(1-vx);
      spd.mean_gcc += rho(irow)/nrow;
      spd.gcc(irow) = sqrt(1-vx);
    }
}


VSAAlgebra::MatrixND VSSpectrumUnfoldingChi2::
calcSmoothingMatrix(const std::string& method, const VSNSpace& krn)
{
  // VSAAlgebra::VecND param(2);
  // param(0) = 1.0;
  // param(1) = 2.0;

  // VSSpectrumFn* spfn = VSSpectrumFn::create("powerlaw",param);

  const unsigned nbin = krn.axis(1).nbin;
  // const double emin = krn.axis(1).lo_bound;
  // const double ebin = krn.axis(1).bin_size;

  // VSAAlgebra::MatrixND sm(nbin,nbin);
  // for(unsigned irow = 0; irow < nbin; irow++)
  //   {
  //     double elo = krn.axis(1).minCoordUnchecked(irow);
  //     double ehi = krn.axis(1).maxCoordUnchecked(irow);
  //     sm(irow,irow) = 
  // 	spfn->integralFlux(emin,emin+ebin)/spfn->integralFlux(elo,ehi);
  //     std::cout << irow << " " << sm(irow,irow) << std::endl;
  //   }

  VSAAlgebra::MatrixND s1;

  if(method == "pol0")
    {
      s1 = VSAAlgebra::MatrixND(nbin,nbin);
      for(unsigned irow = 0; irow < nbin; irow++)
	s1[irow][irow] = 1;
    }
  else if(method == "pol1")
    {
      s1 = VSAAlgebra::MatrixND(nbin-1,nbin);
      for(unsigned irow = 0; irow < nbin-1; irow++)
	{
	  for(unsigned icol = 0; icol < nbin; icol++)
	    {
	      if(irow == icol)
		s1[irow][icol] = -1;
	      else if(icol == irow+1)
		s1[irow][icol] = 1;
	    }
	}
    }
  else if(method == "pol2")
    {
      s1 = VSAAlgebra::MatrixND(nbin-2,nbin);
      for(unsigned irow = 0; irow < nbin-2; irow++)
	{
	  for(unsigned icol = 0; icol < nbin; icol++)
	    {

	      if(irow == icol)
		s1[irow][icol] = -1;
	      else if(icol == irow+1)
		s1[irow][icol] = 2;
	      else if(icol == irow+2)
		s1[irow][icol] = -1;
	    }
	}
    }
  else if(method == "pol3")
    {
      s1 = VSAAlgebra::MatrixND(nbin-3,nbin);
      for(unsigned irow = 0; irow < nbin-3; irow++)
	{
	  for(unsigned icol = 0; icol < nbin; icol++)
	    {

	      if(irow == icol)
		s1[irow][icol] = -1;
	      else if(icol == irow+1)
		s1[irow][icol] = 3;
	      else if(icol == irow+2)
		s1[irow][icol] = -3;
	      else if(icol == irow+3)
		s1[irow][icol] = 1;
	    }
	}
    }
  else
    {
      std::cerr << "Unrecognized smoothing matrix: " 
		<< method << std::endl;
      exit(EXIT_FAILURE);
    }

  //  return sm*s1.getTranspose()*s1*sm;
  return s1.getTranspose()*s1;
}

#if 0

// ============================================================================
// VSSpectrumUnfoldingML
// ============================================================================
VSSpectrumUnfoldingML::
LnLFn::LnLFn():
  m_fn(), m_s(), m_lambda()
{

}

VSSpectrumUnfoldingML::
LnLFn::LnLFn(const std::vector<VSSpectrumFit::Data>& data,
	     VSSpectrumFn* sp_fn,
	     const VSAAlgebra::MatrixND& s,
	     const VSAAlgebra::MatrixND& sm,
	     double lambda):
  m_fn(data,sp_fn), m_s(s), m_sm(sm), m_lambda(lambda)
{

}

double VSSpectrumUnfoldingML::
LnLFn::val(const VSAAlgebra::VecND& a) const
{
  const VSAAlgebra::VecND a2 = m_sm*a;

//   std::cout << "VAL " << m_fn.val(a2) << " " << m_lambda << " "
// 	    << a*m_s*a << std::endl;

  return m_fn.val(a2) + m_lambda*a*m_s*a;
}     

void VSSpectrumUnfoldingML::
LnLFn::val(const VSAAlgebra::VecND& a, 
	   double& chi2,
	   VSAAlgebra::MatrixND& beta,
	   VSAAlgebra::MatrixND& alpha) const
{
  const unsigned ndim = a.ndim();

  for(unsigned ip = 0; ip < ndim; ip++)
    {
      if(a(ip) <= 0)
	{
	  std::ostringstream os;
	  
// 	  os << std::endl
// 	     << "y = " << std::setw(15) << y
// 	     << " ym = " << std::setw(15) << ym
// 	     << " x = " << std::setw(15) << m_data[idata].x
// 	     << std::endl;
	  
	  throw 
	    std::domain_error(std::string(__PRETTY_FUNCTION__)+
			      ": Model value less than or equal to 0."+
			      os.str());
	}
    }

  const VSAAlgebra::VecND a2 = m_sm*a;


  beta.resize(ndim,1);
  beta.set(0.);
  alpha.resize(ndim,ndim);
  alpha.set(0.);
  chi2 = 0;
  m_fn.val(a2,chi2,beta,alpha);

  beta = m_sm*beta;
  alpha = m_sm*alpha*m_sm;

  VSAAlgebra::VecND b = 2*m_s*a;

  chi2 += m_lambda*a*m_s*a;
  //  beta += m_s*a;
  alpha += 2*m_lambda*m_s;

  for(unsigned i = 0; i < beta.ndim(); i++)
    {
//       std::cout << i << " " << beta(i,0) << " " << m_lambda*b(i)
// 		<< std::endl;
      beta(i,0) += m_lambda*b(i);
    }
//   std::cout << "--------------------------------------------" << std::endl;
//   std::cout << a << std::endl << "CHI2 " << chi2 << std::endl;

}

VSSpectrumUnfoldingML::SpectrumFn::SpectrumFn(unsigned nparm,
							    double ebin,
							    double emin):
  VSSpectrumFn(nparm,"fn"), m_binner(ebin,emin), m_ebin(ebin)
{

}
      
double VSSpectrumUnfoldingML::
SpectrumFn::val(const double& log10_egy_tev, 
		const VSAAlgebra::VecND& a) const
{
  double de = std::pow(10,log10_egy_tev + m_ebin/2.) - 
    std::pow(10,log10_egy_tev - m_ebin/2.);

  return a(m_binner.valToBin(log10_egy_tev))/de;
}

double VSSpectrumUnfoldingML::
SpectrumFn::val(const double& log10_egy_tev) const
{
  double de = std::pow(10,log10_egy_tev + m_ebin/2.) - 
    std::pow(10,log10_egy_tev - m_ebin/2.);

  return param(m_binner.valToBin(log10_egy_tev))/de;
}

void VSSpectrumUnfoldingML::
SpectrumFn::dyda(const double& log10_egy_tev, 
		 VSAAlgebra::VecND& dyda) const
{
  dyda = VSAAlgebra::VecND(nparm());
  dyda(m_binner.valToBin(log10_egy_tev)) = 1;  
}

void VSSpectrumUnfoldingML::
SpectrumFn::dyda(const double& log10_egy_tev, 
		 const VSAAlgebra::VecND& a,
		 VSAAlgebra::VecND& dyda) const
{
  dyda = VSAAlgebra::VecND(nparm());
  dyda(m_binner.valToBin(log10_egy_tev)) = 1;  
}
      
double VSSpectrumUnfoldingML::
SpectrumFn::integralFlux(double log10_elo_tev,
			 double log10_ehi_tev,
			 const VSAAlgebra::VecND& a) const
{
  return a(m_binner.valToBin(0.5*(log10_ehi_tev + log10_elo_tev)));
}

double VSSpectrumUnfoldingML::
SpectrumFn::integralFlux(double log10_elo_tev,
			 double log10_ehi_tev) const
{
  return param(m_binner.valToBin(0.5*(log10_ehi_tev + log10_elo_tev)));
}

void VSSpectrumUnfoldingML::
SpectrumFn::integralFluxDerivative(double log10_elo_tev,
				   double log10_ehi_tev,
				   VSAAlgebra::VecND& dyda) const
{
  dyda = VSAAlgebra::VecND(nparm());
  dyda(m_binner.valToBin(0.5*(log10_ehi_tev + log10_elo_tev))) = 1;  
}

VSSpectrumUnfoldingML::
VSSpectrumUnfoldingML(double theta_cut, 
		      const std::string& sp_model,
		      double ebin, double emin, double emax,
		      double log10_enorm,
		      double ul_threshold,
		      const std::string& reg_method,
		      double lambda): 
  VSSpectrumUnfolding(theta_cut,sp_model,ebin,emin,emax),
  m_theta_cut(theta_cut), m_ul_threshold(ul_threshold), 
  m_reg_method(reg_method),  m_lambda(lambda)
{

}

void VSSpectrumUnfoldingML::
reconstruct(const VSAnalysisStage3Data& stage3_data,
	    VSSpectrumData& data) 
{
  VSLimitedErrorsHist<double,double> excess_hist;

  std::vector<VSSpectrumFit::Data> fit_data;

  unsigned nbin = 0;
  double emin = 0;

  const unsigned nrun = stage3_data.nrun();
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      const VSAnalysisStage3Data::RunData& rd = stage3_data.run_data(irun);
      double alpha = 1./rd.sp_noff();
      double offset = (rd.obs_xy()-rd.src_xy()).norm();

      VSSpectrumFit::Data* fd = NULL;

      for(std::vector<VSSpectrumFit::Data>::iterator itr = 
	    fit_data.begin(); itr !=
	    fit_data.end(); ++itr)
	{
	  if(fabs(itr->alpha-alpha) < 1E-3)
	    {
	      fd = &(*itr);
	      break;
	    }	  
	}

      if(fd == NULL)
	{
	  fit_data.resize(fit_data.size()+1);
	  fd = &fit_data.back();
	}

      m_ebin = rd.egy_on_np_hist.binSize();
      nbin = rd.egy_on_np_hist.nBins();
      emin = rd.egy_on_np_hist.loLimit();
      fd->alpha = alpha;

      VSNSpace k = rd.irf_calc().getKernel(rd.egy_on_np_hist.binSize(),
					   rd.egy_on_np_hist.loLimit(),
					   rd.egy_on_np_hist.hiLimit(),
					   offset, m_theta_cut);

      k *= rd.livetime_min*60;

      if(fd->on.ndim() == 0)
	{
	  fd->on = VSSpectrumCalc::toVecND(rd.egy_on_hist.contentsHist());
	  fd->off = VSSpectrumCalc::toVecND(rd.egy_off_hist.contentsHist());
	  fd->krn = k;
	}
      else
	{
	  fd->on += VSSpectrumCalc::toVecND(rd.egy_on_hist.contentsHist());
	  fd->off += 
	    VSSpectrumCalc::toVecND(rd.egy_off_hist.contentsHist());
	  fd->krn += k;
	}

       if(irun == 0) 
 	excess_hist = stage3_data.run_data(irun).egy_on_np_hist;
       else
 	excess_hist += stage3_data.run_data(irun).egy_on_np_hist;

       excess_hist -= stage3_data.run_data(irun).egy_off_np_hist*alpha;      
    }

  VSSpectrumFnPowerLaw pl_fn(2.55,1E-5);

  SpectrumFn* fn = new SpectrumFn(nbin,m_ebin,emin);

  VSAAlgebra::MatrixND s = calcSmoothingMatrix(m_reg_method,m_ebin,nbin,nbin);

  double lambda = 0;

  for(unsigned i = 0; i < fit_data.size(); i++)
    {
      VSAAlgebra::MatrixND k = VSSpectrumCalc::toMatrixND(fit_data[i].krn);
      VSAAlgebra::MatrixND ksq = k*k.getTranspose();
      lambda += (ksq.trace()/s.trace())/fit_data.size();
    }

  // Create scaling matrix ----------------------------------------------------
  VSAAlgebra::MatrixND sm(nbin,nbin);

  for(unsigned ibin = 0; ibin < nbin; ibin++)
    {
      double elo = emin + m_ebin*ibin;
      double ehi = emin + m_ebin*(1.0 + ibin);
      double flux = pl_fn.integralFlux(elo,ehi);
      sm(ibin,ibin) = flux;//std::pow(10,(index+1)*irow*m_ebin);
    }

  LnLFn lnl_reg_fn(fit_data,fn,s,sm,lambda);



  VSSpectrumFit::LnLFn lnl_fn(fit_data,fn);

  VSAAlgebra::VecND param(nbin,1.0);

  for(double gamma = 2; gamma < 4; gamma += 0.05)
    {
      VSAAlgebra::VecND p(nbin);
      VSAAlgebra::VecND p2(nbin);
      pl_fn.setParam(1,gamma);
 
      for(unsigned ibin = 0; ibin < nbin; ibin++)
	{
	  double elo = emin + m_ebin*ibin;
	  double ehi = emin + m_ebin*(1.0 + ibin);
	  double flux = pl_fn.integralFlux(elo,ehi);
	  
	  //      std::cout << elo << " " << ehi << " " << flux << std::endl;
	  p(ibin) = flux/sm(ibin,ibin);
	  p2(ibin) = flux;
	}
 
//       double lnl;
//       VSAAlgebra::MatrixND beta(nbin,1);
//       VSAAlgebra::MatrixND alpha(nbin,nbin);

//       lnl_fn.val(p,lnl,beta,alpha);

      std::cout << "VAL " 
		<< std::setw(15) << gamma 
		<< std::setw(15) << lnl_fn.val(p2) 
		<< std::setw(15) << lnl_reg_fn.val(p) << std::endl;    
      //      std::cout << beta << std::endl;

      if(fabs(gamma-2.5) < 1E-3)
 	param = p;
    }



  VSAMath::NLFitterLM<LnLFn> fitter(lnl_reg_fn);


  double lnl;
  VSAAlgebra::MatrixND beta(nbin,1);
  VSAAlgebra::MatrixND alpha(nbin,nbin);

  lnl_reg_fn.val(param,lnl,beta,alpha);

  for(unsigned ip = 0; ip < param.ndim(); ip++)
    {

      VSAAlgebra::VecND p = param;
      const double dp = 1E-8;

      double y1 = lnl_reg_fn.val(p);
      p(ip) += dp;
      double y2 = lnl_reg_fn.val(p);
      double dy = y2-y1;

      double dyda = dy/dp;
      

      std::cout << "BETA " 
		<< std::setw(5) << ip 
		<< std::setw(15) << param(ip)
		<< std::setw(15) << dy 
		<< std::setw(15) << beta(ip,0) 
		<< std::setw(15) << dyda 
		<< std::setw(15) << beta(ip,0)/dyda
		<< std::endl;
    }

  fitter.initialize(param);
  std::cout << "CHI2 " << fitter.chi2() << std::endl;

  for(unsigned iter = 0; iter < 1; iter++)
    {
      
      for(unsigned i = 0; i < nbin; i++)
	{
	  if(i+20 >= nbin) break;
	  
	  std::cout << iter << " " << i << " " << fitter.chi2();
	  
	  
	  for(unsigned j = 0; j < nbin; j++)
	    {
	      if(j >= i && j < i + 20)
		fitter.free(j);
	      else
		fitter.hold(j);
	    }
	  
	  
	  fitter.initialize(param);
	  fitter.fit();
	  std::cout << " " << fitter.chi2() << std::endl;
	  
	  param = fitter.param();
	}
    }

  fitter.free();
  fitter.initialize(param);
  fitter.fit();

  std::cout << "-------------------" << std::endl;
  for(unsigned i = 0; i < fitter.nparam(); i++)
    std::cout << std::setw(15) << fitter.param(i) 
	      << std::setw(15) << fitter.err(i) 
	      << std::endl;
  std::cout << "CHI2 " << fitter.chi2() << std::endl;

  data.egy_excess_hist = excess_hist;
  data.egy_excess_sp_hist = excess_hist;
  data.flux_hist = data.egy_excess_hist;
  data.flux_bias_hist = data.egy_excess_hist;
  data.flux_rbias_hist = data.egy_excess_hist;

  data.flux_hist.clear();
  data.flux_bias_hist.clear();
  data.flux_rbias_hist.clear();

  for(unsigned ibin = 0; ibin < nbin; ibin++)
    {
      double egy = emin + m_ebin*(0.5+ibin);
      data.flux_hist.accumulate(egy,param(ibin),std::pow(fitter.err(ibin),2));
    }

//   excess = VSSpectrumCalc::toVecND(excess_hist.contentsHist());
//   var = VSSpectrumCalc::toVecND(excess_hist.varianceHist());

//   data.flux_cov_hist = 
//     VSSimple2DHist<double,double>(data.flux_hist.binSize(),
// 				  data.flux_hist.loLimit(),
// 				  data.flux_hist.hiLimit());

//   data.effarea_kernel_nspace = stage3_data.egy_effarea_kernel;
//   reconstruct(krn,excess,var,on,off,alpha,data);

  finalize(data);
}

void VSSpectrumUnfoldingML::
reconstruct(const VSNSpace& kernel_nspace,
	    const VSAAlgebra::VecND& excess,
	    const VSAAlgebra::VecND& var,
	    const VSAAlgebra::VecND& on,
	    const VSAAlgebra::VecND& off,
	    double alpha,
	    VSSpectrumData& data) 
{ 
  VSAAlgebra::MatrixND kernel = VSSpectrumCalc::toMatrixND(kernel_nspace);
  kernel.transpose();

  const unsigned nrow = kernel.nrow();
  const unsigned ncol = kernel.ncol();

  VSAAlgebra::VecND ep = excess;

  // --------------------------------------------------------------------------
  // Normalize kernel and excess

  VSAAlgebra::MatrixND s = 
    calcSmoothingMatrix(m_reg_method,kernel_nspace.axis(0).bin_size,
			nrow,ncol);
  VSAAlgebra::MatrixND k = kernel;

  for(unsigned irow = 0; irow < nrow; irow++)
    {
      if(var[irow] > 0) ep[irow] /= sqrt(var[irow]);

      for(unsigned icol = 0; icol < ncol; icol++)
	{
	  if(var[irow] > 0) k[irow][icol] /= sqrt(var[irow]);
	  //	  else k[irow][icol] = 0;
	}
    }
  
  VSAAlgebra::MatrixND kt = k.getTranspose();
  VSAAlgebra::MatrixND ksq = kt*k;
  
//   ksq = sm.inverse()*ksq;
//   ep = sm.inverse()*ep;

  double lambda = ksq.trace()/s.trace();
 
  double nbin = 0;
  for(unsigned ibin = 0; ibin < excess.ndim(); ibin++)
    if(var(ibin) > 0) nbin++;
  
  Data spd;
  data.lambda_gcc_hist = VSLimitedErrorsHist<double,double>(0.1,-8,0);
  data.lambda_chi2_hist = VSLimitedErrorsHist<double,double>(0.1,-8,0);
  data.lambda_chi2b_hist = VSLimitedErrorsHist<double,double>(0.1,-8,0);
  data.lambda_chi2_ml_hist = VSLimitedErrorsHist<double,double>(0.1,-8,8);
  data.lambda_msbias_hist = VSLimitedErrorsHist<double,double>(0.1,-8,8);
  data.lambda_mvar_hist = VSLimitedErrorsHist<double,double>(0.1,-8,8);
  data.lambda_mse_hist = VSLimitedErrorsHist<double,double>(0.1,-8,8);
  data.lambda_wmse_hist = VSLimitedErrorsHist<double,double>(0.1,-8,8);

  //  std::cout << "NBIN " << nbin << " " << lambda << std::endl;

  for(VSLimitedErrorsHist<double,double>::iterator itr = 
	data.lambda_mse_hist.begin(); itr != data.lambda_mse_hist.end(); ++itr)
    {

      double lam = lambda*std::pow(10,itr->center());

      Data d;
      calcSpectrum(k,s,ep,var,lam,d);
      d.nu = kernel*d.flux;

      VSAAlgebra::VecND mub;
      VSSpectrumFit::LnLFn::mu_bkgnd(on,off,alpha,d.nu,mub);

      unsigned ndf = 0;

      for(unsigned irow = 0; irow < nrow; irow++)
	{
	  if(on(irow)<5) continue;

	  ndf++;

	  double mu = mub(irow);
	  double mu_on = d.nu(irow)+alpha*mu;

	  if(on(irow)) d.chi2_ml += 2*on(irow)*log(on(irow)/mu_on);
	  d.chi2_ml += 2*(mu_on-on(irow));
	}

      std::cout << "LAMBDA " 
		<< std::setw(15) << lam/lambda
		<< std::setw(15) << d.chi2 
		<< std::setw(15) << d.chi2b 
		<< std::setw(15) << d.chi2_ml/ndf
		<< std::setw(15) << d.mse 
		<< std::setw(15) << d.wmse 
		<< std::setw(15) << d.mean_gcc
		<< std::endl;

      data.lambda_gcc_hist.accumulate(itr->center(),d.mean_gcc,0);
      data.lambda_chi2_hist.accumulate(itr->center(),d.chi2,0);
      data.lambda_chi2b_hist.accumulate(itr->center(),d.chi2b,0);
      data.lambda_chi2_ml_hist.accumulate(itr->center(),d.chi2_ml/ndf,0);
      data.lambda_msbias_hist.accumulate(itr->center(),d.msbias,0);
      data.lambda_mvar_hist.accumulate(itr->center(),d.mvar,0);
      data.lambda_mse_hist.accumulate(itr->center(),d.mse,0);
      data.lambda_wmse_hist.accumulate(itr->center(),d.wmse,0);
				      
      data.chi2_mse_graph.addVertex(d.chi2,d.mse);
      data.chi2_ml_mse_graph.addVertex(d.chi2_ml,d.mse);
      data.chi2_wmse_graph.addVertex(d.chi2,d.wmse);
      data.chi2_ml_wmse_graph.addVertex(d.chi2_ml,d.wmse);

      if(itr == data.lambda_mse_hist.begin()) spd = d;
      else if(d.wmse < spd.wmse && d.chi2_ml/ndf < 1)
	{
	  spd = d;
	  spd.lambda = std::pow(10,itr->center());
	}
    }



  if(m_lambda > 0)
    {
      calcSpectrum(k,s,ep,var,m_lambda*lambda,spd);
      spd.nu = kernel*spd.flux;
    }
//   for(unsigned irow = 0; irow < nrow; irow++)
//     std::cout << "NU " << irow << " " << nu(irow)/sqrt(var(irow)) << std::endl;

//     if(var(irow)) nup(irow) /= sqrt(var(irow));



//   VSAAlgebra::MatrixND h = ksq + lambda*s;
//   VSAAlgebra::MatrixND b = h.inverse()*kt;
//   VSAAlgebra::MatrixND bsq = b*b.getTranspose();

//   flux = b*ep;

//   for(unsigned irow = 0; irow < nrow; irow++)
//     flux_var[irow] = bsq[irow][irow];

//   VSAAlgebra::VecND nu = kernel*flux;
//   VSAAlgebra::VecND nup =  nu;
//   for(unsigned irow = 0; irow < nrow; irow++)
//     if(var(irow)) nup(irow) /= sqrt(var(irow));

//   double chi2 = (k*flux-ep)*(k*flux-ep);

  std::cout << "LAMBDA " << spd.lambda << std::endl;
  std::cout << "CHI2   " << spd.chi2 << std::endl;


  data.chi2 = spd.chi2;
  data.lambda = spd.lambda;
  data.egy_excess_sp_hist.clear();
  data.egy_excess_sp_hist.fill(0);

  double chi2 = 0;

  for(unsigned irow = 0; irow < nrow; irow++)
    {
      for(unsigned icol = 0; icol < nrow; icol++)
	data.flux_cov_hist.setBin(irow,icol, spd.flux_cov(irow,icol));

      double sigma = spd.flux[irow]/sqrt(spd.flux_cov(irow,irow));
      double egy = kernel_nspace.axis(0).midCoordUnchecked(irow);

      if(var(irow) > 0) 
	chi2 += std::pow(spd.nu[irow] - excess(irow),2)/var(irow);

//       std::cout << std::setw(20) << egy 
// 		<< std::setw(20) << spd.nu[irow]
// 		<< std::setw(20) << excess(irow)
// 		<< std::setw(20) << sqrt(var(irow))
// 		<< std::setw(20) << chi2 
// 		<< std::endl;


      data.egy_excess_sp_hist.accumulate(egy,spd.nu[irow],0);

      if(sigma > m_ul_threshold)
	data.flux_hist.accumulate(egy, spd.flux[irow], 
				  spd.flux_cov(irow,irow));
      else
	{
	  double ul = 
	    VSAStatistics::heleneUL(spd.flux[irow], 
				    sqrt(spd.flux_cov(irow,irow)), 0.95);
	  data.flux_hist.accumulate(egy,ul,0);
	}


//       std::cout << "FLUX " 
// 		<< std::setw(20) << spd.flux[irow] 
// 		<< std::setw(20) << sqrt(spd.flux_cov(irow,irow)) 
// 		<< std::setw(20) << spd.nu[irow] 
// 		<< std::setw(20) << excess[irow]
// 		<< std::endl;
    }




  for(unsigned irow = 0; irow < nrow;irow++)
    {
      double egy = kernel_nspace.axis(1).midCoordUnchecked(irow);
      data.flux_bias_hist.accumulate(egy, spd.bias[irow], 
				     spd.bias_cov(irow,irow));

      //      double rbias = spd.bias[irow]/spd.flux[irow];
      //      double rbias_var = 

      data.flux_rbias_hist.
	accumulate(egy, spd.bias[irow]/spd.flux[irow], 
		   spd.bias_cov(irow,irow)/std::pow(spd.flux[irow],2));
    }



//   for(unsigned irow = 0; irow < nrow; irow++)
//     std::cout << irow << " " << spd.bias(irow) << " "
// 	      << sqrt(spd.bias_cov(irow,irow)) << std::endl;

}

void VSSpectrumUnfoldingML::
calcSpectrum(const VSAAlgebra::MatrixND& k,
	     const VSAAlgebra::MatrixND& s,
	     const VSAAlgebra::VecND& ep,
	     const VSAAlgebra::VecND& var,
	     double lambda,
	     VSSpectrumUnfolding::Data& spd)
{
  // Create scaling matrix ----------------------------------------------------
  double index = -2.5;

  VSAAlgebra::MatrixND sm(k.nrow(),k.nrow());

  for(unsigned irow = 0; irow < k.nrow(); irow++)
    sm(irow,irow) = std::pow(10,(index+1)*irow*m_ebin);

  // Compute Spectrum ---------------------------------------------------------
  VSAAlgebra::MatrixND kt = k.getTranspose();
  VSAAlgebra::MatrixND ksq = kt*k;

  
  VSAAlgebra::MatrixND h = sm*ksq*sm + lambda*s;
  VSAAlgebra::MatrixND b = h.inverse()*sm*kt;
  VSAAlgebra::MatrixND bsq = b*b.getTranspose();


  spd.flux = sm*b*ep;
  spd.flux_cov = sm*b*b.getTranspose()*sm;

  VSAAlgebra::VecND bp = kt*ep;
  VSAAlgebra::VecND nup = k*spd.flux;
//   VSAAlgebra::VecND nup =  nu;
//   const unsigned nrow = ep.ndim();
//   for(unsigned irow = 0; irow < nrow; irow++)
//     std::cout << "NU " << irow << " " << nup(irow) << std::endl;

//     if(var(irow)) nup(irow) /= sqrt(var(irow));

  spd.chi2 = (k*spd.flux-ep)*(k*spd.flux-ep);
  spd.lambda = lambda;

  // Compute bias -------------------------------------------------------------
  VSAAlgebra::MatrixND ab = -ksq*sm - lambda*s;
  VSAAlgebra::MatrixND bb = k;
  //  VSAAlgebra::MatrixND c = ab.inverse()*bb*(-1);
  VSAAlgebra::MatrixND c = sm*ab.inverse()*bb*(-1);
  VSAAlgebra::MatrixND c2 = sm*(c*k*c-c);

  spd.bias = c*(nup-ep);
  //  spd.bias = sm*c*(nup-ep);
  spd.bias_cov = (c*k*c-c)*(c*k*c-c).getTranspose();
  //  spd.bias_cov = c2*c2.getTranspose();

  VSAAlgebra::MatrixND u = c*c.getTranspose();

  spd.msbias = 0;
  spd.mvar = 0;
  spd.mse = 0;
  spd.wmse = 0;
  spd.mean_gcc = 0;



  VSAAlgebra::MatrixND vx1 = spd.flux_cov;
  VSAAlgebra::SVD svd(vx1);

//   std::cout << "SVD W " << std::endl;
//   std::cout << "SVD THRESHOLD " << svd.thresh() << std::endl;
//   std::cout << svd.w() << std::endl;

  VSAAlgebra::MatrixND vx2 = svd.inverse();

  VSAAlgebra::SVD svd_ksq(ksq);
  
  VSAAlgebra::MatrixND u1 = svd_ksq.u();
  VSAAlgebra::MatrixND u1t = svd_ksq.v().getTranspose();

  VSAAlgebra::MatrixND lam(ksq.nrow(),ksq.nrow());

  for(unsigned irow = 0; irow < ksq.nrow(); irow++)
    lam(irow,irow) = 1/sqrt(svd_ksq.w()(irow));

  VSAAlgebra::SVD svd_s(lam*u1t*s*u1*lam);
  VSAAlgebra::VecND fc1 = lam*u1t*bp;
  VSAAlgebra::VecND fc2 = svd_s.v().getTranspose()*lam*u1t*bp;

  unsigned ndf1 = 0;
  unsigned ndf2 = 0;

  for(unsigned irow = 0; irow < fc1.ndim(); irow++)
    if(fc1(irow) > 2) ndf1++;

  for(unsigned irow = 0; irow < fc2.ndim(); irow++)
    if(fc2(irow) > 2) ndf2++;

//   std::cout << "NDF1 " << ndf1 << std::endl;
//   std::cout << "NDF2 " << ndf2 << std::endl;

  spd.ndf = ndf1;
  double chi2 = 0;

  const unsigned nrow = spd.bias.ndim();
  VSAAlgebra::VecND rho(nrow);

  double nbin = 0;

  spd.gcc_hist = VSLimitedErrorsHist<double,double>(1,0,nrow);

  for(unsigned irow = 0; irow < nrow; irow++)
    {

      //      spd.mse += (1./(double)nrow)*(u(irow,irow) + std::pow(spd.bias(irow),2));
      spd.mse += (1./(double)nrow)*(spd.flux_cov(irow,irow) + 
				    std::pow(spd.bias(irow),2));
 
      double vx = std::min(1.,std::pow(vx1(irow,irow)*vx2(irow,irow),-1));

      rho(irow) = sqrt(1-vx);

      spd.mean_gcc += rho(irow)/nrow;
      spd.gcc_hist.accumulate(irow+0.5,sqrt(1-vx),0);

//       for(unsigned icol = 0; icol < nrow; icol++)
// 	std::cout << icol << " " 
// 		  << std::setw(15) 
// 		  << vx1(irow,icol)/sqrt(vx1(irow,irow)*vx1(icol,icol))
// 		  << std::endl;

//       std::cout << std::setw(5) << irow 
// 		<< std::setw(15) << rho(irow) 
// 		<< std::setw(15) << vx1(irow,irow) 
// 		<< std::setw(15) << vx2(irow,irow) 
// 		<< std::setw(15) << 1./(vx1(irow,irow)*vx2(irow,irow)) 
// 		<< std::setw(15) << vx 
// 		<< std::endl;

      
//       std::cout << std::setw(5) << irow 
// 		<< std::setw(15) << var(irow)
// 		<< std::setw(15) << u(irow,irow)
// 		<< std::setw(15) << spd.flux_cov(irow,irow)
// 		<< std::setw(15) << std::pow(spd.bias(irow),2)
// 		<< std::setw(15) << spd.msbias
// 		<< std::setw(15) << spd.mvar
// 		<< std::setw(15) << spd.mse
// 		<< std::setw(15) << spd.wmse
// 		<< std::setw(15) << spd.chi2b
// 		<< std::endl;

      if(var(irow) > 0)
	{
	  double smc = 1;//std::pow(sm(irow,irow),-2);

	  spd.chi2b += 
	    std::pow(spd.bias(irow),2)/spd.bias_cov(irow,irow);

	  spd.mvar += (1./(double)nrow)*spd.flux_cov(irow,irow)*smc;
	  spd.msbias += (1./(double)nrow)*std::pow(spd.bias(irow),2)*smc;
	  spd.wmse += 
	    (1./(double)nrow)*(spd.flux_cov(irow,irow)*smc + 
			       std::pow(spd.bias(irow),2)*smc)/var(irow);

	  nbin++;
	}
    }

  spd.chi2b /= nbin;

}


#endif

// ============================================================================
// VSSpectrumCalcFactory
// ============================================================================
VSSpectrumCalcFactory::Options::Options():
  np_method("cfm"),
  reg_matrix("pol1"),
  lambda(0),
  scaling_index(1),
  cfm_niter(4),
  spectrum_model("powerlaw"),
  spectrum_log10_enorm(-0.5),
  spectrum_bin_scheme(4.,0.),
  spectrum_bin_width(0.0625),
  theta_cut(0.12),
  upper_limit_threshold(2.5)
{

}

std::auto_ptr<VSSpectrumCalcFactory> VSSpectrumCalcFactory::s_instance;

VSSpectrumCalcFactory::Options VSSpectrumCalcFactory::s_default_options = 
VSSpectrumCalcFactory::Options();

VSSpectrumCalcFactory::VSSpectrumCalcFactory(const Options& opt):
  m_options(opt)
{
  m_egy_bin_width = 1./(double)m_options.spectrum_bin_scheme.first;

  double ecenter = m_options.spectrum_bin_scheme.second;
  const double emin = -1.3;
  const double emax = 2.1;

  m_egy_min = ecenter - m_egy_bin_width/2. -
    lround(fabs(ecenter-emin)/m_egy_bin_width)*m_egy_bin_width;

  m_egy_max = ecenter + m_egy_bin_width/2. +
    lround(fabs(ecenter-emax)/m_egy_bin_width)*m_egy_bin_width;
}

VSSpectrumCalcFactory::~VSSpectrumCalcFactory()
{

}

VSSpectrumCalcFactory* VSSpectrumCalcFactory::getInstance()
{
  if(s_instance.get() == 0)s_instance.reset(new VSSpectrumCalcFactory());
  return s_instance.get();
}

VSSpectrumCalc* VSSpectrumCalcFactory::create()
{
  return new VSSpectrumCalc(m_options.theta_cut,
			    m_options.spectrum_model,
			    m_options.spectrum_log10_enorm);
}

VSSpectrumUnfolding* VSSpectrumCalcFactory::createSpectrumUnfolding()
{
  if(m_options.np_method == "chi2_regularization")
    {
      return new VSSpectrumUnfoldingChi2(m_options.theta_cut,
					 m_options.spectrum_model,
					 m_egy_bin_width,
					 m_egy_min,
					 m_egy_max,
					 m_options.spectrum_log10_enorm,
					 m_options.upper_limit_threshold,
					 m_options.reg_matrix,
					 m_options.lambda,
					 m_options.scaling_index);
    }
  else if(m_options.np_method == "reg_ml")
    {
      return 0;
      // return new VSSpectrumUnfoldingML(m_options.theta_cut,
      // 				       m_options.spectrum_model,
      // 				       m_egy_bin_width,
      // 				       m_egy_min,
      // 				       m_egy_max,
      // 				       m_options.spectrum_log10_enorm,
      // 				       m_options.upper_limit_threshold,
      // 				       m_options.reg_matrix,
      // 				       m_options.lambda);
    }
  else if(m_options.np_method == "cfm")
    {
      if(m_options.cfm_niter <= 0)
	{
	  std::cerr << __PRETTY_FUNCTION__
		    << ": Invalid number of iterations specified. "
		    << std::endl;
	  exit(EXIT_FAILURE);
	}

      return new VSSpectrumUnfoldingCFM(m_options.theta_cut,
					m_options.spectrum_model,
					m_egy_bin_width,
					m_egy_min,
					m_egy_max,
					m_options.spectrum_log10_enorm,
					m_options.upper_limit_threshold,
					m_options.cfm_niter);
    }
  else
    {
      std::cerr 
	<< "Unknown non-parametric spectral reconstruction method: "
	<< m_options.np_method << std::endl;
      exit(EXIT_FAILURE);
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSSpectrumCalcFactory::
configure(VSOptions& options, const std::string& profile, 
	  const std::string& opt_prefix)
{
  options.addCatagory(OPTNAME(opt_prefix,"spectrum"),
		      "Options related to spectral reconstruction.");

  options.findWithValue(OPTNAME(opt_prefix,"spectrum_np_method"),
			s_default_options.np_method,
			"Set the method for non-parametric "
			"spectral reconstruction (unfolding).  Options are: "
			"\"chi2_regularization\" (Chi-Squared Functional w/ "
			"Regularization), "
			"\"ml_regularization\" (Maximum Likelihood Functional "
			"w/ Regularization), "
			"and \"cfm\" (Correction Factors Method)",
			"s3_spectrum");

  options.findWithValue(OPTNAME(opt_prefix,"reg_param"),
			s_default_options.lambda,
			"Set the regularization parameter to be used for "
			"spectral reconstruction.",
			"s3_spectrum");

  options.findWithValue(OPTNAME(opt_prefix,"index"),
			s_default_options.scaling_index,
			"",
			"s3_spectrum");

  options.findWithValue(OPTNAME(opt_prefix,"reg_matrix"),
			s_default_options.reg_matrix,
			"Set the Tikhonov regularization smoothing matrix "
			"used by the unfolding with regularization methods. "
			"Each matrix minimizes a higher successively higher "
			"order derivative of the spectral solution. "
			"Options are: "
			"\"pol0\", \"pol1\", \"pol2\" ",
			"s3_spectrum");

  options.findWithValue(OPTNAME(opt_prefix,"cfm_niter"),
			s_default_options.cfm_niter,
			"Set the number of iterations to be used by the "
			"correction factor method.",
			"s3_spectrum");

  options.findWithValue(OPTNAME(opt_prefix,"spectrum_model"),
			s_default_options.spectrum_model,
			"Set the parameterized spectral model used by the "
			"forward-folding analysis.",
			"s3_spectrum");

  options.findWithValue(OPTNAME(opt_prefix,"spectrum_log10_enorm"),
			s_default_options.spectrum_log10_enorm,
			"Set the normalization energy in log10(E/TeV) of "
			"the parameterized spectral model.",
			"s3_spectrum");
  

  options.findWithValue(OPTNAME(opt_prefix,"spectrum_theta_cut"), 
			s_default_options.theta_cut,
			"Set the aperture size in degrees used for spectral "
			"reconstruction.","s3_spectrum");

  options.findWithValue(OPTNAME(opt_prefix,"upper_limit_threshold"), 
			s_default_options.upper_limit_threshold,
			"Set the threshold in sigma below which the flux "
			"in an energy bin will be "
			"reported as an upper limit.","s3_spectrum");


  options.findWithValue(OPTNAME(opt_prefix,"spectrum_bin_scheme"), 
			s_default_options.spectrum_bin_scheme,
			"Set the spectral bin scheme by specifying the number "
			"of bins per decade and a bin center reference point "
			"in log_10(E [TeV]) -- e.g. 8/0 for 8 bins per decade "
			"centered on 1 TeV.","s3_spectrum");
  
}
