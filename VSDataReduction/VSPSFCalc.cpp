//-*-mode:c++; mode:font-lock;-*-

/*! \file VSPSFCalc.cpp

  Model gamma-ray psf from simulations.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       09/11/2008

  $Id: VSPSFCalc.cpp,v 3.9 2010/06/22 00:00:52 matthew Exp $

*/

#include <VSPSFCalc.hpp>
#include <VSANonlinearFitting.hpp>

using namespace VERITAS;

// ============================================================================
// VSBinnedPSFFn
// ============================================================================
VSBinnedPSFFn::VSBinnedPSFFn(double dthsq):
  VSAFunction::ParamFn<double>(3), m_counts(1), m_dthsq(dthsq)
{
  
}

VSBinnedPSFFn::VSBinnedPSFFn(double dthsq, double counts): 
  VSAFunction::ParamFn<double>(3), m_counts(counts), m_dthsq(dthsq)
{
  
}

double VSBinnedPSFFn::val(const double& x) const 
{ 
  return val(x,param());
}

double VSBinnedPSFFn::val(const double& x, const VSAAlgebra::VecND& a) const 
{ 
  double s1i = std::pow(a(0),-2);
  double s2i = std::pow(a(1),-2);

  double xmax = 0.08;

  double s = 
    m_counts*std::pow(a(2)*(1-exp(-xmax*s1i/2.)) + 
		      (1-a(2))*(1-exp(-xmax*s2i/2.)),-1);

  return s*std::pow(2*M_PI,-1)*
    (a(2)*s1i*exp(-x*s1i/2.) + (1-a(2))*s2i*exp(-x*s2i/2.))*m_dthsq;
}

void VSBinnedPSFFn::dyda(const double& x, VSAAlgebra::VecND& dyda) const 
{
  VSBinnedPSFFn::dyda(x,param(),dyda);
}

void VSBinnedPSFFn::dyda(const double& x, const VSAAlgebra::VecND& a, 
		   VSAAlgebra::VecND& dyda) const 
{
  dyda.resize(3);

  double s1i = std::pow(a(0),-2);
  double s2i = std::pow(a(1),-2);

  double xmax = 0.08;

  double s = 
    m_counts*std::pow(a(2)*(1-exp(-xmax*s1i/2.)) + 
		      (1-a(2))*(1-exp(-xmax*s2i/2.)),-1);

  dyda(0) = s*std::pow(2*M_PI,-1)*
    a(2)*s1i*exp(-x*s1i/2.)*(x-2*a(0)*a(0))/std::pow(a(0),3);
  dyda(1) = s*std::pow(2*M_PI,-1)*
    (1-a(2))*s2i*exp(-x*s2i/2.)*(x-2*a(1)*a(1))/std::pow(a(1),3);
  dyda(2) = s*std::pow(2*M_PI,-1)*
    (s1i*exp(-x*s1i/2.) - s2i*exp(-x*s2i/2.));
    
  dyda *= m_dthsq;
}

double VSBinnedPSFFn::integrate(const VSAAlgebra::VecND& a, 
				double xlo, double xhi)
{
  double s1i = std::pow(a(0),-2);
  double s2i = std::pow(a(1),-2);

  return (exp(-xlo*s2i/2.) - exp(-xhi*s2i/2.) + 
    a(2)*(exp(-xlo*s1i/2.) - exp(-xhi*s1i/2.) - 
	  exp(-xlo*s2i/2.) + exp(-xhi*s2i/2.)));
}

// ============================================================================
// VSPSFEnergyFn
// ============================================================================
VSPSFEnergyFn::VSPSFEnergyFn(double dthsq, unsigned negy, double dloge,
			     double logemin): 
  VSAFunction::ParamFn<VSACoord::CoordND>(7), m_dthsq(dthsq),
  m_ebin(dloge,logemin)
{
  m_fn.resize(negy,VSBinnedPSFFn(dthsq));
}

void VSPSFEnergyFn::setCounts(unsigned iegy, double counts)
{
  vsassert(iegy < m_fn.size());
  m_fn[iegy].setCounts(counts);
}

double VSPSFEnergyFn::val(const VSACoord::CoordND& x) const 
{ 
  return val(x,param());
}

double VSPSFEnergyFn::val(const VSACoord::CoordND& x, 
			  const VSAAlgebra::VecND& a) const 
{ 
  VSAAlgebra::VecND afn(3);

  unsigned iegy = m_ebin.valToBin(x[1]);

  afn(0) = a(0) + x[1]*a(1) + std::pow(x[1],2)*a(2) + std::pow(x[1],3)*a(3);
  afn(1) = a(4) + x[1]*a(5);
  afn(2) = a(6);

  return m_fn[iegy].val(x[0],afn);
}

void VSPSFEnergyFn::dyda(const VSACoord::CoordND& x, VSAAlgebra::VecND& dyda) 
  const 
{
  VSPSFEnergyFn::dyda(x,param(),dyda);
}

void VSPSFEnergyFn::dyda(const VSACoord::CoordND& x, 
			 const VSAAlgebra::VecND& a, 
			 VSAAlgebra::VecND& dyda) const 
{
  dyda.resize(nparm());

  VSAAlgebra::VecND afn(3);

  afn(0) = a(0) + x[1]*a(1) + std::pow(x[1],2)*a(2) + std::pow(x[1],3)*a(3);
  afn(1) = a(4) + x[1]*a(5);
  afn(2) = a(6);

  VSAAlgebra::VecND dfdn;

  unsigned iegy = m_ebin.valToBin(x[1]);

  m_fn[iegy].dyda(x[0],afn,dfdn);

  dyda(0) = dfdn(0);
  dyda(1) = dfdn(0)*x[1];
  dyda(2) = dfdn(0)*std::pow(x[1],2);
  dyda(3) = dfdn(0)*std::pow(x[1],3);
  dyda(4) = dfdn(1);
  dyda(5) = dfdn(1)*x[1];
  dyda(6) = dfdn(2);
}

// ============================================================================
// VSPSFCalc
// ============================================================================
void VSPSFCalc::Data::EnergyData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeScalar("log10_egy_tev",log10_egy_tev);

  thetasq_hist.save(writer->writeStruct("thetasq_hist"));
  thetasq_fit_hist.save(writer->writeStruct("thetasq_fit_hist"));
}

void VSPSFCalc::Data::save(VSOctaveH5WriterStruct* writer) const
{
  egymc_thetasq_hist.save(writer->writeStruct("egymc_thetasq_hist"));
  sigma1_hist.save(writer->writeStruct("sigma1_hist"));
  sigma2_hist.save(writer->writeStruct("sigma2_hist"));
  alpha_hist.save(writer->writeStruct("alpha_hist"));

  sigma1_fit_hist.save(writer->writeStruct("sigma1_fit_hist"));
  sigma2_fit_hist.save(writer->writeStruct("sigma2_fit_hist"));
  alpha_fit_hist.save(writer->writeStruct("alpha_fit_hist"));
  th68_fit_hist.save(writer->writeStruct("th68_fit_hist"));

  const unsigned negy = m_data.size();
  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("egy", negy);
  for(unsigned iegy = 0; iegy < negy; iegy++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(iegy);
      vsassert(ws);  

      m_data[iegy].save(ws);
      delete ws;
    }

  delete wc;
}

VSPSFCalc::VSPSFCalc():
  m_sigma1_offset_hist(),
  m_sigma2_offset_hist(),
  m_alpha_offset_hist(),
  m_th68_offset_hist(),
  m_sigma1_offset_nspace(),
  m_sigma2_offset_nspace(),
  m_alpha_offset_nspace(),
  m_log10_egybin(),
  m_log10_egylo(),
  m_log10_egyhi(),
  m_negy(),
  m_data(), m_data_ptr()
{

}

VSPSFCalc::~VSPSFCalc()
{
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++) delete m_data[idata];
}

void VSPSFCalc::loadSimInfo(VSSimInfoData* sim_info)
{
  vsassert(m_log10_egybin > 0 && (m_log10_egylo != m_log10_egyhi));

  m_data_ptr = NULL;

  for(std::vector<Data*>::iterator itr = m_data.begin();
      itr != m_data.end(); ++itr)
    {
      double doff = fabs(sim_info->wobble_theta_deg-(*itr)->offset_deg);
      
      if(doff < 0.05) 
	{
	  m_data_ptr = *itr;
	  break;
	}
      else if(sim_info->wobble_theta_deg < (*itr)->offset_deg)
	{
	  m_data_ptr = new Data(sim_info->wobble_theta_deg,m_log10_egybin,
				m_log10_egylo,m_log10_egyhi);
	  m_data.insert(itr,m_data_ptr);
	  break;
	}
    }

  if(m_data_ptr == NULL)
    {
      m_data_ptr = new Data(sim_info->wobble_theta_deg,m_log10_egybin,
			    m_log10_egylo,m_log10_egyhi);
      m_data.push_back(m_data_ptr);
    }
}

void VSPSFCalc::loadHeader(const VSHeaderSimulationDatum& sim_header)
{

}

void VSPSFCalc::setEnergyBinning(double egybin, double egylo, double egyhi)
{
  if(m_log10_egybin == 0)
    {
      m_log10_egybin = egybin;
      m_log10_egylo = egylo;
      m_log10_egyhi = egyhi;
      m_negy = lround((egyhi-egylo)/egybin);

      m_sigma1_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
      m_sigma2_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
      m_alpha_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
      m_th68_offset_hist = 
	VSSimple2DHist<double,double>(egybin,egylo,egyhi,0.04,0,1.8);
    }
}

void VSPSFCalc::accumulate(double energy_tev, double thetasq)
{
  m_data_ptr->egymc_thetasq_hist.accumulate(std::log10(energy_tev),
					    thetasq);

  if(thetasq < 0.08)
    m_data_ptr->egymc_hist.accumulate(std::log10(energy_tev));
}

void VSPSFCalc::fit()
{
  const unsigned negy = m_sigma1_offset_hist.nXBins();
  std::vector<VSAFunction::Spline> sigma1_spline(negy);
  std::vector<VSAFunction::Spline> sigma2_spline(negy);
  std::vector<VSAFunction::Spline> alpha_spline(negy);
  std::vector<VSAFunction::Spline> th68_spline(negy);

  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      std::cout << std::string(__PRETTY_FUNCTION__)  
		<< ": Fitting PSF Model " 
		<< " Offset  " << std::setw(15) << m_data[idata]->offset_deg
		<< std::endl;

      fit(m_data[idata]);
      
      for(unsigned iegy = 0; iegy < negy; iegy++)
	{
	  double log10_egy = m_sigma1_offset_hist.xBinToCenter(iegy);
	  double sigma1 = 
	    m_data[idata]->sigma1_fit_hist.countForVal(log10_egy);
	  double sigma2 = 
	    m_data[idata]->sigma2_fit_hist.countForVal(log10_egy);
	  double alpha =
	    m_data[idata]->alpha_fit_hist.countForVal(log10_egy);
	  double th68 =
	    m_data[idata]->th68_fit_hist.countForVal(log10_egy);

	  sigma1_spline[iegy].setPoint(m_data[idata]->offset_deg,sigma1);
	  sigma2_spline[iegy].setPoint(m_data[idata]->offset_deg,sigma2);
	  alpha_spline[iegy].setPoint(m_data[idata]->offset_deg,alpha);
	  th68_spline[iegy].setPoint(m_data[idata]->offset_deg,th68);
	}

    }

  for(unsigned iegy = 0; iegy < negy; iegy++)
    {
      sigma1_spline[iegy].spline();
      sigma2_spline[iegy].spline();
      alpha_spline[iegy].spline();
      th68_spline[iegy].spline();

      for(unsigned iy = 0; iy < m_sigma1_offset_hist.nYBins(); iy++)
	{
	  double offset = m_sigma1_offset_hist.yBinToCenter(iy);
	  m_sigma1_offset_hist.setBin(iegy,iy,sigma1_spline[iegy].val(offset));
	  m_sigma2_offset_hist.setBin(iegy,iy,sigma2_spline[iegy].val(offset));
	  m_alpha_offset_hist.setBin(iegy,iy,alpha_spline[iegy].val(offset));
	  m_th68_offset_hist.setBin(iegy,iy,th68_spline[iegy].val(offset));
	}
    }

  m_sigma1_offset_nspace = VSNSpace(m_sigma1_offset_hist);  
  m_sigma2_offset_nspace = VSNSpace(m_sigma2_offset_hist);  
  m_alpha_offset_nspace = VSNSpace(m_alpha_offset_hist);  
}

void VSPSFCalc::fit(Data* data)
{
  const unsigned negy = data->egymc_thetasq_hist.nXBins();
  const double dthsq = data->egymc_thetasq_hist.yBinSize()*M_PI;
  const double dloge = data->egymc_thetasq_hist.xBinSize();


  VSAMath::Data<VSACoord::CoordND> fit_data;
  VSPSFEnergyFn fn(dthsq,negy,dloge,data->egymc_thetasq_hist.xLoLimit());
  for(unsigned ix = 0; ix < negy; ix++)
    {
      double log10_egy = data->egymc_thetasq_hist.xBinToCenter(ix);
      double count = data->egymc_hist.countForVal(log10_egy);

      fn.setCounts(ix,count);
      Data::EnergyData egy_data(log10_egy);
      egy_data.thetasq_fit_hist.fill(0);
      
      for(unsigned iy = 0; iy < data->egymc_thetasq_hist.nYBins(); iy++)
	{
	  double thsq = data->egymc_thetasq_hist.yBinToCenter(iy);
	  egy_data.thetasq_hist.
	    accumulate(thsq,data->egymc_thetasq_hist.count(ix,iy),
		       data->egymc_thetasq_hist.count(ix,iy));	  
	}

      //      VSAMath::Data<double> fit_data;
      VSLimitedErrorsHist<double,double>& h = egy_data.thetasq_hist;
      for(VSLimitedErrorsHist<double,double>::iterator itr = h.begin();
	  itr != h.end(); ++itr)
	{
	  VSACoord::CoordND c(2);
	      
	  c[0] = itr->center();
	  c[1] = log10_egy;

	  VSAMath::DataPoint<VSACoord::CoordND> p(c,itr->count(),itr->err());
	  if(count > 10) fit_data.insert(p);
	}

      data->m_data.push_back(egy_data);

	  // 	  VSBinnedPSFFn fn(dthsq,count);
// 	  fn.setParam(0,0.08);
// 	  fn.setParam(1,0.15);
// 	  fn.setParam(2,0.5);
	  
// 	  VSAMath::NLFitter< 
// 	  VSAFunction::LnPoissonLikelihood<VSBinnedPSFFn,double> >*
// 	    fitter = VSAMath::NLFitterFactory::createLnL(fit_data,&fn);

// 	  //	  VSAMath::FitLM<VSBinnedPSFFn,double> fitter(fit_data,&fn);

// 	  try
// 	    {
// 	      fitter->initialize();
// 	      fitter->scan(0,0.01,0.2,0.01);
// 	      fitter->scan(1,0.1,0.2,0.01);
// 	      fitter->setLoBound(2,0);
// 	      fitter->setHiBound(2,1.0);
// 	      fitter->fit();
// 	    }
// 	  catch(const std::exception& e)
// 	    {
// 	      std::cout << e.what() << std::endl;
// 	    }
	  
// 	  VSAAlgebra::VecND p = fitter->param();
// 	  VSAAlgebra::MatrixND cov = fitter->cov();

// 	  delete fitter;

// 	  std::cout << std::setw(15) << log10_egy 
// 		    << std::setw(15) << p(0) 
// 		    << std::setw(15) << sqrt(cov(0,0)) 
// 		    << std::setw(15) << p(1) 
// 		    << std::setw(15) << sqrt(cov(1,1))
// 		    << std::setw(15) << p(2) 
// 		    << std::setw(15) << sqrt(cov(2,2))
// 		    << std::endl;


// 	  data->sigma1_hist.accumulate(log10_egy,fabs(p(0)),cov(0,0));
// 	  data->sigma2_hist.accumulate(log10_egy,fabs(p(1)),cov(1,1));
// 	  data->alpha_hist.accumulate(log10_egy,p(2),cov(2,2));

// 	  for(VSLimitedErrorsHist<double,double>::iterator itr = 
// 		egy_data.thetasq_fit_hist.begin();
// 	      itr != egy_data.thetasq_fit_hist.end(); ++itr)
// 	    {
// 	      double v = fn.val(itr->center(),p);
// 	      egy_data.thetasq_fit_hist.accumulate(itr->center(),v,0);
// 	    }
// 	}

      
    }


  fn.setParam(0,0.045);
  fn.setParam(1,-0.027);
  fn.setParam(2,0.007);
  fn.setParam(3,0.);
  fn.setParam(4,0.14);
  fn.setParam(5,0.0);
  fn.setParam(6,0.5);  

  VSAMath::NLFitter< 
  VSAFunction::LnPoissonLikelihood<VSPSFEnergyFn,VSACoord::CoordND> >*
    fitter = VSAMath::NLFitterFactory::createLnL(&fn,fit_data);

  fitter->initialize();
  fitter->fit();

  VSAAlgebra::VecND param = fitter->param();
  VSAAlgebra::MatrixND cov = fitter->cov();

  delete fitter;
  
//   for(unsigned ip = 0; ip < param.ndim(); ip++)
//     {
//       std::cout << std::setw(10) << ip 
// 		<< std::setw(20) << param(ip) 
// 		<< std::setw(20) << sqrt(cov(ip,ip)) 
// 		<< std::endl;
//     }

  for(unsigned ix = 0; ix < negy; ix++)
    {      
      double log10_egy = data->egymc_thetasq_hist.xBinToCenter(ix);

      double sigma1 = param(0) + param(1)*log10_egy + 
	param(2)*std::pow(log10_egy,2) + param(3)*std::pow(log10_egy,3);

      double sigma2 = param(4) + param(5)*log10_egy;
      double alpha = param(6);

      data->sigma1_fit_hist.accumulate(log10_egy,sigma1,0);
      data->sigma2_fit_hist.accumulate(log10_egy,sigma2,0);
      data->alpha_fit_hist.accumulate(log10_egy,alpha,0);

      VSBinnedPSFFn th2_fn(dthsq);
      VSAAlgebra::VecND a(3);
      a(0) = sigma1;
      a(1) = sigma2;
      a(2) = alpha;

      double th68 = 0;
      double s = 0;
      while(s < 0.68)
	{
	  s = th2_fn.integrate(a,0,std::pow(th68,2));
	  if(th68 > 1) break;
	  th68+=0.0002;
	}

      data->th68_fit_hist.accumulate(log10_egy,th68,0);

      for(unsigned iy = 0; iy < data->egymc_thetasq_hist.nYBins(); iy++)
	{
	  double thsq = data->egymc_thetasq_hist.yBinToCenter(iy);

	  VSACoord::CoordND c(2);

	  c[0] = thsq;
	  c[1] = log10_egy;

	  double v = fn.val(c,param);

	  data->m_data[ix].thetasq_fit_hist.accumulate(thsq,v,0);
	}
    }
}

void VSPSFCalc::save(VSOctaveH5WriterStruct* writer) const
{
  const unsigned ndata = m_data.size();

  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("offset", ndata);
  vsassert(wc);

  for(unsigned idata = 0; idata < ndata; idata++)
    {      
      VSOctaveH5WriterStruct* ws = wc->writeStruct(idata);
      vsassert(ws); 
      m_data[idata]->save(ws);
      delete ws;
    }
  delete wc;

  m_sigma1_offset_hist.save(writer->writeStruct("sigma1_offset_hist"));
  m_sigma2_offset_hist.save(writer->writeStruct("sigma2_offset_hist"));
  m_alpha_offset_hist.save(writer->writeStruct("alpha_offset_hist"));
  m_th68_offset_hist.save(writer->writeStruct("th68_offset_hist")); 
}
