#include <VSAnalysisStage3Data.hpp>
#include <VSSpectrumBiasCalc.hpp>
#include <VSSourceInjector.hpp>

using namespace VERITAS;

VSHistRNG::VSHistRNG(const VSLimitedErrorsHist<double,double>& hist)
{
  std::vector<std::pair<double,double> > cdf;

  VSLimitedErrorsHist<double,double> chist = hist.getCumulativeHist();

  double sum = hist.sum();

  for(VSLimitedErrorsHist<double,double>::iterator itr = chist.begin(); 
      itr != chist.end(); ++itr)
    {
      if(itr->count() && cdf.size() == 0)
	cdf.push_back(std::make_pair(itr->val(),0));

      cdf.push_back(std::make_pair(itr->val()+chist.binSize(),
				   itr->count()/sum));
    }
  
  RandomNumbers::GenerateInverseCDF(cdf,10000);
  m_icdf = cdf;
}

double VSHistRNG::rnd(RandomNumbers* rng) const
{
  return rng->InverseCDF(m_icdf);
}

void  VSHistRNG::fill(RandomNumbers* rng,
		      VSLimitedErrorsHist<double,double>& h,
		      double mu) const
{
  const unsigned nevent = rng->Poisson(mu); 
  for(unsigned ievent = 0; ievent < nevent; ievent++)
    h.accumulate(rnd(rng));
}

VSSpectrumBiasCalc::Options VSSpectrumBiasCalc::s_default_options = 
VSSpectrumBiasCalc::Options();

VSSpectrumBiasCalc::VSSpectrumBiasCalc(RandomNumbers* rng,
				       const Options& opt): 
  m_rng(rng), m_options(opt)
{
  VSSpectrumCalcFactory* spcf = VSSpectrumCalcFactory::getInstance();
  m_spectrum_calc = spcf->create();
}

void VSSpectrumBiasCalc::Data::save(VSOctaveH5Writer* writer) const
{
  sp_data.save(writer->writeStruct("sp_data"));
  emc_erec_hist.save(writer->writeStruct("emc_erec_hist"));
  //  sim_src.save(writer->writeStruct("sim_src"));

  bias_hist.save(writer->writeStruct("bias_hist"));
  rbias_hist.save(writer->writeStruct("rbias_hist"));
}

void VSSpectrumBiasCalc::analyze(const std::string& s3_file)
{
  VSAnalysisStage3Data data;

  VSOctaveH5Reader* reader = new VSOctaveH5Reader(s3_file);

  data.load(reader->readStruct("stage3"));

//   for(VSLimitedErrorsHist<double,double>::iterator itr = 
// 	data.egy_off_hist.begin(); itr != data.egy_off_hist.end(); ++itr)
//     {
//       std::cout << itr->center() << " " << itr->count() << std::endl;
//     }

  const std::pair<double,double> ring_radius(0.4,0.5);

  std::vector<VSSourceInjector*> src_injector;

  const unsigned nrun = data.nrun();
  for(unsigned irun = 0; irun < nrun; irun++)
    {
      std::cout << "Initializing run " << irun << std::endl;

      VSSourceInjector* si = 
	new VSSourceInjector(m_rng,m_spectrum_calc->thetaCut(),ring_radius);

      si->initialize(&data.run_data(irun));
      src_injector.push_back(si);
    }



  //    std::cout << data.run_data(irun).noff() << std::endl;


  VSHistRNG hrng(data.egy_off_hist);

  VSLimitedErrorsHist<double,double> mchist = data.egy_off_np_hist;
  mchist.clear();
  mchist.fill(0);

  
  m_sim_dfde_hist = mchist;
  m_sim_edfde_hist = mchist;
  m_sim_e2dfde_hist = mchist;

  m_bias_hist.resize(mchist.nBins(),
		     VSLimitedErrorsHist<double,double>(0.02,-1,1));

  m_bias_stat.resize(mchist.nBins());

  m_mean_bias_hist = mchist;
  m_mean_dev_hist = mchist;
  m_mean_bias_dev_hist = mchist;
  m_ntot_hist = mchist;
  m_ncov_hist = mchist;
  m_coverage_hist = mchist;

  src_injector.front()->fillIntegralHist(mchist);

  m_sim_flux_hist = mchist;

   for(VSLimitedErrorsHist<double,double>::iterator itr =
	 m_sim_flux_hist.begin(); itr != m_sim_flux_hist.end(); ++itr)
     {
       double flux = itr->count();
       double egy = std::pow(10,itr->center());
       double lo = itr->val();
       double hi = itr->val()+m_sim_flux_hist.binSize();
       double de = std::pow(10,hi) - std::pow(10,lo);

       m_sim_dfde_hist.accumulate(itr->center(), flux/de, 0);
       m_sim_edfde_hist.accumulate(itr->center(), flux/de*egy, 0);
       m_sim_e2dfde_hist.accumulate(itr->center(), flux/de*egy*egy, 0);
     }

  m_data.resize(m_options.niter);

  for(unsigned i = 0; i < m_options.niter; i++)
    {
      std::cout << "ITERATION " << i << std::endl;

      VSAnalysisStage3Data d = data;

      for(unsigned irun = 0; irun < nrun; irun++)
	{
	  const double noff = d.run_data(irun).egy_off_np_hist.sum();	    
	  d.run_data(irun).egy_off_np_hist.clear();
	  d.run_data(irun).egy_on_np_hist.clear();

	  d.run_data(irun).egy_off_np_hist.fill(0);
	  d.run_data(irun).egy_on_np_hist.fill(0);

 	  //hrng.fill(m_rng,d.run_data(irun).egy_off_np_hist,(double)noff);
 	  //hrng.fill(m_rng,d.run_data(irun).egy_on_np_hist,(double)0.1*noff);
	  src_injector[irun]->fill(&d.run_data(irun));
	}     

      

//       for(unsigned irun = 0; irun < nrun; irun++)
//     {
//       std::cout << "Filling run " << irun << std::endl;
     
//     }
      VSSpectrumData sp_data;

      m_spectrum_calc->reconstruct(d,sp_data);

      std::cout << "CHI2 " << sp_data.chi2 << std::endl;

      VSLimitedErrorsHist<double,double> bias_hist = mchist;
      VSLimitedErrorsHist<double,double> rbias_hist = mchist;


      bias_hist.clear();
      bias_hist.fill(0);

      rbias_hist.clear();
      rbias_hist.fill(0);

      for(VSLimitedErrorsHist<double,double>::iterator itr =
	    mchist.begin(); itr != mchist.end(); ++itr)
	{
	  double rflux = sp_data.flux_hist.countForVal(itr->center());
	  double rflux_var = sp_data.flux_hist.varForVal(itr->center());
	  double flux = itr->count();

	  if(rflux_var > 0)
	    {
	      if(fabs(rflux-flux) < sqrt(rflux_var))
		m_ncov_hist.accumulate(itr->center(),1,0);	      
	      m_ntot_hist.accumulate(itr->center(),1,0);
	    }

	  //if(rflux_var == 0) continue;

	  bias_hist.setBin(itr->bin(),rflux-flux,rflux_var);

	  if(rflux_var > 0 && fabs((rflux-flux)/flux) < 5)
	    rbias_hist.setBin(itr->bin(),(rflux-flux)/flux,
			      rflux_var/(flux*flux));
	  
	  m_bias_hist[itr->bin()].accumulate((rflux-flux)/flux);
	  m_bias_stat[itr->bin()].accumulate((rflux-flux)/flux);

	  std::cout << std::setw(10) << itr->center() 
		    << std::setw(20) << flux 
		    << std::setw(20) << rflux 
		    << std::setw(20) << rflux-flux 
		    << std::setw(20) << (rflux-flux)/flux 
		    << std::endl;
	}

      m_data[i].sp_data = sp_data;
      m_data[i].bias_hist = bias_hist;
      m_data[i].rbias_hist = rbias_hist;
    }


  for(VSLimitedErrorsHist<double,double>::iterator itr =
	m_mean_bias_hist.begin(); itr != m_mean_bias_hist.end(); ++itr)
    {
      double ntot = m_ntot_hist.countForVal(itr->center());
      double ncov = m_ncov_hist.countForVal(itr->center());

      if(ntot > 0)
	{
	  double cov = ncov/ntot;
	  double cov_err = sqrt(cov*(1-cov)/ntot);
	  m_coverage_hist.accumulate(itr->center(),cov,std::pow(cov_err,2));
	}

      if(m_bias_stat[itr->bin()].count() > 0)
	{
	  std::cout << std::setw(10) << itr->center() 
		    << std::setw(15) << m_bias_stat[itr->bin()].mean()
		    << std::setw(15) << m_bias_stat[itr->bin()].dev()
		    << std::setw(15) << m_bias_stat[itr->bin()].count()
		    << std::endl;

	  double mean_bias = m_bias_stat[itr->bin()].mean();
	  double mean_bias_err = 
	    m_bias_stat[itr->bin()].dev()/
	    sqrt(m_bias_stat[itr->bin()].count());

	  m_mean_bias_hist.setBin(itr->bin(),
				  mean_bias,std::pow(mean_bias_err,2));
	  m_mean_bias_dev_hist.
	    setBin(itr->bin(),
		   mean_bias,std::pow(m_bias_stat[itr->bin()].dev(),2));

	  m_mean_dev_hist.setBin(itr->bin(),
				 m_bias_stat[itr->bin()].dev(),
				 std::pow(m_bias_stat[itr->bin()].dev()/
					  sqrt(2*m_bias_stat[itr->bin()].
					       count()),2));
	}
    }

}

void VSSpectrumBiasCalc::save(VSOctaveH5Writer* writer) const
{
  const unsigned ndata = m_data.size();
  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("data",ndata);
  vsassert(wc);
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(idata);
      m_data[idata].save(ws);
      delete ws;
    }
  delete wc;


  const unsigned nhist = m_bias_hist.size();
  wc = writer->writeCellVector("bias_hist",nhist);
  vsassert(wc);
  for(unsigned ihist = 0; ihist < nhist; ihist++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(ihist);
      m_bias_hist[ihist].save(ws);
      delete ws;
    }
  delete wc;
 
  m_mean_bias_hist.save(writer->writeStruct("mean_bias_hist"));
  m_mean_bias_dev_hist.save(writer->writeStruct("mean_bias_dev_hist"));
  m_mean_dev_hist.save(writer->writeStruct("mean_dev_hist"));
  m_sim_flux_hist.save(writer->writeStruct("sim_flux_hist"));
  m_sim_dfde_hist.save(writer->writeStruct("sim_dfde_hist"));
  m_sim_edfde_hist.save(writer->writeStruct("sim_edfde_hist"));
  m_sim_e2dfde_hist.save(writer->writeStruct("sim_e2dfde_hist"));
  m_coverage_hist.save(writer->writeStruct("coverage_hist"));
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSSpectrumBiasCalc::configure(VSOptions& options, 
				   const std::string& profile,
				   const std::string& opt_prefix)
{
  VSSpectrumCalcFactory::configure(options,profile,opt_prefix);
  VSSourceInjector::configure(options,profile,opt_prefix);

  options.findWithValue(OPTNAME(opt_prefix,"niter"), 
			s_default_options.niter,
			"Select the cut optimization method.  Available "
			"methods are: simple/nspace.");  
}
