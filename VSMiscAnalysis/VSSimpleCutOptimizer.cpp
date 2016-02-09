//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleCutOptimizer.cpp
  

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/17/2005
*/

#include "VSSimpleCutOptimizer.hpp"
#include <VSNSpaceOctaveH5IO.hpp>

using namespace VERITAS;
using namespace SEphem;

void VSSimpleCutOptimizer::OptData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);

  const unsigned ncut = cuts.size();
  for(unsigned icut = 0; icut < ncut; icut++)
    writer->writeScalar(cuts[icut].name,cuts[icut].cut_value);
}

VSSimpleCutOptimizer::VSSimpleCutOptimizer(const Options& opt):
  VSCutOptimizer(opt),
  m_data(),
  m_excess_cumulative(),
  m_excess_var_cumulative(),
  m_excess_rate_cumulative(),
  m_bkgnd_cumulative(),
  m_bkgnd_var_cumulative(),
  m_bkgnd_rate_cumulative(),
  m_bkgnd_counts_cumulative(),
  m_Q_cumulative(),
  m_rate_Q_graph(),
  m_rate_sn_graph(),
  m_fraction_Q_graph(),
  m_fraction_sn_graph()
{

}
 
VSSimpleCutOptimizer::~VSSimpleCutOptimizer()
{

}

void VSSimpleCutOptimizer::optimize(const VSCutOptimizer::Data& data)
{
  m_excess_cumulative       = data.excess;
  m_excess_var_cumulative   = data.excess_var;
  m_excess_rate_cumulative  = data.excess_rate;
  m_bkgnd_counts_cumulative = data.bkgnd_counts;
  m_bkgnd_cumulative        = data.bkgnd;
  m_bkgnd_var_cumulative    = data.bkgnd_var;
  m_bkgnd_rate_cumulative   = data.bkgnd_rate;

  std::cerr << "VSCutOptimizer::optimize(): Integrating counts" << std::endl;
  m_excess_cumulative.partiallyIntegrate();
  m_excess_var_cumulative.partiallyIntegrate();
  m_excess_rate_cumulative.partiallyIntegrate();
  m_bkgnd_cumulative.partiallyIntegrate();
  m_bkgnd_var_cumulative.partiallyIntegrate();
  m_bkgnd_rate_cumulative.partiallyIntegrate();
  m_bkgnd_counts_cumulative.partiallyIntegrate();

  std::cerr << "VSCutOptimizer::optimize(): Calculating signal and background" 
	    << std::endl;
  VSNSpace::Volume vol_good;
  if(m_options.min_events) 
    vol_good = m_bkgnd_counts_cumulative>=double(m_options.min_events);
 
  std::cerr << "VSCutOptimizer::optimize(): Calculating Q " 
	    << std::endl;
  m_Q_cumulative = VSNSpace(m_space);
  for(VSNSpace::Index i = 0; i < m_space.size(); i++)
    {
      VSNSpace::Point p;
      m_space.maxPointOfIndexUnchecked(i,p);

      std::vector<double> x(3);

      x[0] = m_excess_rate_cumulative[i];
      x[1] = m_bkgnd_rate_cumulative[i];
      x[2] = 0.2;

      m_Q_cumulative[i] = (*m_opt_func)(x);
    }

  if(m_options.min_events) m_Q_cumulative.clearOutsideVolume(vol_good);

  VSNSpace::Weight max_rate = m_excess_rate_cumulative.maxWeight();

  std::cout << "MAX RATE: " << max_rate << std::endl;

  for(unsigned ifraction = 0;ifraction < m_options.min_fraction.size();
      ifraction++)
    m_options.min_rate.push_back(m_options.min_fraction[ifraction]*max_rate);

  m_data.resize(m_options.min_rate.size());

  std::cout << std::string(79,'=') << std::endl
	    << "= CUTS" << std::endl
	    << std::string(79,'=') << std::endl;

  std::cout << std::setw(12) << "RATE";
//     	    << std::setw(12) << "Q"
// 	    << std::setw(12) << "SIGNAL"
// 	    << std::setw(12) << "BKGND"
// 	    << std::setw(12) << "TAU [MIN]"
// 	    << std::setw(12) << "R_SIGNAL"
// 	    << std::setw(12) << "R_BKGND";

  for(unsigned idim = 0; idim <  m_space.ndim; idim++)
    std::cout << std::setw(18) << m_space.axes[idim].name;
  std::cout << std::endl;

  for(unsigned irate=0;irate<m_options.min_rate.size();irate++)
    {
      double mr = m_options.min_rate[irate];
      if(mr<0)continue;
      
      VSNSpace::Index max_index = 0;
      VSNSpace::Weight max_Q = 0;

      m_Q_cumulative.maxWeight(m_excess_rate_cumulative>=mr,max_Q,max_index);
  
      double excess = m_excess_cumulative[max_index];
      double excess_var = m_excess_var_cumulative[max_index];
      double excess_rate = m_excess_rate_cumulative[max_index];
      double bkgnd = m_bkgnd_cumulative[max_index];
      double bkgnd_var = m_bkgnd_var_cumulative[max_index];
      double bkgnd_rate = m_bkgnd_rate_cumulative[max_index];

      double Q = max_Q;
      double Q_err = 0;
      double sn = excess/sqrt(bkgnd);
      double sn_err = 0;
      double tau = bkgnd_rate/std::pow(excess_rate,2);
      double tau_err = 0;

      VSNSpace::Point cut_val = m_space.point();
      m_space.maxPointOfIndexUnchecked(max_index, cut_val);
  
      std::cout.setf(std::ios::scientific);

      m_rate_Q_graph.addVertex(excess_rate,0,Q,Q_err);
      m_rate_sn_graph.addVertex(excess_rate,0,sn,sn_err);
      m_fraction_Q_graph.addVertex(excess_rate/max_rate, 0, Q,Q_err);
      m_fraction_sn_graph.addVertex(excess_rate/max_rate,0, sn,sn_err);

      std::cout 
	<< std::setprecision(4)
	<< std::setw(12) << mr;
// 	<< std::setw(12) << max_Q 
// 	<< std::setw(12) << excess
// 	<< std::setw(12) << bkgnd
// 	<< std::setw(12) << tau
// 	<< std::setw(12) << excess_rate 
// 	<< std::setw(12) << bkgnd_rate;

      for(unsigned idim = 0; idim < cut_val.ndim; idim++)
	std::cout << std::setw(18) << cut_val.x[idim];
      std::cout << std::endl;

      m_data[irate].min_rate = mr;
      m_data[irate].Q = Q;
      m_data[irate].Q_err = Q_err;
      m_data[irate].sn = sn;
      m_data[irate].sn_err = sn_err;
      m_data[irate].tau = tau;
      m_data[irate].tau_err = tau_err;
      m_data[irate].excess_rate = excess_rate;
      m_data[irate].excess_rate_err = excess_rate*sqrt(excess_var)/excess;
      m_data[irate].bkgnd_rate = bkgnd_rate;
      m_data[irate].bkgnd_rate_err = bkgnd_rate*sqrt(bkgnd_var)/bkgnd;

      if(m_options.use_theta)
	m_data[irate].theta_cut = sqrt(cut_val.x[0]);

      //  m_data[irate].cuts.push_back(OptCut("theta",sqrt(cut_val.x[0])));
      for(unsigned idim = 0; idim < cut_val.ndim; idim++)
	m_data[irate].cuts.push_back(OptCut(m_space.axes[idim].name,
					    cut_val.x[idim]));
    }

}

void VSSimpleCutOptimizer::test(const VSCutOptimizer::Data& data,
				std::vector<CutResults>& o)
{
  std::cout.setf(std::ios::left);
  std::cout << std::setw(5) << ""
	    << std::setw(25) << "TAU [MIN]"
	    << std::setw(25) << "R_SIGNAL"
	    << std::setw(25) << "R_BKGND"
	    << std::setw(12) << "T (1%)"
	    << std::setw(12) << "T (3%)"
	    << std::setw(12) << "E THRESHOLD"
	    << std::endl;

  const unsigned ndata = m_data.size();

  o.resize(ndata);

  std::vector<double> bkgnd(ndata);
  std::vector<double> bkgnd_var(ndata);
  std::vector<double> signal(ndata);
  std::vector<double> signal_var(ndata);

  for(unsigned idata = 0; idata < ndata; idata++)
    {      
      double th2_cut = m_data[idata].cuts[0].cut_value;

      for(std::vector<VSCutOptimizer::Event>::const_iterator itr = 
	    data.bkgnd_events.begin(); itr != data.bkgnd_events.end(); ++itr)
	{
	  if(!testPoint(m_data[idata].cuts,itr->p,false)) continue;

	  double wt = itr->wt;
	  if(m_options.use_theta) 
	    wt *= th2_cut/std::pow(m_options.theta_cut,2);

	  bkgnd[idata] += wt;
	  bkgnd_var[idata] += wt*wt;
	}

      std::vector<double> ecount(m_sim_count.size());

      for(std::vector<VSCutOptimizer::Event>::const_iterator itr = 
	    data.signal_events.begin(); itr != data.signal_events.end(); ++itr)
	if(testPoint(m_data[idata].cuts,itr->p))
	  {
	    signal[idata] += itr->wt;
	    signal_var[idata] += std::pow(itr->wt,2);

	    ecount[itr->index] += itr->wt;
	  }

      unsigned ipeak = 0;
      double npeak = 0;
      double ethreshold = 0;

      for(unsigned iegy = 0; iegy < ecount.size(); iegy++)
	{
	  o[idata].edrde.addVertex(std::log10(m_sim_energy[iegy]),
			    ecount[iegy]/data.elaptime_min);

	  if(ecount[iegy] > npeak)
	    {
	      ipeak = iegy;
	      npeak = ecount[iegy];
	      ethreshold = m_sim_energy[iegy];
	    }
	}

      for(std::vector<VSCutOptimizer::Event>::const_iterator itr = 
	    data.on_events.begin(); itr != data.on_events.end(); ++itr)
	if(testPoint(m_data[idata].cuts,itr->p))
	  {
	    signal[idata] += itr->wt;
	    signal_var[idata] += std::pow(itr->wt,2);
	  }

      for(std::vector<VSCutOptimizer::Event>::const_iterator itr = 
	    data.off_events.begin(); itr != data.off_events.end(); ++itr)
	if(testPoint(m_data[idata].cuts,itr->p))
	  {
	    signal[idata] -= itr->wt;
	    signal_var[idata] += std::pow(itr->wt,2);
	  }

      double Q = signal[idata]/sqrt(bkgnd[idata]);
      double signal_rate = signal[idata]/data.elaptime_min;
      double signal_rate_err = sqrt(signal_var[idata])/data.elaptime_min;
      double bkgnd_rate = bkgnd[idata]/data.elaptime_min;
      double bkgnd_rate_err = sqrt(bkgnd_var[idata])/data.elaptime_min;
      double tau = bkgnd_rate/std::pow(signal_rate,2);
      double tau_err = 
	tau*sqrt(std::pow(bkgnd_rate_err/bkgnd_rate,2) +
		 2*std::pow(signal_rate_err/signal_rate,2));

      double crab1p_time = 25*1.1*tau/std::pow(0.01,2)/60.;
      double crab3p_time = 25*1.1*tau/std::pow(0.03,2)/60.;

      o[idata].signal_rate = signal_rate;
      o[idata].signal_rate_err = signal_rate_err;
      o[idata].bkgnd_rate = bkgnd_rate;
      o[idata].bkgnd_rate_err = bkgnd_rate_err;
      o[idata].tau = tau;
      o[idata].tau_err = tau_err;

      std::cout << std::setw(5)  << idata 
		<< std::setprecision(3)
		<< std::setw(10) << tau
		<< " +/- " << std::setw(10) << tau_err
		<< std::setw(10) << signal_rate
		<< " +/- " << std::setw(10) << signal_rate_err
		<< std::setw(10) << bkgnd_rate
		<< " +/- " << std::setw(10) << bkgnd_rate_err
		<< std::setw(12) << crab1p_time
		<< std::setw(12) << crab3p_time
		<< std::setw(12) << ethreshold
		<< std::endl;
    }


}

bool VSSimpleCutOptimizer::testPoint(const std::vector<OptCut>& cuts, 
				     const VSNSpace::Point& p, 
				     bool test_th2)
{
  bool passed = true;
  
  const unsigned ncut = cuts.size();
  for(unsigned icut = 0; icut < ncut; icut++)
    {
      if(m_options.use_theta && icut == 0 && !test_th2) continue;

      if(p.x[icut] > cuts[icut].cut_value)
	{
	  passed = false;
	  break;
	}
    }
  
  return passed;
}

void VSSimpleCutOptimizer::save(VSOctaveH5WriterStruct* writer) const
{
  VSCutOptimizer::save(writer);

  VSNSpaceOctaveH5IO* io = new VSNSpaceOctaveH5IO;
      
  io->writeHistogram(writer->writeStruct("s"), m_excess_cumulative);
  io->writeHistogram(writer->writeStruct("b"), m_bkgnd_cumulative);
  io->writeHistogram(writer->writeStruct("Q"), m_Q_cumulative);
      
  delete io;

  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("cuts",m_data.size());
  vsassert(wc);
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      VSOctaveH5WriterStruct* s = wc->writeStruct(idata);
      m_data[idata].save(s);

      VSSimpleCutsCalc* simple_cuts = new VSSimpleCutsCalc;

      const unsigned ncut = m_data[idata].cuts.size();
      for(unsigned icut = 1; icut < ncut; icut++)
	{
	  std::string lo_cut;
	  std::string hi_cut;	      
	  VSDataConverter::toString(hi_cut,
				    m_data[idata].cuts[icut].cut_value);
	  simple_cuts->loadRangedCut(m_data[idata].cuts[icut].name,
				     lo_cut,hi_cut);
	}

      simple_cuts->save(s->writeStruct("simple_cuts"));

      delete s;
    }
  delete wc;
      
  m_rate_Q_graph.save(writer->writeStruct("rate_Q_graph"));
  m_rate_sn_graph.save(writer->writeStruct("rate_sn_graph"));
  m_fraction_Q_graph.save(writer->writeStruct("fraction_Q_graph"));
  m_fraction_sn_graph.save(writer->writeStruct("fraction_sn_graph"));
}
