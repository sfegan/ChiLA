//-*-mode:c++; mode:font-lock;-*-

/*! \file VSNSpaceCutOptimizer.cpp
  

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/17/2005
*/

#include "VSNSpaceCutOptimizer.hpp"
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSNSpaceOptimizer.hpp>

using namespace VERITAS;
using namespace SEphem;

void VSNSpaceCutOptimizer::OptData::save(VSOctaveH5WriterStruct* writer) const
{
  VSNSpaceOctaveH5IO io;
  io.writeVolume(writer->writeStruct("filter"),filter);
  writer->writeCompositeHere(*this);
}

VSNSpaceCutOptimizer::VSNSpaceCutOptimizer(const Options& opt):
  VSCutOptimizer(opt),
  m_data()
{

}
 
VSNSpaceCutOptimizer::~VSNSpaceCutOptimizer()
{

}

void VSNSpaceCutOptimizer::optimize(const VSCutOptimizer::Data& data)
{
  unsigned nthsq = 1;
  if(m_options.use_theta) nthsq = m_space.axes[0].nbin;
  std::vector<Ordering> ordering(nthsq);
  std::cerr << "VSNSpaceCutOptimizer::optimize(): Calculating ordering" 
	    << std::endl;

  VSNSpaceOptimizer* nspace_optimizer = new VSNSpaceOptimizer;

  for(unsigned ibin = 0; ibin < nthsq; ibin++)
    {
      VSNSpace excess      = data.excess;
      VSNSpace excess_var  = data.excess_var;
      VSNSpace excess_rate = data.excess_rate;
      VSNSpace variance    = data.bkgnd_var;
      VSNSpace bkgnd       = data.bkgnd;
      VSNSpace bkgnd_rate  = data.bkgnd_rate;

      if(m_options.use_theta)
	{
	  VSNSpace::Volume marginal_volume(m_space);

	  marginal_volume.setInvert();
	  marginal_volume.clearInsideRange(0,ibin+1,nthsq);
	  
	  excess.marginalize(0,marginal_volume);
	  excess_var.marginalize(0,marginal_volume);
	  variance.marginalize(0,marginal_volume);
	  bkgnd.marginalize(0,marginal_volume);
	  excess_rate.marginalize(0,marginal_volume);
	  bkgnd_rate.marginalize(0,marginal_volume);
	  
	  double thetasq = m_space.axes[0].maxCoordUnchecked(ibin);
	  bkgnd_rate *= thetasq;
	}

      nspace_optimizer->optimize(m_opt_func,excess_rate,bkgnd_rate,0.2);

      ordering[ibin].index = nspace_optimizer->getOrdering();
      std::vector<double> Q = nspace_optimizer->getOrderingQ();
      ordering[ibin].data.resize(ordering[ibin].index.size());

      VSNSpace::Weight sum_excess = 0;
      VSNSpace::Weight sum_excess_var = 0;
      VSNSpace::Weight sum_excess_rate = 0;
      VSNSpace::Weight sum_bkgnd = 0;
      VSNSpace::Weight sum_bkgnd_rate = 0;
      VSNSpace::Weight sum_variance = 0;
      
      const unsigned nindex = ordering[ibin].index.size();
      for(unsigned iindex = 0; iindex<nindex; iindex++)
	{
	  unsigned index = ordering[ibin].index[iindex];
	  sum_excess      += excess.getWeightUnchecked(index);
	  sum_excess_var  += excess_var.getWeightUnchecked(index);
	  sum_excess_rate += excess_rate.getWeightUnchecked(index);
	  sum_bkgnd       += bkgnd.getWeightUnchecked(index);
	  sum_bkgnd_rate  += bkgnd_rate.getWeightUnchecked(index);
	  sum_variance    += variance.getWeightUnchecked(index);
	  
	  Ordering::Data& data = ordering[ibin].data[iindex];

	  data.excess        = sum_excess;
	  data.excess_err    = sqrt(sum_excess_var);
	  data.excess_rate   = sum_excess_rate;
	  data.bkgnd         = sum_bkgnd;
	  data.bkgnd_rate    = sum_bkgnd_rate;
	  data.index         = iindex;
	  data.Q             = Q[iindex];
	}
    }

  delete nspace_optimizer;

  std::cout << std::setw(12) << "RATE"
    	    << std::setw(12) << "Q"
	    << std::setw(12) << "SIGNAL"
	    << std::setw(12) << "BKGND"
	    << std::setw(12) << "TAU [MIN]"
	    << std::setw(12) << "R_SIGNAL"
	    << std::setw(12) << "R_BKGND"
	    << std::setw(18) << "THETA"
	    << std::endl;

  VSSimpleGraph<double,double> rate_Q_graph;
  VSSimpleGraph<double,double> rate_sn_graph;
  VSSimpleGraph<double,double> fraction_Q_graph;
  VSSimpleGraph<double,double> fraction_sn_graph;
 
  VSNSpace::Weight max_rate = data.excess_rate.totalWeight();
  for(unsigned ifraction = 0;ifraction < m_options.min_fraction.size();
      ifraction++)
    m_options.min_rate.push_back(m_options.min_fraction[ifraction]*max_rate);

  if(m_options.min_rate.empty()) m_options.min_rate.push_back(0);

  m_data.resize(m_options.min_rate.size());

  for(unsigned irate=0;irate<m_options.min_rate.size();irate++)
    {
      double mr = m_options.min_rate[irate];
      if(mr<0)continue;

      double thetasq = 0;
      unsigned thsq_index = 0;

      Ordering::Data optimal_ordering;

      const unsigned nthsq = ordering.size();
      for(unsigned ithsq = 0; ithsq < nthsq; ithsq++)
	{
	  const unsigned nindex = ordering[ithsq].index.size();
	  for(unsigned iindex = 0; iindex < nindex; iindex++)
	    {	  
	      if(ordering[ithsq].data[iindex].excess_rate > mr &&
		 ordering[ithsq].data[iindex].Q > optimal_ordering.Q &&
		 std::isfinite(ordering[ithsq].data[iindex].Q))
		{
		  optimal_ordering = ordering[ithsq].data[iindex];
		  if(m_options.use_theta) 
		    thetasq = m_space.axes[0].maxCoordUnchecked(ithsq);
		  thsq_index = ithsq;
		}
	    }
	}

      VSNSpace::Space marginal_space = m_space;
      if(m_options.use_theta) marginal_space.marginalize(0);
      std::vector<unsigned> index = ordering[thsq_index].index;      
      index.erase(index.begin()+optimal_ordering.index,index.end());
      m_data[irate].filter.loadSparse(marginal_space,index);

      double Q = optimal_ordering.Q;
      double Q_err = 0;
      double sn = optimal_ordering.excess/sqrt(optimal_ordering.bkgnd);
      double sn_err = 0;
      double tau = 
	optimal_ordering.bkgnd_rate/std::pow(optimal_ordering.excess_rate,2);
      double tau_err = 0;

      std::cout.setf(std::ios::scientific);

      m_rate_Q_graph.addVertex(optimal_ordering.excess_rate,0,Q,Q_err);
      m_rate_sn_graph.addVertex(optimal_ordering.excess_rate,0,sn,sn_err);
      m_fraction_Q_graph.addVertex(optimal_ordering.excess_rate/max_rate, 
				 0,Q,Q_err);
      m_fraction_sn_graph.addVertex(optimal_ordering.excess_rate/max_rate,
				  0, sn,sn_err);

      std::cout 
	<< std::setprecision(4)
	<< std::setw(12) << mr 
	<< std::setw(12) << optimal_ordering.Q
	<< std::setw(12) << optimal_ordering.excess
	<< std::setw(12) << optimal_ordering.bkgnd
	<< std::setw(12) << tau
	<< std::setw(12) << optimal_ordering.excess_rate
	<< std::setw(12) << optimal_ordering.bkgnd_rate
	<< std::setw(18) << sqrt(thetasq)
	<< std::endl;

      m_data[irate].theta_cut = sqrt(thetasq);
      m_data[irate].min_rate = mr;
      m_data[irate].Q = Q;
      m_data[irate].Q_err = Q_err;
      m_data[irate].sn = sn;
      m_data[irate].sn_err = sn_err;
      m_data[irate].tau = tau;
      m_data[irate].tau_err = tau_err;
      m_data[irate].excess_rate = optimal_ordering.excess_rate;
    }
}

void VSNSpaceCutOptimizer::save(VSOctaveH5WriterStruct* writer) const
{
  VSCutOptimizer::save(writer);

  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("cuts",m_data.size());
  vsassert(wc);
  const unsigned ndata = m_data.size();
  for(unsigned idata = 0; idata < ndata; idata++)
    {
      VSOctaveH5WriterStruct* s = wc->writeStruct(idata);
      m_data[idata].save(s);

      bool is_array_filter = false;
      bool is_scope_filter_default = false;

      const VSNSpace::Volume* array_filter = NULL;
      const VSNSpace::Volume* scope_filter_default = NULL;
      std::vector<const VSNSpace::Volume*> scope_filters;
      
      for(unsigned idim = 0; 
	  idim < m_data[idata].filter.space().ndim; idim++)
	{
	  std::string index = 
	    VSH5DatumElementParser::
	    getIndex(m_data[idata].filter.space().axes[idim].name);
	  
	  if(index == "*")
	    is_scope_filter_default = true;
	  else if(index.empty())
	    is_array_filter = true;      
	}
      
      if(is_scope_filter_default)
	scope_filter_default = &m_data[idata].filter;
      else if(is_array_filter)
	array_filter = &m_data[idata].filter;
      
      VSNSpaceCutsCalc* nspace_cuts =
	new VSNSpaceCutsCalc(array_filter, scope_filter_default, 
			     scope_filters);

      nspace_cuts->save(s->writeStruct("nspace_cuts"));

      delete nspace_cuts;
      delete s;
    }
  delete wc;
      
  m_rate_Q_graph.save(writer->writeStruct("rate_Q_graph"));
  m_rate_sn_graph.save(writer->writeStruct("rate_sn_graph"));
  m_fraction_Q_graph.save(writer->writeStruct("fraction_Q_graph"));
  m_fraction_sn_graph.save(writer->writeStruct("fraction_sn_graph"));
}
