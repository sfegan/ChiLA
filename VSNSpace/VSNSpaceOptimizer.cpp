#include "VSNSpaceOptimizer.hpp"

using namespace VERITAS;

// ============================================================================
// VSNSpaceOptimizer
// ============================================================================
VSNSpaceOptimizer::Options VSNSpaceOptimizer::s_default_options = 
VSNSpaceOptimizer::Options();

VSNSpaceOptimizer::VSNSpaceOptimizer(VSNSpace::Weight min_variance,
				     const Options& opt):
  m_options(opt),
  m_min_variance(min_variance),
  m_ordering()
{

}

VSNSpaceOptimizer::~VSNSpaceOptimizer()
{

}

bool VSNSpaceOptimizer::optimize(VSAMath::Functor< std::vector<double> >* func,
				 const VSNSpace& s, const VSNSpace& b,
				 double alpha)
{
  if(m_options.ordering_method == "cumulative")
    return cumulativeOrdering(func,s,b,alpha);
  else
    return false;
}

bool VSNSpaceOptimizer::
cumulativeOrdering(VSAMath::Functor< std::vector<double> >* func,
		   const VSNSpace& s,const VSNSpace& b, double alpha)
{
  m_ordering.clear();
  m_orderingQ.clear();
  m_ordering.reserve(s.size());

  VSNSpace::Volume vol(s.space());
  for(VSNSpace::Index iindex = 0; iindex < s.size(); iindex++)
    if(s[iindex]) vol.setIndexUnchecked(iindex,true);
  
  VSNSpace::Weight s_sum = 0;
  VSNSpace::Weight b_sum = 0;

  for(VSNSpace::Index iindex = 0; iindex < s.size(); iindex++)
    {
      if(b[iindex]>=m_min_variance && s[iindex] && 
	 !vol.isIndexIsolated(iindex))
	{
	  vol.setIndexUnchecked(iindex,true);
	  m_ordering.push_back(iindex);
	  s_sum += s[iindex];
	  b_sum += b[iindex];
	}
      else vol.setIndexUnchecked(iindex,false);
    }

  VSNSpace::Index iindex = 0;
  std::vector<double> x(3);
  std::vector<VSNSpace::Index>::iterator zindex = m_ordering.begin();
  std::vector<VSNSpace::Index>::iterator nindex = m_ordering.end();
  while(zindex != nindex-1)
    {
      VSNSpace::Weight Q_max = 0;
      std::vector<VSNSpace::Index>::iterator jindex_max = m_ordering.end();

      for(std::vector<VSNSpace::Index>::iterator jindex = zindex; 
	  jindex != nindex; ++jindex)
	{
	  if(!vol.isIndexOnEdge(*jindex)) continue;

	  x[0] = s_sum-s[*jindex];
	  x[1] = b_sum-b[*jindex];
	  x[2] = alpha;

	  VSNSpace::Weight Q = (*func)(x);
	  if(Q>Q_max)Q_max = Q, jindex_max = jindex;
	}

      vsassert(jindex_max != m_ordering.end());

      m_orderingQ.push_back(Q_max);

      vol.setIndexUnchecked(*jindex_max,false);
      s_sum -= s[*jindex_max];
      b_sum -= b[*jindex_max];

      iindex = *zindex;
      *zindex = *jindex_max;
      *jindex_max = iindex;
      zindex++;
    }

  std::reverse(m_ordering.begin(),m_ordering.end());
  std::reverse(m_orderingQ.begin(),m_orderingQ.end());

  return true;
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSNSpaceOptimizer::configure(VSOptions& options,
				  const std::string& profile, 
				  const std::string& opt_prefix)
{
  options.findWithValue(OPTNAME(opt_prefix,"nspace_ordering_method"), 
			s_default_options.ordering_method,
			"Select the ordering method for generating the "
			"NSpace.");  
}

