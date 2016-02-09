//-*-mode:c++; mode:font-lock;-*-

/*! \file VSNSpaceOptimizer.hpp

  Utility class that handles deriving optimized NSpaces.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSNSpaceOptimizer.hpp,v 1.1 2008/12/15 21:07:04 matthew Exp $

*/

#ifndef VSNSPACEOPTIMIZER_HPP
#define VSNSPACEOPTIMIZER_HPP

#include<ostream>
#include<string>
#include<set>
#include<vector>
#include<algorithm>
#include<cmath>
#include<cassert>
#include<stdexcept>
#include<VSNSpace.hpp>
#include<VSOptions.hpp>
#include<VSAFunctor.hpp>

namespace VERITAS 
{
  class VSNSpaceOptimizer
  {
  public:
    struct Options
    {
      Options():
	ordering_method("cumulative")
      { }
      
      std::string                 ordering_method;
    };

    VSNSpaceOptimizer(VSNSpace::Weight min_variance = 0,
		      const Options& opt = s_default_options);
    ~VSNSpaceOptimizer();
    
    bool optimize(VSAMath::Functor< std::vector<double> >* func,
		  const VSNSpace& s, const VSNSpace& b, double alpha);
    bool cumulativeOrdering(VSAMath::Functor< std::vector<double> >* func,
			    const VSNSpace& s, 
			    const VSNSpace& b, double alpha);    

    const std::vector<unsigned>& getOrdering() const { return m_ordering; }
    const std::vector<double>& getOrderingQ() const { return m_orderingQ; }

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");

  private:

    Options                  m_options;

    VSNSpace::Weight         m_min_variance;

    std::vector<unsigned>    m_ordering;
    std::vector<double>      m_orderingQ;
    
    static Options           s_default_options;
  };

}

#endif // defined VSNSPACEOPTIMIZER_HPP
