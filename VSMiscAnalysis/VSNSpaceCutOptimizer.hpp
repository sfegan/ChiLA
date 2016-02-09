//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleCutOptimizer.hpp

  Data structures for array and scope simulation data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSNSpaceCutOptimizer.hpp,v 1.4 2009/12/17 03:39:11 matthew Exp $

*/

#ifndef VSNSPACECUTOPTIMIZER_HPP
#define VSNSPACECUTOPTIMIZER_HPP

#include <VSCutOptimizer.hpp>

namespace VERITAS
{
  class VSNSpaceCutOptimizer : public VSCutOptimizer
  {
  public:

    struct Ordering
    {
    public:

      struct Data
      {
      public:
	Data():
	  index(), 
	  excess(), excess_err(), 
	  excess_rate(), excess_rate_err(),
	  bkgnd(), bkgnd_rate(), variance(), Q()
	{ }
	
	unsigned index;
	double   excess;
	double   excess_err;
	double   excess_rate;
	double   excess_rate_err;
	double   bkgnd;
	double   bkgnd_rate;
	double   variance;
	double   Q;      
      };

      Ordering():
	index(), data() {}

      std::vector<unsigned>       index;
      std::vector<Data>           data;
    };

    struct OptData
    {
      OptData():
	min_rate(),
	Q(),Q_err(),sn(),sn_err(),tau(),tau_err(),
	excess(), excess_err(), excess_rate(), excess_rate_err(),
	bkgnd(), bkgnd_err(), bkgnd_rate(), bkgnd_rate_err(), theta_cut(),
	filter()
      {}

      double              min_rate;
      double              Q;
      double              Q_err;
      double              sn;
      double              sn_err;
      double              tau;
      double              tau_err;
      double              excess;
      double              excess_err;
      double              excess_rate;
      double              excess_rate_err;
      double              bkgnd;
      double              bkgnd_err;
      double              bkgnd_rate;
      double              bkgnd_rate_err;
      double              theta_cut;
      VSNSpace::Volume    filter;

      void save(VSOctaveH5WriterStruct* writer) const;

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,OptData,min_rate);
	H5_ADDMEMBER(c,OptData,Q);
	H5_ADDMEMBER(c,OptData,Q_err);
	H5_ADDMEMBER(c,OptData,sn);
	H5_ADDMEMBER(c,OptData,sn_err);
	H5_ADDMEMBER(c,OptData,tau);
	H5_ADDMEMBER(c,OptData,tau_err);
	H5_ADDMEMBER(c,OptData,excess_rate);
	H5_ADDMEMBER(c,OptData,excess_rate_err);
	H5_ADDMEMBER(c,OptData,bkgnd_rate);
	H5_ADDMEMBER(c,OptData,bkgnd_rate_err);
	H5_ADDMEMBER(c,OptData,theta_cut);
      }
    };
    
    VSNSpaceCutOptimizer(const VSCutOptimizer::Options& opt);
    virtual ~VSNSpaceCutOptimizer();

    virtual void optimize(const VSCutOptimizer::Data& data);
    virtual void test(const VSCutOptimizer::Data& data,
		      std::vector<CutResults>& o) { }

    virtual void save(VSOctaveH5WriterStruct* writer) const;

  private:

    std::vector<OptData> m_data;

    VSSimpleGraph<double,double>             m_rate_Q_graph;
    VSSimpleGraph<double,double>             m_rate_sn_graph;
    VSSimpleGraph<double,double>             m_fraction_Q_graph;
    VSSimpleGraph<double,double>             m_fraction_sn_graph;
  };
}

#endif // VSNSPACECUTOPTIMIZER_HPP
