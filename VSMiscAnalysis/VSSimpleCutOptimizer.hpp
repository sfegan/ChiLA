//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleCutOptimizer.hpp

  Data structures for array and scope simulation data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSSimpleCutOptimizer.hpp,v 1.5 2009/12/17 03:39:11 matthew Exp $

*/

#ifndef VSSIMPLECUTOPTIMIZER_HPP
#define VSSIMPLECUTOPTIMIZER_HPP

#include <VSCutOptimizer.hpp>

namespace VERITAS
{
  class VSSimpleCutOptimizer : public VSCutOptimizer
  {
  public:

    struct OptCut
    {
      OptCut(std::string _name, double _cut_value):
	name(_name), cut_value(_cut_value)
      {}

      std::string      name;
      double           cut_value;
    };

    struct OptData
    {
      OptData():
	min_rate(),
	Q(),Q_err(),sn(),sn_err(),tau(),tau_err(),
	excess_rate(), excess_rate_err(),
	bkgnd_rate(), bkgnd_rate_err(), theta_cut(),
	cuts()
      {}

      double              min_rate;
      double              Q;
      double              Q_err;
      double              sn;
      double              sn_err;
      double              tau;
      double              tau_err;
      double              excess_rate;
      double              excess_rate_err;
      double              bkgnd_rate;
      double              bkgnd_rate_err;
      double              theta_cut;
      std::vector<OptCut> cuts;

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
    
    VSSimpleCutOptimizer(const VSCutOptimizer::Options& opt);
    virtual ~VSSimpleCutOptimizer();

    virtual void optimize(const VSCutOptimizer::Data& data);
    virtual void test(const VSCutOptimizer::Data& data,
		      std::vector<CutResults>& o);

    bool testPoint(const std::vector<OptCut>& cuts, 
		   const VSNSpace::Point& p,
		   bool test_th2 = true);

    virtual void save(VSOctaveH5WriterStruct* writer) const;

  private:

    std::vector<OptData>                     m_data;

    VSNSpace                                 m_excess_cumulative;
    VSNSpace                                 m_excess_var_cumulative;
    VSNSpace                                 m_excess_rate_cumulative;
    VSNSpace                                 m_bkgnd_cumulative;
    VSNSpace                                 m_bkgnd_var_cumulative;
    VSNSpace                                 m_bkgnd_rate_cumulative;
    VSNSpace                                 m_bkgnd_counts_cumulative;
    VSNSpace                                 m_Q_cumulative;

    VSSimpleGraph<double,double>             m_rate_Q_graph;
    VSSimpleGraph<double,double>             m_rate_sn_graph;
    VSSimpleGraph<double,double>             m_fraction_Q_graph;
    VSSimpleGraph<double,double>             m_fraction_sn_graph;
  };
}

#endif // VSSIMPLECUTOPTIMIZER_HPP
