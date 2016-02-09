//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutOptimizer.hpp

  Data structures for array and scope simulation data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSCutOptimizer.hpp,v 1.9 2010/06/20 00:59:19 matthew Exp $

*/

#ifndef VSCUTOPTIMIZER_HPP
#define VSCUTOPTIMIZER_HPP

#include <vector>

#include <SphericalCoords.h>
#include <VSAAlgebra.hpp>
#include <VSEventData.hpp>
#include <VSSimulationData.hpp>
#include <VSMergedCalibrationData.hpp>
#include <VSDiagnosticsData.hpp>
#include <VSAnalysisStage1.hpp>
#include <VSCutsEvaluator.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSSimEnergyWeightCalc.hpp>
#include <VSEventDataVisitor.hpp>
#include <VSH5DatumElement.hpp>
#include <VSAFunctor.hpp>

namespace VERITAS
{
  class VSCutOptimizer : public VSEventDataVisitor
  {
  public:
    class FunctionalQSimple : 
      public VSAMath::Functor< std::vector<double> >
    {
    public:
      FunctionalQSimple():
	VSAMath::Functor< std::vector<double> >() {}
      virtual ~FunctionalQSimple() {}
      
      virtual double operator() (const std::vector< double >& x) const;
//       double val(double s, double b, double alpha) const;
//       double err(double s, double b, double alpha) const;
    };

    class FunctionalLambda : 
      public VSAMath::Functor< std::vector<double> >
    {
    public:
      FunctionalLambda():
	VSAMath::Functor< std::vector<double> >() {}
      virtual ~FunctionalLambda() {}
      
      virtual double operator() (const std::vector< double >& x) const;
    };

    typedef quad<std::string, VSNSpace::Coord, 
		 VSNSpace::Coord, VSNSpace::Coord> AxisDefinition;

    struct Options
    {
      Options();
      
      std::pair<std::string,std::string> optimizer;
      std::vector<AxisDefinition> space;
      double                      theta_cut;
      bool                        use_theta;
      bool                        sim_signal;
      unsigned                    min_events;
      std::vector<unsigned>       itheta_off;
      unsigned                    size_cut_nscope;
      std::vector<double>         min_rate;
      std::vector<double>         min_fraction;
      std::vector<std::string>    bkgnd_selection;
      std::string                 spectrum;
    };

    struct Event
    {
      Event(): p(), wt(1.) {}
      Event(const VSNSpace::Point& _p, double _wt = 1.0, 
	    unsigned _index = 0): 
	p(_p), wt(_wt), index(_index) {}

      VSNSpace::Point p;
      double          wt;
      double          index;

    };

    struct Data
    {
      Data();
      Data(const VSNSpace::Space& space);

      double                            elaptime_min;

      VSNSpace                          excess;
      VSNSpace                          excess_var;
      VSNSpace                          excess_rate;
      VSNSpace                          bkgnd;
      VSNSpace                          bkgnd_var;
      VSNSpace                          bkgnd_rate;
      VSNSpace                          bkgnd_counts;

      std::vector<Event>                bkgnd_events;
      std::vector<Event>                signal_events;
      std::vector<Event>                on_events;
      std::vector<Event>                off_events;
    };

    struct CutResults
    {
      double       tau;
      double       tau_err;
      double       signal_rate;
      double       signal_rate_err;
      double       bkgnd_rate;
      double       bkgnd_rate_err;

      VSSimpleGraph<double,double> edrde;

      void save(VSOctaveH5WriterStruct* writer) const
      {
	writer->writeCompositeHere(*this);
	edrde.save(writer->writeStruct("edrde"));
      }

      static void _compose(VSOctaveH5CompositeDefinition& c)
      {
	H5_ADDMEMBER(c,CutResults,tau);
	H5_ADDMEMBER(c,CutResults,tau_err);
	H5_ADDMEMBER(c,CutResults,signal_rate);
	H5_ADDMEMBER(c,CutResults,signal_rate_err);
	H5_ADDMEMBER(c,CutResults,bkgnd_rate);
	H5_ADDMEMBER(c,CutResults,bkgnd_rate_err);
      }

    };

    VSCutOptimizer(const Options& opt);
    virtual ~VSCutOptimizer();

    virtual void visitRun(const VSAnalysisStage1Data& stage1,
			  const VSTargetTable::Observation& obs,
			  const VSArrayMergedCalibrationData& cal);

    virtual void leaveRun();

    virtual void visitEvent(const VSEventArrayDatum& event);
    virtual void leaveEvent();

    virtual void visitSimEvent(const VSArraySimulationDatum& sim);
    virtual void leaveSimEvent();

    virtual void visitSimHeader(const VSHeaderSimulationDatum& header);
    virtual void leaveSimHeader();
    
    virtual void optimize();
    virtual void optimize(const Data& d) = 0;
    virtual void test(const Data& d, std::vector<CutResults>& o) = 0;

    void accumulate(Data& data, 
		    double begin_fraction, double end_fraction);

    virtual void save(VSOctaveH5WriterStruct* writer) const;

  protected:

    void copyEvents(const std::vector<Event>& events,
		    double begin_fraction, double end_fraction,
		    std::vector<Event>& o);

    Options                           m_options;

    std::string                       m_bkgnd_selection_mode;
    std::vector<double>               m_bkgnd_selection_pars;
    double                            m_bkgnd_theta_cut;
    double                            m_bkgnd_alpha;
    double                            m_alpha;

    VSScaledParameterCalc*            m_sp_calc;
    VSSimEnergyWeightCalc*            m_egywt_calc;
    VSCutsEvaluator*                  m_cuts_calc;
    VSAMath::Functor< std::vector<double> >* m_opt_func;
    VSSpectrumFn*                     m_sim_spectrum;

    bool                              m_is_sim_run;

    VSHeaderSimulationDatum           m_sim_header;
    VSArraySimulationDatum            m_sim;
    VSEventArrayDatum                 m_evt;
    SEphem::SphericalCoords           m_obs_radec_J2000;
    SEphem::SphericalCoords           m_src_radec_J2000;     
    VSAAlgebra::Vec2D                 m_obs_xy;
    VSAAlgebra::Vec2D                 m_src_xy;
    std::vector< VSAAlgebra::Vec2D >  m_off_xy;
    double                            m_livetime;
    double                            m_elaptime;
    double                            m_wt;

    std::vector<VSH5DatumElement<VSEventArrayDatum>*> m_datum_elements;

    VSNSpace::Space                   m_space;
    
    VSNSpace                          m_sim_rate; 
    std::vector<double>               m_sim_count; 
    std::vector<double>               m_sim_energy;

    VSNSpace                          m_sim_on;

    std::vector<Event>                m_bkgnd_events;
    std::vector<Event>                m_sim_events;
    std::vector<Event>                m_on_events;
    std::vector<Event>                m_off_events;

    Data                              m_train_data;
    Data                              m_test_data;

    std::vector<CutResults>           m_train_results;
    std::vector<CutResults>           m_test_results;
  };

  class VSCutOptimizerFactory
  {
  public:
    virtual ~VSCutOptimizerFactory();

    static VSCutOptimizerFactory* getInstance();

    VSCutOptimizer* createCutOptimizer();

    static void configure(VSOptions& options, 
			  const std::string& profile="",
			  const std::string& opt_prefix="");

  private:
    VSCutOptimizerFactory(const VSCutOptimizer::Options& opt = 
			  s_default_options);

    VSCutOptimizer::Options                     m_options;

    static std::auto_ptr<VSCutOptimizerFactory> s_instance;

    static VSCutOptimizer::Options              s_default_options;
  };

}

#endif // VSCUTOPTIMIZER_HPP
