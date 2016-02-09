//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutsEvaluator.hpp

  Class responsible for evaluating whether an event passes cuts.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/08/2007

  $Id: VSCutsEvaluator.hpp,v 3.5 2009/10/08 05:25:15 matthew Exp $

*/

#ifndef VSCUTSEVALUATOR_HPP
#define VSCUTSEVALUATOR_HPP

#include<VSOptions.hpp>
#include<VSCutsData.hpp>
#include<VSCutsCalc.hpp>
#include<VSEventData.hpp>
#include<VSSimpleCutsCalc.hpp>
#include<VSNSpaceCutsCalc.hpp>

namespace VERITAS
{
  class VSCutsEvaluator
  {
  public:
    class Options
    {
    public:
      Options(): 
	simple_cuts_file(), nspace_cuts_file()
      {}
      
      std::vector<std::string>           simple_cuts_file;        
      std::string           nspace_cuts_file;        
    };

    VSCutsEvaluator(const Options& opt = s_default_options);
    ~VSCutsEvaluator();

    bool isSelected(const VSEventArrayDatum& event_data);
    void setDate(const VSTime& date);

    void getScopeParamSet(std::set<std::string>& scope_set);
    void getArrayParamSet(std::set<std::string>& array_set);

    void print();

    void load(VSOctaveH5ReaderStruct* reader);
    void save(VSOctaveH5WriterStruct* writer) const;

    static void configure(VSOptions& options,
			  const std::string& profile = "", 
			  const std::string& opt_prefix = "");

  private:

    Options                                 m_options;

    std::vector<VSSimpleCutsCalc*>   m_simple_cuts_calc;
    VSNSpaceCutsCalc*                m_nspace_cuts_calc;

    std::vector< VSCutsCalc* >       m_cuts_calc;
    std::vector< VSArrayCutsDatum >  m_cuts_results;

    static Options                   s_default_options;
  };
}

#endif // VSCUTSEVALUATOR_HPP
