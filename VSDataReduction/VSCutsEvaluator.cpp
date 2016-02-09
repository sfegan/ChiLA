//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutsEvaluator.cpp

  Class responsible for evaluating whether an event passes cuts.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/09/2007

  $Id: VSCutsEvaluator.cpp,v 3.5 2009/10/08 05:25:15 matthew Exp $

*/

#include<VSFileUtility.hpp>
#include<VSSimpleCutsCalc.hpp>
#include<VSNSpaceCutsCalc.hpp>
#include<VSCutsEvaluator.hpp>

using namespace VERITAS;

VSCutsEvaluator::Options VSCutsEvaluator::s_default_options;

VSCutsEvaluator::VSCutsEvaluator(const Options& opt):
  m_options(opt),
  m_simple_cuts_calc(),
  m_nspace_cuts_calc(),
  m_cuts_calc(),
  m_cuts_results()
{
  for(std::vector<std::string>::iterator itr = 
	m_options.simple_cuts_file.begin(); itr !=
	m_options.simple_cuts_file.end(); ++itr)
    {
      VSSimpleCutsCalc* simple_cuts = new VSSimpleCutsCalc;
      simple_cuts->load(*itr);
      m_simple_cuts_calc.push_back(simple_cuts);
      m_cuts_calc.push_back(simple_cuts);
    }

  m_nspace_cuts_calc = new VSNSpaceCutsCalc;
  m_cuts_calc.push_back(m_nspace_cuts_calc);

  m_cuts_results.resize(m_cuts_calc.size());
}

VSCutsEvaluator::~VSCutsEvaluator()
{

}

void VSCutsEvaluator::getScopeParamSet(std::set<std::string>& scope_set)
{
  const unsigned ncuts_calc = m_cuts_calc.size();
  for(unsigned icuts = 0; icuts < ncuts_calc; icuts++)
    m_cuts_calc[icuts]->getScopeParamSet(scope_set);
}

void VSCutsEvaluator::getArrayParamSet(std::set<std::string>& array_set)
{
  const unsigned ncuts_calc = m_cuts_calc.size();
  for(unsigned icuts = 0; icuts < ncuts_calc; icuts++)
    m_cuts_calc[icuts]->getArrayParamSet(array_set);
}

void VSCutsEvaluator::print()
{
  const unsigned ncuts_calc = m_cuts_calc.size();
  for(unsigned icuts = 0; icuts < ncuts_calc; icuts++)
    m_cuts_calc[icuts]->print();
}

void VSCutsEvaluator::load(VSOctaveH5ReaderStruct* reader)
{
  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("simple_cuts");
  const unsigned nsimple_cuts = wc->dimensions();
  for(unsigned icut = 0; icut < nsimple_cuts; icut++)
    {
      VSOctaveH5ReaderStruct* ws = wc->readStruct(icut);
      VSSimpleCutsCalc* simple_cuts = new VSSimpleCutsCalc;

      simple_cuts->load(ws);
      m_simple_cuts_calc.push_back(simple_cuts);
      m_cuts_calc.push_back(simple_cuts);      
    }
  
  m_nspace_cuts_calc->load(reader->readStruct("nspace_cuts"));
  m_cuts_calc.push_back(m_nspace_cuts_calc);
}

void VSCutsEvaluator::save(VSOctaveH5WriterStruct* writer) const
{
  const unsigned nsimple_cuts = m_simple_cuts_calc.size();
  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("simple_cuts",nsimple_cuts);
  for(unsigned icut = 0; icut < nsimple_cuts; icut++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(icut);
      m_simple_cuts_calc[icut]->save(ws);
      delete ws;
    }
  delete wc;

  m_nspace_cuts_calc->save(writer->writeStruct("nspace_cuts"));
}

bool VSCutsEvaluator::isSelected(const VSEventArrayDatum& event_data)
{
  bool is_selected = true;

  const unsigned ncuts_calc = m_cuts_calc.size();

  if(ncuts_calc == 0)
    return true;

  for(unsigned icuts = 0; icuts < ncuts_calc; icuts++)
    m_cuts_calc[icuts]->getCutResults(m_cuts_results[icuts], event_data);
  
  for(unsigned icuts = 0; icuts < ncuts_calc; icuts++)
    is_selected &= m_cuts_results[icuts].passed_cuts;
  
  return is_selected;
}

void VSCutsEvaluator::setDate(const VSTime& date)
{
  m_cuts_calc.clear();
  for(std::vector<VSSimpleCutsCalc*>::iterator itr = 
	m_simple_cuts_calc.begin(); itr != m_simple_cuts_calc.end(); ++itr)
    {
      if(date >= (*itr)->loDate() && date <= (*itr)->hiDate())
	m_cuts_calc.push_back(*itr);
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSCutsEvaluator::configure(VSOptions& options,
				const std::string& profile, 
				const std::string& opt_prefix)
				
{
  options.addCatagory(OPTNAME(opt_prefix,"cuts"),
 		      "Options for loading and applying cuts.");

  options.findWithValue(OPTNAME(opt_prefix,"simple_cuts"),
			s_default_options.simple_cuts_file,
			"Text file defining a set of cuts used "
			"to select events.",
			OPTNAME(opt_prefix,"cuts"));

  //  VSSimpleCutsCalc::configure(options,profile,opt_prefix);
  VSNSpaceCutsCalc::configure(options,profile,opt_prefix);
}

