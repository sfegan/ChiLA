//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleCutsCalc.cpp

  Class for applying simple box cuts.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/09/2007

  $Id: VSSimpleCutsCalc.cpp,v 3.13 2010/10/20 02:56:15 matthew Exp $

*/

#include <fstream>

#include <VSLineTokenizer.hpp>
#include <VSSimpleCutsCalc.hpp>
#include <VSFileUtility.hpp>

using namespace VERITAS;

VSSimpleCutsCalc::Options VSSimpleCutsCalc::s_default_options;

VSSimpleCutsCalc::VSSimpleCutsCalc(const Options& opt):
  VSCutsCalc(), m_options(opt), 
  m_array_cuts(), m_scope_cuts()
{
  VSFileUtility::expandFilename(m_options.simple_cuts_file);
  if(VSFileUtility::isFile(m_options.simple_cuts_file))
    {
      std::cout << "Loading simple cuts from file: " 
		<< m_options.simple_cuts_file << std::endl;

      load(m_options.simple_cuts_file);
    }
  else if(!m_options.simple_cuts_file.empty())
    {
      std::cerr << m_options.simple_cuts_file << " is not a valid file."
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  for(std::vector<triple<std::string,std::string,std::string> >::iterator
	itr = m_options.simple_ranged_cuts.begin(); 
      itr != m_options.simple_ranged_cuts.end(); ++itr)
    loadRangedCut(itr->first,itr->second,itr->third);

  for(std::vector<std::pair<std::string,std::string> >::iterator
	itr = m_options.simple_pattern_cuts.begin(); 
      itr != m_options.simple_pattern_cuts.end(); ++itr)
    loadPatternCut(itr->first,itr->second);
}

VSSimpleCutsCalc::~VSSimpleCutsCalc()
{

}

void VSSimpleCutsCalc::getCutResults(VSArrayCutsDatum& cut_results, 
				     const VSEventArrayDatum& event)
{
  const unsigned nscope_results = cut_results.scope.size();
  const unsigned nscope = event.scope.size();
  if(nscope_results < nscope)
    {
      cut_results.scope.resize(nscope);
      for(unsigned iscope = nscope_results; iscope < nscope; iscope++)
	if(cut_results.scope[iscope] == NULL)
	  cut_results.scope[iscope] = new VSScopeCutsDatum;
    }

  cut_results.passed_array_cut = evaluateArrayCuts(event);

  unsigned npassed = 0;
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      if(event.scope[iscope])
	{
	  vsassert(cut_results.scope[iscope]);
	  
	  bool passed = event.scope[iscope]->has_image &&
	    event.scope[iscope]->used_in_reconstruction
	    && evaluateScopeCuts(iscope,*event.scope[iscope]);
	  
	  cut_results.scope[iscope]->passed_scope_cut = passed;
	  if(passed)npassed++;
	}
      else if(cut_results.scope[iscope])
	{
	  cut_results.scope[iscope]->passed_scope_cut = false;
	}      
    }
  
  cut_results.npassed_scope_cut = npassed;

  if(cut_results.passed_array_cut && 
     cut_results.npassed_scope_cut >= m_options.nscope)
    cut_results.passed_cuts = true;
  else
    cut_results.passed_cuts = false;
}

bool VSSimpleCutsCalc::load(VSOctaveH5ReaderStruct* reader)
{
  clear();

  VSCutsCalc::load(reader);

  reader->readScalar("nscope_cut",m_options.nscope);
  VSOctaveH5ReaderCellVector* c = reader->readCellVector("array_cuts");
  vsassert(c);
  unsigned narray_cuts = c->dimensions();

  for(unsigned icut = 0; icut < narray_cuts; icut++)
    {
      VSOctaveH5ReaderStruct* s2 = c->readStruct(icut);

      if(VSSimpleCut<VSEventArrayDatum>::isRangedCut(s2))
	{
	  VSSimpleCut<VSEventArrayDatum>* array_cut =
	    new VSSimpleRangedCut<VSEventArrayDatum>;
	  array_cut->load(s2);
	  m_array_cuts.push_back(array_cut);

	}
      else if(VSSimpleCut<VSEventArrayDatum>::isPatternCut(s2))
	{
	  VSSimpleCut<VSEventArrayDatum>* array_cut =
	    new VSSimplePatternCut<VSEventArrayDatum>;
	  array_cut->load(s2);
	  m_array_cuts.push_back(array_cut);
	}

      delete s2;
    }
  delete c;
      
  c = reader->readCellVector("t");
  unsigned nscope = c->dimensions();
  m_scope_cuts.resize(nscope);
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      VSOctaveH5ReaderStruct* ws = c->readStruct(iscope);
      vsassert(ws);
      VSOctaveH5ReaderCellVector* c2 = 
	ws->readCellVector("scope_cuts");
      vsassert(c2);
      unsigned nscope_cuts = c2->dimensions();

      for(unsigned icut = 0; icut < nscope_cuts; icut++)
	{
	  VSOctaveH5ReaderStruct* s2 = c2->readStruct(icut);

	  if(VSSimpleCut<VSEventScopeDatum>::isRangedCut(s2))
	    {
	      VSSimpleCut<VSEventScopeDatum>* scope_cut =
		new VSSimpleRangedCut<VSEventScopeDatum>;	     
	      scope_cut->load(s2);
	      m_scope_cuts[iscope].push_back(scope_cut);
	    }
	  else if(VSSimpleCut<VSEventScopeDatum>::isPatternCut(s2))
	    {
	      VSSimpleCut<VSEventScopeDatum>* scope_cut =
		new VSSimplePatternCut<VSEventScopeDatum>;	     
	      scope_cut->load(s2);
	      m_scope_cuts[iscope].push_back(scope_cut);
	    }

	  delete s2;
	}
      delete c2;
      delete ws;
    }
  delete c;

  return true;
}

bool VSSimpleCutsCalc::load(const std::string& cuts_file)
{
  clear();

  if(VSOctaveH5ReaderStruct::isHDF5(cuts_file))
    {
      VSOctaveH5Reader *reader = new VSOctaveH5Reader(cuts_file);
      bool ret;
      if(reader->isStruct("simple_cuts"))
	{
	  VSOctaveH5ReaderStruct* s = reader->readStruct("simple_cuts");
	  ret = load(s);
	  delete s;
	}
      else ret = load(reader);
      delete reader;
      return ret;
    }
    
  std::ifstream datastream(cuts_file.c_str());
  //  if(!datastream)return false;

  VSLineTokenizer tokenizer;
  VSTokenList tokens;
  std::string line;
  while(getline(datastream,line))
    {
      std::string line_copy = line;
      tokenizer.tokenize(line_copy, tokens);
      if(line_copy.substr(0,1) == "#" || tokens.size() == 0) continue;

      if(tokens[0].string() == "date" && tokens.size() == 3)
	setDateRange(tokens[1].string(),tokens[2].string());
      else if(tokens[0].string() == "date" && tokens.size() == 5)
	{
	  VSTime lo_date;
	  VSTime hi_date;

	  std::string lo_date_string = tokens[1].string();
	  std::string lo_time_string = tokens[2].string();
	  std::string hi_date_string = tokens[3].string();
	  std::string hi_time_string = tokens[4].string();

	  if(!lo_date.setFromString(lo_date_string + " " + lo_time_string))
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": "
			<< "Error parsing timestamp "
			<< lo_date_string + " " + lo_time_string << std::endl;
	      continue;
	    }

	  if(!hi_date.setFromString(hi_date_string + " " + hi_time_string))
	    {
	      std::cerr << std::string(__PRETTY_FUNCTION__) + ": "
			<< "Error parsing timestamp "
			<< hi_date_string + " " + hi_time_string << std::endl;
	      continue;
	    }

	  setDateRange(lo_date,hi_date);
	}
      // Pattern cut ----------------------------------------------------------
      else if(tokens.size() == 2)
	loadPatternCut(tokens[0].string(),tokens[1].string());
      // Ranged cut -----------------------------------------------------------
      else if(tokens.size() == 3)
	loadRangedCut(tokens[0].string(),
		      tokens[1].string(),tokens[2].string());
      else
	{
	  std::cerr << "Error parsing the following line in the simple cuts "
		    << "definition file: " << std::endl << line << std::endl;
	}
    }

  return true;
}

bool VSSimpleCutsCalc::save(VSOctaveH5WriterStruct* writer)
{
  VSCutsCalc::save(writer);
  writer->writeScalar("nscope_cut",m_options.nscope);

  const unsigned narray_cuts = m_array_cuts.size();
  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("array_cuts",narray_cuts);
  vsassert(wc);
  for(unsigned icut = 0; icut < narray_cuts; icut++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(icut);
      m_array_cuts[icut]->save(ws);
      delete ws;
    }
  delete wc;

  const unsigned nscope = m_scope_cuts.size();
  wc = writer->writeCellVector("t",nscope);
  vsassert(wc);
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      const unsigned nscope_cuts = m_scope_cuts[iscope].size();
      VSOctaveH5WriterStruct* ws = wc->writeStruct(iscope);	
      VSOctaveH5WriterCellVector* c2 = 
	ws->writeCellVector("scope_cuts",nscope_cuts);
      for(unsigned icut = 0; icut < nscope_cuts; icut++)
	{
	  VSOctaveH5WriterStruct* s2 = c2->writeStruct(icut);
	  m_scope_cuts[iscope][icut]->save(s2);
	  delete s2;
	}
      delete c2;
      delete ws;
    }
  delete wc;
  return true;
}

bool VSSimpleCutsCalc::evaluateArrayCuts(const VSEventArrayDatum& event)
{
  for(std::vector< ArrayCut* >::iterator itr =
	m_array_cuts.begin(); itr != m_array_cuts.end(); ++itr)
    if(!(*itr)->evaluateCut(event)) return false;

  return true;
}

bool VSSimpleCutsCalc::evaluateScopeCuts(unsigned iscope,
					 const VSEventScopeDatum& scope)
{
  if(iscope >= m_scope_cuts.size()) return true;
  
  for(std::vector< ScopeCut* >::iterator itr =
	m_scope_cuts[iscope].begin(); itr != m_scope_cuts[iscope].end(); ++itr)
    if(!(*itr)->evaluateCut(scope)) return false;

  return true;
}

void VSSimpleCutsCalc::print()
{
  std::cout << std::string(30,'-') 
	    << std::setw(20) << "Simple Cuts"
	    << std::string(30,'-') 
	    << std::endl;

  std::cout << loDate() << " " << hiDate() << std::endl;

  std::cout << "Array Cuts:" << std::endl;

  for(std::vector<ArrayCut*>::iterator itr = m_array_cuts.begin(); itr != 
	m_array_cuts.end(); itr++)
    std::cout << std::setw(5) << "" << (*itr)->print() << std::endl;

  const unsigned nscope = m_scope_cuts.size();
  if(nscope == 0) return;

  std::cout << "Scope Cuts:" << std::endl;

  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      for(std::vector<ScopeCut*>::iterator itr = m_scope_cuts[iscope].begin(); 
	  itr != m_scope_cuts[iscope].end(); itr++)
	std::cout << std::setw(5) << iscope << (*itr)->print() << std::endl;
    }
}

void VSSimpleCutsCalc::getScopeParamSet(std::set< std::string >& scope_set)
{
  for(std::vector< std::vector< ScopeCut* > >::iterator scope_itr = 
	m_scope_cuts.begin(); scope_itr != m_scope_cuts.end();
      ++scope_itr)
    {
      for(std::vector< ScopeCut* >::iterator cut_itr = scope_itr->begin();
	  cut_itr != scope_itr->end(); ++cut_itr)
	scope_set.insert((*cut_itr)->getName());      
    }
  
  scope_set.insert("used_in_reconstruction");
  scope_set.insert("has_image");
}

void VSSimpleCutsCalc::getArrayParamSet(std::set< std::string >& array_set)
{
  for(std::vector< ArrayCut* >::iterator cut_itr = m_array_cuts.begin();
      cut_itr != m_array_cuts.end(); ++cut_itr)
    array_set.insert((*cut_itr)->getName());      
}

void VSSimpleCutsCalc::loadPatternCut(const std::string& datum_element,
				      const std::string& pattern)
{
  if(VSH5DatumElement<VSEventArrayDatum>::hasElement(datum_element))
    {
      VSSimpleCut<VSEventArrayDatum>* array_cut =
	new VSSimplePatternCut<VSEventArrayDatum>(datum_element,
						  pattern);
      m_array_cuts.push_back(array_cut);
    }
  else
    {
      std::cerr << "Invalid datum element name: "
		<< datum_element << std::endl;
      exit(EXIT_FAILURE);
    }  
}

void VSSimpleCutsCalc::loadRangedCut(const std::string& datum_element,
				     const std::string& lo_cut,
				     const std::string& hi_cut)
{
  if(VSH5DatumElement<VSEventArrayDatum>::hasElement(datum_element))
    {
      VSSimpleCut<VSEventArrayDatum>* array_cut =
	new VSSimpleRangedCut<VSEventArrayDatum>(datum_element,
						 lo_cut,hi_cut);
      m_array_cuts.push_back(array_cut);
    }
  else if(VSH5DatumElement<VSEventScopeDatum>::
	  hasElement(datum_element))
    {
      std::string scope_index = 
	VSH5DatumElementParser::getIndex(datum_element);
      
      std::vector<unsigned> scope_id_vec;
      
      if(scope_index.empty())
	{
	  std::cerr 
	    << "In  VSSimpleCutsCalc::load(): "
	    << "Error parsing scope index." << std::endl;
	  exit(EXIT_FAILURE);
	}
      else if(scope_index == "*")
	{
	  scope_id_vec.push_back(0);
	  scope_id_vec.push_back(1);		  
	  scope_id_vec.push_back(2);
	  scope_id_vec.push_back(3);
	}
      else
	{
	  VSDatumConverter< std::vector<unsigned> >::
	    fromString(scope_id_vec,scope_index.c_str());
	}
      
      for(std::vector<unsigned>::iterator itr = scope_id_vec.begin();
	  itr != scope_id_vec.end(); ++itr)
	{
	  if(*itr >= m_scope_cuts.size())
	    m_scope_cuts.resize(*itr+1);
	  
	  VSSimpleCut<VSEventScopeDatum>* scope_cut =
	    new VSSimpleRangedCut<VSEventScopeDatum>
	    (datum_element,lo_cut,hi_cut);
	  m_scope_cuts[*itr].push_back(scope_cut);
	}
    }
  else
    {
      std::cerr << "Invalid datum element name: "
		<< datum_element << std::endl;
      exit(EXIT_FAILURE);
    } 
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSSimpleCutsCalc::configure(VSOptions& options,
				 const std::string& profile, 
				 const std::string& opt_prefix)
				
{
  options.findWithValue(OPTNAME(opt_prefix,"simple_cuts"),
			s_default_options.simple_cuts_file,
			"Text file defining a set of cuts used "
			"to select events.",
			OPTNAME(opt_prefix,"cuts"));

  options.findWithValue(OPTNAME(opt_prefix,"simple_ranged_cuts"),
			s_default_options.simple_ranged_cuts,
			"Define rectangular cuts on one or more parameters "
			"as a comma-separated list.  Cuts on each parameter "
			"should be defined using the syntax: "
			"<parameter>/<locut>/<hicut>.  To apply only a low or "
			"high cut, leave the other cut blank, for "
			"example: "
			"msc_width//1.  Cuts defined on the command-line will "
			"be applied in addition to those defined in the cuts "
			"file."
			"",
			OPTNAME(opt_prefix,"cuts"));

  options.findWithValue(OPTNAME(opt_prefix,"simple_pattern_cuts"),
			s_default_options.simple_pattern_cuts,
			"Define pattern cuts on one or more parameters as "
			"as a comma-separated list.  Each cut should be "
			"defined using the syntax: "
			"<parameter>/<pattern>.  Pattern cuts can either be "
			"exclusive or inclusive.  By default the cut is "
			"satisfied when the parameter equals the pattern.  If "
			"the pattern is defined with a ! preceding it, then "
			"the cut will be satisified when the parameter does "
			"not equal the pattern.",
			OPTNAME(opt_prefix,"cuts"));

  options.findWithValue(OPTNAME(opt_prefix,"simple_cuts_nscope"),
			s_default_options.nscope,
			"Set the minimum number of telescopes that must "
			"pass simple cuts.",
			OPTNAME(opt_prefix,"cuts"));
}
