//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimpleCutsCalc.hpp

  Class for applying simple box cuts.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/08/2007

  $Id: VSSimpleCutsCalc.hpp,v 3.9 2008/12/04 02:34:57 matthew Exp $

*/

#ifndef VSSIMPLECUTSCALC_HPP
#define VSSIMPLECUTSCALC_HPP

#include<VSCutsData.hpp>
#include<VSEventData.hpp>
#include<VSCutsCalc.hpp>
#include<VSOptions.hpp>
#include<VSH5DatumElement.hpp>

namespace VERITAS
{
  // ==========================================================================
  // VSSimpleCut
  // ==========================================================================
  template<typename T>
  class VSSimpleCut
  {
  public:
    VSSimpleCut();
    VSSimpleCut(const std::string& _datum_element_name);
    virtual ~VSSimpleCut();
    
    virtual bool evaluateCut(const T& datum) = 0;
    virtual std::string print() = 0;

    std::string getName() { return m_datum_element->getName(); }

    void initialize(const std::string& datum_element_name)
    {
      delete m_datum_element;
      m_datum_element = 
	VSH5DatumElement<T>::createDatumElement(datum_element_name);
    }

    static bool isRangedCut(VSOctaveH5ReaderStruct* reader)
    {
      std::string cut_type;
      reader->readString("cut_type",cut_type);
      if(cut_type == "ranged")
	return true;
      else
	return false;
    }

    static bool isPatternCut(VSOctaveH5ReaderStruct* reader)
    {
      std::string cut_type;
      reader->readString("cut_type",cut_type);
      if(cut_type == "pattern")
	return true;
      else
	return false;
    }

    virtual bool save(VSOctaveH5WriterStruct* writer) = 0;
    virtual bool load(VSOctaveH5ReaderStruct* reader) = 0;

    VSH5DatumElement<T>* m_datum_element;
  };

  template< typename T >
  VSSimpleCut<T>::VSSimpleCut():
    m_datum_element(NULL)
  {
    
  }
  
  template< typename T >
  VSSimpleCut<T>::VSSimpleCut(const std::string& _datum_element_name):
    m_datum_element(NULL)
  {
    m_datum_element = 
      VSH5DatumElement<T>::createDatumElement(_datum_element_name);
  }

  template< typename T >
  VSSimpleCut<T>::~VSSimpleCut()
  {
    delete m_datum_element;
  }

  // ==========================================================================
  // VSSimpleRangedCut
  // ==========================================================================
  template<typename T>
  class VSSimpleRangedCut : public VSSimpleCut<T>
  {
  public:    
    VSSimpleRangedCut();
    VSSimpleRangedCut(const std::string& _datum_element_name,
		      const std::string& _lo_cut, const std::string& _hi_cut);
    VSSimpleRangedCut(const std::string& _datum_element_name,
		      double _lo_cut, bool _has_lo_cut,
		      double _hi_cut, bool _has_hi_cut);
    virtual ~VSSimpleRangedCut() { }

    bool evaluateCut(const T& datum);

    bool save(VSOctaveH5WriterStruct* writer);
    bool load(VSOctaveH5ReaderStruct* reader);

    std::string print()
    {
      std::ostringstream os;
      os << std::setw(30) 
	 << VSH5DatumElementParser::getElement(datum_element_name);

      if(has_lo_cut)
	os << std::setw(20) << lo_cut;
      else
	os << std::setw(20) << "-";

      if(has_hi_cut)
	os << std::setw(20) << hi_cut;
      else
	os << std::setw(20) << "-";
	  
      return os.str();
    }

    double getLoCut() { return lo_cut; }
    std::string getLoCutString() 
    {
      if(has_lo_cut)
	return VSDataConverter::toString(lo_cut,true);
      else
	return "-";
    }

    double getHiCut() { return hi_cut; }
    std::string getHiCutString() 
    {
      if(has_hi_cut)
	return VSDataConverter::toString(hi_cut,true);
      else
	return "-";
    }

    bool hasLoCut() { return has_lo_cut;} 
    bool hasHiCut() { return has_hi_cut; }

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSSimpleRangedCut,datum_element_name);
      H5_ADDMEMBER(c,VSSimpleRangedCut,lo_cut);
      H5_ADDMEMBER(c,VSSimpleRangedCut,hi_cut);
      H5_ADDMEMBER(c,VSSimpleRangedCut,has_lo_cut);
      H5_ADDMEMBER(c,VSSimpleRangedCut,has_hi_cut);
    }

  private:
    VSSimpleRangedCut(const VSSimpleRangedCut& cut);
    VSSimpleRangedCut& operator=(const VSSimpleRangedCut& cut);

    std::string datum_element_name;
    double      lo_cut;
    double      hi_cut;
    bool        has_lo_cut;
    bool        has_hi_cut;
  };

  template< typename T >
  VSSimpleRangedCut<T>::VSSimpleRangedCut():
    VSSimpleCut<T>(), lo_cut(), hi_cut(), has_lo_cut(), has_hi_cut()
  {

  }

  template< typename T >
  VSSimpleRangedCut<T>::
  VSSimpleRangedCut(const std::string& _datum_element_name,
		    double _lo_cut, bool _has_lo_cut,
		    double _hi_cut, bool _has_hi_cut):
    VSSimpleCut<T>(_datum_element_name),
    datum_element_name(_datum_element_name),
    lo_cut(_lo_cut),hi_cut(_hi_cut),
    has_lo_cut(_has_lo_cut),has_hi_cut(_has_hi_cut)
  {

  }
  
  template< typename T >
  VSSimpleRangedCut<T>::
  VSSimpleRangedCut(const std::string& _datum_element_name,
		    const std::string& _lo_cut,
		    const std::string& _hi_cut):
    VSSimpleCut<T>(_datum_element_name),
    datum_element_name(_datum_element_name)
  {
    if(_lo_cut == "-" || _lo_cut.empty())
      has_lo_cut = false;
    else
      {
	has_lo_cut = true;
	VSDatumConverter<double>::fromString(lo_cut,_lo_cut.c_str());
      }

    if(_hi_cut == "-" || _hi_cut.empty())
      has_hi_cut = false;
    else
      {
	has_hi_cut = true;
	VSDatumConverter<double>::fromString(hi_cut,_hi_cut.c_str());
      }
  }
  
  template< typename T >
  bool VSSimpleRangedCut<T>::evaluateCut(const T& datum)
  {
    double x;
    if(!VSSimpleCut<T>::m_datum_element->getValue(datum,x))
      return false;

    if((has_lo_cut && x < lo_cut) || (has_hi_cut && x >hi_cut))
      return false;      
    else
      return true;
  }
  
  template< typename T >
  bool VSSimpleRangedCut<T>::save(VSOctaveH5WriterStruct* writer)
  {
    writer->writeString("cut_type","ranged");
    if(!writer->writeCompositeHere(*this))
      return false;
    else
      return true;
  }

  template< typename T >
  bool VSSimpleRangedCut<T>::load(VSOctaveH5ReaderStruct* reader)
  {
    std::string cut_type;
    reader->readString("cut_type",cut_type);
    vsassert(cut_type == "ranged");
    reader->readCompositeHere(*this);
    VSSimpleCut<T>::initialize(datum_element_name);
    return true;
  }

  // ==========================================================================
  // VSSimplePatternCut
  // ==========================================================================
  template<typename T>
  class VSSimplePatternCut : public VSSimpleCut<T>
  {
  public:    
    VSSimplePatternCut();
    VSSimplePatternCut(const std::string& _datum_element_name,
		       const std::string& _patterns);
    virtual ~VSSimplePatternCut() { }

    bool evaluateCut(const T& datum);

    bool save(VSOctaveH5WriterStruct* writer);
    bool load(VSOctaveH5ReaderStruct* reader);

    std::string print()
    {
       std::ostringstream os;
       os << std::setw(30) 
	  << VSH5DatumElementParser::getElement(datum_element_name);

       std::string patterns_string;
       if(exclusive)
	 patterns_string += "!";

       patterns_string += VSDataConverter::toString(patterns);

       os << std::setw(20) << patterns_string;
       return os.str();
    }

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSSimplePatternCut,datum_element_name);
      H5_ADDMEMBER(c,VSSimplePatternCut,exclusive);
    }

  private:
    VSSimplePatternCut(const VSSimplePatternCut& cut);
    VSSimplePatternCut& operator=(const VSSimplePatternCut& cut);

    std::string             datum_element_name;
    bool                    exclusive;
    std::vector< unsigned > patterns;
  };

  template< typename T >
  VSSimplePatternCut<T>::VSSimplePatternCut():
    VSSimpleCut<T>(), exclusive(), patterns()
  {

  }

  template< typename T >
  VSSimplePatternCut<T>::
  VSSimplePatternCut(const std::string& _datum_element_name,
		     const std::string& _patterns):
    VSSimpleCut<T>(_datum_element_name),
    datum_element_name(_datum_element_name),
    exclusive(),
    patterns()
  {
    if(_patterns.substr(0,1) == "!")
      {
	exclusive = true;
	VSDatumConverter< std::vector<unsigned> >::
	  fromString(patterns,_patterns.substr(1,_patterns.size()-1).c_str());
      }
    else
      {
	exclusive = false;
	VSDatumConverter< std::vector<unsigned> >::
	  fromString(patterns,_patterns.c_str());
      }
  }
  
  template< typename T >
  bool VSSimplePatternCut<T>::evaluateCut(const T& datum)
  {
    double x;
    if(!VSSimpleCut<T>::m_datum_element->getValue(datum,x)) return false;

    for(std::vector< unsigned >::iterator itr = patterns.begin();
	itr != patterns.end(); ++itr)
      {
	if(exclusive && (unsigned)lround(x) == *itr) return false;
	else if(!exclusive && (unsigned)lround(x) == *itr) return true;
      }

    if(exclusive) return true;
    else return false;
  }
  
  template< typename T >
  bool VSSimplePatternCut<T>::save(VSOctaveH5WriterStruct* writer)
  {
    writer->writeString("cut_type","pattern");
    writer->writeVector("patterns",patterns);
    if(!writer->writeCompositeHere(*this))
      return false;
    else
      return true;
  }

  template< typename T >
  bool VSSimplePatternCut<T>::load(VSOctaveH5ReaderStruct* reader)
  {
    std::string cut_type;
    reader->readString("cut_type",cut_type);
    vsassert(cut_type == "pattern");
    reader->readVector("patterns",patterns);
    reader->readCompositeHere(*this);
    VSSimpleCut<T>::initialize(datum_element_name);
    return true;
  }
  
  // ==========================================================================
  // VSSimpleCutsCalc
  // ==========================================================================
  class VSSimpleCutsCalc : public VSCutsCalc
  {
  public:
    class Options
    {
    public:
      typedef triple<std::string,std::string,std::string> Triple;
      typedef std::pair<std::string,std::string> Pair;
      Options(): 
	simple_ranged_cuts(), simple_pattern_cuts(), simple_cuts_file(),
	nscope(2)
      {}
      
      std::vector< Triple > simple_ranged_cuts;
      std::vector< Pair >   simple_pattern_cuts;
      std::string           simple_cuts_file;  
      unsigned              nscope;
    };

    typedef VSSimpleCut<VSEventArrayDatum> ArrayCut;
    typedef VSSimpleCut<VSEventScopeDatum> ScopeCut;

    VSSimpleCutsCalc(const Options& opt = s_default_options);
    virtual ~VSSimpleCutsCalc();

    void getCutResults(VSArrayCutsDatum& cut_results, 
		       const VSEventArrayDatum& event);
			
    virtual void clear()
    {
      m_array_cuts.clear();
      m_scope_cuts.clear();
    }

    virtual bool load(VSOctaveH5ReaderStruct* reader);
    virtual bool load(const std::string& text_file);
    virtual bool save(VSOctaveH5WriterStruct* writer);
    
    virtual void getScopeParamSet(std::set< std::string >& scope_set);
    virtual void getArrayParamSet(std::set< std::string >& array_set);

    bool evaluateArrayCuts(const VSEventArrayDatum& event);
    bool evaluateScopeCuts(unsigned iscope, const VSEventScopeDatum& scope);

    void print();
    void dump(std::ostream& stream) const;

    void loadRangedCut(const std::string& datum_element,
		       const std::string& lo_cut,
		       const std::string& hi_cut);

    void loadPatternCut(const std::string& datum_element,
			const std::string& pattern);

    static void configure(VSOptions& options,
			  const std::string& profile = "", 
			  const std::string& opt_prefix = "");
      
  private:

    Options                                 m_options;

    std::vector< ArrayCut* >                m_array_cuts;
    std::vector< std::vector< ScopeCut* > > m_scope_cuts;

    static Options                          s_default_options;
  };
}

#endif // defined VSSIMPLECUTSCALC_HPP
