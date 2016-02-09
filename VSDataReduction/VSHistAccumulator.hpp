//-*-mode:c++; mode:font-lock;-*-

/*! \file VSHistAccumulator.hpp

  Class for filling histograms.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    1.0
  \date       07/29/2006

  $Id: VSHistAccumulator.hpp,v 3.4 2008/04/05 03:39:41 matthew Exp $

*/

#ifndef VSHISTACCUMULATOR_HPP
#define VSHISTACCUMULATOR_HPP

#include <VSOctaveIO.hpp>
#include <VSEventData.hpp>
#include <VSDatumElementExtractor.hpp>
#include <VSLineTokenizer.hpp>
#include <VSSimpleHist.hpp>

namespace VERITAS
{
  
  template <typename HIST>
  class VSHist1DAccumulator
  {
  public:
    VSHist1DAccumulator(const std::string& hist_def):
      m_hist(), m_datum_element()
    {
      std::string tmp = hist_def;

      VSLineTokenizer tokenizer;
      VSTokenList tokens;
      tokenizer.tokenize(tmp, tokens);

      double lo, hi, binsize;

      std::string element_name = tokens[0].string();
      VSDatumConverter< double >::
	fromString(binsize,tokens[1].string().c_str());
      VSDatumConverter< double >::
	fromString(lo,tokens[2].string().c_str());
      VSDatumConverter< double >::
	fromString(hi,tokens[3].string().c_str());

      vsassert(lo < hi);
      vsassert(binsize < hi-lo);

      std::string hist_name = VSH5DatumElementParser::getElement(element_name);

      if(VSH5DatumElement<VSEventArrayDatum>::hasElement(element_name))      
	m_datum_element = new VSEventArrayDatumElement(element_name);
      else if(VSH5DatumElement<VSEventScopeDatum>::hasElement(element_name))
	m_datum_element = new VSEventScopeDatumElement(element_name);

      m_hist = HIST(binsize,lo,hi,hist_name);
    }

    ~VSHist1DAccumulator()
    {
      delete m_datum_element;
    }

    HIST& getHist() { return m_hist; }

    void accumulate(VSEventArrayDatum& datum) 
    {
      double value;
      if(m_datum_element->getValue(datum,value)) m_hist.accumulate(value);
    }

  private:

    HIST                  m_hist;
    VSEventDatumElement*  m_datum_element;
  };

  template <typename HIST>
  class VSHist2DAccumulator
  {
  public:
    VSHist2DAccumulator(const std::string& hist_def):
      m_hist(), m_datum_element()
    {
      std::string tmp = hist_def;

      VSLineTokenizer tokenizer;
      VSTokenList tokens;
      tokenizer.tokenize(tmp, tokens);

      double xlo, xhi, xbinsize;
      double ylo, yhi, ybinsize;

      std::string element_name1 = tokens[0].string();
      std::string element_name2 = tokens[1].string();
      std::string hist_name1 = 
	VSH5DatumElementParser::getElement(element_name1);
      std::string hist_name2 = 
	VSH5DatumElementParser::getElement(element_name2);
      std::string hist_name = hist_name1 + ":" + hist_name2;

      if(tokens.size() == 5)
	{
	  VSDatumConverter< double >::
	    fromString(xbinsize,tokens[2].string().c_str());
	  VSDatumConverter< double >::
	    fromString(xlo,tokens[3].string().c_str());
	  VSDatumConverter< double >::
	    fromString(xhi,tokens[4].string().c_str());

	  m_hist = HIST(xbinsize,xlo,xhi,hist_name);
	}
      else if(tokens.size() == 8)
	{
	  VSDatumConverter< double >::
	    fromString(xbinsize,tokens[2].string().c_str());
	  VSDatumConverter< double >::
	    fromString(xlo,tokens[3].string().c_str());
	  VSDatumConverter< double >::
	    fromString(xhi,tokens[4].string().c_str());
	  VSDatumConverter< double >::
	    fromString(ybinsize,tokens[5].string().c_str());
	  VSDatumConverter< double >::
	    fromString(ylo,tokens[6].string().c_str());
	  VSDatumConverter< double >::
	    fromString(yhi,tokens[7].string().c_str());

	  m_hist = HIST(xbinsize,xlo,xhi,ybinsize,ylo,yhi,hist_name);
	}
      else
	{
	  std::cerr << "Error parsing histogram definition: "
		    << hist_def << std::endl;
	  exit(EXIT_FAILURE);
	}

      if(VSH5DatumElement<VSEventArrayDatum>::hasElement(element_name1))      
	m_datum_element.first = new VSEventArrayDatumElement(element_name1);
      else if(VSH5DatumElement<VSEventScopeDatum>::hasElement(element_name1))
	m_datum_element.first = new VSEventScopeDatumElement(element_name1);

      if(VSH5DatumElement<VSEventArrayDatum>::hasElement(element_name2))      
	m_datum_element.second = new VSEventArrayDatumElement(element_name2);
      else if(VSH5DatumElement<VSEventScopeDatum>::hasElement(element_name2))
	m_datum_element.second = new VSEventScopeDatumElement(element_name2);
    }

    ~VSHist2DAccumulator()
    {
      delete m_datum_element.first;
      delete m_datum_element.second;
    }

    HIST& getHist() { return m_hist; }

    void accumulate(VSEventArrayDatum& datum) 
    {
      double value1, value2;
      if(m_datum_element.first->getValue(datum,value1) &&
	 m_datum_element.second->getValue(datum,value2))
	m_hist.accumulate(value1,value2);
    }

  private:

    HIST                  m_hist;
    std::pair<VSEventDatumElement*,VSEventDatumElement*>  m_datum_element;
  };

} // namespace VERITAS

#endif // VSHISTACCUMULATOR_HPP
