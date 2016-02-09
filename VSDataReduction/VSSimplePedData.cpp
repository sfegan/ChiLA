//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimplePedData.hpp
  Simple pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/18/2005

  $Id: VSSimplePedData.cpp,v 3.1 2007/12/04 18:05:03 sfegan Exp $

*/

#include<fstream>
#include<memory>
#include<vsassert>

#include <VSLineTokenizer.hpp>

#include "VSSimpleStat.hpp"
#include "VSSimplePedData.hpp"

using namespace VERITAS;

bool VSSimplePedData::load(const std::string& filename)
{
  std::ifstream datastream(filename.c_str());
  if(!datastream)return false;
  
  m_data.clear();

  std::string line;
  VSLineTokenizer tokenizer;
  VSTokenList tokens;
  unsigned iline=0;
  while(getline(datastream,line))
    {
      iline++;
      std::string line_copy = line;
      tokenizer.tokenize(line_copy, tokens);
      if(tokens.size() == 0)continue;
      if(tokens.size() < 4)
        {
          std::cerr << filename << ": line " << iline
                    << ": " << tokens.size() 
		    << " columns found (>=4 expected)" << std::endl
                    << "Line: " << line << std::endl;
          continue;
        }
      
      unsigned iscope = 0;
      unsigned ichan = 0;
      double ped = 0;
      
      tokens[0].to(iscope);
      tokens[1].to(ichan);
      tokens[2].to(ped);

      if(iscope >= m_data.size())m_data.resize(iscope+1);
      if(ichan >= m_data[iscope].size())m_data[iscope].resize(ichan+1);
      m_data[iscope][ichan].ped = ped;
      unsigned nsample = tokens.size()-3;
      m_data[iscope][ichan].dev.resize(nsample);
      for(unsigned isample=0;isample<nsample;isample++)
	tokens[3+isample].to(m_data[iscope][ichan].dev[isample]);
      m_data[iscope][ichan].suppressed = false;
    }

  return true;
}

bool VSSimplePedData::save(const std::string& filename)
{
  std::ostream* filestream = 0;
  std::ostream* datastream = &std::cout;
  if(!filename.empty())
    {
      filestream = new std::ofstream(filename.c_str());
      if(!filestream->good())return false;
      datastream = filestream;
    }

  unsigned nscope = m_data.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nchan = m_data[iscope].size();
      for(unsigned ichan=0;ichan<nchan;ichan++)
	{
	  unsigned nsample = m_data[iscope][ichan].dev.size();
	  (*datastream) << iscope << ' ' << ichan << ' '
		     << m_data[iscope][ichan].ped;
	  for(unsigned isample=0;isample<nsample;isample++)
	    (*datastream) << ' ' << m_data[iscope][ichan].dev[isample];
	  (*datastream) << std::endl;
	}
    }
  delete filestream;
  return true;
}

void VSSimplePedData::suppress(const double lo, const double hi, 
			       const unsigned window)
{
  unsigned nscope = m_data.size();
  unsigned isample = window;

  // This loop to find window size to use for suppression (if necessary!)
  if(isample==0)
    {
      for(unsigned iscope=0;isample==0&&iscope<nscope;iscope++)
	{
	  unsigned nchan = m_data[iscope].size();
	  for(unsigned ichan=0;isample==0&&ichan<nchan;ichan++)
	    {
	      unsigned nsample = m_data[iscope][ichan].dev.size();
	      if(nsample)isample=nsample;
	    }
	}
      vsassert(isample!=0);
    }
  
  isample--;

  // This loop to suppress
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nchan = m_data[iscope].size();
      if(nchan==0)continue;
      std::vector<std::pair<bool,double> > dev_list(nchan);
      for(unsigned ichan=0;ichan<nchan;ichan++)
	dev_list[ichan].first=!m_data[iscope][ichan].suppressed,
	  dev_list[ichan].second=m_data[iscope][ichan].dev[isample];
      double meddev = median(dev_list);
      double locut = meddev*lo;
      double hicut = meddev*hi;

      for(unsigned ichan=0;ichan<nchan;ichan++)
	if((m_data[iscope][ichan].dev[isample]<locut)
	   ||(m_data[iscope][ichan].dev[isample]>hicut))
	  m_data[iscope][ichan].suppressed=true;
    }
}
