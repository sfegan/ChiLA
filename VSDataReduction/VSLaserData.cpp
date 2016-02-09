//-*-mode:c++; mode:font-lock;-*-

/*! \file VSLaserData.hpp
  Simple pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/18/2005

  $Id: VSLaserData.cpp,v 3.10 2009/11/24 18:25:11 matthew Exp $

*/

#include<fstream>
#include<memory>
#include<vsassert>
#include<cctype>

#include <VSLineTokenizer.hpp>
#include <VSFileUtility.hpp>
#include <VSSimpleStat.hpp>

#include "VSLaserData.hpp"

using namespace VERITAS;

bool VSLaserData::load(const std::string& filename)
{
  std::vector<std::string> laser_scope;
  VSDataConverter::fromString(laser_scope,filename);

  std::string default_laser_file;
  std::map<unsigned,std::string> laser_file;

  for(unsigned i = 0; i < laser_scope.size(); i++)
    {
      std::string::size_type pos = laser_scope[i].find_first_of('/');

      if(VSFileUtility::isFile(laser_scope[i]))
	default_laser_file = laser_scope[i];
      else if(pos != std::string::npos)
	{
	  std::string file = 
	    laser_scope[i].substr(pos,laser_scope[i].size()-pos);
	  std::string scope_string = laser_scope[i].substr(0,pos);
	  unsigned scope_id = 0;

	  if(!VSFileUtility::isFile(file))
	    {
	      std::cerr << "Laser file does not exist: " 
			<< file << std::endl;
	      return false;
	    }

	  VSDataConverter::fromString(scope_id,scope_string);

	  laser_file[scope_id] = file;
	}
      else
	{
	  std::cerr << "Laser file does not exist: " 
		    << laser_scope[i] << std::endl;
	  return false;
	}
    }

  if(!default_laser_file.empty())
    {
      VSFileUtility::expandFilename(default_laser_file);
      VSOctaveH5Reader* reader = new VSOctaveH5Reader(default_laser_file);
      vsassert(reader);

      if(!reader->isStruct("stage1.laser"))
	{
	  std::cerr << "No laser data found in file: " 
		    << default_laser_file << std::endl;
	  return false;
	}

      VSOctaveH5ReaderStruct* s = reader->readStruct("stage1.laser");
      vsassert(s);
      load(s);
      delete s;
      delete reader;
    }

  for(std::map<unsigned,std::string>::iterator itr = laser_file.begin();
      itr != laser_file.end(); ++itr)
    {
      VSFileUtility::expandFilename(itr->second);
      VSOctaveH5Reader* reader = new VSOctaveH5Reader(itr->second);
      vsassert(reader);

      if(!reader->isStruct("stage1.laser"))
	{
	  std::cerr << "No laser data found in file: " 
		    << itr->second << std::endl;
	  return false;
	}

      VSOctaveH5ReaderStruct* s = reader->readStruct("stage1.laser");
      vsassert(s);

      load(itr->first,s);

      delete s;
      delete reader;
    }

  return true;
}

void VSLaserData::suppress(const double lo, const double hi)
{
  unsigned nscope = scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nchan = scope[iscope].chan.size();
      for(unsigned ichan=0;ichan<nchan;ichan++)
	if((!scope[iscope].chan[ichan].uncalibrated)
	   &&((scope[iscope].chan[ichan].gain<lo)
	      ||(scope[iscope].chan[ichan].gain>hi)))
	  scope[iscope].chan[ichan].uncalibrated = true;
    }
}

void VSLaserData::clear()
{
  m_runno=0;
  m_threshold_nchan=0;
  m_threshold_dc=0;
  scope.clear();
}

void VSLaserData::load(VSOctaveH5ReaderStruct* reader)
{
  clear();
  reader->readCompositeHere(*this);
  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("scope");
  vsassert(wc);
  unsigned nscope = wc->dimensions();

  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      if(wc->isEmpty(iscope))
	continue;
      
      VSOctaveH5ReaderStruct* wstel = wc->readStruct(iscope);	
      vsassert(wstel);
      wstel->readCompositeHere(scope[iscope]);
      wstel->readCompositeVector("crate",scope[iscope].crate);
      scope[iscope].signal_hist.load(wstel->readStruct("signal_hist"));
      scope[iscope].nchan_flash_hist.load(wstel->readStruct("nchan_hist"));
      scope[iscope].nchan_logain_hist.
	load(wstel->readStruct("nchan_logain_hist"));
      scope[iscope].nchan_logain_hist.
	load(wstel->readStruct("nchan_higain_hist"));

      if(scope[iscope].runno == 0) scope[iscope].runno = m_runno;

      VSOctaveH5ReaderStruct* wschan = wstel->readStruct("chan");  
      vsassert(wschan);
      wschan->readCompositeVectorHere(scope[iscope].chan);
      VSOctaveH5ReaderCellVector* wchist = 
	wschan->readCellVector("signal_hist");
      vsassert(wchist);
      unsigned nchan = wchist->dimensions();

      for(unsigned ichan=0;ichan<nchan;ichan++) 
	{
	  VSOctaveH5ReaderStruct* s = wchist->readStruct(ichan);
	  scope[iscope].chan[ichan].signal_hist.load(s);
	  delete s;
	}

      delete wchist;
      delete wschan;
      delete wstel;
    }

  delete wc;
}

bool VSLaserData::load(unsigned iscope, VSOctaveH5ReaderStruct* reader)
{
  VSLaserData::ScopeData sd;

  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("scope");
  vsassert(wc);
  unsigned nscope = wc->dimensions();

  if(nscope <= iscope || wc->isEmpty(iscope))
    {
      delete wc;
      return false;
    }

  VSOctaveH5ReaderStruct* wstel = wc->readStruct(iscope);	
  vsassert(wstel);
  wstel->readCompositeHere(sd);
  wstel->readCompositeVector("crate",sd.crate);
  sd.signal_hist.load(wstel->readStruct("signal_hist"));
  sd.nchan_flash_hist.load(wstel->readStruct("nchan_hist"));
  sd.nchan_logain_hist.load(wstel->readStruct("nchan_logain_hist"));
  sd.nchan_logain_hist.load(wstel->readStruct("nchan_higain_hist"));
  
  VSOctaveH5ReaderStruct* wschan = wstel->readStruct("chan");  
  vsassert(wschan);
  wschan->readCompositeVectorHere(sd.chan);
  VSOctaveH5ReaderCellVector* wchist = wschan->readCellVector("signal_hist");
  vsassert(wchist);
  unsigned nchan = wchist->dimensions();
  
  for(unsigned ichan=0;ichan<nchan;ichan++) 
    {
      VSOctaveH5ReaderStruct* s = wchist->readStruct(ichan);
      sd.chan[ichan].signal_hist.load(s);
      delete s;
    }
  
  delete wchist;
  delete wschan;
  delete wstel;

  if(sd.runno == 0) reader->readScalar("m_runno",sd.runno);

  if(iscope >= scope.size()) scope.resize(iscope+1);
  scope[iscope] = sd;

  delete wc;
  return true;
}

void VSLaserData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);
  const unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("scope", nscope);

  for(unsigned iscope=0;iscope<nscope;iscope++) 
    {
      if(!scope[iscope].nchan)
	continue;

      VSOctaveH5WriterStruct* wstel = wc->writeStruct(iscope);
      vsassert(wstel);
      wstel->writeCompositeHere(scope[iscope]);
      wstel->writeCompositeVector("crate",scope[iscope].crate);
      scope[iscope].signal_hist.save(wstel->writeStruct("signal_hist"));
      scope[iscope].nchan_flash_hist.save(wstel->writeStruct("nchan_hist"));
      scope[iscope].nchan_logain_hist.
	save(wstel->writeStruct("nchan_logain_hist"));
      scope[iscope].nchan_higain_hist.
	save(wstel->writeStruct("nchan_higain_hist"));

      VSOctaveH5WriterStruct* wschan = wstel->writeStruct("chan");
      vsassert(wschan);
      wschan->writeCompositeVectorHere(scope[iscope].chan);      
      unsigned nchan = scope[iscope].chan.size();
      VSOctaveH5WriterCellVector* wchist = 
	wschan->writeCellVector("signal_hist",  nchan);
      vsassert(wchist);

      for(unsigned ichan = 0;ichan<nchan;ichan++)
	{
	  VSOctaveH5WriterStruct* s = wchist->writeStruct(ichan);
	  scope[iscope].chan[ichan].signal_hist.save(s);
	  delete s;
	}

      delete wchist;
      delete wschan;
      delete wstel;	
    }

  delete wc;
}
