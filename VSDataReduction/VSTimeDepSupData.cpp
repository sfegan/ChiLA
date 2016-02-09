//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTimeDepSupData.cpp

  Time dependant suppressed data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       10/11/2006

  $Id: VSTimeDepSupData.cpp,v 3.2 2007/12/04 18:05:04 sfegan Exp $

*/

#include<VSTimeDepSupData.hpp>
#include<VSCentralizedDBAccess.hpp>

using namespace VERITAS;

unsigned VSTimeDepSupData::
getSuppressedCount(unsigned iscope, unsigned ichan) const  
{
  unsigned count = 0;
  for(unsigned islice=0;islice<m_nslice;islice++)
    if(isSuppressed(islice,iscope,ichan))count++;
  return count;
}

void VSTimeDepSupData::
suppressFromIMon(const std::vector<VBFRunInfo::Slice>& slice)
{
  if(!m_nslice)return;

  vsassert(m_nslice == slice.size());
  for(unsigned islice=0;islice<m_nslice;islice++)
    vsassert(slice[islice].event_num_hi-slice[islice].event_num_lo
	   == m_events_per_slice);

  const VSTime lo_time = slice[0].event_time_lo;
  const VSTime hi_time = slice[m_nslice-1].event_time_hi;

  const unsigned nscope = m_scope.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const unsigned nchan = m_scope[iscope].size();
      if(nchan==0)continue;

      std::vector<VSCentralizedDBAccess::PowerDatum> data;
      VSCentralizedDBAccess::getInstance()->
	getPower(iscope, lo_time, hi_time, data);

      for(unsigned idatum=0;idatum<data.size();idatum++)
	{
	  unsigned ichan = data[idatum].ichan;
	  if(ichan >= nchan)continue;

	  // Never unsuppress a slice 
	  if(data[idatum].power_on)continue;

	  VSTime& st(data[idatum].db_start_time);

	  unsigned islice=0;
	  while((islice<m_nslice)
		&&(st>=slice[islice].event_time_lo))islice++;
	  if((islice==m_nslice)&&(st>=hi_time))continue;
	  islice--;
	      
	  VSTime& et(data[idatum].db_end_time);
	  while((islice<m_nslice)&&(et>slice[islice].event_time_lo))
	    setSuppressed(islice++, iscope, ichan);
	}	  
    }
}

void VSTimeDepSupData::clear()
{
  m_events_per_slice = 0;
  m_nslice           = 0;
  m_scope.clear();
}

void VSTimeDepSupData::load(VSOctaveH5ReaderStruct* reader)
{
  m_scope.clear();

  reader->readScalar("events_per_slice", m_events_per_slice);
  reader->readScalar("nslice",m_nslice);

  VSOctaveH5ReaderCellVector* c = reader->readCellVector("scope");
  unsigned nscope = c->dimensions();
  m_scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
#if 0
    if(c->isCellVector(iscope,0))
      {
	VSOctaveH5ReaderCellVector* c2 = c->readCellVector(iscope);
	unsigned nchan = c2->dimensions();
	m_scope[iscope].resize(nchan);
	for(unsigned ichan=0;ichan<nchan;ichan++)
	  c2->readVector(ichan,m_scope[iscope][ichan]);
	delete c2;
      }
#else
    if(!c->isEmpty(iscope))
      {
	unsigned nchan;
	unsigned nword;
	c->dimensions(iscope,nword,nchan);
	uint32_t* data = new uint32_t[nchan*nword];
	c->readMatrix(iscope,data);
	m_scope[iscope].resize(nchan);
	for(unsigned ichan=0;ichan<nchan;ichan++)
	  {
	    m_scope[iscope][ichan].resize(nword);
	    for(unsigned iword=0;iword<nword;iword++)
	      m_scope[iscope][ichan][iword] = data[iword*nchan+ichan];
	  }
	delete[] data;
      }
#endif
    delete c;
}

void VSTimeDepSupData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeScalar("events_per_slice", m_events_per_slice);
  writer->writeScalar("nslice",m_nslice);
  unsigned nscope = m_scope.size();
  VSOctaveH5WriterCellVector* c = writer->writeCellVector("scope",nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nchan = m_scope[iscope].size();
#if 0
      if(nchan)
	{
	  VSOctaveH5WriterCellVector* c2 = c->writeCellVector(iscope,nchan);
	  for(unsigned ichan=0;ichan<nchan;ichan++)
	    c2->writeVector(ichan,m_scope[iscope][ichan]);
	  delete c2;
	}
#else
      unsigned nword = m_nslice?((m_nslice-1)>>5)+1:0;
      unsigned ndata = nchan*nword;
      if(ndata)
	{
	  uint32_t* data = new uint32_t[ndata];
	  for(unsigned ichan=0;ichan<nchan;ichan++)
	    for(unsigned iword=0;iword<nword;iword++)
	      data[iword*nchan+ichan] = m_scope[iscope][ichan][iword];
	  c->writeMatrix(iscope,nword,nchan,data);
	  delete[] data;
	}
#endif
    }
  delete c;
}
