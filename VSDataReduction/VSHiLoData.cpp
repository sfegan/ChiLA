//-*-mode:c++; mode:font-lock;-*-

/*! \file VSHiLoData.hpp
  Simple pedestal data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/18/2005

  $Id: VSHiLoData.cpp,v 3.3 2008/02/23 23:28:35 sfegan Exp $

*/

#include<vsassert>

#include "VSHiLoData.hpp"

using namespace VERITAS;

void VSHiLoData::clear()
{
  runno                = 0;
  has_hi_lo_gain_ratio = false;
  has_lo_gain_ped      = false;
  has_lo_gain_switch   = false;
  scope.clear();
}

void VSHiLoData::load(VSOctaveH5ReaderStruct* reader)
{
  clear();
  reader->readCompositeHere(*this);
  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("scope");
  vsassert(wc);

  unsigned nscope = wc->dimensions();
  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      if(wc->isEmpty(iscope))continue;
      
      VSOctaveH5ReaderStruct* wstel = wc->readStruct(iscope);	
      vsassert(wstel);
      wstel->readCompositeHere(scope[iscope]);
      wstel->readCompositeVector("chan",scope[iscope].chan);
      delete wstel;
    }

  delete wc;
}

void VSHiLoData::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);

  const unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("scope", nscope);
  vsassert(wc);

  for(unsigned iscope=0;iscope<nscope;iscope++) 
    {
      if(scope[iscope].chan.empty())continue;

      VSOctaveH5WriterStruct* wstel = wc->writeStruct(iscope);
      vsassert(wstel);
      wstel->writeCompositeHere(scope[iscope]);
      wstel->writeCompositeVector("chan",scope[iscope].chan);      
      delete wstel;	
    }

  delete wc;
}
