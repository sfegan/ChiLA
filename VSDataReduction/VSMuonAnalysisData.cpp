//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMuonAnalysisData.cpp

  Data structures for muon analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/13/2007

  $Id: VSMuonAnalysisData.cpp,v 3.4 2007/12/04 18:05:03 sfegan Exp $

*/

#include<VSMuonAnalysisData.hpp>

using namespace VERITAS;

void VSMuonAnalysisData::load(VSOctaveH5ReaderStruct* reader)
{
  clear();

  VSOctaveH5ReaderCellVector* c = reader->readCellVector("scope");
  vsassert(c);
  unsigned nscope = c->dimensions();

  scope.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!c->isEmpty(iscope))
      c->readCompositeVector(iscope,scope[iscope]);

  delete c;
}

void VSMuonAnalysisData::save(VSOctaveH5WriterStruct* writer) const
{
  unsigned nscope = scope.size();
  VSOctaveH5WriterCellVector* c = writer->writeCellVector("scope",nscope);

  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!scope[iscope].empty())
      c->writeCompositeVector(iscope,scope[iscope]);

  delete c;
}
