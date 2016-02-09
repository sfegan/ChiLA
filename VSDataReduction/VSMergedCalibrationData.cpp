//-*-mode:c++; mode:font-lock;-*-

/*! \file VSMergedCalibrationData.cpp

  Time independent merged calibration data (laser, ped, suppression)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/15/2006

  $Id: VSMergedCalibrationData.cpp,v 3.7 2008/11/18 02:58:39 matthew Exp $

*/

#include<VSMergedCalibrationData.hpp>

using namespace VERITAS;

void VSScopeMergedCalibrationBase::load(VSOctaveH5ReaderStruct* s)
{
  s->readCompositeHere(*this);
  median_dev_hist.load(s->readStruct("median_dev_hist"));
}

void VSScopeMergedCalibrationBase::save(VSOctaveH5WriterStruct* s) const
{
  s->writeCompositeHere(*this);
  median_dev_hist.save(s->writeStruct("median_dev_hist"));
}

void VSArrayMergedCalibrationData::load(VSOctaveH5ReaderStruct* s)
{
  clear();  
  s->readComposite("array_info",*this);  
  if(s->isCellVector("scope_info"))
    {
      VSOctaveH5ReaderCellVector* wc_scope = s->readCellVector("scope_info");
      vsassert(wc_scope);
      unsigned nscope = wc_scope->dimensions();
      scope.resize(nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	if(!wc_scope->isEmpty(iscope))
	  {
	    VSOctaveH5ReaderStruct* ws = wc_scope->readStruct(iscope);	
	    vsassert(ws);
	    scope[iscope].load(ws);
	    delete ws;
	  }
      delete wc_scope;
    }
  else s->readCompositeVector("scope_info",scope);

  VSOctaveH5ReaderCellVector* wc_chan = s->readCellVector("channel_info");
  vsassert(wc_chan);
  unsigned nscope = wc_chan->dimensions();
  vsassert(scope.size()==nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!wc_chan->isEmpty(iscope))
      {
	VSOctaveH5ReaderStruct* ws = wc_chan->readStruct(iscope);	
	vsassert(ws);
	ws->readComposite("suppressed_count",scope[iscope].suppressed_count);
	ws->readCompositeVectorHere(scope[iscope].channel);
	delete ws;
      }
  delete wc_chan;
}

void VSArrayMergedCalibrationData::save(VSOctaveH5WriterStruct* s) const
{
  s->writeComposite("array_info",*this);
  const unsigned nscope = scope.size();  
  VSOctaveH5WriterCellVector* wc_scope = 
    s->writeCellVector("scope_info", nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope].nchan)
      {
	VSOctaveH5WriterStruct* ws = wc_scope->writeStruct(iscope);
	vsassert(ws);
	scope[iscope].save(ws);
	delete ws;
      }
  delete wc_scope;

  VSOctaveH5WriterCellVector* wc_chan = 
    s->writeCellVector("channel_info", nscope);
  vsassert(wc_chan);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(scope[iscope].nchan)
      {
	VSOctaveH5WriterStruct* ws = wc_chan->writeStruct(iscope);	
	vsassert(ws);
	ws->writeComposite("suppressed_count",scope[iscope].suppressed_count);
	ws->writeCompositeVectorHere(scope[iscope].channel);
	delete ws;
      }
  delete wc_chan;
}
