//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimulationData.cpp

  Data structures for array and scope simulation data

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/14/2007

  $Id: VSSimulationData.cpp,v 3.6 2008/01/11 04:57:04 sfegan Exp $

*/

#include<VSSimulationData.hpp>

using namespace VERITAS;

void VSHeaderSimulationDatum::clear()
{
  run_number         = 0;
  date_of_sims       = VSTime();
  simulation_package = 0;
  simulator_name     = "";
  date_of_array      = VSTime();
  corsika_atm_model  = 0;
  obs_altitude_m     = 0;
  sim_config         = "";
  database_name      = "";
  tables.clear();
}

void VSHeaderSimulationDatum::load(VSOctaveH5ReaderStruct* reader) 
{
  clear();
  reader->readCompositeHere(*this);
  reader->readCompositeVector("table",tables);
}

void VSHeaderSimulationDatum::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeCompositeHere(*this);
  writer->writeCompositeVector("table",tables);
}

VSArraySimulationWriter::
VSArraySimulationWriter(VSOctaveH5WriterStruct* s, unsigned ntel = 0):
  m_array_writer()
#if 0
, m_scope_writer()
#endif
{
  m_array_writer = 
    s->writeCompositeExpandableVectorHere<VSArraySimulationDatum>();

#if 0
  VSOctaveH5WriterCellVector* c = s->writeCellVector("t",ntel);
  for(unsigned iscope=0;iscope<ntel;iscope++)
    m_scope_writer[iscope] =
      c->writeCompositeExpandableVector<VSScopeSimulationDatum>(iscope);
#endif
}

VSArraySimulationWriter::~VSArraySimulationWriter()
{
  delete m_array_writer;
#if 0
  unsigned nscope = m_scope_writer.size();
  for(unsigned iscope=nscope;iscope<m_scope_writer.size();iscope++)
    delete m_scope_writer[iscope];
#endif
}

bool VSArraySimulationWriter::append(const VSArraySimulationDatum& x)
{
  bool status = m_array_writer->append(x);

#if 0
  for(unsigned iscope=0;iscope<m_scope_writer.size();iscope++)
    if(m_scope_writer[iscope])
      status &= m_scope_writer[iscope]->append(*x.scope[iscope]);
#endif
  return status;
}

VSArraySimulationReader::
VSArraySimulationReader(VSOctaveH5ReaderStruct* s,
			const MemberSubset& array_subset):
  m_array_reader()
#if 0
, m_scope_reader()
#endif
{
  m_array_reader = 
    s->readCompositeVectorElementsHere<VSArraySimulationDatum>(std::string(),
							       array_subset);

#if 0
  VSOctaveH5ReaderCellVector* c = s->readCellVector("t");
  vsassert(c);
  nscope = c->dimensions();

  m_scope_reader.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(!c->isEmpty(iscope))
      m_scope_reader[iscope] =
	c->readCompositeVectorElements<VSScopeSimulationDatum>(iscope);
#endif
}

VSArraySimulationReader::~VSArraySimulationReader()
{
  delete m_array_reader;
#if 0
  unsigned nscope = m_scope_reader.size();
  for(unsigned iscope=nscope;iscope<m_scope_reader.size();iscope++)
    delete m_scope_reader[iscope];
#endif
}

bool VSArraySimulationReader::element(VSArraySimulationDatum& x, 
				      unsigned index)
{
  bool status = m_array_reader->element(x, index);

#if 0
  unsigned nscope = m_scope_reader.size();
  if(x.scope.size() > nscope)
    for(unsigned iscope=nscope;iscope<m_scope_reader.size();iscope++)
      delete x.scope[iscope];
  if(x.scope.size() != nscope)x.scope.resize(nscope);

  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(m_scope_reader[iscope])
      {
	if(x.scope[iscope]==0)x.scope[iscope] = new VSScopeSimulationDatum;
	status &= m_scope_reader[iscope]->element(*x.scope[iscope], index);
      }
    else if(x.scope[iscope]!=0)
      {
	delete x.scope[iscope];
	x.scope[iscope] = 0;
      }
#endif

  return status;
}

bool VSArraySimulationReader::
loadAllEvents(VSOctaveH5ReaderStruct* s, 
	      std::vector<VSArraySimulationDatum>& x)
{
  bool status = true;

  x.clear();
 
  VSArraySimulationReader* reader = new VSArraySimulationReader(s);
  unsigned nevent = reader->rows();
  x.resize(nevent);

  for(unsigned ievent=0;ievent<nevent;ievent++)
    status &= reader->element(x[ievent],ievent);

  return status;
}
