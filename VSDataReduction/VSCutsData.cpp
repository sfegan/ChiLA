//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutsData.cpp

  Data structures for array and scope cuts

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/08/2007

  $Id: VSCutsData.cpp,v 3.5 2007/12/04 18:05:03 sfegan Exp $

*/

#include<VSCutsData.hpp>

using namespace VERITAS;

VSArrayCutsWriter::
VSArrayCutsWriter(VSOctaveH5WriterStruct* s,
		  const std::vector<unsigned>& nchan):
  m_array_writer(), m_scope_writer(nchan.size())
{
  m_array_writer = s->writeCompositeExpandableVectorHere<VSArrayCutsDatum>();

  if(nchan.size())
    {
      VSOctaveH5WriterCellVector* c = s->writeCellVector("t",nchan.size());
      for(unsigned iscope=0;iscope<nchan.size();iscope++)
	if(nchan[iscope])
	  m_scope_writer[iscope] =
	    c->writeCompositeExpandableVector<VSScopeCutsDatum>(iscope);
    }
}

VSArrayCutsWriter::~VSArrayCutsWriter()
{
  delete m_array_writer;
  unsigned nscope = m_scope_writer.size();
  for(unsigned iscope=nscope;iscope<m_scope_writer.size();iscope++)
    delete m_scope_writer[iscope];
}

bool VSArrayCutsWriter::append(const VSArrayCutsDatum& x)
{
  bool status = m_array_writer->append(x);

  vsassert(m_scope_writer.size() == x.scope.size());

  for(unsigned iscope=0;iscope<m_scope_writer.size();iscope++)
    if(m_scope_writer[iscope])
      status &= m_scope_writer[iscope]->append(*x.scope[iscope]);

  return status;
}

VSArrayCutsReader::VSArrayCutsReader():
  m_array_reader(), m_scope_reader()
{

}

VSArrayCutsReader::
VSArrayCutsReader(VSOctaveH5ReaderStruct* s):
  m_array_reader(), m_scope_reader()
{
  m_array_reader = s->readCompositeVectorElementsHere<VSArrayCutsDatum>();

  if(s->isCellVector("t"))
    {
      VSOctaveH5ReaderCellVector* c = s->readCellVector("t");
      vsassert(c);
      unsigned nscope = c->dimensions();

      m_scope_reader.resize(nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	if(!c->isEmpty(iscope))
	  m_scope_reader[iscope] =
	    c->readCompositeVectorElements<VSScopeCutsDatum>(iscope);
    }
}

VSArrayCutsReader::~VSArrayCutsReader()
{
  delete m_array_reader;
  unsigned nscope = m_scope_reader.size();
  for(unsigned iscope=nscope;iscope<m_scope_reader.size();iscope++)
    delete m_scope_reader[iscope];
}

bool VSArrayCutsReader::element(VSArrayCutsDatum& x, unsigned index)
{
  bool status = m_array_reader->element(x, index);

  unsigned nscope = m_scope_reader.size();
  if(x.scope.size() > nscope)
    for(unsigned iscope=nscope;iscope<m_scope_reader.size();iscope++)
      delete x.scope[iscope];
  if(x.scope.size() != nscope)x.scope.resize(nscope);

  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(m_scope_reader[iscope])
      {
	if(x.scope[iscope]==0)x.scope[iscope] = new VSScopeCutsDatum;
	status &= m_scope_reader[iscope]->element(*x.scope[iscope], index);
      }
    else if(x.scope[iscope]!=0)
      {
	delete x.scope[iscope];
	x.scope[iscope] = 0;
      }

  return status;
}

bool VSArrayCutsReader::loadAllArrayCuts(VSOctaveH5ReaderStruct* s, 
					 std::vector<VSArrayCutsDatum>& x)
{
  bool status = true;

  x.clear();
 
  VSArrayCutsReader* reader = new VSArrayCutsReader(s);
  unsigned nevent = reader->rows();
  x.resize(nevent);

  for(unsigned ievent=0;ievent<nevent;ievent++)
    status &= reader->element(x[ievent],ievent);

  return status;
}
