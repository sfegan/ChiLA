//-*-mode:c++; mode:font-lock;-*-

/*! \file VSNSpaceCutsCalc.cpp

  Calculate whether event is in the array and scope filter volume

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/09/2007

  $Id: VSNSpaceCutsCalc.cpp,v 3.14 2009/10/08 05:26:34 matthew Exp $

*/

#include<VSNSpaceCutsCalc.hpp>
#include<VSNSpaceOctaveH5IO.hpp>
#include<VSFileUtility.hpp>

using namespace VERITAS;

// ============================================================================
// VSNSpaceFilterDatum<VSEventArrayDatum>
// ============================================================================
VSNSpaceFilterDatum<VSEventArrayDatum>::
VSNSpaceFilterDatum(const VSNSpace::Volume& vol):
  m_volume(vol), m_array_datums(), m_scope_datums(), m_is_scope_datum(),
  m_ndatum(0)
{
  for(std::vector<VSNSpace::Axis>::const_iterator itr = 
	m_volume.space().axes.begin();
      itr != m_volume.space().axes.end(); ++itr) 
    {
      std::string scope_index = 
	VSH5DatumElementParser::getIndex(itr->name);

      if(scope_index == "")
	{
	  m_array_datums.push_back(VSH5DatumElement<VSEventArrayDatum>::
				   createDatumElement(itr->name));
	  m_is_scope_datum.push_back(false);
	}
      else
	{
	  m_scope_datums.push_back(VSH5DatumElement<VSEventScopeDatum>::
				   createDatumElement(itr->name));
	  m_is_scope_datum.push_back(true);
	}
      m_ndatum++;
    }
}

VSNSpaceFilterDatum<VSEventArrayDatum>::~VSNSpaceFilterDatum()
{
  for(std::vector<VSH5DatumElement<VSEventArrayDatum>*>::iterator 
	itr = m_array_datums.begin(); itr != m_array_datums.end(); ++itr)
    delete (*itr);

  for(std::vector<VSH5DatumElement<VSEventScopeDatum>*>::iterator 
	itr = m_scope_datums.begin(); itr != m_scope_datums.end(); ++itr)
    delete (*itr);
}

bool VSNSpaceFilterDatum<VSEventArrayDatum>::
isDatumInVolume(const VSEventArrayDatum& x) const
{
  VERITAS::VSNSpace::Point point(m_ndatum); 
  for(unsigned i = 0; i < m_ndatum; i++)
    {
      vsassert(!m_is_scope_datum[i]);
      if(!m_array_datums[i]->getValue(x,point.x[i])) return false;
    }      

  return m_volume.isPointInVolume(point);
}

bool VSNSpaceFilterDatum<VSEventArrayDatum>::
isDatumInVolume(const VSEventArrayDatum& x, unsigned index) const
{
  VERITAS::VSNSpace::Point point(m_ndatum);   
  unsigned iscope = 0;
  unsigned iarray = 0;
  for(unsigned i = 0; i < m_ndatum; i++)
    {
      if(m_is_scope_datum[i])	
	{
	  if(!m_scope_datums[iscope]->getValue(*(x.scope[index]),point.x[i]))
	    return false;

	  iscope++;
	}
      else
	{
	  if(!m_array_datums[iarray]->getValue(x,point.x[i]))
	    return false;

	  iarray++;
	}
    }      

  return m_volume.isPointInVolume(point);
}

const VSNSpace::Volume& VSNSpaceFilterDatum<VSEventArrayDatum>::
getVolume() const
{
  return m_volume;
}

// ============================================================================
// VSNSpaceCutsCalc
// ============================================================================

VSNSpaceCutsCalc::Options VSNSpaceCutsCalc::s_default_options;

VSNSpaceCutsCalc::VSNSpaceCutsCalc(const Options& opt):
  VSCutsCalc(), m_options(opt),
  m_array_filter(), m_scope_filter_default(), m_scope_filter()
{
  if(!m_options.nspace_cuts_file.empty()) 
    VSFileUtility::expandFilename(m_options.nspace_cuts_file);

  if(VSFileUtility::isFile(m_options.nspace_cuts_file) && 
     !m_options.nspace_cuts_file.empty() &&
     VSOctaveH5ReaderStruct::isHDF5(m_options.nspace_cuts_file) )
    {
      std::cout << "Loading nspace cuts from file: " 
		<< m_options.nspace_cuts_file << std::endl;

      VSOctaveH5Reader *reader = 
	new VSOctaveH5Reader(m_options.nspace_cuts_file);

      if(reader->isStruct("nspace_cuts"))
	{
	  VSOctaveH5ReaderStruct* s = reader->readStruct("nspace_cuts");
	  load(s);
	  delete s;
	}
      else load(reader);
      delete reader;
    }
  else if(!m_options.nspace_cuts_file.empty())
    {
      std::cerr << m_options.nspace_cuts_file << " is not a valid file."
		<< std::endl;
      exit(EXIT_FAILURE);
    }
}

VSNSpaceCutsCalc::
VSNSpaceCutsCalc(const VSNSpace::Volume* array_vol,
		 const VSNSpace::Volume* scope_vol_default,
		 const std::vector< const VSNSpace::Volume* >& scope_vol):
  VSCutsCalc(), m_array_filter(), m_scope_filter_default(), m_scope_filter()
{
  if(array_vol)m_array_filter = new ArrayFilter(*array_vol);
  if(scope_vol_default)
    m_scope_filter_default = new ArrayFilter(*scope_vol_default);
  m_scope_filter.resize(scope_vol.size());
  for(unsigned iscope=0;iscope<m_scope_filter.size();iscope++)
    if(scope_vol[iscope])
      m_scope_filter[iscope] = new ScopeFilter(*scope_vol[iscope]);
}

VSNSpaceCutsCalc::~VSNSpaceCutsCalc()
{
  delete m_array_filter;
  delete m_scope_filter_default;
  for(unsigned iscope=0;iscope<m_scope_filter.size();iscope++)
    delete m_scope_filter[iscope];
}

void VSNSpaceCutsCalc::
getCutResults(VSArrayCutsDatum& cut_results, const VSEventArrayDatum& event)
{
  const unsigned nscope_results = cut_results.scope.size();
  const unsigned nscope = event.scope.size();
  if(nscope_results < nscope)
    {
      cut_results.scope.resize(nscope);
      for(unsigned iscope = nscope_results; iscope < nscope; iscope++)
	if(cut_results.scope[iscope] == NULL)
	  cut_results.scope[iscope] = new VSScopeCutsDatum;
    }

  if(m_array_filter)
    cut_results.passed_array_cut = m_array_filter->isDatumInVolume(event);
  else
    cut_results.passed_array_cut = true;

  unsigned npassed = 0;
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      if(event.scope[iscope] == NULL)
	continue;
      else if(iscope < m_scope_filter.size() && m_scope_filter[iscope])
	{
	  bool passed = event.scope[iscope]->has_image
	    && m_scope_filter[iscope]->isDatumInVolume(*event.scope[iscope]);
	  
	  cut_results.scope[iscope]->passed_scope_cut = passed;
	  if(passed)npassed++;
	}
      else if(m_scope_filter_default)
	{
	  // bool passed = event.scope[iscope]->has_image
// 	    && m_scope_filter_default->isDatumInVolume(*event.scope[iscope]);
	  bool passed = event.scope[iscope]->has_image
 	    && m_scope_filter_default->isDatumInVolume(event,iscope);

	  cut_results.scope[iscope]->passed_scope_cut = passed;
	  if(passed)npassed++;
	}
      else
	{
	  cut_results.scope[iscope]->passed_scope_cut = true;
	  npassed++;
	}
    }

  cut_results.npassed_scope_cut = npassed;

  if(cut_results.passed_array_cut && 
     cut_results.npassed_scope_cut >= m_options.nscope)
    cut_results.passed_cuts = true;
  else
    cut_results.passed_cuts = false;
}

void VSNSpaceCutsCalc::clear()
{
  delete m_array_filter, m_array_filter = NULL;
  delete m_scope_filter_default, m_scope_filter_default = NULL;

  const unsigned nscope = m_scope_filter.size();
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    delete m_scope_filter[iscope];

  m_scope_filter.clear();
}

bool VSNSpaceCutsCalc::load(VSOctaveH5ReaderStruct* reader)
{
  reader->readScalar("nscope_cut",m_options.nscope);

  VSNSpaceOctaveH5IO io;
  if(reader->isStruct("array"))
    {
      VSOctaveH5ReaderStruct* s2 = reader->readStruct("array");
      VSNSpace::Volume* array_vol = new VSNSpace::Volume;
      io.readVolume(s2,*array_vol);
      m_array_filter = new ArrayFilter(*array_vol);
      delete array_vol;
      delete s2;
    }

  if(reader->isStruct("scope_default"))
    {
      VSOctaveH5ReaderStruct* s2 = 
	reader->readStruct("scope_default");
      VSNSpace::Volume* scope_vol_default = new VSNSpace::Volume;
      io.readVolume(s2,*scope_vol_default);
      m_scope_filter_default = new ArrayFilter(*scope_vol_default);
      delete scope_vol_default;
      delete s2;
    }

  if(reader->isCellVector("scope"))
    {
      VSOctaveH5ReaderCellVector* c = reader->readCellVector("scope");
      unsigned nscope = c->dimensions();
      m_scope_filter.resize(nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	if(c->isStruct(iscope))
	  {
	    VSOctaveH5ReaderStruct* s2 = c->readStruct(iscope);
	    VSNSpace::Volume* scope_vol = new VSNSpace::Volume;
	    io.readVolume(s2,*scope_vol);
	    m_scope_filter[iscope] = new ScopeFilter(*scope_vol);
	    delete scope_vol;
	    delete s2;
	  }
      delete c;
    }

  return true;
}

bool VSNSpaceCutsCalc::save(VSOctaveH5WriterStruct* writer)
{
  writer->writeScalar("nscope_cut",m_options.nscope);

  VSNSpaceOctaveH5IO io;

  if(m_array_filter != NULL)
    {
      VSOctaveH5WriterStruct* ws = writer->writeStruct("array");
      io.writeVolume(ws,m_array_filter->getVolume());
    }

  if(m_scope_filter_default != NULL)
    {
      VSOctaveH5WriterStruct* ws = writer->writeStruct("scope_default");
      io.writeVolume(ws,m_scope_filter_default->getVolume());
    }

  if(m_scope_filter.size() > 0)
    {
      const unsigned nscope = m_scope_filter.size();
      VSOctaveH5WriterCellVector* wc = writer->writeCellVector("scope",nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  if(m_scope_filter[iscope] == NULL) continue;
	  
	  VSOctaveH5WriterStruct* ws = wc->writeStruct(iscope);
	  io.writeVolume(ws,m_scope_filter[iscope]->getVolume());
	}
    }

  return true;
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSNSpaceCutsCalc::configure(VSOptions& options,
				 const std::string& profile, 
				 const std::string& opt_prefix)
				
{
  options.findWithValue(OPTNAME(opt_prefix,"nspace_cuts"),
			s_default_options.nspace_cuts_file,
			"HDF5 file defining a set of nspace volumes "
			"used to select events.",
			OPTNAME(opt_prefix,"cuts"));

  options.findWithValue(OPTNAME(opt_prefix,"nspace_cuts_nscope"), 
			s_default_options.nscope,
			"Minimum number of telescopes which must pass all "
			"scope-level cuts.",
			OPTNAME(opt_prefix,"cuts"));
}
