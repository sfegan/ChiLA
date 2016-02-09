//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventData.cpp

  Data structures for event analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/13/2007

  $Id: VSEventData.cpp,v 3.11 2008/07/29 22:34:22 matthew Exp $

*/

#include<VSEventData.hpp>

using namespace VERITAS;

VSEventArrayDatum::VSEventArrayDatum(const VSEventArrayDatum& o)
{
  *this = o;
}

VSEventArrayDatum& VSEventArrayDatum::operator=(const VSEventArrayDatum& o)
{
  if(this == &o) return *this;

  for(unsigned iscope = 0; iscope < scope.size(); iscope++)
    delete scope[iscope];

  scope.resize(o.scope.size());

  for(unsigned iscope = 0; iscope < o.scope.size(); iscope++)
    {
      if(o.scope[iscope])
	{
	  scope[iscope] = new VSEventScopeDatum;
	  *scope[iscope] = *o.scope[iscope];
	}
      else scope[iscope] = NULL;
    }

  event_num      = o.event_num;
  abs_event_time = o.abs_event_time;
  mjd            = o.mjd;
  event_time     = o.event_time;
  elapsed_ticks  = o.elapsed_ticks;
  live_ticks     = o.live_ticks;
  trigger_mask   = o.trigger_mask;
  l3_sent_mask   = o.l3_sent_mask;
  has_datum_mask = o.has_datum_mask;
  has_image_mask = o.has_image_mask;
  used_in_reconstruction_mask = o.used_in_reconstruction_mask;
  mean_array_zn  = o.mean_array_zn;
  mean_array_az  = o.mean_array_az;
  nscope_image   = o.nscope_image;
  chi2e          = o.chi2e;
  chi2R          = o.chi2R;
  mean_fov_x     = o.mean_fov_x;
  mean_fov_y     = o.mean_fov_y;
  mean_derotated_fov_x = o.mean_derotated_fov_x;
  mean_derotated_fov_y = o.mean_derotated_fov_y;
  zn             = o.zn;
  az             = o.az;
  ra             = o.ra;
  dec            = o.dec;
  ra_J2000       = o.ra_J2000;
  dec_J2000      = o.dec_J2000;
  R              = o.R;
  Rx             = o.Rx;
  Ry             = o.Ry;
  deltael        = o.deltael;
  deltaew        = o.deltaew;
  deltaRl        = o.deltaRl;
  deltaRw        = o.deltaRw;
  theta0         = o.theta0;
  theta1         = o.theta1;
  N2             = o.N2;
  msc_width      = o.msc_width;
  msc_length     = o.msc_length;
  msc_disp       = o.msc_disp;
  mlt_log10_energy = o.mlt_log10_energy;
  mlt_log10_energy_chi2 = o.mlt_log10_energy_chi2;
  theta          = o.theta;

  return *this;
}

VSEventDataWriter::
VSEventDataWriter(VSOctaveH5WriterStruct* s,
		  std::vector<unsigned> nchan, unsigned ntarget):
  m_array_writer(), m_theta_writer(ntarget), m_scope_writer(nchan.size())
{
  m_array_writer = 
    s->writeCompositeExpandableVectorHere<VSEventArrayDatum>();

  VSOctaveH5WriterCellVector* ct = s->writeCellVector("theta",ntarget);
  for(unsigned itarget=0;itarget<ntarget;itarget++)
    m_theta_writer[itarget] = ct->writeExpandableVector<double>(itarget);

  VSOctaveH5WriterCellVector* c = s->writeCellVector("t",nchan.size());
  for(unsigned iscope=0;iscope<nchan.size();iscope++)
    if(nchan[iscope])
      m_scope_writer[iscope] =
	c->writeCompositeExpandableVector<VSEventScopeDatum>(iscope);
}

VSEventDataWriter::~VSEventDataWriter()
{
  delete m_array_writer;

  unsigned ntarget = m_theta_writer.size();
  for(unsigned itarget=0;itarget<ntarget;itarget++)
    delete m_theta_writer[itarget];

  unsigned nscope = m_scope_writer.size();
  for(unsigned iscope=0;iscope<nscope;iscope++)
    delete m_scope_writer[iscope];
}

bool VSEventDataWriter::append(const VSEventArrayDatum& x)
{
  bool status = m_array_writer->append(x);

  for(unsigned itarget=0;itarget<m_theta_writer.size();itarget++)
    status &= m_theta_writer[itarget]->append(x.theta[itarget]);

  for(unsigned iscope=0;iscope<m_scope_writer.size();iscope++)
    if(m_scope_writer[iscope])
      status &= m_scope_writer[iscope]->append(*x.scope[iscope]);

  return status;
}

VSEventDataReader::
VSEventDataReader(VSOctaveH5ReaderStruct* s,
		  const MemberSubset& array_subset,
		  const MemberSubset& scope_subset):
  m_array_reader(), m_theta_reader(), m_scope_reader(), 
  m_no_has_intrinsic()
{
  m_array_reader = 
    s->readCompositeVectorElementsHere<VSEventArrayDatum>(std::string(),
							  array_subset);

  if(array_subset.empty() || array_subset.find("theta")!=array_subset.end())
    {
      VSOctaveH5ReaderCellVector* ct = s->readCellVector("theta");
      vsassert(ct);
      unsigned ntarget = ct->dimensions();

      m_theta_reader.resize(ntarget);
      for(unsigned itarget=0;itarget<ntarget;itarget++)
	m_theta_reader[itarget] = ct->readVectorElements<double>(itarget);
    }

  if(array_subset.empty() || array_subset.find("scope")!=array_subset.end())
    {
      VSOctaveH5ReaderCellVector* c = s->readCellVector("t");
      vsassert(c);
      unsigned nscope = c->dimensions();

      m_scope_reader.resize(nscope);
      for(unsigned iscope=0;iscope<nscope;iscope++)
	if(!c->isEmpty(iscope))
	  {
	    m_scope_reader[iscope] =
	      c->readCompositeVectorElements<VSEventScopeDatum>(iscope,
								scope_subset);
	    if(!m_scope_reader[iscope]->hasElement("intrinsic_width"))
	      m_no_has_intrinsic = true;
	  }
    }
}

VSEventDataReader::~VSEventDataReader()
{
  delete m_array_reader;
  unsigned ntarget = m_theta_reader.size();
  for(unsigned itarget=0;itarget<ntarget;itarget++)
    delete m_theta_reader[itarget];
  unsigned nscope = m_scope_reader.size();
  for(unsigned iscope=nscope;iscope<m_scope_reader.size();iscope++)
    delete m_scope_reader[iscope];
}

bool VSEventDataReader::element(VSEventArrayDatum& x, unsigned index)
{
  bool status = m_array_reader->element(x, index);

  unsigned ntarget = m_theta_reader.size();
  if(x.theta.size() != ntarget)x.theta.resize(ntarget);
  for(unsigned itarget=0;itarget<ntarget;itarget++)    
    status &= m_theta_reader[itarget]->element(x.theta[itarget], index);

  unsigned nscope = m_scope_reader.size();
  for(unsigned iscope=nscope;iscope<x.scope.size();iscope++)
    delete x.scope[iscope];
  if(x.scope.size() != nscope)x.scope.resize(nscope);

  for(unsigned iscope=0;iscope<nscope;iscope++)
    if(m_scope_reader[iscope])
      {
	if(x.scope[iscope]==0)x.scope[iscope] = new VSEventScopeDatum;
	status &= m_scope_reader[iscope]->element(*x.scope[iscope], index);
      }
    else if(x.scope[iscope]!=0)
      {
	delete x.scope[iscope];
	x.scope[iscope] = 0;
      }

  if(m_no_has_intrinsic)
    for(unsigned iscope=0;iscope<nscope;iscope++)
      {
	VSEventScopeDatum* xs = x.scope[iscope];
	if(xs)
	  xs->intrinsic_length = xs->fp_length,
	    xs->intrinsic_width = xs->fp_width;
      }

  return status;
}

bool VSEventDataReader::loadAllEvents(VSOctaveH5ReaderStruct* s, 
				      std::vector<VSEventArrayDatum>& x)
{
  bool status = true;

  x.clear();
 
  VSEventDataReader* reader = new VSEventDataReader(s);
  unsigned nevent = reader->rows();
  x.resize(nevent);

  for(unsigned ievent=0;ievent<nevent;ievent++)
    status &= reader->element(x[ievent],ievent);

  return status;
}
