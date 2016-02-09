//-*-mode:c++; mode:font-lock;-*-

/*! \file VSScaledParameterCalc.hpp

  Class which calculates scaled parameters

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/09/2008

  $Id: VSScaledParameterCalc.cpp,v 3.10 2010/06/22 00:00:32 matthew Exp $

*/

#include <fstream>

#include<VSScaledParameterCalc.hpp>
#include<VSFileUtility.hpp>
#include<VSLineTokenizer.hpp>
#include<VSScaledParameterLibrary.hpp>
#include<VSNSpaceOctaveH5IO.hpp>

//static double NEGINF = -std::numeric_limits<double>::infinity();
static double POSINF = std::numeric_limits<double>::infinity();

using namespace VERITAS;

VSScaledParameterCalc::Options 
VSScaledParameterCalc::s_default_options = VSScaledParameterCalc::Options();

VSScaledParameterCalc::VSScaledParameterCalc(const Options& opt):
  VSLTLibraryCalc(), m_options(opt), m_lt_width(), m_lt_length(), m_lt_disp()
{

}

VSScaledParameterCalc::
VSScaledParameterCalc(const std::vector<unsigned>& nchan,
		      double zn_deg, double az_deg, double ped_rms,
		      std::ostream* stream, const Options& opt):
  VSLTLibraryCalc(), m_options(opt), m_lt_width(), m_lt_length(), m_lt_disp()
{
  load(nchan,zn_deg,az_deg,ped_rms);
}

VSScaledParameterCalc::~VSScaledParameterCalc()
{
  clear();
}

bool VSScaledParameterCalc::
getScopeSP(State& state, unsigned iscope, 
	   double R, double N, 
	   double fp_width_deg, double fp_length_deg, double fp_disp_deg, 
	   double fp_dist_deg,
	   double& scope_sc_width, double& scope_sc_length, 
	   double& scope_sc_disp)
{
  if(iscope<m_lt_width.size())
    {
      const SCParameterPair& width_lt(m_lt_width[iscope]);
      const SCParameterPair& length_lt(m_lt_length[iscope]);
      const SCParameterPair& disp_lt(m_lt_disp[iscope]);

      unsigned ndim = VSScaledParameterLibraryWriter::defaultNDim();
      if(width_lt.exp)
	ndim = width_lt.exp->ndim();
      else if(length_lt.exp)
	ndim = length_lt.exp->ndim();
      else if(disp_lt.exp)
	ndim = disp_lt.exp->ndim();

      VSNSpace::Point p;
      VSScaledParameterLibraryWriter::point(p, R, N, fp_disp_deg, ndim);

      if(!width_lt.mask->space().isPointCompatible(p))
	{
	  scope_sc_width            = 0;
	  scope_sc_length           = 0;
	  scope_sc_disp             = 0;
	  return false;
	}
      else if(!(*width_lt.mask)[p] || !(*length_lt.mask)[p] || 
	      !(*disp_lt.mask)[p])
	{
	  scope_sc_width            = 0;
	  scope_sc_length           = 0;
	  scope_sc_disp             = 0;
	  return false;
	}

      bool ret = true;

      ret &= setScaledParameter(width_lt, p, fp_width_deg,
				scope_sc_width, 
				state.sum_width_w, state.sum_width);

      ret &= setScaledParameter(length_lt, p, fp_length_deg,
				scope_sc_length, 
				state.sum_length_w, state.sum_length);

      ret &= setScaledParameter(disp_lt, p, fp_disp_deg,
				scope_sc_disp, 
				state.sum_disp_w, state.sum_disp);

      return ret;
    }
  else
    {
      scope_sc_width            = 0;
      scope_sc_length           = 0;
      scope_sc_disp             = 0;
      return false;
    }
}

bool VSScaledParameterCalc::
getArraySP(State& state, double& array_msc_width, double& array_msc_length, 
	   double& array_msc_disp)
{
  bool ret = true;
  if(state.sum_width_w) 
    array_msc_width = state.sum_width/state.sum_width_w;
  else 
    array_msc_width = POSINF, ret = false;
  
  if(state.sum_length_w)
    array_msc_length = state.sum_length/state.sum_length_w;
  else 
    array_msc_length = POSINF, ret = false;
  
  if(state.sum_disp_w)
    array_msc_disp = state.sum_disp/state.sum_disp_w;
  else 
    array_msc_disp = POSINF, ret = false;

  return ret;
}

void VSScaledParameterCalc::clear()
{
  for(unsigned iscope=0; iscope<m_lt_width.size(); iscope++)
    {
      delete m_lt_width[iscope].exp;
      delete m_lt_width[iscope].rms;
      delete m_lt_width[iscope].mask;
      delete m_lt_length[iscope].exp;
      delete m_lt_length[iscope].rms;
      delete m_lt_length[iscope].mask;
      delete m_lt_disp[iscope].exp;
      delete m_lt_disp[iscope].rms;
      delete m_lt_disp[iscope].mask;
    }

  m_lt_width.clear();
  m_lt_length.clear();
  m_lt_disp.clear();
}

bool VSScaledParameterCalc::load(const std::vector<unsigned>& nchan,
				 double zn_deg, double az_deg, double ped_rms)
{
  m_has_lt = false;
  clear();
  if(m_options.sp_lookup_file.empty()) return false;
  else return load(m_options.sp_lookup_file,nchan,zn_deg,az_deg,ped_rms);
}

bool VSScaledParameterCalc::load(const VSTime& date,
				 const std::vector<unsigned>& nchan,
				 double zn_deg, double az_deg, double ped_rms)
{
  m_has_lt = false;
  clear();

  return VSLTLibraryCalc::load(m_options.sp_lookup_file,date,
			       nchan,zn_deg,az_deg,ped_rms);
}

bool VSScaledParameterCalc::load(const std::string& sp_lookup_file,
				 const std::vector<unsigned>& nchan,
				 double zn_deg, double az_deg, double ped_rms)
{
  clear();

  std::string name = sp_lookup_file;
  VSFileUtility::expandFilename(name);

  if(!VSFileUtility::isFile(name))
    {
      std::cerr << std::string(__PRETTY_FUNCTION__) + ": " << std::endl
		<< "Failed to load lookup table: "
		<< name << std::endl;
      std::cerr << "File does not exist." << std::endl;
      return false;
    }
  else if(!VSOctaveH5ReaderStruct::isHDF5(name))
    {
      std::cerr << std::string(__PRETTY_FUNCTION__) + ": " << std::endl
		<< "Failed to load lookup table: "
		<< name << std::endl;
      std::cerr << "File is not HDF5." << std::endl;
      return false;
    }

  unsigned nscope = nchan.size();
  m_lt_width.resize(nscope);
  m_lt_length.resize(nscope);
  m_lt_disp.resize(nscope);

  try
    {
      std::cout << "Loading Scaled Parameter Lookup: "
		<< name << std::endl;

      typedef VSScaledParameterLibraryWriter SPL;
      VSOctaveH5Reader* file = new VSOctaveH5Reader(name);
      VSScaledParameterLibraryReader library(file);
      
      for(unsigned iscope=0; iscope<nscope; iscope++)
	if(nchan[iscope])
	  {
	    VSNSpace* sp_exp = new VSNSpace;
	    VSNSpace* sp_rms = 0;
	    VSNSpace* sp_mask = new VSNSpace;

	    library.readAndInterpolate(*sp_exp,*sp_mask,
				       SPL::SPS_WIDTH_EXPECTED, 
				       zn_deg, az_deg, iscope,
				       ped_rms);
	    sp_rms = library.readAndInterpolate(SPL::SPS_WIDTH_RMS, 
						zn_deg, az_deg, iscope,
						ped_rms);

	    if(sp_exp && sp_rms)
	      {
		m_lt_width[iscope].exp = sp_exp;
		m_lt_width[iscope].rms = sp_rms;
		m_lt_width[iscope].mask = sp_mask;
	      }
	    else	      
	      {
		delete sp_exp;
		delete sp_rms;
		delete sp_mask;
	      }

	    sp_exp = new VSNSpace;
	    sp_mask = new VSNSpace;

	    library.readAndInterpolate(*sp_exp,*sp_mask,
				       SPL::SPS_LENGTH_EXPECTED, 
				       zn_deg, az_deg,iscope,
				       ped_rms);
	    sp_rms = library.readAndInterpolate(SPL::SPS_LENGTH_RMS, 
						zn_deg, az_deg, iscope,
						ped_rms);

	    if(sp_exp && sp_rms)
	      {
		m_lt_length[iscope].exp = sp_exp;
		m_lt_length[iscope].rms = sp_rms;
		m_lt_length[iscope].mask = sp_mask;
	      }
	    else
	      {
		delete sp_exp;
		delete sp_rms;
		delete sp_mask;
	      }

	    sp_exp = new VSNSpace;
	    sp_mask = new VSNSpace;

	    library.readAndInterpolate(*sp_exp,*sp_mask,
				       SPL::SPS_DISP_EXPECTED, 
				       zn_deg, az_deg, iscope,
				       ped_rms);
	    sp_rms = library.readAndInterpolate(SPL::SPS_DISP_RMS, 
						zn_deg, az_deg, iscope,
						ped_rms);

	    if(sp_exp && sp_rms)
	      {
		m_lt_disp[iscope].exp = sp_exp, 
		m_lt_disp[iscope].rms = sp_rms;
		m_lt_disp[iscope].mask = sp_mask;
	      }
	    else
	      {
		delete sp_exp;
		delete sp_rms;
		delete sp_mask;
	      }
	  }
	  
      delete file;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr
	<< "Caught instance of VSOctaveH5Exception loading" << std::endl
	<< "scaled paremeter lookup tables from:" << name << std::endl
	<< e.message() << std::endl;
      exit(EXIT_FAILURE);
    }

  m_has_lt = true;
  return true;
}

bool VSScaledParameterCalc::load(VSOctaveH5ReaderStruct* reader)
{
  m_has_lt = false;
  clear();

  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("t");

  if(!wc) return false;

  const unsigned nscope = wc->dimensions();

  VSNSpaceOctaveH5IO io;

  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(!wc->isStruct(iscope)) continue;

      VSOctaveH5ReaderStruct* ws = wc->readStruct(iscope);
      vsassert(ws);

      VSNSpace* sp_width_exp = new VSNSpace;
      VSNSpace* sp_width_rms = new VSNSpace;
      VSNSpace* sp_width_mask = new VSNSpace;

      VSNSpace* sp_length_exp = new VSNSpace;
      VSNSpace* sp_length_rms = new VSNSpace;
      VSNSpace* sp_length_mask = new VSNSpace;

      VSNSpace* sp_disp_exp = new VSNSpace;
      VSNSpace* sp_disp_rms = new VSNSpace;
      VSNSpace* sp_disp_mask = new VSNSpace;

      bool status = true;

      status &= 
	io.readHistogram(ws->readStruct("sp_width_exp"),*sp_width_exp);
      status &= 
	io.readHistogram(ws->readStruct("sp_width_rms"),*sp_width_rms);

      if(ws->isStruct("sp_width_mask"))
	sp_width_mask->load(ws->readStruct("sp_width_mask"));
      else
	{
	  *sp_width_mask = *sp_width_exp;
	  sp_width_mask->clear(1.0);
	}

      status &= 
	io.readHistogram(ws->readStruct("sp_length_exp"),*sp_length_exp);
      status &= 
	io.readHistogram(ws->readStruct("sp_length_rms"),*sp_length_rms);

      if(ws->isStruct("sp_length_mask"))
	sp_length_mask->load(ws->readStruct("sp_length_mask"));
      else
	{
	  *sp_length_mask = *sp_length_exp;
	  sp_length_mask->clear(1.0);
	}

      status &= 
	io.readHistogram(ws->readStruct("sp_disp_exp"),*sp_disp_exp);
      status &= 
	io.readHistogram(ws->readStruct("sp_disp_rms"),*sp_disp_rms);

      if(ws->isStruct("sp_disp_mask"))
	sp_disp_mask->load(ws->readStruct("sp_disp_mask"));
      else
	{
	  *sp_disp_mask = *sp_disp_exp;
	  sp_disp_mask->clear(1.0);
	}

      if(!status)
	{
	  delete sp_width_exp;
	  delete sp_width_rms;
	  delete sp_width_mask;
	  delete sp_length_exp;
	  delete sp_length_rms;
	  delete sp_length_mask;
	  delete sp_disp_exp;
	  delete sp_disp_rms;
	  delete sp_disp_mask;
	  delete ws;
	  delete wc;

	  clear();
	  return false;
	}

      SCParameterPair sp_width(sp_width_exp,sp_width_rms,sp_width_mask);
      SCParameterPair sp_length(sp_length_exp,sp_length_rms,sp_length_mask);
      SCParameterPair sp_disp(sp_disp_exp,sp_disp_rms,sp_disp_mask);
					 
      m_lt_width.push_back(sp_width);
      m_lt_length.push_back(sp_length);
      m_lt_disp.push_back(sp_disp);

      delete ws;
    }

  delete wc;

  m_has_lt = true;
  return true;
}

void VSScaledParameterCalc::save(VSOctaveH5WriterStruct* writer) const
{
  const unsigned nscope = m_lt_width.size();

  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("t",nscope);

  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(!m_lt_width[iscope].exp || !m_lt_width[iscope].rms ||
	 !m_lt_length[iscope].exp || !m_lt_length[iscope].rms ||
	 !m_lt_disp[iscope].exp || !m_lt_disp[iscope].rms)
	continue;

      VSOctaveH5WriterStruct* ws = wc->writeStruct(iscope);

      m_lt_width[iscope].exp->save(ws->writeStruct("sp_width_exp"));
      m_lt_width[iscope].rms->save(ws->writeStruct("sp_width_rms"));
      m_lt_width[iscope].mask->save(ws->writeStruct("sp_width_mask"));

      m_lt_length[iscope].exp->save(ws->writeStruct("sp_length_exp"));
      m_lt_length[iscope].rms->save(ws->writeStruct("sp_length_rms"));
      m_lt_length[iscope].mask->save(ws->writeStruct("sp_length_mask"));

      m_lt_disp[iscope].exp->save(ws->writeStruct("sp_disp_exp"));
      m_lt_disp[iscope].rms->save(ws->writeStruct("sp_disp_rms"));
      m_lt_disp[iscope].mask->save(ws->writeStruct("sp_disp_mask"));

      delete ws;
    }

  delete wc;
}

VSScaledParameterCalc::VSScaledParameterCalc(const VSScaledParameterCalc& o)
{
  m_options = o.m_options;
    
  m_has_lt = o.m_has_lt;
  
  const unsigned nscope = o.m_lt_width.size();
  m_lt_width.resize(nscope);
  m_lt_length.resize(nscope);
  m_lt_disp.resize(nscope);
  
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(o.m_lt_width[iscope].exp && o.m_lt_width[iscope].rms)
	{
	  m_lt_width[iscope].exp = 
	    new VSNSpace(*o.m_lt_width[iscope].exp);
	  m_lt_width[iscope].rms = 
	    new VSNSpace(*o.m_lt_width[iscope].rms);
	  m_lt_width[iscope].mask = 
	    new VSNSpace(*o.m_lt_width[iscope].mask);
	}
      else
	{
	  m_lt_width[iscope].exp = NULL;
	  m_lt_width[iscope].rms = NULL;
	  m_lt_width[iscope].mask = NULL;
	}

      if(o.m_lt_length[iscope].exp && o.m_lt_length[iscope].rms)
	{
	  m_lt_length[iscope].exp = 
	    new VSNSpace(*o.m_lt_length[iscope].exp);
	  m_lt_length[iscope].rms = 
	    new VSNSpace(*o.m_lt_length[iscope].rms);
	  m_lt_length[iscope].mask = 
	    new VSNSpace(*o.m_lt_length[iscope].mask);
	}
      else
	{
	  m_lt_length[iscope].exp = NULL;
	  m_lt_length[iscope].rms = NULL;
	  m_lt_length[iscope].mask = NULL;
	}

      if(o.m_lt_disp[iscope].exp && o.m_lt_disp[iscope].rms)
	{
	  m_lt_disp[iscope].exp = 
	    new VSNSpace(*o.m_lt_disp[iscope].exp);
	  m_lt_disp[iscope].rms = 
	    new VSNSpace(*o.m_lt_disp[iscope].rms);
	  m_lt_disp[iscope].mask = 
	    new VSNSpace(*o.m_lt_disp[iscope].mask);
	}
      else
	{
	  m_lt_disp[iscope].exp = NULL;
	  m_lt_disp[iscope].rms = NULL;
	  m_lt_disp[iscope].mask = NULL;
	}
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSScaledParameterCalc::
configure(VSOptions& options, const std::string& profile, 
	  const std::string& opt_prefix)
{
  options.findWithValue(OPTNAME(opt_prefix,"sc_parameter_lookup"),
			s_default_options.sp_lookup_file,
			"Set the file name from which to load the scaled "
			"parameter look-up tables, from which to calculate "
			"mean scaled width, length and disp.",
			OPTNAME(opt_prefix,"shape"));

  options.findWithValue(OPTNAME(opt_prefix,"msc_weight_power"),
			s_default_options.sp_weight_power,
			"Set the weighting for combining telescope scaled "
			"parameter values into the mean scaled value. The "
			"weighing is given by the 1/sigma to the power of "
			"this value. A value of zero weights all telescopes "
			"equally, a value of \"inf\" means only the value "
			"with the smallest sigma is used. Suggested values "
			"are 0,1,2 or inf.",
			OPTNAME(opt_prefix,"shape"));
}
