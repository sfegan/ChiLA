//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEnergyCalcLT.cpp

  Class which calculates reconstructed energy.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \author     Matthew Wood                \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/09/2008

  $Id: VSEnergyCalcLT.cpp,v 3.12 2010/06/22 00:00:15 matthew Exp $

*/

#include<VSEnergyCalcLT.hpp>
#include<VSFileUtility.hpp>
#include<VSScaledParameterLibrary.hpp>
#include<VSNSpaceOctaveH5IO.hpp>
#include<VSLTLibraryCalc.hpp>

static double NEGINF = -std::numeric_limits<double>::infinity();

using namespace VERITAS;

VSEnergyCalcLT::Options 
VSEnergyCalcLT::s_default_options = VSEnergyCalcLT::Options();

VSEnergyCalcLT::VSEnergyCalcLT(const Options& opt):
  VSLTLibraryCalc(), m_options(opt), m_lt()
{

}

VSEnergyCalcLT::VSEnergyCalcLT(const std::vector<unsigned>& nchan,
			       double zn_deg, double az_deg, double ped_rms,
			       std::ostream* stream, const Options& opt):
  VSLTLibraryCalc(), m_options(opt), m_lt(), m_egy_estimator()
{
  load(nchan,zn_deg,az_deg,ped_rms);
}

VSEnergyCalcLT::~VSEnergyCalcLT()
{
  clear();
}

bool VSEnergyCalcLT::calcEnergy(VSEventArrayDatum& event_data)
{
  State state;
  unsigned nscope = event_data.scope.size();
  for(unsigned iscope=0;iscope<nscope; iscope++)
    {
      VSEventScopeDatum* scope_data = event_data.scope[iscope];
      
      bool used_in_reconstruction = 
	((0x1<<iscope)&event_data.used_in_reconstruction_mask);
      
		       // scope_data->R, scope_data->N, scope_data->fp_disp,
		       // scope_data->fp_dist,

      if(scope_data && used_in_reconstruction)
	getScopeEnergy(state, iscope, scope_data,
		       scope_data->lt_log10_energy, 
		       scope_data->lt_log10_energy_err);
    }
  return getArrayEnergy(state, event_data.mlt_log10_energy,
			event_data.mlt_log10_energy_chi2);
}

bool VSEnergyCalcLT::
getScopeEnergy(State& state, 
	       unsigned iscope, double R, double N, double fp_disp_deg, 
	       double fp_dist_deg,
	       double& scope_log10_e, double& scope_log10_e_rms)
{

  return true;
}

bool VSEnergyCalcLT::
getScopeEnergy(State& state, 
	       unsigned iscope, VSEventScopeDatum* sd,
	       double& scope_log10_e, double& scope_log10_e_rms)
{
  if(iscope<m_lt.size()&&m_lt[iscope].exp&&m_lt[iscope].rms)
    {
      const SCParameterPair& lt(m_lt[iscope]);

      unsigned ndim = lt.exp->ndim();
      VSNSpace::Point p(ndim);
      point(p,sd,m_egy_estimator);
		
      scope_log10_e     = NEGINF;
      scope_log10_e_rms = 0;

      if(!lt.mask->space().isPointCompatible(p)) return false;
      else if(!(*lt.mask)[p]) return false;

      if((lt.exp->interpolateWeight(p, scope_log10_e))
	 &&(lt.rms->interpolateWeight(p, scope_log10_e_rms))
	 &&(scope_log10_e_rms>0.0)
	 &&((m_options.energy_rms_cut<=0.0)
	    ||(scope_log10_e_rms<m_options.energy_rms_cut))
	 && sd->fp_dist <= m_options.energy_quality_cuts.first &&
	 sd->R <= m_options.energy_quality_cuts.second)
	{
	  double x_si2 = 1.0/(scope_log10_e_rms*scope_log10_e_rms);

	  state.sum_1_sigma2  += x_si2;  x_si2 *= scope_log10_e;
	  state.sum_e_sigma2  += x_si2;  x_si2 *= scope_log10_e;
	  state.sum_e2_sigma2 += x_si2;
	  
	  if(std::isfinite(m_options.energy_weight_power))
	    {
	      double w = 1.0;
	      if(m_options.energy_weight_power > 0)
		w = pow(scope_log10_e_rms,
			-m_options.energy_weight_power);
	      state.sum_w += w;
	      state.sum_e += w * scope_log10_e;
	    }
	  else if((state.sum_w == 0)||(scope_log10_e_rms < state.sum_w))
	    {
	      double w = scope_log10_e_rms;
	      state.sum_w = w;
	      state.sum_e = w * scope_log10_e;
	    }
	}
      else return false;
    }
  else
    {
      scope_log10_e     = 0;
      scope_log10_e_rms = 0;
      return false;
    }

  return true;
}

bool VSEnergyCalcLT::
getArrayEnergy(State& state, 
	       double& array_log10_e, double& array_log10_e_chi2)
{
  if(state.sum_w)
    {
      array_log10_e = state.sum_e/state.sum_w;
      array_log10_e_chi2 = 
	state.sum_e2_sigma2
	- 2.0*array_log10_e*state.sum_e_sigma2
	+ array_log10_e*array_log10_e*state.sum_1_sigma2;
      // Handle case of negative chi2 due to FP rounding here
      if(array_log10_e_chi2<0)array_log10_e_chi2=0;
    }
  else 
    {
      array_log10_e = NEGINF;
      array_log10_e_chi2 = NEGINF;
      return false;
    }
  return true;
}

void VSEnergyCalcLT::point(VSNSpace::Point& p, VSEventScopeDatum* sd,
			   EnergyEstimator ee)
{
  
  if(ee == EE_N_R)
    {
      p.resize(2);
      p.x[0] = sd->R;
      p.x[1] = std::log10(sd->N);
    }
  else if(ee == EE_N_R_DISP)
    {
      p.resize(3);
      p.x[0] = sd->R;
      p.x[1] = std::log10(sd->N);
      p.x[2] = sd->fp_disp;
    }
  else if(ee == EE_LAMBDAD_R)
    {
      p.resize(2);
      p.x[0] = sd->R;
      p.x[1] = std::log10(sd->lambdad);
    }
  else if(ee == EE_LAMBDAD_G)
    {
      p.resize(2);
      p.x[0] = sd->G;
      p.x[1] = std::log10(sd->lambdad);
    }
  else if(ee == EE_LAMBDAD_ETA)
    {
      double eta = -std::log(tan(sd->theta1/2));

      p.resize(2);
      p.x[0] = eta;
      p.x[1] = std::log10(sd->lambdad);
    }
  else if(ee == EE_LAMBDAD_G_ETA)
    {
      double eta = -std::log(tan(sd->theta1/2));

      p.resize(3);
      p.x[0] = sd->G;
      p.x[1] = std::log10(sd->lambdad);
      p.x[2] = eta;
    }
  else
    {
      std::cerr << "Unrecognized energy estimator: " << ee << std::endl;
      exit(EXIT_FAILURE);
    }
}

void VSEnergyCalcLT::clear()
{
  for(unsigned iscope=0; iscope<m_lt.size(); iscope++)
    {
      delete m_lt[iscope].exp;
      delete m_lt[iscope].rms;
      delete m_lt[iscope].mask;
    }
  
  m_lt.clear();
  m_egy_estimator = EE_N_R;
}

bool VSEnergyCalcLT::load(const std::vector<unsigned>& nchan,
			  double zn_deg, double az_deg, double ped_rms)
{
  m_has_lt = false;
  clear();

  if(m_options.energy_lookup_file.empty()) return false;
  else return load(m_options.energy_lookup_file,nchan,zn_deg,az_deg,ped_rms);
}

bool VSEnergyCalcLT::load(const VSTime& date,
			  const std::vector<unsigned>& nchan,
			  double zn_deg, double az_deg, double ped_rms)
{
  m_has_lt = false;
  clear();

  return VSLTLibraryCalc::load(m_options.energy_lookup_file,date,
			       nchan,zn_deg,az_deg,ped_rms);
}

bool VSEnergyCalcLT::load(const std::string& egy_lookup_file,
			  const std::vector<unsigned>& nchan,
			  double zn_deg, double az_deg, double ped_rms)
{
  m_has_lt = false;
  clear();

  std::string name = egy_lookup_file;
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
  m_lt.resize(nscope);

  try
    {
      std::cout << "Loading Energy Estimator Lookup: "
		<< name << std::endl;

      typedef VSScaledParameterLibraryWriter SPL;
      VSOctaveH5Reader* file = new VSOctaveH5Reader(name);
      VSScaledParameterLibraryReader library(file);
      
      for(unsigned iscope=0; iscope<nscope; iscope++)
	if(nchan[iscope])
	  {
	    VSNSpace* lt_exp = new VSNSpace;
	    VSNSpace* lt_rms = 0;
	    VSNSpace* lt_mask = new VSNSpace;	    

	    library.
	      readAndInterpolate(*lt_exp,*lt_mask,
				 SPL::SPS_ENERGY_EXPECTED, 
				 zn_deg, az_deg, iscope, ped_rms);
	    lt_rms = library.
	      readAndInterpolate(SPL::SPS_ENERGY_RMS, 
				 zn_deg, az_deg, iscope, ped_rms);
	   

	    if(lt_exp->ndim() == 2) 
	      {
		if(lt_exp->axis(0).name == "R" && lt_exp->axis(1).name == "N")
		  m_egy_estimator = EE_N_R;
		else if(lt_exp->axis(0).name == "R" && 
			lt_exp->axis(1).name == "lambdad")
		  m_egy_estimator = EE_LAMBDAD_R;
		else if(lt_exp->axis(0).name == "G" && 
			lt_exp->axis(1).name == "lambdad")
		  m_egy_estimator = EE_LAMBDAD_G;
		else if(lt_exp->axis(0).name == "eta" && 
			lt_exp->axis(1).name == "lambdad")
		   m_egy_estimator = EE_LAMBDAD_ETA;
		else
		  m_egy_estimator = EE_N_R;
	      }
	    else if(lt_exp->ndim() == 3) 
	      {
		if(lt_exp->axis(0).name == "R" && 
		   lt_exp->axis(1).name == "N" &&
		   lt_exp->axis(2).name == "fp_disp")
		  m_egy_estimator = EE_N_R_DISP;
		else if(lt_exp->axis(0).name == "G" && 
			lt_exp->axis(1).name == "lambdad" &&
			lt_exp->axis(2).name == "eta")
		  m_egy_estimator =  EE_LAMBDAD_G_ETA;
		else
		  m_egy_estimator = EE_N_R_DISP;
	      }
	    else
	      {
		std::cerr << "Invalid energy lookup table dimension: "
			  << lt_exp->ndim() << std::endl;
		exit(EXIT_FAILURE);
	      }

	    if(lt_exp && lt_rms)
	      {
		m_lt[iscope].exp = lt_exp; 
		m_lt[iscope].rms = lt_rms;
		m_lt[iscope].mask = lt_mask;
	      }
	    else
	      {
		m_egy_estimator = EE_N_R;
		delete lt_exp; 
		delete lt_rms;
		delete lt_mask;
	      }
	  }
      
      delete file;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr
	<< "Caught instance of VSOctaveH5Exception loading" << std::endl
	<< "energy lookup tables from: " << name << std::endl
	<< e.message() << std::endl;
      exit(EXIT_FAILURE);
    }

  m_has_lt = true;
  return true;
}

bool VSEnergyCalcLT::load(VSOctaveH5ReaderStruct* reader)
{
  clear();
  m_has_lt = false;

  VSOctaveH5ReaderCellVector* wc = reader->readCellVector("t");

  if(!wc) return false;

  const unsigned nscope = wc->dimensions();

  VSNSpaceOctaveH5IO io;

  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      VSOctaveH5ReaderStruct* ws = wc->readStruct(iscope);
      vsassert(ws);

      VSNSpace* egy_lookup_exp = new VSNSpace;
      VSNSpace* egy_lookup_rms = new VSNSpace;
      VSNSpace* egy_lookup_mask = new VSNSpace;

      bool status = true;

      status &= 
	io.readHistogram(ws->readStruct("egy_lookup_exp"),*egy_lookup_exp);
      status &= 
	io.readHistogram(ws->readStruct("egy_lookup_rms"),*egy_lookup_rms);

      if(ws->isStruct("egy_lookup_mask"))
	egy_lookup_mask->load(ws->readStruct("egy_lookup_mask"));
      else
	{
	  *egy_lookup_mask = *egy_lookup_exp;
	  egy_lookup_mask->clear(1.0);
	}

      if(!status)
	{
	  delete egy_lookup_exp;
	  delete egy_lookup_rms;
	  delete egy_lookup_mask;
	  delete ws;
	  delete wc;
	  return false;
	}

      SCParameterPair p(egy_lookup_exp,egy_lookup_rms,egy_lookup_mask);
      m_lt.push_back(p);

      delete ws;
    }

  delete wc;

  m_has_lt = true;
  return true;
}

void VSEnergyCalcLT::save(VSOctaveH5WriterStruct* writer) const
{
  const unsigned nscope = m_lt.size();

  VSOctaveH5WriterCellVector* wc = writer->writeCellVector("t",nscope);
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      VSOctaveH5WriterStruct* ws = wc->writeStruct(iscope);
      m_lt[iscope].exp->save(ws->writeStruct("egy_lookup_exp"));
      m_lt[iscope].rms->save(ws->writeStruct("egy_lookup_rms"));
      m_lt[iscope].mask->save(ws->writeStruct("egy_lookup_mask"));
      delete ws;
    }

  delete wc;
}

VSEnergyCalcLT::VSEnergyCalcLT(const VSEnergyCalcLT& o)
{
  m_options = o.m_options;
  m_lt.resize(o.m_lt.size());
  const unsigned nscope = o.m_lt.size();
  for(unsigned iscope = 0; iscope < nscope; iscope++)
    {
      if(o.m_lt[iscope].exp && o.m_lt[iscope].rms) 
	{
	  m_lt[iscope].exp = new VSNSpace(*o.m_lt[iscope].exp);
	  m_lt[iscope].rms = new VSNSpace(*o.m_lt[iscope].rms);
	  m_lt[iscope].mask = new VSNSpace(*o.m_lt[iscope].mask);
	}
      else 
	{
	  m_lt[iscope].exp = NULL;
	  m_lt[iscope].rms = NULL;
	  m_lt[iscope].mask = NULL;
	}
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSEnergyCalcLT::
configure(VSOptions& options, const std::string& profile, 
	  const std::string& opt_prefix)
{
  options.findWithValue(OPTNAME(opt_prefix,"energy_lookup"),
			s_default_options.energy_lookup_file,
			"Set the file name from which to load the energy "
			"look-up tables.",
			OPTNAME(opt_prefix,"energy"));

  options.findWithValue(OPTNAME(opt_prefix,"energy_weight_power"),
			s_default_options.energy_weight_power,
			"Set the weighting for combining individual telescope "
			"energy etimates into the array energy estimate. The "
			"weighing is given by the 1/sigma to the power of "
			"this value. A value of zero weights all telescopes "
			"equally, a value of \"inf\" means only the value "
			"with the smallest sigma is used. Suggested values "
			"are 0,1,2 or inf.",
			OPTNAME(opt_prefix,"energy"));

  options.findWithValue(OPTNAME(opt_prefix,"energy_rms_cut"),
			s_default_options.energy_rms_cut,
			"Set threshold on RMS of per telescope (log10) energy "
			"estimate above which the telescope event will not "
			"be used in array energy estimate. A value of zero "
			"means the cut will not be applied. For example, "
			"a value of 0.097 ~= log10(1.25) specifies a cut "
			"at an RMS error of 25%.",
			OPTNAME(opt_prefix,"energy"));

  options.findWithValue(OPTNAME(opt_prefix,"energy_quality_cuts"),
			s_default_options.energy_quality_cuts,
			"Set the per-telescope quality cuts: maximum distance "
                        "from the center of the camera to the centroid in "
			"degrees and maximum distance of the shower core from "
			"the telescope in meters. "
			"Cuts should be specified as: dist/R",
			OPTNAME(opt_prefix,"energy"));
}
