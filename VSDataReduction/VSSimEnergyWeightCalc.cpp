//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimEnergyWeightCalc.cpp

  Class that reweights simulated energies.

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \author     Tim Arlen                   \n
              UCLA                        \n
              arlen@astro.ucla.edu        \n

  \version    1.0
  \date       03/09/2008

  $Id: VSSimEnergyWeightCalc.cpp,v 3.10 2010/10/20 17:51:19 matthew Exp $

*/

#include<VSSimEnergyWeightCalc.hpp>
#include<VSAQuadrature.hpp>

using namespace VERITAS;

/////////////////////////////////////////////////////////////////////
// To Choose which energy weighting to apply, specify with command //
// line option: -sim_energy_weight                                 //
//   ex: -sim_energy_weight=powerlaw,2.5                          //
//       -sim_energy_weight=powerlawexp_cutoff,2.5,0.3,2.0        //
//         (second would give hard exp. cutoff at 300 GeV).        //
/////////////////////////////////////////////////////////////////////

// ============================================================================
// VSSpectrumFn
// ============================================================================
VSSpectrumFn::VSSpectrumFn(unsigned nparm, const std::string& name): 
  VSAFunction::ParamFn<double>(nparm),
  m_name(name),
  m_log10_norm_egy()
{

}

VSSpectrumFn* VSSpectrumFn::create(const std::string& sp)
{
  std::vector<std::string> spv;
  VSDataConverter::fromString(spv,sp); 
  return create(spv);
}

VSSpectrumFn* VSSpectrumFn::create(const std::vector<std::string>& sp)
{
  std::vector<std::string> args = sp;
  vsassert(sp.size()>0);
  std::string sp_type = sp[0];

  args.erase(args.begin());

  if(sp_type == "powerlaw")
    return new VSSpectrumFnPowerLaw(args);
  else if (sp_type == "powerlawexp")
    return new VSSpectrumFnPowerLawExp(args);
  else if (sp_type == "powerlawexp_cutoff")
    {
      // if size < 4, program exits
      vsassert(sp.size() >= 4);
      double index;
      double cutoff;
      double alpha;
      VSDataConverter::fromString(index,sp[1]);
      VSDataConverter::fromString(cutoff,sp[2]);
      VSDataConverter::fromString(alpha,sp[3]);
      return new VSSpectrumFnPowerLawExpCutoff(index,cutoff,alpha);
    }
  else
    {
      std::cerr << "Unknown spectrum type: " << sp_type
		<< std::endl;
      exit(EXIT_FAILURE);
    }
}

VSSpectrumFn* VSSpectrumFn::create(const std::string& sp_type,
				   const VSAAlgebra::VecND& param)
{

  if(sp_type == "powerlaw")
    return new VSSpectrumFnPowerLaw(param(1),param(0));
  else if (sp_type == "powerlawexp")
    return new VSSpectrumFnPowerLawExp(param(1),param(2),param(0));
  else if (sp_type == "powerlawexp_cutoff")
    {
      exit(EXIT_FAILURE);
    }
  else
    {
      std::cerr << "Unknown spectrum type: " << sp_type
		<< std::endl;
      exit(EXIT_FAILURE);
    }
}

// ============================================================================
// VSSpectrumFnPowerLaw
// ============================================================================
VSSpectrumFnPowerLaw::
VSSpectrumFnPowerLaw(const std::vector<std::string>& args):
  VSSpectrumFn(2,"powerlaw"), m_index(2.5), m_flux_constant(3.2E-7)
{
  if(args.size() == 1)
    VSDataConverter::fromString(m_index,args[0]);      
  else if(args.size() >= 2)
    {
      VSDataConverter::fromString(m_index,args[0]);  
      VSDataConverter::fromString(m_flux_constant,args[1]);  
    }

  setParam(0,m_flux_constant);
  setParam(1,m_index);
}

VSSpectrumFnPowerLaw::VSSpectrumFnPowerLaw(double index, 
					   double flux_constant):
  VSSpectrumFn(2,"powerlaw"), m_index(index), m_flux_constant(flux_constant)
{
  setParam(0,flux_constant);
  setParam(1,index);
}

double VSSpectrumFnPowerLaw::integralFlux(double log10_elo_tev,
					  double log10_ehi_tev,
					  const VSAAlgebra::VecND& a) const
{
  const double lam = -a(1);
  
  if(fabs(lam+1) < 1E-3) 
    {
      const double x = lam+1;
      const double a1 = (log10_elo_tev - m_log10_norm_egy)*log(10);
      const double a2 = (log10_ehi_tev - m_log10_norm_egy)*log(10);

      return a(0)*((a2-a1) + 0.5*(a2*a2+a1*a1)*x + 
		   (1./6.)*(a2*a2*a2+a1*a1*a1)*x*x +
		   (1./24.)*(a2*a2*a2*a2+a1*a1*a1*a1)*x*x*x);
    }
  else return std::pow(10,-lam*m_log10_norm_egy)*
    a(0)*std::pow(1+lam,-1)*
    (std::pow(10,log10_ehi_tev*(1+lam)) - 
     std::pow(10,log10_elo_tev*(1+lam)));
}

void VSSpectrumFnPowerLaw::
integralFluxDerivative(double log10_elo_tev,
		       double log10_ehi_tev,
		       VSAAlgebra::VecND& dyda) const
{
  const double lam = -param(1);
  double elo_tev = std::pow(10,log10_elo_tev);
  double ehi_tev = std::pow(10,log10_ehi_tev);

  dyda.resize(2);

  const double y1 = 
    elo_tev*std::pow(10,(log10_elo_tev - m_log10_norm_egy)*lam);
  const double y2 = 
    ehi_tev*std::pow(10,(log10_ehi_tev - m_log10_norm_egy)*lam);
  
  if(fabs(lam+1) < 1E-3) 
    {
      const double x = lam+1;
      const double a1 = (log10_elo_tev - m_log10_norm_egy)*log(10);
      const double a2 = (log10_ehi_tev - m_log10_norm_egy)*log(10);
      
      dyda(0) = 
	(a2-a1) + 
	0.5*(a2*a2+a1*a1)*x + 
	(1./6.)*(a2*a2*a2+a1*a1*a1)*x*x +
	(1./24.)*(a2*a2*a2*a2+a1*a1*a1*a1)*x*x*x;	
  
      dyda(1) = 
	param(0)*(0.5*(a2*a2+a1*a1) +
		  2*(1./6.)*(a2*a2*a2+a1*a1*a1)*x +
		  3*(1./24.)*(a2*a2*a2*a2+a1*a1*a1*a1)*x*x);
    }
  else
    {
      dyda(0) = std::pow(1+lam,-1)*(y2 - y1);
      dyda(1) = param(0)*(y2-y1)*std::pow(1+lam,-2) -
	param(0)*std::pow(1+lam,-1)*log(10.)*
	(y2*(log10_ehi_tev - m_log10_norm_egy) - 
	 y1*(log10_elo_tev - m_log10_norm_egy));
    }
}

double VSSpectrumFnPowerLaw::val(const double& log10_egy_tev, 
				 const VSAAlgebra::VecND& a) const
{
  return a(0)*std::pow(10,-(log10_egy_tev-m_log10_norm_egy)*a(1));
}

void VSSpectrumFnPowerLaw::dyda(const double& log10_egy_tev, 
				VSAAlgebra::VecND& dyda) const
{
  VSSpectrumFnPowerLaw::dyda(log10_egy_tev,param(),dyda);
}

void VSSpectrumFnPowerLaw::dyda(const double& log10_egy_tev, 
				const VSAAlgebra::VecND& a,
				VSAAlgebra::VecND& dyda) const
{
  dyda.resize(2);

  dyda(0) = std::pow(10,-(log10_egy_tev-m_log10_norm_egy)*a(1));
  dyda(1) = -(log10_egy_tev-m_log10_norm_egy)*log(10.)*
    a(0)*std::pow(10,-(log10_egy_tev-m_log10_norm_egy)*a(1));
}

// ============================================================================
// VSSpectrumFnPowerLawExp
// ============================================================================
VSSpectrumFnPowerLawExp::
VSSpectrumFnPowerLawExp(const std::vector<std::string>& args):
  VSSpectrumFn(3,"powerlawexp"), m_index(2.5), 
  m_log10_ecut_tev(0), 
  m_flux_constant(3.2E-7)
{
  if(args.size() == 1)
    VSDataConverter::fromString(m_index,args[0]);      
  else if(args.size() == 2)
    {
      VSDataConverter::fromString(m_index,args[0]);  
      VSDataConverter::fromString(m_log10_ecut_tev,args[1]);  
    }
  else if(args.size() >= 3)
    {
      VSDataConverter::fromString(m_index,args[0]);  
      VSDataConverter::fromString(m_log10_ecut_tev,args[1]); 
      VSDataConverter::fromString(m_flux_constant,args[2]); 
    }

  setParam(0,m_flux_constant);
  setParam(1,m_index);
  setParam(2,m_log10_ecut_tev);
}

VSSpectrumFnPowerLawExp::VSSpectrumFnPowerLawExp(double index, 
						 double log10_ecut_tev, 
						 double flux_constant):
  VSSpectrumFn(3,"powerlawexp"),
  m_index(index), m_log10_ecut_tev(log10_ecut_tev), 
  m_flux_constant(flux_constant)
{
  setParam(0,flux_constant);
  setParam(1,index);
  setParam(2,log10_ecut_tev);
}

double VSSpectrumFnPowerLawExp::integralFlux(double log10_elo_tev,
					     double log10_ehi_tev,
					     const VSAAlgebra::VecND& a) const
{
  typedef VSAFunction::ParamMemberFn<VSSpectrumFnPowerLawExp,double> Fn;

  Fn fn(this,&VSSpectrumFn::dfdlog10e,a);
  return VSAMath::Quadrature::integrate(fn,log10_elo_tev,log10_ehi_tev);
}

void VSSpectrumFnPowerLawExp::
integralFluxDerivative(double log10_elo_tev,
		       double log10_ehi_tev,
		       VSAAlgebra::VecND& dyda) const
{
  vsassert(1);
}

double VSSpectrumFnPowerLawExp::val(const double& log10_egy_tev, 
				 const VSAAlgebra::VecND& a) const
{
  double egy_tev = std::pow(10,log10_egy_tev);
  double enorm = std::pow(10,m_log10_norm_egy);
  double ecut = std::pow(10,a(2));
  //  double a = exp(-std::pow(10,m_log10_norm_egy)/ecut);

  return a(0)*std::pow(egy_tev/enorm,-a(1))*exp(-egy_tev/ecut)/
    exp(-std::pow(10,m_log10_norm_egy)/ecut);
}

void VSSpectrumFnPowerLawExp::dyda(const double& log10_egy_tev, 
				   VSAAlgebra::VecND& dyda) const
{
  vsassert(1);
}

void VSSpectrumFnPowerLawExp::dyda(const double& log10_egy_tev, 
				   const VSAAlgebra::VecND& a,
				   VSAAlgebra::VecND& dyda) const
{
  vsassert(1);
}

// ============================================================================
// VSSpectrumFnPowerLawExpCutoff
// ============================================================================
VSSpectrumFnPowerLawExpCutoff::
VSSpectrumFnPowerLawExpCutoff(double index, 
			      double energy_cutoff,
			      double alpha,
			      double flux_constant):
  VSSpectrumFn(4,"powerlawexpcutoff"),
  m_index(index), m_energy_cutoff(energy_cutoff), m_alpha(alpha),
  m_flux_constant(flux_constant)
{
  /*
    Simulates a source with spectrum: 
       F(E) = flux_constant*E^(index)*exp(alpha*[1-*E/energy_cutoff])
         (if E > energy_cutoff)
    Useful for ebl studies of distant sources.
  */

  setParam(0,flux_constant);
  setParam(1,index);
  setParam(2,energy_cutoff);
  setParam(3,alpha);
}

double VSSpectrumFnPowerLawExpCutoff::val(const double& log10_egy_tev, 
					  const VSAAlgebra::VecND& a) const
{
  double egy_tev = std::pow(10,log10_egy_tev);
  if(egy_tev < a(2))
    {
      // Simple Power law:
      return a(0)*std::pow(egy_tev,a(1));
    }
  else
    {
      // Power law with a hard exponential cutoff:
      return a(0)*std::pow(egy_tev,a(1))*exp(a(3)*(1.0 - egy_tev/a(2)));
    }
}

void VSSpectrumFnPowerLawExpCutoff::dyda(const double& log10_egy_tev, 
					 VSAAlgebra::VecND& dyda) const
{
  vsassert(1);
}

void VSSpectrumFnPowerLawExpCutoff::dyda(const double& log10_egy_tev, 
					 const VSAAlgebra::VecND& a,
					 VSAAlgebra::VecND& dyda) const
{
  vsassert(1);
}

double VSSpectrumFnPowerLawExpCutoff::integralFlux(double log10_elo_tev,
						   double log10_ehi_tev,
						   const VSAAlgebra::VecND& a) 
  const
{
  vsassert(1);
  return 0;
}

void VSSpectrumFnPowerLawExpCutoff::
integralFluxDerivative(double log10_elo_tev,
		       double log10_ehi_tev,
		       VSAAlgebra::VecND& dyda) const
{
  vsassert(1);
}

// ============================================================================
// VSSimEnergyWeightCalc
// ============================================================================

VSSimEnergyWeightCalc::Options 
VSSimEnergyWeightCalc::s_default_options = VSSimEnergyWeightCalc::Options();

VSSimEnergyWeightCalc::VSSimEnergyWeightCalc(const Options& opt):
  m_spectrum(),
  m_energy_weight(),
  m_energy_norm(),
  m_apply_weighting(),
  m_options(opt)
{
  if(!m_options.sim_energy_weight.empty())
    {
      m_apply_weighting = true;
      m_spectrum = VSSpectrumFn::create(m_options.sim_energy_weight);
    }
  else
    m_spectrum = new VSSpectrumFnPowerLaw;
}

VSSimEnergyWeightCalc::~VSSimEnergyWeightCalc()
{
  delete m_spectrum;
}

void VSSimEnergyWeightCalc::calcWeighting(const VSHeaderSimulationDatum& 
					  sim_header)
{
  m_energy_weight.clear();

  double wsum = 0;
  double nevent = 0;
  for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
	sim_header.tables.begin(); itr != sim_header.tables.end(); ++itr)
    {
      double area = std::pow(itr->sampling_radius_m,2);
      double log10_egy_tev = std::log10(itr->energy_tev);
      double wegy = itr->energy_tev*m_spectrum->diffFlux(log10_egy_tev);
      
      nevent += itr->event_count;
      wsum += area*wegy;	  
      m_energy_weight[itr->energy_tev] = area*wegy/itr->event_count;
    }
  
  m_energy_norm = nevent/wsum;
}

double VSSimEnergyWeightCalc::calcWeight(double egy_tev)
{
  if(!m_apply_weighting || m_energy_weight.size() <= 1) return 1;
  
  std::map<double,double>::iterator itr = m_energy_weight.upper_bound(egy_tev);
  if(itr == m_energy_weight.end()) --itr;

  double egy2 = itr->first;
  double w2 = itr->second*m_energy_norm;
  --itr;
  double egy1 = itr->first;
  double w1 = itr->second*m_energy_norm;

  double w = w1 + (egy_tev-egy1)*(w2-w1)/(egy2-egy1);
  return w;
}

VSSpectrumFn* VSSimEnergyWeightCalc::getSpectrum()
{ 
  return m_spectrum; 
}


#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSSimEnergyWeightCalc::configure(VSOptions& options, 
				      const std::string& profile, 
				      const std::string& opt_prefix)
{
  options.findWithValue(OPTNAME(opt_prefix,"sim_energy_weight"),
			s_default_options.sim_energy_weight,
			"Set the parameters used to reweight the energy "
			"spectrum of a simulated gamma-ray file.");
}

#ifdef TEST_MAIN

// g++ -O3 -g -Wall -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DUSEALLOCA -D H5_USE_16_API -fPIC -I/usr/include/mysql -I../VSUtility -I../VSSimDB -I../Physics -I../VSShower -I../VSOptics -I../VSAnalysis -I../VSNSpace -I../VSDataReduction -I../SEphem -I. -I../VSCommon -I/usr/local/veritas/include -I$HDF5DIR/include -DTEST_MAIN -o test VSSimEnergyWeightCalc.cpp -L../VSAnalysis -L../VSCommon -lVSAnalysis -lVSCommon -L../VSUtility -lVSUtility  -lpthread -L$HDF5DIR/lib -lhdf5

int main()
{
  VSSpectrumFnPowerLawExp fn(2.5,0,1E-7);

  std::cout << fn.val(0) << " " << fn.val(-1) << std::endl;

  std::cout << fn.integralFlux(-1,0) << std::endl;
}

#endif
