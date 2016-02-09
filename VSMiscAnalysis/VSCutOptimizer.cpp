//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCutOptimizer.cpp
  

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/17/2005
*/

#include <cmath>

#include <Astro.h>

#include "VSCutOptimizer.hpp"
#include "VSSimpleCutOptimizer.hpp"
#include "VSNSpaceCutOptimizer.hpp"
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSDataModel.hpp>

using namespace VERITAS;
using namespace SEphem;

// ============================================================================
// VSCutOptimizer::FunctionalQSimple
// ============================================================================
double VSCutOptimizer::FunctionalQSimple::operator() 
  (const std::vector<double>& x) const
{
  return x[0]/sqrt(x[1]);
}

// ============================================================================
// VSCutOptimizer::FunctionalLambda
// ============================================================================
double VSCutOptimizer::FunctionalLambda::operator() 
  (const std::vector<double>& x) const
{
  return x[0]/sqrt(x[1]);
}

// ============================================================================
// VSCutOptimizer
// ============================================================================

VSCutOptimizer::Options::Options():
  optimizer("simple","Q_simple"),
  space(),
  theta_cut(0.2),
  use_theta(false),
  min_events(5), 
  size_cut_nscope(2),
  min_rate(),
  min_fraction(),
  bkgnd_selection(),
  spectrum("powerlaw,2.5,3.2E-7")
{ 
  bkgnd_selection.push_back("theta_off");
  bkgnd_selection.push_back("0.2");

  space.push_back(AxisDefinition("msc_width",0.0,1.4,0.05));
  space.push_back(AxisDefinition("msc_length",0.5,2.0,0.1));
  space.push_back(AxisDefinition("neglog:events.N2",-3.5,-2.0,0.02));
}

inline double mymax(double x, double y) 
{
  if(std::isnan(y))return y;
  return std::max(x,y);
}

VSCutOptimizer::Data::Data():
  excess(),
  excess_var(),
  bkgnd(),
  bkgnd_var(),
  bkgnd_counts(),
  bkgnd_events(),
  signal_events(),
  on_events(),
  off_events()
{

}

VSCutOptimizer::Data::Data(const VSNSpace::Space& space):
  excess(space),
  excess_var(space),
  excess_rate(space),
  bkgnd(space),
  bkgnd_var(space),
  bkgnd_rate(space),
  bkgnd_counts(space),
  bkgnd_events(),
  signal_events(),
  on_events(),
  off_events()
{

}

VSCutOptimizer::VSCutOptimizer(const Options& opt):
  VSEventDataVisitor(),
  m_options(opt),
  m_bkgnd_selection_mode(),
  m_bkgnd_selection_pars(),
  m_bkgnd_theta_cut(0.2),
  m_bkgnd_alpha(1.0),
  m_sp_calc(),
  m_egywt_calc(),
  m_cuts_calc(),
  m_opt_func(),
  m_is_sim_run(),
  m_sim_header(),
  m_sim(),
  m_evt(), 
  m_obs_radec_J2000(),
  m_src_radec_J2000(),
  m_obs_xy(),
  m_src_xy(),
  m_off_xy(),
  m_livetime(), m_elaptime(), m_wt(1.0),
  m_space()
{
  vsassert(m_options.bkgnd_selection.size());
  m_bkgnd_selection_mode = m_options.bkgnd_selection[0];
  m_bkgnd_selection_pars.resize(m_options.bkgnd_selection.size()-1);

  const unsigned npar = m_options.bkgnd_selection.size();
  for(unsigned ipar = 1; ipar < npar; ipar++)
    VSDataConverter::fromString(m_bkgnd_selection_pars[ipar-1],
				m_options.bkgnd_selection[ipar]);

  if(m_bkgnd_selection_mode == "on") 
    vsassert(m_bkgnd_selection_pars.size() == 1);
  else if(m_bkgnd_selection_mode == "theta_off") 
    {
      vsassert(m_bkgnd_selection_pars.size() == 1);
      m_bkgnd_theta_cut = m_bkgnd_selection_pars[0];
    }
  else if(m_bkgnd_selection_mode == "ring") 
    vsassert(m_bkgnd_selection_pars.size() == 2);
  else
    {
      std::cerr << "Unknown background selection mode: "
		<< m_bkgnd_selection_mode << std::endl;
      exit(EXIT_FAILURE);
    }

  m_cuts_calc = new VSCutsEvaluator;
  m_egywt_calc = new VSSimEnergyWeightCalc;
  m_sim_spectrum = VSSpectrumFn::create(m_options.spectrum);

  if(m_options.optimizer.second == "Q_simple")
    m_opt_func = new FunctionalQSimple;
  else if(m_options.optimizer.second == "lambda")
    m_opt_func = new FunctionalLambda;
  else
    {
      std::cerr << "Unknown optimization functional: "
		<< m_options.optimizer.second << std::endl;
      exit(EXIT_FAILURE);
    }

  // DISP 0.5 2.5 0.1

  const double min_th2 = 0.00;
  const double max_th2 = 0.05;
  const double dx_th2 = 0.0005;

  if(m_options.use_theta)
    m_options.space.
      insert(m_options.space.begin(),
	     AxisDefinition("sq:theta1",min_th2, max_th2, dx_th2));
  
  const unsigned ndim = m_options.space.size();
  m_space.resize(ndim);
  for(unsigned idim = 0; idim < ndim; idim++)
    {
      std::cout << idim << " " 
		<< std::setw(20) << m_options.space[idim].first
		<< std::setw(15) << m_options.space[idim].second
		<< std::setw(15) << m_options.space[idim].third
		<< std::setw(15) << m_options.space[idim].fourth
		<< std::endl;

      m_space.axes[idim] = VSNSpace::Axis(m_options.space[idim].second, 
					  m_options.space[idim].third, 
					  m_options.space[idim].fourth, 
					  0, m_options.space[idim].first);

      VSH5DatumElement<VSEventArrayDatum>* element = 
	VSH5DatumElement<VSEventArrayDatum>::
	createDatumElement(m_options.space[idim].first);
      m_datum_elements.push_back(element);
    }

  m_train_data = Data(m_space);
  m_test_data = Data(m_space);

  m_sim_rate = VSNSpace(m_space);
}

VSCutOptimizer::~VSCutOptimizer()
{
  delete m_egywt_calc;
  delete m_cuts_calc;
  delete m_sim_spectrum;
}

void VSCutOptimizer::visitRun(const VSAnalysisStage1Data& stage1,
			      const VSTargetTable::Observation& obs,
			      const VSArrayMergedCalibrationData& cal)
{
  // --------------------------------------------------------------------------
  // Load lookup tables
  // --------------------------------------------------------------------------
  m_sp_calc = new VSScaledParameterCalc;
  m_sp_calc->load(stage1.run_info.nchan,stage1.run_info.zn_mean_deg,
		  stage1.run_info.az_mean_deg,cal.mean_scaled_dev);

  m_src_radec_J2000 = obs.src_radec_J2000;	


//   std::pair< Angle, Angle > src_xy;
//   Astro::raDecToXY(m_src_radec_J2000,origin_radec_J2000,src_xy);

  if(m_is_sim_run)
    {
      m_obs_radec_J2000 = 
	SphericalCoords::makeLatLongDeg(stage1.sim_info->obs_dec_deg,
					stage1.sim_info->obs_ra_deg);
    }
  else
    {
      m_obs_radec_J2000 = obs.obs_radec_J2000;	

    }

  // Determine off regions
  std::pair<Angle,Angle> src_xy;
  std::pair<Angle,Angle> obs_xy;

  Astro::raDecToXY(m_src_radec_J2000,m_src_radec_J2000,src_xy);
  Astro::raDecToXY(m_obs_radec_J2000,m_src_radec_J2000,obs_xy);  

//   std::cout << "SRC RADEC "
// 	    << std::setw(15) << m_src_radec_J2000.thetaDeg()
// 	    << std::setw(15) << m_src_radec_J2000.phiDeg()
// 	    << std::endl;

//   std::cout << "OBS RADEC "
// 	    << std::setw(15) << m_obs_radec_J2000.thetaDeg()
// 	    << std::setw(15) << m_obs_radec_J2000.phiDeg()
// 	    << std::endl;

  m_src_xy = VSAAlgebra::Vec2D(src_xy.first.degPM(),src_xy.second.degPM());
  m_obs_xy = VSAAlgebra::Vec2D(obs_xy.first.degPM(),obs_xy.second.degPM());

  m_off_xy.clear();

  VSOffRegion::getRegions(m_src_xy,m_obs_xy,m_options.theta_cut,0,m_off_xy);

//   std::cout << "SRC XY " 
// 	    << std::setw(15) << m_src_xy.x()
// 	    << std::setw(15) << m_src_xy.y()
// 	    << std::endl;

//   std::cout << "OBS XY " 
// 	    << std::setw(15) << m_obs_xy.x()
// 	    << std::setw(15) << m_obs_xy.y()
// 	    << std::endl;

//   for(unsigned ioff = 0; ioff < m_off_xy.size(); ioff++)
//     {
//       std::cout << std::setw(6) << ioff 
// 		<< std::setw(15) << m_off_xy[ioff].x() 
// 		<< std::setw(15) << m_off_xy[ioff].y() 
// 		<< std::setw(15) << m_off_xy[ioff].d(m_obs_xy)
// 		<< std::setw(15) << m_off_xy[ioff].d(m_src_xy)
// 		<< std::endl;
//     }

  m_wt = 1;
  m_alpha = 1./m_off_xy.size();
  m_bkgnd_alpha = 1./m_off_xy.size();
}

void VSCutOptimizer::leaveRun()
{
  if(!m_is_sim_run)
    {
      m_livetime += m_evt.elapsed_ticks/1.0E7;
      m_elaptime += m_evt.live_ticks/1.0E7;
    }

//   // Sim Signal ---------------------------------------------------------------
//   m_excess += m_sim_on;
//   m_excess_var += m_sim_on;

//   // Sim Bkgnd ----------------------------------------------------------------
//   m_bkgnd += m_sim_off;
//   m_bkgnd_var += m_sim_off;
//   m_bkgnd_counts += m_sim_off;

  // Data ---------------------------------------------------------------------
  //  double alpha_bkgnd = 1;

//   if(m_bkgnd_selection_mode == "theta_off")
//     alpha = 1./double(m_theta_off.size());
//   std::cout << "ALPHA " << alpha << std::endl;

//   m_bkgnd_counts += m_off;

//   VSNSpace bkgnd = m_off;
//   bkgnd *= alpha;
//   m_bkgnd += bkgnd;

//   VSNSpace bkgnd_var = m_off;
//   bkgnd_var *= alpha*alpha;
//   m_bkgnd_var += bkgnd_var;

//   if(!m_options.sim_signal)
//     {
//       double alpha = 1./double(m_off_xy.size());

//       VSNSpace excess = m_on;
//       excess -= m_off*alpha;
//       m_excess += excess;

//       VSNSpace excess_var = m_on;
//       excess_var += m_off*(alpha*alpha);
//       m_excess_var += excess_var;
//     }      

  delete m_sp_calc;
  m_sp_calc = NULL;
}

void VSCutOptimizer::visitEvent(const VSEventArrayDatum& evt)
{
  m_evt = evt;
  m_sp_calc->calcSP(m_evt);

  if(!m_cuts_calc->isSelected(m_evt)) return;

  SEphem::SphericalCoords event_radec = 
    SEphem::SphericalCoords::makeLatLongDeg(m_evt.dec_J2000, m_evt.ra_J2000);


  std::pair<Angle,Angle> xy;
  Astro::raDecToXY(event_radec,m_src_radec_J2000,xy);
  VSAAlgebra::Vec2D event_xy(xy.first.degPM(),xy.second.degPM());

  VSNSpace::Point p = m_space.point();

  const unsigned nelement = m_datum_elements.size();
  for(unsigned iel = 0; iel < nelement; iel++)
    {
      double val;
      bool ret = m_datum_elements[iel]->getValue(m_evt,val);
      if(!ret) return;

      if(m_options.optimizer.first == "simple")
	p.x[iel] = mymax(m_space.axes[iel].lo_bound + 
			 m_space.axes[iel].bin_size/2., val);
      else
	p.x[iel] = val;
    }

  // --------------------------------------------------------------------------
  // Sim gamma-ray event
  // --------------------------------------------------------------------------
  if(m_is_sim_run && m_sim.corsika_particle_id == 1)   
    {
      if(m_evt.theta1 < m_options.theta_cut)
	m_sim_events.push_back(Event(p,m_wt,m_sim.table_index));
    }
  // --------------------------------------------------------------------------
  // Sim background event
  // --------------------------------------------------------------------------
  else if(m_is_sim_run && m_sim.corsika_particle_id > 1) 
    {
      if(m_bkgnd_selection_mode == "on")
	{
	  double theta = m_evt.theta0;

	  if(m_options.use_theta) p.x[0] = 0;
	  if(theta < m_bkgnd_selection_pars[0])
	    {
	      m_bkgnd_events.push_back(Event(p,m_bkgnd_alpha));
	      //	      m_sim_off.accumulate(p);     
	    }
	}
    }
  // --------------------------------------------------------------------------
  // Data event
  // --------------------------------------------------------------------------
  else
    {
      if(!m_options.sim_signal)
	{
	  if(m_evt.theta1 < m_options.theta_cut)
	    {
	      m_on_events.push_back(Event(p));
	    }

	  // Accumulate off counts --------------------------------------------
	  for(std::vector< VSAAlgebra::Vec2D >::const_iterator itr =
		m_off_xy.begin(); itr != m_off_xy.end(); ++itr)
	    {
	      double dr_off = itr->d(event_xy);
	      if(m_options.use_theta) p.x[0] = std::pow(dr_off,2);

	      if(dr_off < m_options.theta_cut) 
		m_off_events.push_back(Event(p,m_alpha));
	    }
	}

      if(m_options.use_theta) p.x[0] = 0;

      if(m_bkgnd_selection_mode == "on")
	{
	  double theta = m_evt.theta0;
	  if(theta < m_bkgnd_selection_pars[0] && m_evt.theta1 > 0.25) 
	    m_bkgnd_events.push_back(Event(p,m_bkgnd_alpha));
	}
      else if(m_bkgnd_selection_mode == "theta_off")
	{
	  // Accumulate off counts --------------------------------------------
	  for(std::vector< VSAAlgebra::Vec2D >::const_iterator itr =
		m_off_xy.begin(); itr != m_off_xy.end(); ++itr)
	    {
	      double dr_off = itr->d(event_xy);

	      if(dr_off < m_bkgnd_theta_cut)
		m_bkgnd_events.push_back(Event(p,m_bkgnd_alpha));
	    }
	}
      else if(m_bkgnd_selection_mode == "ring")
	{
	  if(m_evt.theta0 > m_bkgnd_selection_pars[0] &&
	     m_evt.theta0 < m_bkgnd_selection_pars[1])
	    m_bkgnd_events.push_back(Event(p,m_bkgnd_alpha));
	}
      else
	{
	  std::cerr << "Unknown background selection mode: "
		    << m_bkgnd_selection_mode << std::endl;
	  exit(EXIT_FAILURE);
	}
    }
}

void VSCutOptimizer::leaveEvent()
{

}
	
void VSCutOptimizer::visitSimEvent(const VSArraySimulationDatum& sim)
{
  m_sim = sim;

  unsigned index = m_sim.table_index;

  double area = M_PI*std::pow(m_sim_header.tables[index].sampling_radius_m,2);

  m_wt = area/m_sim_count[index]*m_sim_energy[index]*
    m_sim_spectrum->diffFlux(log10(m_sim_energy[index]))*log(10)*60*0.0625;

    //m_egywt_calc->calcWeight(sim.energy_tev);
}

void VSCutOptimizer::leaveSimEvent()
{
  
}

void VSCutOptimizer::
visitSimHeader(const VSHeaderSimulationDatum& sim_header)
{
  m_sim_header = sim_header;
  m_is_sim_run = true;
  m_egywt_calc->calcWeighting(sim_header);
  m_sim_count.resize(sim_header.tables.size(), 0);
  m_sim_energy.resize(sim_header.tables.size(), 0);
  
  for(std::vector< VSTableSimulationDatum >::const_iterator itr = 
	sim_header.tables.begin(); itr != sim_header.tables.end(); ++itr)    
    {		 
      if(itr->primary_id == 1)
	{
	  m_sim_count[itr->table_index] += itr->event_count;
	  m_sim_energy[itr->table_index] = itr->energy_tev;
	}
    }
}

void VSCutOptimizer::leaveSimHeader()
{
  m_is_sim_run = false;
}

void VSCutOptimizer::optimize()
{
  for(std::vector<Event>::iterator itr = 
 	m_sim_events.begin(); itr != m_sim_events.end(); ++itr)
    {
      itr->wt *= m_livetime/60.;
    }

  accumulate(m_train_data,0.,0.5);
  accumulate(m_test_data,0.5,1.0);
  
  optimize(m_train_data);

  std::cout 
    << std::string(79,'=') << std::endl
    << "= RESULTS (TRAIN SAMPLE)" << std::endl
    << std::string(79,'=') << std::endl;

  test(m_train_data,m_train_results);

  std::cout 
    << std::string(79,'=') << std::endl
    << "= RESULTS (TEST SAMPLE)" << std::endl
    << std::string(79,'=') << std::endl;

  test(m_test_data,m_test_results);
}

void VSCutOptimizer::copyEvents(const std::vector<Event>& events,
				double begin_fraction, double end_fraction,
				std::vector<Event>& o)
{
  int lo = std::max((int)0,(int)lround(begin_fraction*events.size()));
  int hi = std::min((int)events.size(),
		    (int)lround(end_fraction*events.size()));  

  if(lo == hi) return;
  o.resize(hi-lo);  
  std::copy( events.begin() + lo, events.begin() + hi, o.begin());    
}

void VSCutOptimizer::accumulate(Data& data, 
				double begin_fraction, double end_fraction)
{
  copyEvents( m_bkgnd_events, begin_fraction, end_fraction, data.bkgnd_events);
  copyEvents( m_sim_events, begin_fraction, end_fraction, 
	      data.signal_events);
  copyEvents( m_on_events, begin_fraction, end_fraction, data.on_events);
  copyEvents( m_off_events, begin_fraction, end_fraction, data.off_events);

  for(std::vector<Event>::iterator itr = data.bkgnd_events.begin(); 
      itr != data.bkgnd_events.end(); ++itr)
    {
      data.bkgnd.accumulate(itr->p,itr->wt);
      data.bkgnd_var.accumulate(itr->p,std::pow(itr->wt,2));
      data.bkgnd_counts.accumulate(itr->p);
    }

   if(m_options.use_theta) 
     {
       // Use assumption of flat background distribution ----------------------
       for(VSNSpace::Index i = 0; i < m_space.size(); i++)
	 {
	   VSNSpace::Cell c;
	   m_space.cellOfIndexUnchecked(i,c);
	   
	   if(c.i[0] == 0)
	     {
	       double bkgnd = data.bkgnd[c];
	       double bkgnd_var = data.bkgnd_var[c];

	       for(c.i[0] = 0; c.i[0] < m_space.axes[0].nbin; c.i[0]++)
		 {
		   double wt = 
		     m_space.axes[0].bin_size/std::pow(m_options.theta_cut,2);

		   data.bkgnd[c] = bkgnd*wt;
		   data.bkgnd_var[c] = bkgnd_var*wt*wt;
		 }
	     }
	 }
     }

   for(std::vector<Event>::iterator itr = 
	data.signal_events.begin(); itr != data.signal_events.end(); ++itr)
    {
      data.excess.accumulate(itr->p,itr->wt);
      data.excess_var.accumulate(itr->p,std::pow(itr->wt,2));
    }

   for(std::vector<Event>::iterator itr = 
	data.on_events.begin(); itr != data.on_events.end(); ++itr)
    {
      data.excess.accumulate(itr->p,itr->wt);
      data.excess_var.accumulate(itr->p,std::pow(itr->wt,2));
    }

  for(std::vector<Event>::iterator itr = 
	data.off_events.begin(); itr != data.off_events.end(); ++itr)
    {
      data.excess.accumulate(itr->p,-itr->wt);
      data.excess_var.accumulate(itr->p,std::pow(itr->wt,2));
    }

  double elaptime = m_elaptime*(end_fraction-begin_fraction);

//   if(m_options.sim_signal) data.excess_rate = m_sim_rate;
//   else 
//     {
  data.excess_rate = data.excess;
  data.excess_rate *= 60/elaptime;
      //    }

  data.bkgnd_rate = data.bkgnd;
  data.bkgnd_rate *= 60/elaptime;
  data.elaptime_min = elaptime/60.;
}

void VSCutOptimizer::save(VSOctaveH5WriterStruct* writer) const
{
  VSNSpaceOctaveH5IO* io = new VSNSpaceOctaveH5IO;

  // Write directly accumulated data
  io->writeHistogram(writer->writeStruct("excess"), m_train_data.excess);
  io->writeHistogram(writer->writeStruct("bkgnd"), m_train_data.bkgnd);
  io->writeHistogram(writer->writeStruct("bkgnd_counts"), 
		     m_train_data.bkgnd_counts);

  writer->writeScalar("livetime",m_livetime);
  writer->writeScalar("elaptime",m_elaptime);

  delete io;

  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("train_results",m_train_results.size());
  vsassert(wc);
  for(unsigned idata = 0; idata < m_train_results.size(); idata++)
    {
      VSOctaveH5WriterStruct* s = wc->writeStruct(idata);
      m_train_results[idata].save(s);
      delete s;
    }
  delete wc;

  wc = writer->writeCellVector("test_results",m_test_results.size());
  vsassert(wc);
  for(unsigned idata = 0; idata < m_test_results.size(); idata++)
    {
      VSOctaveH5WriterStruct* s = wc->writeStruct(idata);
      m_test_results[idata].save(s);
      delete s;
    }
  delete wc;
}

// ============================================================================
// VSCutOptimizerFactory
// ============================================================================

std::auto_ptr<VSCutOptimizerFactory> VSCutOptimizerFactory::s_instance;

VSCutOptimizer::Options VSCutOptimizerFactory::s_default_options = 
VSCutOptimizer::Options();

VSCutOptimizerFactory::VSCutOptimizerFactory(const 
					     VSCutOptimizer::Options& opt):
  m_options(opt)
{

}

VSCutOptimizerFactory::~VSCutOptimizerFactory()
{

}

VSCutOptimizerFactory* VSCutOptimizerFactory::getInstance()
{
  if(s_instance.get() == 0)s_instance.reset(new VSCutOptimizerFactory());
  return s_instance.get();
}

VSCutOptimizer* VSCutOptimizerFactory::createCutOptimizer()
{
  std::cerr 
    << "CUT OPTIMIZER:  " << m_options.optimizer.first << std::endl
    << "OPTIMIZER_MODE: " << m_options.optimizer.second << std::endl;

  if(m_options.optimizer.first == "simple")
    {
      return new VSSimpleCutOptimizer(m_options);
    }
  else if(m_options.optimizer.first == "nspace")
    {
      return new VSNSpaceCutOptimizer(m_options);
    }
  else
    {
      std::cerr << "Unknown optimizer type: "
		<< m_options.optimizer.first << std::endl;
      exit(EXIT_FAILURE);
    }
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSCutOptimizerFactory::configure(VSOptions& options,
				      const std::string& profile, 
				      const std::string& opt_prefix)
{
  VSScaledParameterCalc::configure(options,profile,opt_prefix);
  VSSimEnergyWeightCalc::configure(options,profile,opt_prefix);

  options.findWithValue(OPTNAME(opt_prefix,"optimizer"), 
			s_default_options.optimizer,
			"Select the cut optimization method.  Available "
			"methods are: simple/nspace.");  

  options.findWithValue(OPTNAME(opt_prefix,"space"), 
			s_default_options.space,
			"Define the parameter space on which the optimizer "
			"will operate.");  

  options.findWithValue(OPTNAME(opt_prefix,"theta_cut"), 
			s_default_options.theta_cut, 
			""); 

  options.findBoolValue(OPTNAME(opt_prefix,"use_theta"), 
			s_default_options.use_theta, true,
			"Include theta as an independent parameter in the "
			"optimization.");  

  options.findBoolValue(OPTNAME(opt_prefix,"sim_signal"), 
			s_default_options.sim_signal, true,
			"Use simulated gamma-rays for the signal sample.  "
			"If this option is set to false, the signal sample "
			"will be constructed from the gamma-ray excess "
			"present in the provided data files.");

  options.findWithValue(OPTNAME(opt_prefix,"min_events"), 
			s_default_options.min_events,
			"Minimum number of events per cell in the background "
			"histogram for cell to be considered in "
			"optimization.");  

  options.findWithValue(OPTNAME(opt_prefix,"size_cut_nscope"), 
			s_default_options.size_cut_nscope, 
			"Number of telescope images that must pass size "
			"cut. Value of zero implies that all images present "
			"must pass.");

  options.findWithValue(OPTNAME(opt_prefix,"min_rate"), 
			s_default_options.min_rate, 
			"Set the minimum acceptable gamma-ray rate, in "
			"gamma/min after cuts. Multiple values can be "
			"given."); 

  options.findWithValue(OPTNAME(opt_prefix,"min_fraction"), 
			s_default_options.min_fraction, 
			"Set the minimum acceptable gamma-ray rate as a "
			"fraction of the total gamma-ray rate before cuts.  "
			"Multiple values can be given.  This can be used in "
			"conjunction with the min_rate option."); 

  options.findWithValue(OPTNAME(opt_prefix,"bkgnd_selection"), 
			s_default_options.bkgnd_selection, 
			"Set the region of the FoV from which "
			"to draw the background sample.  The option should "
			"specify a comma-separated list starting with the "
			"method name followed by the parameters used to "
			"configure that method.  Examples: on,<rmax> "
			"theta_off,<theta>."); 

  options.findWithValue(OPTNAME(opt_prefix,"spectrum"), 
			s_default_options.spectrum, 
			"");

//   bool use_msc_disp = false;
//   options.findBoolValue("use_msc_disp", use_msc_disp, true,
// 			"Use mean scaled displacement.");

// #if 0  
//   double size_cut=0;
//   options.findWithValue("size_cut", size_cut, 
// 			"Size cut to apply in optimization.");
// #endif
}
