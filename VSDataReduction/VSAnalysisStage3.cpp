//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage3.cpp

  Stage 3 analysis

  \author     Matthew Wood                \n
              UCLA                        \n
              mdwood@astro.ucla.edu       \n

  \version    $Revision: 3.72 $
  \date       07/29/2007

  $Id: VSAnalysisStage3.cpp,v 3.72 2010/10/20 03:57:55 matthew Exp $

*/

#include <fstream>

#include <VSFileLock.hpp>
#include <VSLineTokenizer.hpp>
#include <VSFileUtility.hpp>
#include <VSEventData.hpp>
#include <VSRunInfoData.hpp>
#include <VSSimulationData.hpp>
#include <VSAnalysisStage3.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSPointing.hpp>
#include <VSAReconstruction.hpp>
#include <VSStarCatalog.hpp>
#include <VSInstrumentResponseCalc.hpp>
#include <Angle.h>

using namespace VERITAS;
using namespace SEphem;

static const char*const VERSION = 
  "$Id: VSAnalysisStage3.cpp,v 3.72 2010/10/20 03:57:55 matthew Exp $";

static const char*const REVISION = 
  "$Revision: 3.72 $";


VSAnalysisStage3::Options VSAnalysisStage3::s_default_options;

VSAnalysisStage3::Options::Options():
  max_offset_cut(1.7),
  min_energy_gev(10.),
  max_energy_gev(10E4),
  sky_bin_width_deg(0.02),
  source_exclusion_radius(0.3),
  star_exclusion_vmag_limit(6.0),
  source_pos("xy_deg",0.,0.),
  no_sim(false),
  no_run_results(false),
  write_run_hists(false),
  nevents(),
  coord_system("radec"),
  rng_seed(1)
{

}

VSAnalysisStage3::VSAnalysisStage3():
  m_cuts_evaluator(), 
  m_src_radec_J2000(),
  m_options(s_default_options)
{
  // Histogram Definitions ----------------------------------------------------
  // std::vector< std::string > hist_defs;

  // VSFileUtility::expandFilename(m_options.hist_file);
  // VSFileUtility::expandFilename(m_options.sim_library);

  // if(VSFileUtility::isFile(m_options.hist_file) && 
  //    !m_options.hist_file.empty())
  //   {
  //     std::ifstream hists(m_options.hist_file.c_str());
  //     std::string line;
  //     while(getline(hists,line))
  // 	hist_defs.push_back(line);
  //   }
  // else if(!m_options.hist_file.empty())
  //   {
  //     std::cerr 
  // 	<< std::string(__PRETTY_FUNCTION__) << ": "
  // 	<< m_options.hist_file << " is not a valid file."
  // 	<< std::endl;
  //     exit(EXIT_FAILURE);
  //   }

  // Add the list of histograms -----------------------------------------------
//   for(std::vector< std::string >::iterator itr = hist_defs.begin();
//       itr != hist_defs.end(); ++itr)
//     addHist(*itr);

  // Load the effective area library ------------------------------------------
  // if(VSFileUtility::isFile(m_options.sim_library))
  //   {
  //     VSOctaveH5Reader* sim_library = 
  // 	new VSOctaveH5Reader(m_options.sim_library);

  //     std::cout << "Loading simulations library: " 
  // 		<< m_options.sim_library << std::endl;

  //     m_sim_library_reader = 
  // 	new VSScaledParameterLibraryReader(sim_library);
  //     vsassert(m_sim_library_reader);
  //   }
  
//   if(m_options.mlm_source_model == "pointsource" && 
//      m_options.results_calc_method == "mlm" &&
//      m_sim_library_reader == NULL)
//     {
//       std::cerr 
// 	<< "Error: Effective/area psf library must be specified when using "
// 	<< "'pointsource' mlm source model." << std::endl;
//       exit(EXIT_FAILURE);
//     }


  m_array_subset.insert("ra_J2000");
  m_array_subset.insert("dec_J2000");
  m_array_subset.insert("mean_derotated_fov_x");
  m_array_subset.insert("mean_derotated_fov_y");
  m_array_subset.insert("mean_fov_x");
  m_array_subset.insert("mean_fov_y");
  m_array_subset.insert("dec_J2000");
  m_array_subset.insert("theta0");
  m_array_subset.insert("theta1");
  m_array_subset.insert("scope");
  m_array_subset.insert("theta");
  m_array_subset.insert("mlt_log10_energy");  
  m_array_subset.insert("mlt_log10_energy_chi2");  
  m_array_subset.insert("R");  
  m_array_subset.insert("msc_width");    
  m_array_subset.insert("msc_length");  
  m_array_subset.insert("msc_disp");  
  m_array_subset.insert("N2");  
  m_array_subset.insert("used_in_reconstruction_mask");  
  m_array_subset.insert("has_image_mask");  
  m_array_subset.insert("trigger_mask");  

  m_scope_subset.insert("N");
  m_scope_subset.insert("R");
  m_scope_subset.insert("intrinsic_width");
  m_scope_subset.insert("intrinsic_length");
  m_scope_subset.insert("fp_dist");
  m_scope_subset.insert("fp_disp");
  m_scope_subset.insert("fp_width");
  m_scope_subset.insert("fp_length");
  m_scope_subset.insert("fp_xc");
  m_scope_subset.insert("fp_yc");
  m_scope_subset.insert("sc_width");
  m_scope_subset.insert("sc_length");
  m_scope_subset.insert("sc_disp");
  m_scope_subset.insert("G");
  m_scope_subset.insert("delta2l");
  m_scope_subset.insert("delta2m");
  m_scope_subset.insert("lambdac");
  m_scope_subset.insert("lambdad");
  m_scope_subset.insert("nimage");
  m_scope_subset.insert("used_in_reconstruction");

  m_cuts_evaluator = new VSCutsEvaluator;
  m_cuts_evaluator->getScopeParamSet(m_scope_subset);
  m_cuts_evaluator->getArrayParamSet(m_array_subset);
  m_cuts_evaluator->print();
}

VSAnalysisStage3::~VSAnalysisStage3()
{
  delete m_cuts_evaluator;
}

void VSAnalysisStage3::runStage3(const std::list<std::string>& files,
				 VSOctaveH5Writer* writer)
{
  VSAnalysisData analysis_data(REVISION,VERSION);

  VSEventDataDispatcher* dispatcher = new VSEventDataDispatcher;
  dispatcher->loadList(files);

  std::vector< std::string > stg2_files = dispatcher->getFiles();

  // --------------------------------------------------------------------------
  // Determine coordinates of source 
  // --------------------------------------------------------------------------
  std::cout 
    << std::string(30,'-')
    << std::setw(19) << "Run List     "
    << std::string(30,'-')
    << std::endl;

  std::cout << std::setw(9) << "RUN"
	    << std::setw(20) << "SOURCE"
	    << std::setw(18) << "MODE"
	    << std::setw(13) << "OBS RA"
	    << std::setw(13) << "OBS DEC"
	    << std::endl << std::endl;

  VSAAlgebra::Vec3D src_radec;
  std::vector< VSAAlgebra::Vec3D > obs_radec;
  std::vector<unsigned> nrun_obs;
  std::map< unsigned, VSAnalysisStage3Visitor::Run > run_list;


  double max_dist_deg = 0;
  bool is_sim = false;
  std::set< std::string > src_names;
  std::string src_name;

  for(std::vector<std::string>::iterator itr = stg2_files.begin();
      itr != stg2_files.end(); ++itr)
    {
      VSFileLockBSD fl(*itr);

      VSOctaveH5Reader *reader = new VSOctaveH5Reader(*itr);
      VSRunInfoData run_info;
      VSOctaveH5ReaderStruct* stg1_struct = reader->readStruct("stage1");
      vsassert(stg1_struct);
      run_info.load(stg1_struct->readStruct("run_info"));

      VSTargetTable::Observation obs;

      obs.load(reader->readStruct("observation"));

      if(reader->readStruct("sim_event") != NULL) is_sim = true;

      std::cout 
	<< std::setw(9) << run_info.run_number
	<< std::setw(20) << obs.name
	<< std::setw(18) << obs.mode_string
	<< std::setw(13) << obs.obs_radec_J2000.longitude().hmsString(1) 
	<< std::setw(13) << obs.obs_radec_J2000.latitude().dmsString(1) 
	<< std::endl;

      if(obs.mode == VSTargetTable::Observation::OM_DRIFT ||
	 obs.mode == VSTargetTable::Observation::OM_UNKNOWN)
	{
	  std::cerr << std::string(79,'=') << std::endl;
	  std::cerr << "Skipping Run " << run_info.run_number
		    << " for invalid observation mode." << std::endl;
	  std::cerr << std::string(79,'=') << std::endl;
	  itr = stg2_files.erase(itr);
	  itr--;
	  delete reader;
	  continue;
	}

      double src_dec_J2000 = obs.src_radec_J2000.latitudeRad();
      double src_ra_J2000  = obs.src_radec_J2000.longitudeRad();
      double obs_dec_J2000 = obs.obs_radec_J2000.latitudeRad();
      double obs_ra_J2000  = obs.obs_radec_J2000.longitudeRad();

//       reader->readScalar("observation.src_dec_J2000_rad",src_dec_J2000);
//       reader->readScalar("observation.src_ra_J2000_rad",src_ra_J2000);
//       reader->readScalar("observation.obs_dec_J2000_rad",obs_dec_J2000);
//       reader->readScalar("observation.obs_ra_J2000_rad",obs_ra_J2000);

      const unsigned nobs = obs_radec.size();
      bool found_obs = false;
      for(unsigned iobs = 0; iobs < nobs; iobs++)
	{
	  SphericalCoords radec(obs_radec[iobs].theta(),obs_radec[iobs].phi());

	  if(radec.separation(obs.obs_radec_J2000).deg() < 0.02)
	    {
	      obs_radec[iobs] *= (double)nrun_obs[iobs];
	      obs_radec[iobs] += 
		VSAAlgebra::Vec3D::makePolar(M_PI/2.-obs_dec_J2000,
					     obs_ra_J2000);
	      nrun_obs[iobs]++;
	      obs_radec[iobs] /= (double)nrun_obs[iobs];
	      found_obs = true;
	      run_list[run_info.run_number].iptg = iobs;
	    }
	}

      if(!found_obs)
	{
	  run_list[run_info.run_number].iptg = obs_radec.size();
	  obs_radec.
	    push_back(VSAAlgebra::Vec3D::makePolar(M_PI/2.-obs_dec_J2000,
						   obs_ra_J2000));
	  nrun_obs.push_back(1);	  
	}

      double dist_deg = 
	obs.src_radec_J2000.separation(obs.obs_radec_J2000).deg();

      if(dist_deg > max_dist_deg) max_dist_deg = dist_deg;

      src_names.insert(obs.name);

      src_radec += 
	VSAAlgebra::Vec3D::makePolar(M_PI/2.-src_dec_J2000,src_ra_J2000);

      delete reader;
    }
  
  double nfiles = stg2_files.size();
  src_radec *= (1./nfiles);
  m_src_radec_J2000 = SphericalCoords(src_radec.theta(),src_radec.phi());

  SphericalCoords origin_radec_J2000 = m_src_radec_J2000;
  
  if(m_options.source_pos.first == "xy_deg")
    {
      Astro::xyToRaDec(m_options.source_pos.second,
		       m_options.source_pos.third,
		       origin_radec_J2000,
		       m_src_radec_J2000);
    }
  else if(m_options.source_pos.first == "radec_deg")
    {
      m_src_radec_J2000.setLatLongDeg(m_options.source_pos.third,
				      m_options.source_pos.second);
    }
  else if(!m_options.source_pos.first.empty())
    {
      std::cerr 
	<< std::string(__PRETTY_FUNCTION__) << ": "
	<< "Unrecognized coordinate type in 'source_pos': "
	<< m_options.source_pos.first << std::endl;
      exit(EXIT_FAILURE);
    }

  src_name = *src_names.begin();

  std::cout << std::endl
	    << "Source Position:   " 
	    << std::setw(13) 
	    << m_src_radec_J2000.longitude().hmsString(1) 
	    << std::setw(13)
	    << m_src_radec_J2000.latitude().dmsString(1) 
	    << std::endl;

  std::cout << std::endl
	    << "Observation Coordinates:   " << std::endl;

  std::cout << std::setw(5) << "IOBS"
	    << std::setw(5) << "NRUN"
	    << std::setw(13) << "RA" 
	    << std::setw(13) << "DEC"
	    << std::setw(17) << "OFFSET [DEG]"
	    << std::setw(17) << "PHI [DEG]"
	    << std::endl << std::endl;

  const unsigned nobs = obs_radec.size();
  for(unsigned iobs = 0; iobs < nobs; iobs++)
    {
      m_ptg_radec_J2000.push_back(SphericalCoords(obs_radec[iobs].theta(),
						  obs_radec[iobs].phi()));

      std::cout 
	<< std::setw(5) << iobs
	<< std::setw(5) << nrun_obs[iobs]
	<< std::setw(13) << m_ptg_radec_J2000[iobs].longitude().hmsString(1)
	<< std::setw(13) << m_ptg_radec_J2000[iobs].latitude().dmsString(1) 
	<< std::setw(17) 
	<< m_ptg_radec_J2000[iobs].separation(m_src_radec_J2000).deg()
	<< std::setw(17) 
	<< m_ptg_radec_J2000[iobs].directionTo(m_src_radec_J2000).deg()
	<< std::endl;
    }

  for(std::map< unsigned, VSAnalysisStage3Visitor::Run >::iterator
	itr = run_list.begin(); itr != run_list.end(); ++itr)
    itr->second.ptg_radec = m_ptg_radec_J2000[itr->second.iptg];

  // --------------------------------------------------------------------------
  // Define exclusion regions
  // --------------------------------------------------------------------------
  std::vector< VSStarCatalog::Star > stars;
  Angle max_distance((max_dist_deg+m_options.max_offset_cut+0.25)*M_PI/180.);

  if(!is_sim)
    VSStarCatalog::getStars(stars,origin_radec_J2000,max_distance,
			    m_options.star_exclusion_vmag_limit);
  
  if(stars.size() > 0)
    std::cout 
      << std::setw(35) << "Star Exclusion Regions "
      << std::string(44,'-') << std::endl
      << std::setw(15) << "V MAG"
      << std::setw(13) << "RA" << std::setw(13) << "DEC"
      << std::setw(13) << "X" << std::setw(13) << "Y"
      << std::endl << std::endl;

  VSExclusionRegion exclusion_region;

  for(std::vector< VSStarCatalog::Star >::iterator itr = 
	stars.begin(); itr != stars.end(); ++itr)
    {
      std::pair< Angle, Angle > xy;
      Astro::raDecToXY(itr->radec_j2000,origin_radec_J2000,xy);
      exclusion_region.addStar(xy.first.degPM(),xy.second.degPM(),
			       itr->radec_j2000,0.25,itr->v_mag);

      std::cout 
	<< std::setw(15) << itr->v_mag
	<< std::setw(13) << itr->radec_j2000.longitude().hmsString(1) 
	<< std::setw(13) << itr->radec_j2000.latitude().dmsString(1) 
	<< std::setw(13) << xy.first.degPM()
	<< std::setw(13) << xy.second.degPM()
	<< std::endl;
   }

  std::pair< Angle, Angle > src_xy;
  Astro::raDecToXY(m_src_radec_J2000,origin_radec_J2000,src_xy);
  if(m_options.source_exclusion_radius > 0)
    exclusion_region.addSource(src_xy.first.degPM(),src_xy.second.degPM(),
			       m_src_radec_J2000,	       
			       m_options.source_exclusion_radius);

  // --------------------------------------------------------------------------
  // Construct the stage3 analysis visitor
  // --------------------------------------------------------------------------
  m_visitor = new VSAnalysisStage3Visitor(m_cuts_evaluator,
					  writer,
					  m_options.sky_bin_width_deg,
					  m_options.max_offset_cut,
					  m_options.coord_system,
					  src_name,
					  origin_radec_J2000,
					  m_src_radec_J2000,
					  m_ptg_radec_J2000,
					  exclusion_region,
					  run_list,
					  m_options.no_run_results,
					  m_options.write_run_hists,
					  m_options.rng_seed);

  dispatcher->setVisitor(m_visitor);

  // --------------------------------------------------------------------------
  // Loop over stage2 files
  // --------------------------------------------------------------------------
  dispatcher->processFiles(m_options.nevents,true);  
  
  // --------------------------------------------------------------------------
  // Clean up
  // --------------------------------------------------------------------------
  delete m_visitor;

  exclusion_region.save(writer);

  VSOctaveH5WriterCellVector* wc = 
    writer->writeCellVector("files", stg2_files.size());
  for(unsigned ifile = 0; ifile < stg2_files.size(); ifile++)
    wc->writeString(ifile,stg2_files[ifile]);
  delete wc;

  VSOctaveH5WriterStruct* analyze_struct = writer->writeStruct("analyze");  
  analysis_data.stop();
  analysis_data.save(analyze_struct);
  delete analyze_struct;
  analysis_data.printTiming(std::cout);
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSAnalysisStage3::configure(VSOptions& options,
				 const std::string& profile, 
				 const std::string& opt_prefix)
{
  VSIntegralAnalysisFactory::configure(options,profile,opt_prefix);
  VSCutsEvaluator::configure(options,profile,opt_prefix);
  VSScaledParameterCalc::configure(options,profile,opt_prefix);
  VSEnergyCalcLT::configure(options,profile,opt_prefix);
  VSInstrumentResponseCalc::configure(options,profile,opt_prefix);
  VSSimEnergyWeightCalc::configure(options,profile,opt_prefix);
  VSSourceInjector::configure(options,profile,opt_prefix);
  VSSpectrumCalcFactory::configure(options,profile,opt_prefix);


  options.addCatagory(OPTNAME(opt_prefix,"sim"),
		      "Options related to processing simulation "
		      "files.");

  options.addCatagory(OPTNAME(opt_prefix,"energy"),
		      "Options related to energy reconstruction.");

  options.findBoolValue(OPTNAME(opt_prefix,"no_sim"), 
			s_default_options.no_sim, true,
			"Skip processing of simulation data.","s3_sim");

  options.findBoolValue(OPTNAME(opt_prefix,"no_run_results"),
			s_default_options.no_run_results, true,
			"Disable computation of separate analysis results "
			"for each run.");

  options.findBoolValue(OPTNAME(opt_prefix,"write_run_hists"), 
			s_default_options.write_run_hists, true,
			"Write results histograms for individual runs.");

  options.findWithValue(OPTNAME(opt_prefix,"nevents"), 
			s_default_options.nevents,
			"Set the maximum number of events to process.");

  options.findWithValue(OPTNAME(opt_prefix,"max_offset_cut"), 
			s_default_options.max_offset_cut,
			"Set the cut on the maximum offset in degrees of the "
			"reconstructed direction of an event relative to the "
			"center of the FoV.");

  options.findWithValue(OPTNAME(opt_prefix,"source_exclusion_radius"), 
			s_default_options.source_exclusion_radius,
			"Radius of excluded region in degrees used when "
			"the 'exclude_source' option is enabled.  Set this "
			"to zero to disable source exclusion.");

  options.findWithValue(OPTNAME(opt_prefix,"star_exclusion_vmag_limit"), 
			s_default_options.star_exclusion_vmag_limit,
			"Limiting V-band magnitude which will be used "
			"to select stars for exclusion.");

  options.findWithValue(OPTNAME(opt_prefix,"source_pos"), 
			s_default_options.source_pos,
			"Set the source position coordinates for which "
			"results (significance,rate,flux,etc.) "
			"will be computed.  If no "
			"position is defined, then the "
			"source position will be determined from the "
			"average over all runs.  Coordinates are "
			"specified using the syntax: coord_type/phi/theta.  "
			"Available coordinate types are: xy_deg, radec_deg, "
			"gal_deg."
			"The coordinate type 'xy_deg' can be used to set the "
			"source position in terms of a relative offset from "
			"the nominal source position in the chosen coordinate "
			"system.");  
  
  options.findWithValue(OPTNAME(opt_prefix,"sky_bin_width"), 
			s_default_options.sky_bin_width_deg,
			"Set the sky map bin width in degrees.");

  options.findWithValue(OPTNAME(opt_prefix,"coord_system"), 
			s_default_options.coord_system,
			"Coordinate system used for generating sky maps.  "
			"Options are radec/galactic.");

  options.findWithValue(OPTNAME(opt_prefix,"rng_seed"),
			s_default_options.rng_seed,
			"Set the random number generator seed.");



}

