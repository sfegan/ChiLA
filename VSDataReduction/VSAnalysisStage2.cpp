//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAnalysisStage2.cpp

  Driver for stage 2 analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 3.113 $
  \date       05/10/2006

  $Id: VSAnalysisStage2.cpp,v 3.113 2010/03/23 00:47:47 matthew Exp $

*/

#include <sstream>
#include <iomanip>

#include <VSFileUtility.hpp>
#include <VSAnalysisStage2.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSEnergyCalcLT.hpp>

using namespace VERITAS;
using namespace SEphem;

static const char*const VERSION = 
  "$Id: VSAnalysisStage2.cpp,v 3.113 2010/03/23 00:47:47 matthew Exp $";

static const char*const REVISION = 
  "$Revision: 3.113 $";

static const char*const LOGO = 
  "  ____      ________    _ __    ___     ____    \n"
  " / / /     / ____/ /_  (_) /   /   |    \\ \\ \\   \n"
  "/ / /     / /   / __ \\/ / /   / /| |     \\ \\ \\   __            _  \n"
  "\\ \\ \\    / /___/ / / / / /___/ ___ |     / / /  (__|_ _. _  _   ) \n"
  " \\_\\_\\   \\____/_/ /_/_/_____/_/  |_|    /_/_/   __)|_(_|(_|(/_ /_ \n"
  "                                                        ._|       \n"
  "\n"
  "$Id: VSAnalysisStage2.cpp,v 3.113 2010/03/23 00:47:47 matthew Exp $\n";

// ----------------------------------------------------------------------------
// __   _____   _             _         _    ___ _                 ___
// \ \ / / __| /_\  _ _  __ _| |_  _ __(_)__/ __| |_ __ _ __ _ ___|_  )
//  \ V /\__ \/ _ \| ' \/ _` | | || (_-< (_-<__ \  _/ _` / _` / -_)/ /
//   \_/ |___/_/ \_\_||_\__,_|_|\_, /__/_/__/___/\__\__,_\__, \___/___|
//                              |__/                     |___/
// ----------------------------------------------------------------------------

VSAnalysisStage2::Options::Options():
  integration_lo_zero_sample(4,6),
  integration_hi_zero_sample(4,0),
  integration_threshold_frac(0.5),
  integration_window_start(0.5),
  integration_window_width(8),
  integration_apply_increase(true),
  integration_threshold_charge(100.0),
  integration_window_start_increase(0.5),
  integration_window_width_increase(3),
  ped_suppress_mode("interval"),
  ped_suppress_hi(1.5),
  ped_suppress_lo(1/1.5),
  ped_suppress_fraction(0.85),
  ped_suppress_scale(2.0),
  ped_event_count_min(100),
  gain_suppress_hi(1.3),
  gain_suppress_lo(1/1.3),
  no_laser(false),
  permissive_laser(false),
  unity_gains(false),
  primary_cleaning_args(),
  secondary_cleaning_args(),
  cleaning_scale("pedrms"),
  all_events(false),
  nimage_cut(1),
  nscope_cut(2),
  l2_trigger_threshold(3),
  l3_trigger_threshold(2),
  scope_gain(4,1.0),
  scope_dev_scaling(),
  scope_hi_lo_gain_ratio(4,6.0),
  scope_suppress(),
  auto_suppress_l2_chan(true),
  pad_zero_suppressed_chan(false),
  chan_no_pmt(),
  chan_suppress(),
  scope_pos(),
  no_set_scope_pos_from_simulations(false),
  camera("veritas/499"),
  cam_rotation(4,0.0),
  atmosphere(),
  method(1),
  weighting("ellipticity"),
  theta_cut(0),
  no_diagnostics(false),
  no_slow_diagnostics(false),
  print_frequency(10000),
  pointing_string("positioner"),
  tracking_recorrections_file(),
  tracking_recorrections_date(),
  tracking_target(""),
  no_commanded_target(false),
  wobble_region_radius(0.18),
  enable_l2_corrections(true),
  theta_targets(),
  no_reorder(),
#ifndef NOTHREADS
  nthreads(),
  no_vbf_reader_thread(),
#endif
  npackets(),
  no_muon_analysis(false),
  muon_raw_ring_radius(0.15),
  muon_nimage_cut(50),       //0),
  muon_radius_min_cut(0.7),  //0),
  muon_radius_max_cut(1.4),  //std::numeric_limits<double>::infinity()),
  muon_rms_max_cut(0.2),     //std::numeric_limits<double>::infinity()),
  muon_ring_edge_dist_max_cut(1.8), //std::numeric_limits<double>::infinity()),
  muon_centroid_radius_ratio_max_cut(0.75), //std::numeric_limits<double>::infinity())
#ifdef MUON_TEST_NON_UNIFORMITY
  muon_non_uniformity_beta(0),
#endif
  primary_qc(0,0,std::numeric_limits<double>::infinity(),2),
  secondary_qc(0,0,std::numeric_limits<double>::infinity(),2),
  nspace_cuts(),
  software_trigger_masks(),
  qc_masks(),
  simple_trigger_masks(0,"0123"),
  limited_dt_range(0.001,0.025),
  sc_parameter_lookup_file(),
  msc_weight_power(1.0),
  psf_poly_radial(),
  psf_poly_tangential()
{
  primary_cleaning_args.push_back("regional");
  primary_cleaning_args.push_back("4");
  primary_cleaning_args.push_back("3");
  primary_cleaning_args.push_back("9");

#if 0
  unsigned suppress_all_scope[] = { }; //128, 249, 259, 498, 499 };
  for(unsigned iscope=0;iscope<scope_pos.size();iscope++)
    for(unsigned isuppress=0;
	isuppress<sizeof(suppress_all_scope)/sizeof(*suppress_all_scope);
	isuppress++)
      {
	unsigned ichan = suppress_all_scope[isuppress];
	chan_no_pmt.push_back(std::make_pair<unsigned,unsigned>(iscope,ichan));
      }
#endif

  theta_targets.push_back(target_type("auto","auto"));

  scope_dev_scaling.push_back(LinearCoefficients(0.801,0.265));
  scope_dev_scaling.push_back(LinearCoefficients(0.837,0.294));
  scope_dev_scaling.push_back(LinearCoefficients(0.900,0.145));
  scope_dev_scaling.push_back(LinearCoefficients(1.000,0.000));
}

VSAnalysisStage2::Options VSAnalysisStage2::s_default_options;

VSAnalysisStage2::VSAnalysisStage2(): m_options(s_default_options)
{
  // nothing to see here
}

VSAnalysisStage2::~VSAnalysisStage2()
{
  // nothing to see here
}

void calculateMasks(std::vector<unsigned>& bit_masks,
		    const std::vector<std::string>& str_masks)
{
  std::set<unsigned> masks;
  for(unsigned itrig=0;itrig<str_masks.size();itrig++)
    {
      std::string str_mask(str_masks[itrig]);
      std::set<unsigned> mask;
      for(unsigned ichar=0;ichar<str_mask.size();ichar++)
	if(isdigit(str_mask[ichar]))mask.insert(str_mask[ichar]-'0');
      unsigned m=0;
      for(std::set<unsigned>::const_iterator imask=mask.begin();
	  imask!=mask.end();imask++)
	m |= 1<<(*imask);
      masks.insert(m);
    }
  for(std::set<unsigned>::const_iterator imask=masks.begin();
      imask!=masks.end();imask++)bit_masks.push_back(*imask);
}

void VSAnalysisStage2::
runStage2(const std::string& vbf_filename, 
	  VSAnalysisStage1Data* stage1, 
	  VSOctaveH5WriterStruct* writer,
	  const SEphem::SphericalCoords& earth_position,
	  double earth_elevation,
	  VSAnalysisStage1Data* pad_stage1,
	  bool no_db, bool no_l3, bool no_threads, bool no_verbose,
	  bool do_use_overflow)
{
  if(!no_verbose)vstream << LOGO << std::endl;
  vsassert(stage1);

  VSAnalysisData analysis_data(REVISION,VERSION);

  if((stage1->run_info.got_run_number)&&(!no_verbose))
    vstream << "Run number:  " << stage1->run_info.run_number
	    << std::endl << std::endl;

  // --------------------------------------------------------------------------
  // Suppress channels based on pedestal variances
  // --------------------------------------------------------------------------
  
  if(m_options.ped_suppress_mode == "median")
    {
      if((m_options.ped_suppress_lo>=0)&&(m_options.ped_suppress_hi>=0))
	{
	  std::vector<std::pair<double,double> > cuts;
	  stage1->pedestals.
	    suppressByMedianFraction(m_options.ped_suppress_lo, 
				     m_options.ped_suppress_hi,
				     stage1->suppress, cuts,
				     m_options.ped_event_count_min);
	  if(!no_verbose)
	    {
	      vstream << "Pedestal RMS cuts:" << std::endl;
	      unsigned nscope = cuts.size();
	      for(unsigned iscope=0;iscope<nscope;iscope++)
		{
		  vstream << "  T" << iscope+1 << ": "
			  << cuts[iscope].first << ' '
			  << cuts[iscope].second << std::endl;
		}
 	      vstream << std::endl;
	    }
	}
    }
  else if(m_options.ped_suppress_mode == "interval")
    {
      if((m_options.ped_suppress_fraction>=0.0)
	 &&(m_options.ped_suppress_fraction<=1.0)
	 &&(m_options.ped_suppress_scale>=0.0))
	{
	  std::vector<std::pair<double,double> > cuts;
	  stage1->pedestals.
	    suppressByContainmentFraction(m_options.ped_suppress_fraction,
					  m_options.ped_suppress_scale,
					  stage1->suppress, cuts,
					  m_options.ped_event_count_min);
	  if(!no_verbose)
	    {
	      vstream << "Pedestal RMS cuts:" << std::endl;
	      unsigned nscope = cuts.size();
	      for(unsigned iscope=0;iscope<nscope;iscope++)
		{
		  vstream << "  T" << iscope+1 << ": "
			  << cuts[iscope].first << ' '
			  << cuts[iscope].second << std::endl;
		}
 	      vstream << std::endl;
	    }
	}
    }
  else if(m_options.ped_suppress_mode != "none")
    {
      vstream << "Unknown pedestal suppression mode: " 
	      << m_options.ped_suppress_mode << std::endl;
      vsassert(0);
    }

  // --------------------------------------------------------------------------
  // Suppress channels based on LASER gain
  // --------------------------------------------------------------------------

  if((stage1->laser)&&(m_options.no_laser))
    {
      delete stage1->laser;
      stage1->laser=0;
    }
  else if((stage1->laser==0)&&(!m_options.no_laser))
    {
      vstream 
	<< "No LASER data supplied but option to not require LASER" <<std::endl
	<< "data was not selected. Either supply laser data or set" <<std::endl
	<< "option." << std::endl;
      vsassert(0);
    }

  if((stage1->laser)
     &&(m_options.gain_suppress_lo>=0)&&(m_options.gain_suppress_hi>=0))
    {
      stage1->laser->suppress(m_options.gain_suppress_lo, 
			      m_options.gain_suppress_hi);
    }

  // --------------------------------------------------------------------------
  // Load compiled in targets if none were supplied in the Stage 1 file
  // --------------------------------------------------------------------------

  if(!stage1->target_table)stage1->target_table = new VSTargetTableData;
  if(stage1->target_table->empty())
    *(stage1->target_table) = VSTargetTable::getCompiledTargets();

  // --------------------------------------------------------------------------
  // Recorrect the pointing from the stage1 file if requested
  // --------------------------------------------------------------------------

  if((!m_options.tracking_recorrections_file.empty())
     ||(!m_options.tracking_recorrections_date.empty()))
    {
      vstream << "Re-appling tracking corrections\n";
      
      if((m_options.pointing_string == "positioner")
	 &&(stage1->db_pointing))
	{
	  std::vector<CorrectionParameters*> 
	    scope_cp(stage1->db_pointing->scope.size());
	  
	  std::vector<std::string> 
	    scope_cp_str(stage1->db_pointing->scope.size());

	  unsigned nfile = m_options.tracking_recorrections_file.size();
	  for(unsigned ifile=0; ifile<nfile; ifile++)
	    {
	      unsigned iscope = 
		m_options.tracking_recorrections_file[ifile].first;

	      if((iscope>=stage1->db_pointing->scope.size())
		 ||(stage1->db_pointing->scope[iscope].empty()))
		{
		  vstream << "WARN: Tracking re-correction requested for T"
			  << iscope+1 << " which is either invalid\n"
			  << "or has no tracking data" << std::endl;
		  continue;
		}
	      
	      std::string name = 
		m_options.tracking_recorrections_file[ifile].second;
	      VSFileUtility::expandFilename(name);

	      CorrectionParameters cp;
	      if(!cp.load(name.c_str()))
		{
		  vstream << "WARN: Tracking re-correction for T"
			  << iscope+1 << " could not be loaded\n"
			  << "from file: " << name  << std::endl;
		  continue;
		}
	    
	      delete scope_cp[iscope];
	      scope_cp[iscope] = new CorrectionParameters(cp);
	      scope_cp_str[iscope] = std::string("file: ")+name;
	    }

	  unsigned ndate = m_options.tracking_recorrections_date.size();
	  for(unsigned idate=0; idate<ndate; idate++)
	    {
	      std::string scopes = 
		m_options.tracking_recorrections_date[idate].first;

	      std::string sdate = 
		m_options.tracking_recorrections_date[idate].second;

	      VSTime date;
	      bool next_please = false;

	      if(sdate == "next")
		{
		  next_please = true;
		  date = stage1->run_info.lo_event_time;
		}
	      else
		{
		  date.setFromString(sdate);
		  if(!date.isOK())
		    {
		      vstream << "WARN: Tracking re-correction requested for "
			      << "date \"" << sdate << "\" which is not valid"
			      << std::endl;
		      continue;
		    }
		}
		  
	      if(scopes=="*")
		{
		  scopes = "";
		  for(unsigned iscope=0;
		      iscope<stage1->db_pointing->scope.size();iscope++)
		    if(!stage1->db_pointing->scope[iscope].empty())
		      {
			char c[] = "0\0";
			c[0] = char(unsigned('0')+iscope);
			scopes += std::string(c);
		      }
		}

	      for(unsigned ichar=0;ichar<scopes.size();ichar++)
		{
		  std::string sscope = scopes.substr(ichar,1);
		  if(!isdigit(sscope[0]))
		    {
		      vstream << "WARN: Tracking re-correction requested for "
			      << "telescope id \"" << sscope << "\" which\n"
			      << "is not a valid number" << std::endl;
		      continue;
		    }

		  unsigned iscope = unsigned(sscope[0] - '0');
		  if((iscope>=stage1->db_pointing->scope.size())
		     ||(stage1->db_pointing->scope[iscope].empty())
		     ||(stage1->misc_db == 0)
		     ||(iscope>=stage1->misc_db->scope.size()))
		    {
		      vstream << "WARN: Tracking re-correction requested for T"
			      << iscope+1 << " which is either invalid\n"
			      << "or has no tracking data" << std::endl;
		      continue;
		    }
		  
		  VSMiscellaneousDBData::CorrectionParametersDatum cpd;
		  bool found = false;
		  
		  if(next_please)
		    found = stage1->misc_db->scope[iscope].
		      getNextCorrectionParameters(date,cpd);
		  else
		    found = stage1->misc_db->scope[iscope].
		      getCorrectionParameters(date,cpd);
		  
		  if(!found)
		    {
		      vstream << "WARN: Tracking re-correction for T"
			      << iscope+1 << " could not be found in DB\n"
			      << "from date: " << date.toString() << std::endl;
		      continue;
		    }
	    
		  delete scope_cp[iscope];
		  scope_cp[iscope] = new CorrectionParameters(cpd);
		  scope_cp_str[iscope] = 
		    std::string("db: ")+cpd.db_start_time.toString();
		}
	    }

	  for(unsigned iscope=0;iscope<scope_cp.size();iscope++)
	    if(scope_cp[iscope])
	      {
		double off_x = 0;
		double off_y = 0;
		stage1->db_pointing->recorrect(*scope_cp[iscope],
					       iscope,&off_x,&off_y);
		delete scope_cp[iscope];

		vstream << 'T' << iscope+1 << " corrections: " 
			<< scope_cp_str[iscope] << '\n'
			<< 'T' << iscope+1 
			<< " mean offset of new coordinates WRT old: "
			<< std::fixed
			<< std::setprecision(6) << Angle::toDeg(off_x) << ',' 
			<< std::setprecision(6) << Angle::toDeg(off_y) 
			<< '\n';
	      }
	  vstream << '\n';
	}
      else
	{
	  vstream 
	    << "WARN: Tracking re-correction requested but pointing mode is "
	    << "not \"positioner\"\n"
	    << "or pointing information is not available from stage 1\n"
	    << std::endl;
	}
    }

  // --------------------------------------------------------------------------
  // Write stage 1 information
  // --------------------------------------------------------------------------

  VSOctaveH5WriterStruct* s = writer->writeStruct("stage1");
  stage1->save(s);
  delete s;

  if(pad_stage1)
    {
      VSOctaveH5WriterStruct* s = writer->writeStruct("pad_stage1");
      delete pad_stage1->laser;
      pad_stage1->laser=0;
      pad_stage1->save(s);
      delete s;
    }

  // --------------------------------------------------------------------------
  // Set up the camera and array definition
  // --------------------------------------------------------------------------
  VSChannelMap channel_map(stage1->run_info.lo_event_time);

  unsigned nscope = stage1->run_info.nchan.size();
  VSAReconstruction::ArrayInfo array_info(nscope);
  array_info.elevation = earth_elevation;

  unsigned camera_nchan = 0;
  const float* camera_xcoord = 0;
  const float* camera_ycoord = 0;
  const int (*camera_neighbors)[NUM_NEIGHBORS] = 0;

  if(m_options.camera == "veritas/499")
    {
      camera_nchan = sizeof(VC499GroundYcoord)/sizeof(*VC499GroundYcoord);
      camera_xcoord = VC499GroundXcoord;
      camera_ycoord = VC499GroundYcoord;
      camera_neighbors = VC499Neighbors;
    }
  else
    {
      std::cerr << "Unknown camera \"" << m_options.camera
		<< "\" .. aborting" << std::endl;
      exit(EXIT_FAILURE);
    }

  if((stage1->sim_info)&&(!m_options.no_set_scope_pos_from_simulations))
    {
      double sum_z = 0;
      m_options.scope_pos = stage1->sim_info->scope_positions;
      for(unsigned iscope=0;iscope<m_options.scope_pos.size();iscope++)
	sum_z += m_options.scope_pos[iscope].third;
      double mean_z = sum_z/double(m_options.scope_pos.size());
      for(unsigned iscope=0;iscope<m_options.scope_pos.size();iscope++)
	m_options.scope_pos[iscope].third -= mean_z;
      array_info.elevation = mean_z;
    }

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const double t = Angle::frDeg(m_options.cam_rotation[iscope]);
      const double c = cos(t);
      const double s = sin(t);

      if(iscope < m_options.scope_pos.size())
	{
	  array_info.scopes[iscope].ri = 
	    VSAAlgebra::Vec3D(m_options.scope_pos[iscope].first,
			      m_options.scope_pos[iscope].second,
			      m_options.scope_pos[iscope].third);
	}
      else
	{
	  pos_type scope_pos = channel_map.posForScope(iscope);
	  array_info.scopes[iscope].ri = 
	    VSAAlgebra::Vec3D(scope_pos.first,
			      scope_pos.second,
			      scope_pos.third);
	}

      array_info.scopes[iscope].pij.resize(camera_nchan);
      for(unsigned ichan=0;ichan<camera_nchan;ichan++)
	{
	  double x = Angle::frDeg(camera_xcoord[ichan]);
	  double y = Angle::frDeg(camera_ycoord[ichan]);
	  if((iscope<m_options.cam_rotation.size())
	     &&(m_options.cam_rotation[iscope] != 0))
	    {
	      const double _x = c*x - s*y;
	      const double _y = s*x + c*y;
	      x = _x;
	      y = _y;
	    }
	  array_info.scopes[iscope].pij[ichan] = VSAAlgebra::Vec2D(x,y);
	}

      array_info.scopes[iscope].Di = 12;
    }

  // --------------------------------------------------------------------------
  // Create the atmospheric profile
  // --------------------------------------------------------------------------

  VSAAtmosphere atmo;
  if(!m_options.atmosphere.empty())atmo = VSAAtmosphere(m_options.atmosphere);
  if(!atmo.good())atmo = VSAAtmosphere::usStandard();

  // --------------------------------------------------------------------------
  // Create the reconstruction
  // --------------------------------------------------------------------------

  VSAReconstruction::ArrayQualityCuts qc(nscope);
  qc.nuse_in_reconstruction_min = s_default_options.primary_qc.fourth;
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      qc.scopes[iscope].nimage_min = s_default_options.primary_qc.first;
      qc.scopes[iscope].N_min      = s_default_options.primary_qc.second;
      qc.scopes[iscope].dist2_max = 
	Angle::frDeg(s_default_options.primary_qc.third)
	*Angle::frDeg(s_default_options.primary_qc.third);
    }

  VSAReconstruction::ArrayQualityCuts sqc(nscope);
  sqc.nuse_in_reconstruction_min = s_default_options.secondary_qc.fourth;
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      sqc.scopes[iscope].nimage_min = s_default_options.secondary_qc.first;
      sqc.scopes[iscope].N_min      = s_default_options.secondary_qc.second;
      sqc.scopes[iscope].dist2_max = 
	Angle::frDeg(s_default_options.secondary_qc.third)
	*Angle::frDeg(s_default_options.secondary_qc.third);
    }

  VSAReconstruction recon(array_info, atmo, qc);

  VSAReconstruction::Method method = VSAReconstruction::M_2;
  if(m_options.method==1)method = VSAReconstruction::M_1;
  else if(m_options.method==2)method = VSAReconstruction::M_2;
  else if(m_options.method==3)method = VSAReconstruction::M_3;

  VSAReconstruction::ScopeWeighting weighting = 
    VSAReconstruction::SW_ELLIPTICITY_MEMO;
  if(m_options.weighting=="one")weighting =
    VSAReconstruction::SW_ONE;
  else if(m_options.weighting=="size")weighting =
    VSAReconstruction::SW_SIZE;
  else if(m_options.weighting=="ellipticity")weighting =
    VSAReconstruction::SW_ELLIPTICITY_MEMO;
  else if(m_options.weighting=="size_ellipticity")weighting =
    VSAReconstruction::SW_SIZE_ELLIPTICITY_MEMO;
  else if(m_options.weighting=="ellipticity_traditional")weighting =
    VSAReconstruction::SW_ELLIPTICITY_TRAD;
  else if(m_options.weighting=="width")weighting =
    VSAReconstruction::SW_WIDTH;
  else 
    {
      vstream 
	<< "Weighting scheme \"" << m_options.weighting << "\" unknown."
	<< std::endl;
      vsassert(0);
    }

  // --------------------------------------------------------------------------
  // Create the muon analysis
  // --------------------------------------------------------------------------

  std::vector<std::vector<double> > all_gain(nscope);
  std::vector<VSMuonAnalysis*> muon_analysis;
  if(!m_options.no_muon_analysis)
    {
      muon_analysis.resize(nscope);
      for(unsigned iscope=0; iscope<nscope; iscope++)
  	muon_analysis[iscope] = new 
	  VSMuonAnalysis(array_info.scopes[iscope],
			 Angle::frDeg(m_options.muon_raw_ring_radius),
			 m_options.muon_nimage_cut,
			 Angle::frDeg(m_options.muon_radius_min_cut),
			 Angle::frDeg(m_options.muon_radius_max_cut),
			 Angle::frDeg(m_options.muon_rms_max_cut),
			 Angle::frDeg(m_options.muon_ring_edge_dist_max_cut),
			 m_options.muon_centroid_radius_ratio_max_cut
#ifdef MUON_TEST_NON_UNIFORMITY
			 , m_options.muon_non_uniformity_beta
#endif
			 );
    }

  // --------------------------------------------------------------------------
  // Create the cuts calculator
  // --------------------------------------------------------------------------

  VSCutsCalc* cuts_calc = 0;

  if(!m_options.nspace_cuts.empty())
    {
      cuts_calc = new VSNSpaceCutsCalc;

      std::string name  = m_options.nspace_cuts;
      VSFileUtility::expandFilename(name);      
      VSOctaveH5Reader *reader = new VSOctaveH5Reader(name);

      cuts_calc->load(reader);

      if(!no_verbose)
	vstream << "Loaded nspace cuts: " << name << std::endl;      

      delete reader;
    }

  // --------------------------------------------------------------------------
  // Construct the helpers
  // --------------------------------------------------------------------------

  VSCleaner* clean = 
    VSCleanerFactory::getCleaner(m_options.primary_cleaning_args,
				 camera_nchan, camera_neighbors);

  VSCleaner* sclean = 0;
  if(!m_options.secondary_cleaning_args.empty())
    sclean = VSCleanerFactory::getCleaner(m_options.secondary_cleaning_args,
					  camera_nchan, camera_neighbors);

  VBFAnalysisStage2::CleaningScale cleaning_scale;
  if((m_options.cleaning_scale == "pedrms")
     ||(m_options.cleaning_scale == "pedvar"))
    cleaning_scale=VBFAnalysisStage2::CS_PEDRMS;
  else if(m_options.cleaning_scale == "gain")
    cleaning_scale=VBFAnalysisStage2::CS_GAIN;
  else if(m_options.cleaning_scale == "unity")
    cleaning_scale=VBFAnalysisStage2::CS_UNITY;
  else
    {
      std::cerr << "Unknown cleaning scaling: " << m_options.cleaning_scale
		<< std::endl;
      exit(EXIT_FAILURE);
    }

  VSPointing* pointing = 0;
  VSPointing* my_pointing = 0;
  VSTargetPointing* target_pointing = 0;

  if(stage1->misc_db)
    {
      unsigned ntpd = stage1->misc_db->scope.size();
      std::vector<VSTargetPointing::TargetData> tpd(ntpd);
      for(unsigned itpd=0;itpd<ntpd;itpd++)
	tpd[itpd] = stage1->misc_db->scope[itpd].tracking_targets;
      target_pointing = new VSTargetPointing(tpd,earth_position);
    }

  if(m_options.pointing_string == "positioner")
    if(stage1->db_pointing)
      my_pointing = pointing = new VSDBPointing(*stage1->db_pointing);
    else if(!no_db)
      my_pointing = pointing = new VSDBPointing;
    else 
      {
	std::cerr
	  << "WARN: Pointing records not present in stage 1 file and DB is disabled.\n"
	  << "WARN: Attempting to use interpolated L3 pointing.\n"
	  << std::endl;
	m_options.pointing_string = "l3";
      }
  
  if(m_options.pointing_string == "l3")
    if(stage1->l3_pointing)
      my_pointing = pointing = new VSL3Pointing(*stage1->l3_pointing);
    else 
      {
	std::cerr
	  << "WARN: L3 pointing records are  not present in stage 1 file.\n"
	  << "WARN: Attempting to use L3 records directly for pointing.\n"
	  << std::endl;
	m_options.pointing_string = "l3direct";
      }

  if(m_options.pointing_string == "l3direct")
    if(!no_l3)
      {
	pointing = VSDirectL3Pointing::getInstance();

	if(stage1->l3_pointing)
	  my_pointing = new VSL3Pointing(*stage1->l3_pointing);
	else 
	  {
	    std::cerr
	    << "WARN: L3 pointing records are  not present in stage 1 file.\n"
	    << "WARN: Mean RA/Dec and Az/El values will be wrong.\n"
	    << std::endl;
	  }
      }
    else
      {
	std::cerr
	  << "FATAL: Direct L3 pointing requested but L3 processing has been\n"
	  << "explicityly disabled." << std::endl;
	exit(EXIT_FAILURE);
      }

  if(m_options.pointing_string == "target")
    if(target_pointing)
      pointing = target_pointing;
    else
      {
	std::cerr
	  << "FATAL: Target pointing requested but target is not known." 
	  << std::endl;
	exit(EXIT_FAILURE);
      }

  if(pointing == 0)
    {
      std::cerr << "Unknown pointing information source: "
		<< m_options.pointing_string << std::endl;
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // Get observation and targets
  // --------------------------------------------------------------------------

  SEphem::SphericalCoords mean_radec;
  double rms_ra_rad = 0;
  double rms_dec_rad = 0;
  if(my_pointing)
    my_pointing->getMeanTarget(mean_radec, rms_ra_rad, rms_dec_rad,
			       stage1->run_info.lo_event_time,
			       stage1->run_info.hi_event_time, earth_position);

  SEphem::SphericalCoords mean_azzn;
  double rms_az_rad = 0;
  double rms_zn_rad = 0;
  double mean_zn_rad = 0;
  if(my_pointing)
    my_pointing->getMeanAzZn(mean_azzn, rms_az_rad, rms_zn_rad, mean_zn_rad,
			     stage1->run_info.lo_event_time,
			     stage1->run_info.hi_event_time);

  VSTargetTable target_table(stage1->target_table);
  if(!target_table.empty())target_table.mergeCompiledTargets();

  std::string demand_source_name = "";
  
  for(std::vector<target_type>::const_iterator itarget = 
	m_options.theta_targets.begin(); 
      itarget != m_options.theta_targets.end(); itarget++)
    if((itarget->second=="auto")&&(itarget->first!="auto"))
      demand_source_name = itarget->first;
  
  VSTargetTable::Observation obs_mp;
    
  if(stage1->sim_info == NULL)
    obs_mp = target_table.getObservation(mean_radec, rms_ra_rad, rms_dec_rad,
					 mean_azzn, rms_az_rad, rms_zn_rad,
					 demand_source_name,
					 stage1->run_info.lo_event_time);
  else
    {
      obs_mp.src_radec.setLatLongDeg(stage1->sim_info->src_dec_deg,
				     stage1->sim_info->src_ra_deg);
      obs_mp.obs_radec.setLatLongDeg(stage1->sim_info->obs_dec_deg,
				     stage1->sim_info->obs_ra_deg);
      obs_mp.src_radec_J2000.setLatLongDeg(stage1->sim_info->src_dec_deg,
					   stage1->sim_info->src_ra_deg);
      obs_mp.obs_radec_J2000.setLatLongDeg(stage1->sim_info->obs_dec_deg,
					   stage1->sim_info->obs_ra_deg);
      obs_mp.wobble_theta_rad = stage1->sim_info->wobble_theta_deg*M_PI/180.;
      obs_mp.wobble_phi_rad   = stage1->sim_info->wobble_phi_deg*M_PI/180.;
 
      if(stage1->sim_info->wobble_theta_deg < 0.05)
	{
	  obs_mp.mode = VSTargetTable::Observation::OM_SIM_ON;
	  obs_mp.name = "sim_on";
	  obs_mp.mode_string = "sim_on";
	}
      else
	{
	  obs_mp.mode = VSTargetTable::Observation::OM_SIM_WOBBLE;
	  obs_mp.name = "sim_wobble";
	  obs_mp.mode_string = "sim_wobble";
	}
    }

  VSTargetTable::Observation obs;
  bool obs_verified = false;
  if(target_pointing)
    obs_verified = 
      target_pointing->getObservation(obs, 
				      stage1->run_info.lo_event_time,
				      stage1->run_info.hi_event_time);

  if(m_options.no_commanded_target)obs_verified = false;
  if(!obs_verified)obs = obs_mp;
  
  VSOctaveH5WriterStruct* obs_s = writer->writeStruct("observation");
  obs.save(obs_s);
  delete obs_s;

  VSOctaveH5WriterStruct* obs_mp_s = 
    writer->writeStruct("observation_mean_pointing");
  obs_mp.save(obs_mp_s);
  delete obs_mp_s;

  if(!no_verbose)
    {
      vstream << "Observation: " << obs.name << '/' << obs.mode_string;
      if(obs_verified)vstream << " (verified)";
      vstream << std::endl << std::endl;
    }

  std::vector<VBFAnalysisStage2::Coords> coords = 
    getTargets(m_options.theta_targets, obs, target_table, 
	       stage1->run_info.lo_event_time, m_options.wobble_region_radius);

  if((!coords.empty())&&(!no_verbose))
    {
      vstream << "Targets:" << std::endl;
      for(unsigned itarget=0;itarget<coords.size();itarget++)
	{
	  vstream << "  " << coords[itarget].name << " ";
	  if(coords[itarget].name.length()<35)
	    vstream << std::string(35-coords[itarget].name.length(),' ');
	  vstream << coords[itarget].coord.longitude().hmsString(1) << ' '
		  << coords[itarget].coord.latitude().dmsString(1)
		  << std::endl;
	}
      vstream << std::endl;
    }

  // --------------------------------------------------------------------------
  // Mean scaled parameter and energy lookup tables
  // --------------------------------------------------------------------------
  VSScaledParameterCalc sp_calc;
  VSEnergyCalc energy_calc;
  
  // --------------------------------------------------------------------------
  // Settings
  // --------------------------------------------------------------------------

  VBFAnalysisStage2::Settings settings;

  settings.earth_position               = earth_position;
  settings.process_all_events           = m_options.all_events;
  settings.enable_l2_corrections        = m_options.enable_l2_corrections;
  settings.method                       = method;
  settings.weighting                    = weighting;
  settings.theta_targets                = coords;
  settings.theta_cut                    = m_options.theta_cut;
  settings.nimage_cut                   = m_options.nimage_cut;
  settings.nscope_cut                   = m_options.nscope_cut;
  settings.do_not_write_diagnostics     = m_options.no_diagnostics;
  settings.no_slow_diagnostics          = m_options.no_slow_diagnostics;
  settings.l2_trigger_threshold         = m_options.l2_trigger_threshold;
  settings.l3_trigger_threshold         = m_options.l3_trigger_threshold;
  settings.print_frequency              = m_options.print_frequency;
  settings.cleaning_scale               = cleaning_scale;
  settings.ped_event_count_min          = m_options.ped_event_count_min;
  settings.permissive_laser             = m_options.permissive_laser;
  settings.unity_gain                   = m_options.unity_gains;
  settings.auto_suppress_l2_chan        = m_options.auto_suppress_l2_chan;
  settings.pad_zero_suppressed_chan     = m_options.pad_zero_suppressed_chan;
  settings.chan_suppress                = m_options.chan_suppress;
  settings.chan_has_no_pmt              = m_options.chan_no_pmt;
  settings.scope_suppress               .resize(nscope,false);
  for(unsigned isuppress=0;isuppress<m_options.scope_suppress.size();
      isuppress++)
    if(m_options.scope_suppress[isuppress]<nscope)
      settings.scope_suppress[m_options.scope_suppress[isuppress]] = true;
  settings.scope_gain                   = m_options.scope_gain;
  settings.scope_dev_scaling            = m_options.scope_dev_scaling;
  settings.scope_hi_lo_gain_ratio       = m_options.scope_hi_lo_gain_ratio;
  settings.integration_lo_zero_sample   = m_options.integration_lo_zero_sample;
  settings.integration_hi_zero_sample   = m_options.integration_hi_zero_sample;
  settings.integration_threshold_frac   = m_options.integration_threshold_frac;
  settings.integration_window_start     = m_options.integration_window_start;
  settings.integration_window_width     = m_options.integration_window_width;
  settings.integration_apply_increase   = m_options.integration_apply_increase;
  settings.integration_threshold_charge 
                            = m_options.integration_threshold_charge;
  settings.integration_window_start_increase 
                            = m_options.integration_window_start_increase;
  settings.integration_window_width_increase
                            = m_options.integration_window_width_increase;
  settings.neighbor_nchan               = camera_nchan;
  settings.neighbors                    = camera_neighbors;
  settings.verbose                      = !no_verbose;
  settings.software_trigger_masks       .clear();
  settings.qc_masks                     .clear();
  settings.limited_dt_range             = m_options.limited_dt_range;
    
  settings.psf_poly_radial              .set(m_options.psf_poly_radial);
  settings.psf_poly_tangential          .set(m_options.psf_poly_tangential);
  
  calculateMasks(settings.software_trigger_masks,
		 m_options.software_trigger_masks);
  calculateMasks(settings.qc_masks, m_options.qc_masks);

  // --------------------------------------------------------------------------
  // Construct the analysis class
  // --------------------------------------------------------------------------

  VBFAnalysisStage2::SecondaryCleaning secondary_clean(sclean, &sqc);

  uint32_t dispatcher_flags = VSSimpleVBFDispatcher::DISPATCH_NONE;

#ifndef NOTHREADS
  if((!no_threads)&&(!m_options.no_vbf_reader_thread))
    dispatcher_flags |= VSSimpleVBFDispatcher::DISPATCH_THREADED;
  if((!no_threads)&&(m_options.nthreads))dispatcher_flags |= 
    VSSimpleVBFDispatcher::DISPATCH_NTHREAD_SET(m_options.nthreads);
#endif
  if(do_use_overflow)
    dispatcher_flags |= VSSimpleVBFDispatcher::DISPATCH_INTEGRATE_OVERFLOW;

  VBFAnalysisStage2* analysis = 
    new VBFAnalysisStage2(clean, secondary_clean,
			  pointing, &recon, muon_analysis, 
			  &sp_calc, &energy_calc, cuts_calc,
			  *stage1, 
			  VSChannelMap(stage1->run_info.lo_event_time),
			  settings, writer,
			  pad_stage1);

  unsigned nevent = 0;
  try 
    {
      VBFStandardRationalizedClock clock(no_l3, 
					 stage1->run_info.lo_event_time,
					 stage1->run_info.hi_event_time);
      VSSimpleVBFDispatcher dispatcher(analysis);
      dispatcher.resetClock(&clock);
      dispatcher.openFile(vbf_filename.c_str());

      vstream << "Stage2 event loop starting - " << dispatcher.numPackets() 
	      << " packets" << std::endl;

      VSSimpleVBFDispatcher::catchSignalAndStopProcessingFile(SIGINT);
      nevent = dispatcher.dispatchAllPackets(m_options.npackets,
					     dispatcher_flags);

      vstream << "Stage2 event loop complete - " << nevent << " packets" 
	      << std::endl << std::endl;
    }
  catch (const std::exception& e)
    {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
    }

  // --------------------------------------------------------------------------
  // Clean up
  // --------------------------------------------------------------------------

  delete analysis;
  delete clean;
  delete my_pointing;
  delete target_pointing;
  delete cuts_calc;
  for(unsigned iscope=0;iscope<muon_analysis.size();iscope++)
    delete muon_analysis[iscope];

  // --------------------------------------------------------------------------
  // Stop the time counter and print out the execution time
  // --------------------------------------------------------------------------

  writer->flush();
  VSOctaveH5WriterStruct* analyze_struct = writer->writeStruct("analyze");  
  analysis_data.stop();
  analysis_data.save(analyze_struct);
  delete analyze_struct;
  analysis_data.printTiming(vstream,vbf_filename,nevent);
}

#define OPTNAME(x,y) std::string(x)+std::string(y)

void VSAnalysisStage2::configure(VSOptions& options,
				 const std::string& profile,
				 const std::string& opt_prefix)
{
  if(profile == "analysis")
    {

    }
  else if(profile == "diagnostics")
    {
      s_default_options.permissive_laser                   = true;
      s_default_options.nscope_cut                         = 1;
      //s_default_options.ped_suppress_mode                  = "none";
    }
  else if(profile == "simulations")
    {
      s_default_options.no_muon_analysis                   = true;
      s_default_options.pointing_string                    = "l3direct";
    }

  options.addCatagory("s2_int",
		      "Stage 2 trace integration options. Configure the "
		      "parameters of the trace integration algorithm, such "
		      "as window width, start point and growth.");

  options.addCatagory("s2_sup","Stage 2 channel suppression options.");

  options.addCatagory("s2_gain","Stage 2 gain options.");

  options.addCatagory("s2_misc", "Stage 2 miscellaneous options.");

  options.addCatagory("s2_diag", "Stage 2 diagnostics options.");

  options.addCatagory("s2_clean", 
		      "Stage 2 event processing, event cleaning and quality "
		      "cuts.");

  options.addCatagory("s2_recon",
		      "Stage 2 array reconstruction options");

  options.addCatagory("s2_tar",
		      "Stage 2 options controlling selection of targets");

  options.addCatagory("s2_point",
		      "Stage 2 options selecting and configuring pointing "
		      "determination algorithm.");

  options.addCatagory("s2_muon",
		      "Stage 2 muon analysis options.");

  options.addCatagory("s2_shape",
		      "Stage 2 gamma hadron separation.");

  options.addCatagory("s2_energy",
		      "Stage 2 energy reconstruction.");

  options.findWithValue(OPTNAME(opt_prefix,"hi_gain_zero_sample"),
			s_default_options.integration_hi_zero_sample,
			"Set the minimum sample number on each telescope from "
			"which HIGH gain channels can be integrated. This can "
			"usually be set at ZERO.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"lo_gain_zero_sample"),
			s_default_options.integration_lo_zero_sample,
			"Set the minimum sample number on each telescope from "
			"from which LOW gain channels can be integrated. This "
			"should usually be set to a value which excludes the "
			"remnants of the HIGH gain signal which may remain at "
			"the edge of the integration window.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_start"),
			s_default_options.integration_window_start,
			"Sample number to integrate from. The window start "
			"is defined relative to T-zero, the point at which "
			"the trace reaches some fraction of the peak height. "
			"Positive values means the integration starts before "
			"the T-zero point.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_width"),
			s_default_options.integration_window_width,
			"Number of samples to integrate when calculating "
			"charge in each channel.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"peak_threshold_fraction"),
			s_default_options.integration_threshold_frac,
			"Fraction of peak of trace at which T0 is defined.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"apply_window_increase"),
			s_default_options.integration_apply_increase,
			"Apply increase in window size for large pulses.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_increase_threshold"),
			s_default_options.integration_threshold_charge,
			"Total charge in channel above which the integration "
			"window is increased. The window size is increased "
			"and shifted earlier in the trace logarithmically "
			"with the charge above this threshold.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_start_increase"),
			s_default_options.integration_window_start_increase,
			"Amount of advance to apply to the start of the "
			"integration window (in samples) per decade that the "
			"total charge exceeds the threshold.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"window_width_increase"),
			s_default_options.integration_window_width_increase,
			"Amount to increase the width of the integration "
			"window (in samples) per decade that the total charge "
			"exceeds the threshold.",
			"s2_int");

  options.findWithValue(OPTNAME(opt_prefix,"ped_suppress_mode"),
			s_default_options.ped_suppress_mode,
			"Set the pedestal suppression mode, should be either "
			"\"median\", \"interval\" or \"none\".",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"ped_suppress_lo"),
			s_default_options.ped_suppress_lo,
			"Fraction of median pedestal variance to below which "
			"channel is suppressed. This option applies to "
			"\"median\" mode only.",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"ped_suppress_hi"), 
			s_default_options.ped_suppress_hi,
			"Fraction of median pedestal variance to above which "
			"channel is suppressed. This option applies to "
			"\"median\" mode only.",
			"s2_sup");
  
  options.findWithValue(OPTNAME(opt_prefix,"ped_suppress_fraction"),
			s_default_options.ped_suppress_fraction,
			"Fraction of channels to keep unsuppressed before "
			"scaling of interval is applied. This option applies "
			"to \"interval\" mode only\".",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"ped_suppress_scale"), 
			s_default_options.ped_suppress_scale,
			"Scaling to apply to interval which contains desired "
			"fraction of channels. A value larger than 1.0 makes "
			"the interval larger, retaining more channels. A "
			"value less than 1.0 decreases the interval, "
			"suppressing more channels. This option applies "
			"to \"interval\" mode only\".",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"ped_event_count_min"), 
			s_default_options.ped_event_count_min,
			"Set the minimum number of pedestal events required "
			"per time slice. Any channels which don't have enough "
			"events in a slice are suppressed.",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"gain_suppress_lo"),
			s_default_options.gain_suppress_lo,
			"Inverse channel gain, below which channel is "
			"suppressed.",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"gain_suppress_hi"),
			s_default_options.gain_suppress_hi,
			"Inverse channel gain, above which channel is "
			"suppressed.",
			"s2_sup");

  options.findBoolValue(OPTNAME(opt_prefix,"no_laser"), 
			s_default_options.no_laser, true,
			"Do not use the LASER data for L2 timing, channel "
			"suppression or relative gain equalization.",
			"s2_gain");
  
  options.findBoolValue(OPTNAME(opt_prefix,"permissive_laser"), 
			s_default_options.permissive_laser, true,
			"Do not die if a telescope is missing LASER data, "
			"instead use unity for all LASER quantities.",
			"s2_gain");

  options.findBoolValue(OPTNAME(opt_prefix,"unity_gains"), 
			s_default_options.unity_gains, true,
			"Override gains from laser file with 1.0 for every "
			"channel. Keep timing information.",
			"s2_gain");
  
  options.findWithValue(OPTNAME(opt_prefix,"scope_gain_multiplier"),
			s_default_options.scope_gain,
			"Set the gain multiplier for each telescope.",
			"s2_gain");

  options.findWithValue(OPTNAME(opt_prefix,"scope_hi_lo_gain_ratio"),
			s_default_options.scope_hi_lo_gain_ratio,
			"Set the hi/lo gain ratio for each telescope.  This "
			"will be the default ratio assumed for any channel "
			"that does not have a defined gain ratio in the hilo "
			"calibration file.",
			"s2_gain");

  options.findWithValue(OPTNAME(opt_prefix,"scope_dev_scaling"),
			s_default_options.scope_dev_scaling,
			"Set the coefficients of the linear scaling between "
			"telescope pedestal RMS and the mean scaled RMS used "
			"to estimate the NSB",
			"s2_misc");

  options.findBoolValue(OPTNAME(opt_prefix,"no_l2_corrections"), 
			s_default_options.enable_l2_corrections, false,
			"Disable crate timing corrections using the "
			"L2 pulses.",
			"s2_misc");

  options.findBoolValue(OPTNAME(opt_prefix,"pad_zero_suppressed_chan"), 
			s_default_options.pad_zero_suppressed_chan, false,
			"Inject a random signal into zero suppressed channels "
			"using the pedestal rms.  This option should be used "
			"when running on simulations in which pixels "
			"with no cherenkov photons are zero suppressed.  It "
			"is not recommended to use this option with "
			"real zero-suppressed data.",
			"s2_misc");

  options.findWithValue(OPTNAME(opt_prefix,"cleaning_scale"),
			s_default_options.cleaning_scale,
			"Set the scale of the cleaning algorithm. Valid "
			"settings are \"pedrms\", \"gain\", and \"unity\".",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix,"cleaning"),
			s_default_options.primary_cleaning_args,
			"Configure the primary cleaning algorithm. Valid "
			"forms are \"picbnd,piclevel,bndlevel\" or "
			"\"regional,rlevel,rsize,ilevel\".",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix,"secondary_cleaning"),
			s_default_options.secondary_cleaning_args,
			"Configure the secondary cleaning algorithm. "
			"Secondary cleaning is applied to calculation of "
			"gamma/hadron separation parameters. Valid choices "
			"are as in previous option.",
			"s2_clean");

  options.findBoolValue(OPTNAME(opt_prefix,"all_events"), 
			s_default_options.all_events, true,
			"Do not restrict processing to only L2 events.",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix,"nimage_cut"),
			s_default_options.nimage_cut,
			"Set the minimum number of channels required after "
			"cleaning for a telescope image to be considered in "
			"the reconstruction.",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix,"nscope_cut"), 
			s_default_options.nscope_cut,
			"Set the minimum number of telescope images required "
			"to proceed with reconstruction.",
			"s2_clean");

  options.findBoolValue(OPTNAME(opt_prefix,"no_diagnostics"), 
			s_default_options.no_diagnostics, true,
			"Disable writing of diagnostic information to output "
			"file. Not recommended if you plan to run stage3 "
			"on the output file.",
			"s2_diag");

  options.findBoolValue(OPTNAME(opt_prefix,"no_slow_diagnostics"), 
			s_default_options.no_slow_diagnostics, true,
			"Disable calculation of the \"slow\" diagnostics. "
			"This does not disable the calculation of diagnostic "
			"values needed by stage3.",
			"s2_diag");
  
  options.findWithValue(OPTNAME(opt_prefix,"l2_trigger_threshold"),
			s_default_options.l2_trigger_threshold,
			"Set number of channels required for L2 trigger. "
			"This is used to flag channels which seem to cause "
			"L2 triggers with too few channels.",
			"s2_diag");

  options.findWithValue(OPTNAME(opt_prefix,"l3_trigger_threshold"),
			s_default_options.l3_trigger_threshold,
			"Set number of telescope required for L3 trigger. "
			"This is used to flag events which have too few "
			"telescope events. These events are counted against "
			"the dead time.",
			"s2_diag");

  options.findWithValue(OPTNAME(opt_prefix,"dt_fit_range"),
			s_default_options.limited_dt_range,
			"Set range of Delta-T used to estimate live time. "
			"Minimum and maximum time differences should be given "
			"in seconds.",
			"s2_diag");

  options.findWithValue(OPTNAME(opt_prefix,"scope_suppress"),
			s_default_options.scope_suppress,
			"Suppress telescope for all events in the run.",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"chan_suppress"),
			s_default_options.chan_suppress,
			"Suppress telescope channel for all events in the "
			"run. Specify as \"scope/channel\".",
			"s2_sup");
  
  options.findBoolValue(OPTNAME(opt_prefix,"no_auto_suppress_l2_chan"),
			s_default_options.auto_suppress_l2_chan, false,
			"Do not automatically treat L2 channels as channels "
			"which do not have a valid PMT signal.",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"chan_no_pmt"),
			s_default_options.chan_no_pmt,
			"Specify that telescope channel does not have a PMT. "
			"This option suppresses the channels for all events "
			"in the run, but differs from \"chan_suppress\" by "
			"how the channels are handled in the diagnostics. "
			"Specify as \"scope/channel\".",
			"s2_sup");

  options.findWithValue(OPTNAME(opt_prefix,"method"),
			s_default_options.method,
			"Set the reconstruction method, should be 1,2 or 3.",
			"s2_recon");
  
  options.findWithValue(OPTNAME(opt_prefix,"weighting"),
			s_default_options.weighting,
			"Set the reconstruction weighting, should be \"one\", "
			"\"size\", \"ellipticity\", \"size_ellipticity\", "
			"\"ellipticity_traditional\", or \"width\".",
			"s2_recon");

  options.findWithValue(OPTNAME(opt_prefix,"scope_pos"),
			s_default_options.scope_pos,
			"Set the position for the telescopes relative to the "
			"array. A comma separated list of coordinate triples "
			"should be given, specifying the EW/NS/UD position of "
			"each telescope in meters.  Any telescope positions "
			"left undefined will be set to the known telescope "
			"position on the basis of the run start date.",
			"s2_recon");

  options.findBoolValue(OPTNAME(opt_prefix,"no_set_scope_pos_from_sims"),
			s_default_options.no_set_scope_pos_from_simulations,
			true,
			"Do not set the telescope positions from the values "
			"stored in the simulated VBF file.",
			"s2_recon");

  options.findWithValue(OPTNAME(opt_prefix,"cam_rotation"),
			s_default_options.cam_rotation,
			"Set the rotation of each camera plane. Must be "
			"specified as a comma separated list of angles in "
			"degrees.",
			"s2_recon");

  options.findWithValue(OPTNAME(opt_prefix,"atmosphere"),
			s_default_options.atmosphere,
			"Set file name of CORSIKA atmosphere model to load. "
			"If the file name is empty then the US standard model "
			"is used.",
			"s2_recon");

  options.findWithValue(OPTNAME(opt_prefix,"theta_targets"),
			s_default_options.theta_targets,
			"Specify list of targets for which to calculate theta "
			"with respect to. Targets should be specified as a "
			"comma separated list of string pairs. The first "
			"element of the pair specifies the target name or "
			"coordinates, the second specifies the offset. "
			"Use the \"-list_targets\" option to get a list of "
			"target names and coordinate formats. The offset "
			"parameter can have the value \"on\", \"off\" or "
			"be in the form \"theta@phi\" to specify wobbles "
			"from the target position. Be careful specifying "
			"wobbles, if an observation was taken in wobble mode "
			"of 0.3 deg North, then you should specify 0.6@0.0 "
			"to find the correct off-source location.",
			"s2_tar");

  options.findWithValue(OPTNAME(opt_prefix,"wobble_off_region_radius"),
			s_default_options.wobble_region_radius,
			"Set the size of the off region required when "
			"calculating the maximum number of off regions than "
			"fit a given wobble offset.",
			"s2_tar");

  options.findWithValue(OPTNAME(opt_prefix,"theta_cut"),
			s_default_options.theta_cut,
			"Specify a cut on theta, above which the event is not "
			"saved to the file. Useful for cutting down on file "
			"size while optimizing parameters. Specifying a value "
			"<=0 disables the cut.",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix,"pointing"),
			s_default_options.pointing_string,
			"Set the pointing information source. Valid options "
			"are (1) \"positioner\", use the database entries "
			"written by "
			"the positioner code, (2) \"l3\", interpolate from "
			"the tracking entries in VBF written by L3, (4) "
			"\"l3_direct\", use the VBF L3 entries directly, "
			"without interpolation, or (4) \"target\", assume the "
			"telescopes are accurately tracking an astronomical "
			"target.",
			"s2_point");

  options.findWithValue(OPTNAME(opt_prefix,"tracking_recorrections_file"),
			s_default_options.tracking_recorrections_file,
			"Reapply a tracking model to the raw tracking "
			"records. This option allows a new set of corrections "
			"to be applied. The corrections are loaded from the "
			"given file. Note, for this option to have any "
			"effect, you must select the \"positioner\" pointing "
			"source and the tracking records must have been saved "
			"in the stage 1 file. This option should specify "
			"a comma separated list of telescopes to correct and "
			"the files which hold the corrections for each, for "
			"example: \"0/corr_T1.dat,1/corr_T2.dat\".",
			"s2_point");

  options.findWithValue(OPTNAME(opt_prefix,"tracking_recorrections_date"),
			s_default_options.tracking_recorrections_date,
			"Reapply a tracking model to the raw tracking "
			"records. This option allows a new set of corrections "
			"to be applied. The corrections are loaded from the "
			"database by the given date. Note, for this option "
			"to have any effect, you must select the "
			"\"positioner\" pointing source and the tracking "
			"records must have been saved in the stage 1 file. "
			"This option should specify a comma separated list "
			"of (multiple) telescopes to correct and the date for "
			"the corrections for each. For example: "
			"\"01/2007-10-10,23/2007-11-01\" would load T1,T2 "
			"corrections for the Oct 10, 2007 and T3,T4 from "
			"Nov 1. As a simplification \"*/2007-10-10\" will "
			"load all telescopes for the given date, and "
			"specifying the date as \"next\" will load the next"
			"set of corrections stored DB from the date of the "
			"run.",
			"s2_point");

  options.findWithValue(OPTNAME(opt_prefix,"demand_tracking_target"),
			s_default_options.tracking_target,
			"Set the name of the target to assume when the "
			"pointing mode is set to \"target\". If not set "
			"the tracking target records set by the positioners "
			"are used.",
			"s2_tar");

  options.findBoolValue(OPTNAME(opt_prefix,"no_commanded_target"),
			s_default_options.no_commanded_target, true,
			"Do not use the commanded target from the positioner "
			"to determine the observation.",
			"s2_tar");

  options.findBoolValue(OPTNAME(opt_prefix,"no_muon_analysis"),
			s_default_options.no_muon_analysis, true,
			"Disable muon analysis.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix,"muon_raw_ring_radius"),
			s_default_options.muon_raw_ring_radius,
			"Set the half-width of the ring used in the raw "
			"muon analysis. All channels within this distance of "
			"the ring will be included in the raw analysis. A "
			"value of zero disables the raw analysis.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix,"muon_nimage_cut"),
			s_default_options.muon_nimage_cut,
			"Set the minimum number of channels in each image "
			"before attempting muon ring fitting.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix,"muon_radius_min_cut"),
			s_default_options.muon_radius_min_cut,
			"Set the lower radius a fitted ring can have "
			"for its information will to be saved.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix,"muon_radius_max_cut"),
			s_default_options.muon_radius_max_cut,
			"Set the higher radius a fitted ring can have "
			"for its information will to be saved.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix,"muon_width_max_cut"),
			s_default_options.muon_rms_max_cut,
			"Set the maximum RMS width a ring can have "
			"for its information will to be saved.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix,"muon_ring_edge_dist_max_cut"),
			s_default_options.muon_ring_edge_dist_max_cut,
			"Set the maximum distance between the center of the "
			"camera to the outer edge of the fitted muon ring.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix,"muon_non_uniformity"),
			s_default_options.muon_non_uniformity_beta,
			"Set the non uniformity beta parameter.",
			"s2_muon");
  
  options.findWithValue(OPTNAME(opt_prefix,
				"muon_centroid_radius_ratio_max_cut"),
			s_default_options.muon_centroid_radius_ratio_max_cut,
			"Set the maximum ratio of the image centroid with "
			"respect to the ring center to the ring radius.",
			"s2_muon");

  options.findWithValue(OPTNAME(opt_prefix, "qc"),
			s_default_options.primary_qc,
			"Primary quality cuts: set the minimum number of "
			"image channels, minimum size, maximum distance from "
			"the center of the camera to the centroid, and "
			"minimum number of telescopes passing cuts. "
			"Cuts should be specified as: "
			"\"nimage/size/dist/nscope\", for example "
			"\"4/200/1.5/2\".",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix, "qc_masks"),
			s_default_options.qc_masks,
			"Set a mask of telescopes which must pass the "
			"primary quality cuts for before the event is written "
			"to the output file. For "
			"example 01,02,12,13,23 will suppress events in which "
			"only the T1/T4 pair is triggered, allowing all other "
			"combinations of 2 telescopes and any 3 or 4 "
			"telescope trigger.",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix, "secondary_qc"),
			s_default_options.secondary_qc,
			"Secondary quality cuts, used only in conjunction "
			"with secondary cleaning for calculation of "
			"gamma/hadron parameters. Has no effect if secondary "
			"cleaning is not defined. See previous option for "
			"details of how to define the cuts.",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix,"nspace_cuts"),
			s_default_options.nspace_cuts,
			"Specify the name of the file containing the set of "
			"array and/or scope-level nspace filters which will "
			"be applied to the data.  Events are marked as either "
			"passing or failing this cut if a space is provided.",
			"s2_shape");

  options.findWithValue(OPTNAME(opt_prefix,"print_frequency"),
			s_default_options.print_frequency,
			"Set the frequency of printing event numbers.",
			"s2_misc");

  options.findBoolValue(OPTNAME(opt_prefix,"no_vbf_reader_thread"),
			s_default_options.no_vbf_reader_thread, true,
			"Do not use a separate thread to read and decode VBF "
			"packets.",
			"s2_misc");

#ifndef NOTHREADS
  options.findWithValue(OPTNAME(opt_prefix,"nthreads"),
			s_default_options.nthreads, 
			"Set the number of threads to use in the analysis.",
			"s2_misc");
#endif
  
  options.findWithValue(OPTNAME(opt_prefix,"npackets"),
			s_default_options.npackets, 
			"Set the maximum number of packets to process.",
			"s2_misc");

  options.findWithValue(OPTNAME(opt_prefix,"software_trigger_masks"),
			s_default_options.software_trigger_masks, 
			"Set software trigger mask apply to each event. For "
			"example 01,02,12,13,23 will suppress events in which "
			"only the T1/T4 pair is triggered, allowing all other "
			"combinations of 2 telescopes and any 3 or 4 "
			"telescope trigger.",
			"s2_clean");

  options.findWithValue(OPTNAME(opt_prefix,"psf_poly_radial"),
			s_default_options.psf_poly_radial,
			"Set the polynomial of the radial PSF function. The "
			"polynomial gives the PSF as a function of off-axis "
			"distance squared.",
			"s2_shape");

  options.findWithValue(OPTNAME(opt_prefix,"psf_poly_tangential"),
			s_default_options.psf_poly_tangential,
			"Set the polynomial of the tangential PSF function. "
			"The polynomial gives the PSF as a function of "
			"off-axis distance squared.",
			"s2_shape");

  VSEnergyCalcLT::configure(options,profile,opt_prefix);
  VSScaledParameterCalc::configure(options,profile,opt_prefix);
}

std::vector<VBFAnalysisStage2::Coords> VSAnalysisStage2::
getTargets(std::vector<target_type> targets, 
	   const VSTargetTable::Observation& observation,
	   const VSTargetTable& target_table, const VSTime& approx_time,
	   double wobble_region_radius)
{
  std::vector<VBFAnalysisStage2::Coords> coords;
  
  bool has_auto = false;

  // Go through the list of inputs and compute targets ------------------------

  for(std::vector<target_type>::const_iterator itarget = targets.begin();
      itarget != targets.end(); itarget++)
    if(itarget->second=="auto")
      {
	if(!has_auto)
	  {
	    std::vector<VSTargetTable::Target> auto_targets =
	      target_table.getTargetForObservation(observation, 
						   wobble_region_radius);
	
	    for(std::vector<VSTargetTable::Target>::iterator itarget = 
		  auto_targets.begin(); itarget!=auto_targets.end(); itarget++)
	      coords.push_back(*itarget);

	    has_auto = true;
	  }
      }
    else
      {
	VBFAnalysisStage2::Coords c;
	c.name = itarget->first;
	if(!itarget->second.empty())
	  {
	    c.name += std::string("/");
	    c.name += itarget->second;
	  }
	
	c.coord = target_table.getTarget(itarget->first, itarget->second, 
					 approx_time);

	coords.push_back(c);
      }

  return coords;
}
