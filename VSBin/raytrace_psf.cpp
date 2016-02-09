//-*-mode:c++; mode:font-lock;-*-

/*! \file raytrace_psf.cpp

  Calculate the PSF of an optics configuration

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \date    07/09/2005
  \version 1.0
  \note
*/

#include <iostream>
#include <string>
#include <cstdlib>

#include <Vec3D.hpp>
#include <Vec4D.hpp>

#include <VSOptions.hpp>
#include <VSDBFactory.hpp>
#include <VSDBParameterTable.hpp>
#include <VSOTelescopeArray.hpp>
#include <VSORayTracer.hpp>
#include <VSSimpleStat.hpp>

using namespace VERITAS;
using namespace Physics;

enum DataSource { DS_NONE, DS_INIFILE, DS_DUMP, DS_DB };

int main(int argc, char ** argv)
{
  std::string progname(*argv);
  VSOptions options(argc,argv,true);

  // --------------------------------------------------------------------------
  // Options
  // --------------------------------------------------------------------------

  int exit_value = EXIT_SUCCESS;
  bool print_usage = false;
  if(options.find("h","Print this help message")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this help message")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  bool dump_rays = false;
  options.findBoolValue("dump_rays",dump_rays,true,
			"Dump the coordinates of each ray to screen "
			"[rad_x, rad_y, deg_x, deg_y].");

  double direction_theta = 0;
  double direction_phi = 0;
  options.findWithValue("theta",direction_theta,"Set the off-axis angle");
  options.findWithValue("phi",direction_phi,
			"Set the off-axis azimuthal angle");

  DataSource ds_from = DS_NONE;
  std::string ds_from_name = "array.ini";
  if(options.findWithValue("from_ini", ds_from_name,
			   "Create random optics instance from specified "
			   "parameter (INI) file. If the parameter file "
			   "name is given as an empty string, the "
			   "canonical values are used.") != 
     VSOptions::FS_NOT_FOUND)ds_from = DS_INIFILE;

  if(options.findWithValue("from_dump", ds_from_name,
			   "Load optics instance from dump file") != 
     VSOptions::FS_NOT_FOUND)ds_from = DS_DUMP;

  if(options.findWithValue("from_db", ds_from_name,
			   "Load optics instance from database. Argument "
			   "should be in the form of database[/optics_id]") != 
     VSOptions::FS_NOT_FOUND)ds_from = DS_DB;

  VSOArrayParameters::setCanonicalValuesFromOptions(options);
  VSDBFactory::configure(&options);

  std::string rngstatefile(RandomNumbers::defaultFilename());
  options.findWithValue("rng_state_file",rngstatefile,
			"Set random number generator state file");
  
  unsigned num_scope = 1;
  unsigned num_mirror = 0;
  unsigned num_pixel = 0;

  double emission_start    = -2.0;
  double emission_end      = -2.0;
  double emission_radius_i = 0.0;
  double emission_radius_o = 0.55;
  double emission_angle_i  = 0.0;
  double emission_angle_o  = 0.0;

  unsigned num_photons     = 100; // x 1000 photons

  double elevation = 90.0002;
  double azimuth = 0;

  options.findWithValue("scope",num_scope);
  options.findWithValue("mirror",num_mirror);
  options.findWithValue("pixel",num_pixel);

  options.findWithValue("n",num_photons,
			"Number of rays that will be traced in thousands.");
  options.findWithValue("emission_start",emission_start);
  options.findWithValue("emission_end",emission_end);
  options.findWithValue("emission_radius_i",emission_radius_i);
  options.findWithValue("emission_radius_o",emission_radius_o);
  options.findWithValue("emission_angle_i",emission_angle_i);
  options.findWithValue("emission_angle_o",emission_angle_o);

  options.findWithValue("elevation",elevation);
  options.findWithValue("azimuth",azimuth);

  // --------------------------------------------------------------------------
  // Connect to the database if required
  // --------------------------------------------------------------------------
  
  VSDatabase* db = 0;
  if(ds_from == DS_DB)
    db = VSDBFactory::getInstance()->createVSDB();

  // --------------------------------------------------------------------------
  // Create random number generator
  // --------------------------------------------------------------------------

  RandomNumbers rng(rngstatefile.c_str());

  // --------------------------------------------------------------------------
  // Read/Create the array and INI parameters
  // --------------------------------------------------------------------------

  VSOArrayParameters param;
  VSOTelescopeArray array;
  bool have_param = false;

  if(ds_from == DS_INIFILE)
    {
      if(ds_from_name.empty()) 
	{
	  param.reset(true); 
	  have_param=true; 
	}
      else 
	{
	  have_param=param.readFromArrayINIFile(ds_from_name);
	  if(!have_param)
	    {
	      std::cerr << "Could not load array parameters (INI) file " 
			<< ds_from_name << std::endl;
	      return EXIT_FAILURE;
	    }
	}
      array.generateFromArrayParameters(param, rng);
    }
  else if(ds_from == DS_DUMP)
    {
      array.readFromShortDump(ds_from_name);
    }
  else if(ds_from == DS_DB)
    {
      unsigned optics_id = 0;
      std::string database = ds_from_name;
      std::string::size_type ichar = ds_from_name.find('/');
      if((ichar != std::string::npos)
	 &&(ichar < ds_from_name.size()-1))
	{
	  database = ds_from_name.substr(0,ichar);
	  VSDataConverter::fromString(optics_id,ds_from_name.substr(ichar+1));
	}
									
      db->useDatabase(database);
      have_param = param.readFromDatabase(db,optics_id);
      array.readFromDatabase(db,optics_id);
    }

  // --------------------------------------------------------------------------
  // Set up some variables
  // --------------------------------------------------------------------------
  Vec3D beam_offset;
  options.findWithValue("beam_offset",beam_offset);

  VSOTelescope* scope = array.telescopeByHexID(num_scope);
  if(scope == 0)
    {
      std::cerr << "Array does not contain telescope " << num_scope 
		<< std::endl;
      std::cerr << "Usage: " << progname 
		<< " [options]" << std::endl
		<< std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return EXIT_FAILURE;
    }

  Vec3D beam_base;
  if((num_mirror != 0)&&(scope->mirrorByHexID(num_mirror)!=0))
    beam_base = scope->mirrorByHexID(num_mirror)->pos();

  Vec3D beam_direction(0,-1,0);
  if((num_pixel != 0)&&(scope->pixelByHexID(num_pixel)!=0))
    {
      beam_direction = 
	(scope->pixelByHexID(num_pixel)->pos()+scope->focalPlanePosition());
      beam_direction /= beam_direction.Norm();
      beam_direction.y = -beam_direction.y;
    }

  scope->pointTelescopeAzEl(azimuth * M_PI/180., elevation * M_PI/180.);

  // --------------------------------------------------------------------------
  // ALL KNOWN OPTIONS HAVE BEEN PROCESSED
  // --------------------------------------------------------------------------

  argc--,argv++;

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": Unknown command line options:";
      for(int iopt=0;iopt<argc;iopt++)
	if(*argv[iopt]=='-')std::cerr << ' ' << argv[iopt];
      std::cerr << std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;
    }
  else if(argc!=0)
    {
      std::cerr << progname 
		<< ": unwanted arguments required (found " << argc << ")"
		<< std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;
    }  

  if(ds_from == DS_NONE)
    {
      std::cerr << progname 
		<< ": must specify one of \"-from_ini\", \"-from_dump\" and "
		<< "\"-from_db\"" << std::endl << std::endl;
      print_usage=true;
      exit_value=EXIT_FAILURE;      
    }

  if(print_usage)
    {
      std::cerr << "Usage: " << progname 
		<< " [options]" << std::endl
		<< std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return exit_value;
    }

  // --------------------------------------------------------------------------
  // Rescale variables
  // --------------------------------------------------------------------------

  num_photons       *= 1000;

  emission_start    *= scope->curvatureRadius();
  emission_end      *= scope->curvatureRadius();
  emission_radius_i *= scope->aperture();
  emission_radius_o *= scope->aperture();
  emission_angle_i  *= M_PI/180.0;
  emission_angle_o  *= M_PI/180.0;

  direction_theta   *= M_PI/180.0;
  direction_phi     *= M_PI/180.0;

  beam_offset.x *= scope->aperture();
  beam_offset.y *= scope->curvatureRadius();
  beam_offset.z *= scope->aperture();
  beam_offset += beam_base;

  Vec3D beam_center(scope->position());
  beam_offset.Rotate(Vec3D(1,0,0)*scope->elevation());
  beam_offset.Rotate(Vec3D(0,0,-1)*scope->azimuth());
  beam_center += beam_offset;

#if 0
  options.findWithValue("beam_center",beam_center);
#endif

  beam_direction.Rotate(Vec3D(1,0,0)*direction_theta);
  beam_direction.Rotate(Vec3D(0,1,0)*direction_phi);

  elevation         *= M_PI/180.0;
  azimuth           *= M_PI/180.0;

  beam_direction.Rotate(Vec3D(1,0,0)*elevation);
  beam_direction.Rotate(Vec3D(0,0,-1)*azimuth);

  std::cerr << "Beam center:       " << beam_center << std::endl
	    << "Beam_direction:    " << beam_direction << std::endl
	    << "Scope_elevation:   " 
	    << scope->elevation()*180./M_PI << std::endl
	    << "Scope_azimuth:     " << scope->azimuth()*180./M_PI << std::endl
	    << "Emission_start:    " << emission_start << std::endl
	    << "Emission_end:      " << emission_end << std::endl
	    << "Emission_radius_i: " << emission_radius_i << std::endl
	    << "Emission_radius_o: " << emission_radius_o << std::endl
	    << "Emission_angle_i:  " << emission_angle_i << std::endl
	    << "Emission_angle_o:  " << emission_angle_o << std::endl;

  // --------------------------------------------------------------------------
  // Raytracer
  // --------------------------------------------------------------------------

  double tocm = scope->pixelSpacing();
  double todeg = 180.0/M_PI/scope->focalPlanePosition().y;

  VSORayTracer tracer(array,rng);

  unsigned nhit = 0;
  VSSimpleStat2<double> x;
  VSSimpleStat2<double> y;
  VSSimpleStat2<double> xx;
  VSSimpleStat2<double> xy;
  VSSimpleStat2<double> yy;

  std::vector< std::pair< double, double > > ray_coord;

  for(unsigned i=0;i<num_photons;i++)
    {
      Particle ray;
      tracer.beam(ray,
		  beam_center, beam_direction, 
		  emission_start, emission_end, 
		  emission_radius_i, emission_radius_o,
		  emission_angle_i, emission_angle_o);
      
      VSORayTracer::TraceInfo info;
      tracer.trace(ray,info,scope);

      if(info.pixel_hexid != 0)
	{
	  ray_coord.
	    push_back(std::make_pair(info.fplane_x*tocm,info.fplane_z*tocm));

	  if(dump_rays)
	    {
	      std::cout << info.fplane_x*tocm << ' ' 
			<< info.fplane_z*tocm << ' '
			<< info.fplane_x*tocm*todeg << ' ' 
			<< info.fplane_z*tocm*todeg << '\n';
	    }
	  else
	    {
	      nhit++;
	      x.accumulate(info.fplane_x);
	      y.accumulate(info.fplane_z);
	      xx.accumulate(info.fplane_x*info.fplane_x);
	      xy.accumulate(info.fplane_x*info.fplane_z);
	      yy.accumulate(info.fplane_z*info.fplane_z);
	    }
	}
    }

  if(dump_rays)
    return EXIT_SUCCESS;

  // --------------------------------------------------------------------------
  // Calculate RMS
  // --------------------------------------------------------------------------

  double mx = x.mean()*tocm;
  double my = y.mean()*tocm;
  double mxx = xx.mean()*tocm*tocm-mx*mx;
  double mxy = xy.mean()*tocm*tocm-mx*my;
  double myy = yy.mean()*tocm*tocm-my*my;

  std::vector< double > theta;

  for(unsigned i = 0; i < ray_coord.size(); i++)
    {
      double th = sqrt(std::pow(mx-ray_coord[i].first,2) +
		       std::pow(my-ray_coord[i].second,2));
      theta.push_back(th);
    }

  double th68, th68_err;
  double th80, th80_err;

  containmentIntervalUpper(theta,0.68,th68,th68_err);
  containmentIntervalUpper(theta,0.80,th80,th80_err);
  
  double l = sqrt((mxx+myy+sqrt(4*mxy*mxy+(mxx-myy)*(mxx-myy)))/2.0);
  double w = sqrt((mxx+myy-sqrt(4*mxy*mxy+(mxx-myy)*(mxx-myy)))/2.0);

  std::cout << "Nray hits: " << nhit << std::endl
	    << "X-mean  :  " << mx << "  [" << mx*todeg << ']' << std::endl
	    << "Y-mean:    " << my << "  [" << my*todeg << ']' << std::endl
	    << "X-rms:     " << sqrt(mxx) 
	    << "  [" << sqrt(mxx)*todeg << ']' << std::endl
	    << "Y-rms:     " << sqrt(myy) 
	    << "  [" << sqrt(myy)*todeg << ']' << std::endl
	    << "Eigen-L:   " << l << "  [" << l*todeg << ']' << std::endl
	    << "Eigen-W:   " << w << "  [" << w*todeg << ']' << std::endl
	    << "R-rms:     " << sqrt(mxx+myy) 
	    << "  [" << sqrt(mxx+myy)*todeg << ']' << std::endl
	    << "FWHM:      " << 1.67*sqrt(mxx+myy) 
	    << "  [" << 1.67*sqrt(mxx+myy)*todeg << ']' << std::endl
    	    << "68%:       " << 2*th68 << " +/- " << 2*th68_err
	    << "  [" << 2*th68*todeg << " +/- " << 2*th68_err*todeg 
	    << "]" << std::endl
       	    << "80%:       " << 2*th80 << " +/- " << 2*th80_err
	    << "  [" << 2*th80*todeg << " +/- " << 2*th80_err*todeg 
	    << "]" << std::endl;

  // --------------------------------------------------------------------------
  // Done
  // --------------------------------------------------------------------------
  return EXIT_SUCCESS;
}
