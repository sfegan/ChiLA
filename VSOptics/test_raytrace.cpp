//-*-mode:c++; mode:font-lock;-*-

/*! \file test_raytrace.cpp

  Program to test ray tracing code

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \date    12/07/2004
  \version 0.2
  \note
*/

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include <Vec3D.hpp>
#include <Vec4D.hpp>

#include <VSOptions.hpp>
#include <VSDatabase.hpp>
#include <VSDBFactory.hpp>

#include "VSOArrayParameters.hpp"
#include "VSOTelescopeArray.hpp"
#include "VSORayTracer.hpp"

using namespace Physics;
using namespace VERITAS;

const char* ts_name[] = { "TS_NONE",
			  "TS_DOES_INTERSECT_GROUND",
			  "TS_NO_SCOPE",
			  "TS_OUTSIDE_REFLECTOR_IP",
			  "TS_MISSED_REFLECTOR_SPHERE",
			  "TS_NO_MIRROR",
			  "TS_MIRROR_REMOVED",
			  "TS_MISSED_MIRROR_SPHERE",
			  "TS_MISSED_MIRROR_EDGE",
			  "TS_ABSORBED_AT_MIRROR",
			  "TS_TRAVELLING_AWAY_FROM_FOCAL_PLANE",
			  "TS_NO_PIXEL",
			  "TS_ABSORBED_AT_CONCENTRATOR",
			  "TS_PE_GENERATED" };

void usage(std::ostream& stream, const std::string& program)
{
  stream << program << " [options]" << std::endl;
  stream << std::endl;
  stream << "Array configuation source (one must be selected)" << std::endl;
  stream << "  -inifile FILENAME                    Generate random array from INI file" << std::endl;
  stream << "                                         If FILENAME is empty, use default" << std::endl;
  stream << "                                         array parameters." << std::endl;
  stream << "  -db DATABASE                         Read array configuration from database" << std::endl;
  stream << "  -optics_id id                        Use configuration #id from database" << std::endl;
  stream << std::endl;
  stream << "Random array generation options" << std::endl;
  VSOArrayParameters::printOptions(stream);
  stream << std::endl;
  stream << "Photon test beam setup:" << std::endl;
  stream << "  -scope N                             Base beam on telescope N" << std::endl;
  stream << "  -mirror N                            Base beam on mirror N" << std::endl;
  stream << "  -pixel N                             Direct beam into pixel N" << std::endl;
  stream << "  -beam_offset (X,Y,Z)                 Center beam at X,Y,Z in reflector coords" << std::endl;
  stream << "  -beam_center (X,Y,Z)                 Center beam at X,Y,Z in array coords" << std::endl;
  stream << "  -emission_start DIST                 Start photon emission region at D" << std::endl;
  stream << "  -emission_end DIST                   Stop photon emission region at D" << std::endl;
  stream << "  -emission_radius_i RADIUS            Inner radius of emission region" << std::endl;
  stream << "  -emission_radius_o RADIUS            Outer radius of emission region" << std::endl;
  stream << "  -emission_angle_i ANGLE              Inner radius of emission region" << std::endl;
  stream << "  -emission_angle_o ANGLE              Outer radius of emission region" << std::endl;
  stream << "  -direction_theta ANGLE               Offset emission direction by ANGLE" << std::endl;
  stream << "  -direction_phi ANGLE                 Set direction of offset angle" << std::endl;
  stream << "  -elevation ANGLE                     Set elevation of beam direction" << std::endl;
  stream << "  -azimuth ANGLE                       Set azimuth of beam direction" << std::endl;
  stream << std::endl;
  stream << "Simulation options:" << std::endl;
  stream << "  -n N                                 Number of photons to sample" << std::endl;
  stream << "  -dump_array FILENAME                 Dump array to FILENAME" << std::endl;
}

int main(int argc, char ** argv)
{
  RandomNumbers rng("random.seed");

  // -------------------------------------------------------------------------
  // Command Line Options
  // --------------------------------------------------------------------------
  VSOptions options(argc,argv);
  VSOArrayParameters::setCanonicalValuesFromOptions(options);
  VSDBFactory::configure(&options);

  // --------------------------------------------------------------------------
  // Generate the array
  // --------------------------------------------------------------------------
  VSOTelescopeArray array;
  std::string name;
  if(options.findWithValue("inifile",name) == VSOptions::FS_FOUND)
    {
      VSOArrayParameters param;
      if(name.empty())param.reset(true);
      else param.readFromArrayINIFile(name);
      array.generateFromArrayParameters(param, rng);
    }
  else if(options.findWithValue("db",name) == VSOptions::FS_FOUND)
    {
      unsigned optics_id = 0;
      options.findWithValue("optics_id",optics_id);

      VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
      db->useDatabase(name);

      // ----------------------------------------------------------------------
      // Retrieve an array configuration
      // ----------------------------------------------------------------------
      array.readFromDatabase(db,optics_id);
      delete db;
    }
  else
    {
      usage(std::cerr, options.arg0());
      return EXIT_FAILURE;
    }

  std::string dumpfile;
  if(options.findWithValue("dump_array",dumpfile))
    {
      std::ofstream stream(dumpfile.c_str());
      array.dumpShort(stream);
    }

  bool quiet = false;
  if(options.find("q") != VSOptions::FS_NOT_FOUND)quiet=true;

  unsigned num_scope = 1;
  unsigned num_mirror = 0;
  unsigned num_pixel = 0;

  options.findWithValue("scope",num_scope);
  options.findWithValue("mirror",num_mirror);
  options.findWithValue("pixel",num_pixel);

  VSOTelescope* scope = array.telescopeByHexID(num_scope);
  if(scope == 0)
    {
      std::cerr << "Array does not contain telescope " << num_scope 
		<< std::endl;
      return EXIT_FAILURE;
    }

  double emission_start    = -2.0;
  double emission_end      = -2.0;
  double emission_radius_i = 0.0;
  double emission_radius_o = 0.55;
  double emission_angle_i  = 0.0;
  double emission_angle_o  = 0.0;

  unsigned num_photons     = 100; // x 1000 photons

  Vec3D beam_base;
  if((num_mirror != 0)&&(scope->mirrorByHexID(num_mirror)!=0))
    beam_base = scope->mirrorByHexID(num_mirror)->pos();

  Vec3D beam_offset;

  double direction_theta   = 0.0;
  double direction_phi     = 0.0;

  Vec3D beam_direction(0,-1,0);
  if((num_pixel != 0)&&(scope->pixelByHexID(num_pixel)!=0))
    {
      beam_direction = 
	(scope->pixelByHexID(num_pixel)->pos()+scope->focalPlanePosition());
      beam_direction /= beam_direction.Norm();
      beam_direction.y = -beam_direction.y;
    }

  double elevation = scope->elevation() / M_PI * 180.0;
  double azimuth = scope->azimuth() / M_PI * 180.0;  

  // --------------------------------------------------------------------------
  // Set variables from Options
  // --------------------------------------------------------------------------
  options.findWithValue("n",num_photons);
  options.findWithValue("emission_start",emission_start);
  options.findWithValue("emission_end",emission_end);
  options.findWithValue("emission_radius_i",emission_radius_i);
  options.findWithValue("emission_radius_o",emission_radius_o);
  options.findWithValue("emission_angle_i",emission_angle_i);
  options.findWithValue("emission_angle_o",emission_angle_o);
  options.findWithValue("beam_offset",beam_offset);
  options.findWithValue("direction_theta",direction_theta);
  options.findWithValue("direction_phi",direction_phi);
  options.findWithValue("elevation",elevation);
  options.findWithValue("azimuth",azimuth);

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

  options.findWithValue("beam_center",beam_center);

  beam_direction.Rotate(Vec3D(1,0,0)*direction_theta);
  beam_direction.Rotate(Vec3D(0,1,0)*direction_phi);

  elevation         *= M_PI/180.0;
  azimuth           *= M_PI/180.0;

  beam_direction.Rotate(Vec3D(1,0,0)*elevation);
  beam_direction.Rotate(Vec3D(0,0,-1)*azimuth);

  if(!quiet)
    std::cerr << "Beam center:       " << beam_center << std::endl
	      << "Beam_direction:    " << beam_direction << std::endl
	      << "Emission_start:    " << emission_start << std::endl
	      << "Emission_end:      " << emission_end << std::endl
	      << "Emission_radius_i: " << emission_radius_i << std::endl
	      << "Emission_radius_o: " << emission_radius_o << std::endl
	      << "Emission_angle_i:  " << emission_angle_i << std::endl
	      << "Emission_angle_o:  " << emission_angle_o << std::endl;

  // --------------------------------------------------------------------------
  // 
  // --------------------------------------------------------------------------
  if(argc>1)
    {
      std::cerr << "Unrecognised options: ";
      for(int i=1;i<argc;i++)
	{
	  if(i!=1)std::cerr << ' ';
	  std::cerr << argv[i]; 
	}
      std::cerr << std::endl;
      return EXIT_FAILURE;
    }

  // --------------------------------------------------------------------------
  // Raytracer
  // --------------------------------------------------------------------------
  VSORayTracer tracer(array,rng);

  int c[sizeof(ts_name)/sizeof(*ts_name)];
  for(unsigned i=0;i<sizeof(c)/sizeof(*c);i++)c[i]=0;

  for(unsigned i=0;i<num_photons;i++)
    {
      Particle ray;
      tracer.beam(ray,
		  beam_center, beam_direction, 
		  emission_start, emission_end, 
		  emission_radius_i, emission_radius_o,
		  emission_angle_i, emission_angle_o);
      
      VSORayTracer::TraceInfo info;
      tracer.trace(ray,info);

      unsigned s = (unsigned)info.status;
      c[s]++;

      info.write(std::cout,true); 
      std::cout << std::endl;
    }

  if(!quiet)
    for(unsigned i=0;i<sizeof(c)/sizeof(*c);i++)
      std::cerr 
	<< std::left << std::setw(35) << ts_name[i] << ' ' << c[i] << std::endl;
}
