//-*-mode:c++; mode:font-lock;-*-

/*! \file simple_reconstruct.cpp
  Simple reconstruction program to illustrate use of various DB components

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       09/16/2005
*/

//#define NO_OUTPUT

#include <list>
#include <iostream>
#include <iterator>

#include <RandomNumbers.hpp>

#include <VSOptions.hpp>
#include <VSDatabase.hpp>
#include <VSDBFactory.hpp>
#include <VSDBParameterTable.hpp>
#include <VSSimDB.hpp>
#include <VSSimDBTables.hpp>
#include <VSOTelescopeArray.hpp>
#include <VSTestTimer.hpp>
#include <VSAAtmosphere.hpp>
#include <VSAReconstruction.hpp>

//using namespace Physics;
using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

struct ChannelInfo
{
  ChannelInfo(): valid(false), vec() { } 
  bool valid;
  VSAAlgebra::Vec2D vec;
};

struct Datum
{
  Datum(): pixel_id(), pe_count(), pe_mean_time_ns() { }
  unsigned pixel_id;
  unsigned pe_count;
  float pe_mean_time_ns;
};

int main(int argc, char** argv)
{
  std::string progname(*argv);

  VSOptions options(argc,argv,true);

  // --------------------------------------------------------------------------
  // Options
  // --------------------------------------------------------------------------

  VSDBFactory::configure(&options);

#if 0
  std::ostringstream stream;
  stream << "/tmp/rng_state_uid" << getuid() << ".dat";
  std::string rngstatefile(stream.str());
  options.findWithValue("rng_state_file",rngstatefile,
			"Set random number generator state file");
#endif

  bool print_usage = false;
  int exit_code = EXIT_SUCCESS;

  if(options.find("h","Print this message") != VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help", "Print this message") != VSOptions::FS_NOT_FOUND)
    print_usage=true;

  bool single_event_mode = false;
  unsigned event_no = 0;
  if(options.findWithValue("ev",event_no,
			   "Process only one event, writing details of chi^2 "
			   "grid search. Number of event to process should be "
			   "given as an argument to this option") 
     != VSOptions::FS_NOT_FOUND)single_event_mode = true;

  unsigned max_events = 0;
  options.findWithValue("max_events",max_events,
			"Maximum number of events to process, or zero to "
			"process all events");
  
  unsigned imethod = 2;
  options.findWithValue("method",imethod,
			"Set the reconstruction method, should be 1,2 or 3");

  double plate_scale = 1.0;
  options.findWithValue("platescale",plate_scale,
			"Set the plate scale correction "
			"(1.0 gives no correction)");

  std::vector<unsigned> scope_suppress_id;
  options.findWithValue("scope_suppress",scope_suppress_id,
			"Comma separated list of telescope numbers to "
			"suppress in analysis");

  double fov = 3.5;
  options.findWithValue("fov",fov,
			"Set the camera field of view");

  unsigned l1_pe = 6;
  options.findWithValue("l1_pe",l1_pe,
			"Set the L1 trigger level in PE");

  unsigned l2_mult = 3;
  options.findWithValue("l2_mult",l2_mult,
			"Set the L2 trigger multiplicity in channels");

  bool l2_force = false;
  if(options.find("l2_force",
		  "Force L2 trigger on all telescopes when L3 "
		  "trigger occurs") != VSOptions::FS_NOT_FOUND)
    l2_force = true;

  unsigned l3_mult = 2;
  options.findWithValue("l3_mult",l3_mult,
			"Set the L3 trigger multiplicity in telescopes");

  unsigned cleaning_pe = unsigned(floor(0.1878*10*3.0));
  options.findWithValue("cleaning_pe",cleaning_pe,
			"Set the cleaning level in PE");

  // --------------------------------------------------------------------------
  // Get the database and table names from the command line arguments
  // --------------------------------------------------------------------------

  argc--,argv++;

  if (argc!=2)
    {
      print_usage=true;
      exit_code=EXIT_FAILURE;
    }

  if(print_usage)
    {
      std::cerr << "Usage: " << progname 
		<< " [options] database table_name" 
		<< std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return exit_code;
    }

  std::string database(*argv);
  argc--; argv++;

  std::string tablename(*argv);
  argc--; argv++;

  // --------------------------------------------------------------------------
  // Create the database connection
  // --------------------------------------------------------------------------

  VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
  if(db->useDatabase(database) < 0)
    {
      std::cerr << "Could not connect to database " << database << std::endl;
      return EXIT_FAILURE;
    };

  VSSimDB db_sim(db);
  VSDBParameterTable db_param(db);

  // --------------------------------------------------------------------------
  // Get the table parameters from the DB
  // --------------------------------------------------------------------------

  VSSimDBTableParam* param;
  param = db_sim.getDataTableByName(tablename);
  if(!param)
    {
      std::cerr << "Could not get information for table " << tablename 
		<< std::endl;
      return EXIT_FAILURE;      
    }

  // --------------------------------------------------------------------------
  // Get the definition of the array from the database
  // --------------------------------------------------------------------------

  VSOTelescopeArray array;
  array.readFromDatabase(db, param->fOpticsID);

  std::vector<bool> scope_suppress;
  std::vector<std::vector<ChannelInfo> > channel_info;
  std::vector<unsigned> max_channel;

  const unsigned nscope = array.numTelescopes();

  scope_suppress.resize(nscope);
  channel_info.resize(nscope);
  max_channel.resize(nscope);
  
  for(std::vector<unsigned>::const_iterator iid=scope_suppress_id.begin();
      iid!=scope_suppress_id.end(); iid++)
    if(*iid<scope_suppress.size())scope_suppress[*iid]=true;

  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      unsigned nvalid=0;
      unsigned npix=array.telescope(iscope)->numPixels();

      max_channel[iscope]=0;

      channel_info[iscope].resize(npix);
      for(unsigned ipix=0;ipix<npix; ipix++)
	{
	  Physics::Vec3D r = 
	    array.telescope(iscope)->pixel(ipix)->
	    incomingSkyVectorAtZenith(plate_scale);

	  channel_info[iscope][ipix].vec.set(r.x, -r.y);
	  
	  double theta = acos(r * Physics::Vec3D(0,0,-1));
	  
	  if(theta <= fov/2/180*M_PI)
	    {
	      max_channel[iscope]=ipix;
	      channel_info[iscope][ipix].valid=true;
	      nvalid++;
	    }
	}
      
      if(npix==0)
	{
	  max_channel[iscope] = 0xFFFFFFFFU;
	}

      std::cerr << progname << ": scope " << iscope << " has " << nvalid 
		<< " pixels (max id:" << max_channel[iscope] << ')'
		<< std::endl;
    }

  // --------------------------------------------------------------------------
  // Set up the reconstruction
  // --------------------------------------------------------------------------

  VSAReconstruction::ArrayInfo array_info(nscope);
  array_info.elevation = array.altitude()/100.0;
  for(unsigned iscope=0;iscope<nscope;iscope++)
    {
      const VSOTelescope* scope = array.telescope(iscope);

      unsigned npix = scope->numPixels();
      VSAReconstruction::ScopeInfo& scope_info(array_info.scopes[iscope]);
      scope_info.pij.resize(npix);

      scope_info.ri.set(scope->pos().x/100, scope->pos().y/100,
			scope->pos().z/100-array_info.elevation);
      scope_info.Di = scope->aperture();

      for(unsigned ipix=0;ipix<npix; ipix++)
	scope_info.pij[ipix] = channel_info[iscope][ipix].vec;
    }  
  
  VSAAtmosphere atmo = VSAAtmosphere::usStandard();
  VSAReconstruction reconstruction(array_info, atmo);

  // --------------------------------------------------------------------------
  // Get a list of all (or single, if requested) complete events
  // --------------------------------------------------------------------------
 
  std::string tablename_ev = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    +tablename
    + std::string(VSIMDB_TABLE_POSTFIX_EVENTS);

  std::string tablename_tel = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    +tablename
    + std::string(VSIMDB_TABLE_POSTFIX_SCOPES);

  std::string tablename_pe = 
    std::string(VSIMDB_TABLE_PREFIX_DATA)
    +tablename
    + std::string(VSIMDB_TABLE_POSTFIX_PES);

  std::string ev_condition = "EventComplete=1";
  if(single_event_mode)
    ev_condition += " AND EventID=" + VSDataConverter::toString(event_no);
  else if(max_events!=0)
    ev_condition += " AND EventID<" + VSDataConverter::toString(max_events);    
  VSDBStatement* stmt = 
    db->createSelectQuery(tablename_ev, ev_condition, "*",
			  VSDatabase::FLAG_NO_BUFFER);

  if(!stmt)
    {
      std::cerr << progname
		<< ": could no query table: " << tablename_ev << std::endl;
      return EXIT_FAILURE;
    }
  
  VSSimDBEventData event;
  std::list<VSSimDBEventData> events;

  stmt->bindToResult(event.fEventID);
  stmt->bindToResult(event.fTargetZenithRad);
  stmt->bindToResult(event.fTargetAzimuthRad);
  stmt->bindToResult(event.fPrimaryZenithRad);
  stmt->bindToResult(event.fPrimaryAzimuthRad);
  stmt->bindToResult(event.fPrimaryCoreEastM);
  stmt->bindToResult(event.fPrimaryCoreNorthM);
  stmt->bindToResult(event.fPrimaryCoreUpASLM);
  stmt->bindToResult(event.fNumHitScopes);
  stmt->bindToResult(event.fEventComplete);

  if(stmt->execute() < 0)
    {
      std::cerr << progname
		<< ": error executing query on: " << tablename_ev << std::endl;
      return EXIT_FAILURE;
    }

  VSTestTimer<VSTimerCoreIA32> timer_query;
  VSTestTimer<VSTimerCoreIA32> timer_reconstruction;
  //VSTestTimer<VSTimerCoreUNIX> timer_query;
  //VSTestTimer<VSTimerCoreUNIX> timer_reconstruction;

  timer_query.start();
  while(stmt->retrieveNextRow())events.push_back(event);
  timer_query.stop();
  delete stmt;

  std::cerr << progname << ": loaded " << events.size() 
	    << " event headers in "  << timer_query << " seconds" << std::endl;

  // --------------------------------------------------------------------------
  // Create the database queries and bind variables
  // --------------------------------------------------------------------------

  VSDBStatement* stmt_tel = 
    db->createSelectQuery(tablename_tel, "EventID=?", 
			  "ScopeID,ScopeZenithRad,ScopeAzimuthRad");

  if(!stmt_tel)
    {
      std::cerr << progname
		<< ": could no query table: " << tablename_tel << std::endl;
      return EXIT_FAILURE;
    }

  VSDBStatement* stmt_pix = 
    db->createQuery(std::string("SELECT PixelID,COUNT(PixelID),"
				"AVG(PixelTimeNS) FROM ")
		    + tablename_pe
		    + std::string(" WHERE EventID=? AND ScopeID=? "
				  "AND PixelID<=? GROUP BY PixelID"));

  if(!stmt_pix)
    {
      std::cerr << progname
		<< ": could no query table: " << tablename_pe << std::endl;
      return EXIT_FAILURE;
    }

  uint32_t event_id;
  uint16_t scope_id;

  float scope_zn;
  float scope_az;

  uint32_t max_pixel_id;
  Datum datum;

  stmt_tel->bindToParam(event_id);
  stmt_tel->bindToResult(scope_id);
  stmt_tel->bindToResult(scope_zn);
  stmt_tel->bindToResult(scope_az);

  stmt_pix->bindToParam(event_id);
  stmt_pix->bindToParam(scope_id);
  stmt_pix->bindToParam(max_pixel_id);
  stmt_pix->bindToResult(datum.pixel_id);
  stmt_pix->bindToResult(datum.pe_count);
  stmt_pix->bindToResult(datum.pe_mean_time_ns);
  
  // --------------------------------------------------------------------------
  // Query database and reconstruct events
  // --------------------------------------------------------------------------

  std::vector<std::vector<Datum> > data;
  data.resize(nscope);
  for(unsigned iscope=0;iscope<nscope;iscope++)
    data[iscope].reserve(max_channel[iscope]+1);

  std::vector<double> scopes_zn;
  std::vector<double> scopes_az;
  scopes_az.resize(nscope);
  scopes_zn.resize(nscope);

  timer_query.reset();
  unsigned ntrig_ev = 0;
  for(std::list<VSSimDBEventData>::const_iterator ievent = events.begin();
      ievent!=events.end(); ievent++)
    {
      unsigned ntrig_l3 = 0;

      std::vector<bool> scopes_l2;
      scopes_l2.resize(nscope);

      event_id = ievent->fEventID;

      for(unsigned iscope=0;iscope<nscope;iscope++)
	data[iscope].clear();

      timer_query.start();
      if(stmt_tel->execute()<0)
	{
	  std::cerr << progname
		    << ": error executing query on: " << tablename_tel 
		    << std::endl;
	  return EXIT_FAILURE;
	}
      
      while(stmt_tel->retrieveNextRow())
	{
	  if(scope_suppress[scope_id])continue;

	  uint32_t ntrig_l2 = 0;
	  max_pixel_id = max_channel[scope_id];

	  scopes_zn[scope_id] = scope_zn;
	  scopes_az[scope_id] = scope_az;
	  
	  if(stmt_pix->execute()<0)
	    {
	      std::cerr << progname
			<< ": error executing query on: " << tablename_pe 
			<< std::endl;
	      return EXIT_FAILURE;
	    }

	  while(stmt_pix->retrieveNextRow())
	    {
	      if((datum.pixel_id<channel_info[scope_id].size())&&
		 (!channel_info[scope_id][datum.pixel_id].valid))continue;
	      if(datum.pe_count>=l1_pe)ntrig_l2++;
	      if(datum.pe_count>=cleaning_pe)
		data[scope_id].push_back(datum);
	    }

	  if(ntrig_l2>=l2_mult)
	    {
	      ntrig_l3++;
	      scopes_l2[scope_id]=true;
	    }
	}

      timer_query.stop();

#ifndef NO_OUTPUT
      if(single_event_mode)
	std::cerr 
	  << "event:   " << ievent->fEventID << std::endl
	  << "tar zn:  " << ievent->fTargetZenithRad/M_PI*180 << std::endl
	  << "tar az:  " << ievent->fTargetAzimuthRad/M_PI*180 << std::endl
 	  << "pri zn:  " << ievent->fPrimaryZenithRad/M_PI*180 << std::endl
	  << "pri az:  " << ievent->fPrimaryAzimuthRad/M_PI*180 << std::endl
	  << "core ew: " << ievent->fPrimaryCoreEastM << std::endl
	  << "core ns: " << ievent->fPrimaryCoreNorthM << std::endl
	  << "core ud: " << ievent->fPrimaryCoreUpASLM << std::endl
	  << "L3:      " << ntrig_l3 << std::endl;
      else
	std::cout 
	  << ievent->fEventID
	  << ' ' << ievent->fTargetZenithRad
	  << ' ' << ievent->fTargetAzimuthRad
	  << ' ' << ievent->fPrimaryZenithRad
	  << ' ' << ievent->fPrimaryAzimuthRad
	  << ' ' << ievent->fPrimaryCoreEastM
	  << ' ' << ievent->fPrimaryCoreNorthM
	  << ' ' << ievent->fPrimaryCoreUpASLM
	  << ' ' << ntrig_l3;
#endif
	    
      if(ntrig_l3>=l3_mult)
	{
#ifndef NO_OUTPUT
	  if(!single_event_mode)
	    std::cout << ' ' << 1;
#endif

	  // ******************************************************************
	  // R E C O N S T R U C T I O N
	  // ******************************************************************

	  timer_reconstruction.start();

	  std::vector<VSAReconstruction::ScopeImage> images(nscope);
	  unsigned nray = 0;
	  for(unsigned iscope=0;iscope<nscope;iscope++)
	    {
	      images[iscope].zenith = ievent->fTargetZenithRad;
	      images[iscope].azimuth = ievent->fTargetAzimuthRad;
	      images[iscope].camera_rotation = 0;
	      
	      unsigned npix = data[iscope].size();
	      if((npix==0)||((l2_force==false)&&(scopes_l2[iscope]==false)))
		 continue;

	      images[iscope].pixels.reserve(npix);
	      for(unsigned ipix=0;ipix<npix;ipix++)
		{
		  VSAReconstruction::PixelImage pix_image;
		  pix_image.j = data[iscope][ipix].pixel_id;
		  pix_image.nij = data[iscope][ipix].pe_count;
		  pix_image.tij = data[iscope][ipix].pe_mean_time_ns;
		  images[iscope].pixels.push_back(pix_image);
		}
	    }

	  VSAReconstruction::Reconstruction recon;
	  VSAReconstruction::Method method = VSAReconstruction::M_2;
	  if(imethod==1)method = VSAReconstruction::M_1;
	  else if(imethod==2)method = VSAReconstruction::M_2;
	  else if(imethod==3)method = VSAReconstruction::M_3;
	  reconstruction.reconstruct(recon, method, images);

	  VSAReconstruction::ArrayParameters param;
	  reconstruction.calculateArrayParameters(param, images, recon);

	  timer_reconstruction.stop();

	  Vec3D R0;
	  if(recon.e.z() >= -DBL_EPSILON)
	    R0 = Vec3D(1e8, 1e8, 0);
	  else
	    R0 = recon.R - recon.e*(recon.R.z()/recon.e.z());

	  Vec3D DR = R0-Vec3D(ievent->fPrimaryCoreEastM,
			      ievent->fPrimaryCoreNorthM,
			      ievent->fPrimaryCoreUpASLM/100
			      -array_info.elevation);
	  Vec3D epri = Vec3D(sin(ievent->fTargetZenithRad)*
			     sin(ievent->fTargetAzimuthRad),
			     sin(ievent->fTargetZenithRad)*
			     cos(ievent->fTargetAzimuthRad),
			     cos(ievent->fTargetZenithRad));
	  
#ifndef NO_OUTPUT
	  if(single_event_mode)
	    std::cerr
	      << "e:       " << recon.e << std::endl
	      << "theta:   " << acos(-recon.e*epri)/M_PI*180 << std::endl
	      << "R:       " << recon.R << std::endl
	      << "R0:      " << R0 << std::endl
	      << "DR:      " << DR << std::endl
	      << "chi2e:   " << recon.chi2e << std::endl
	      << "chi2r:   " << recon.chi2R << std::endl
	      << "deltael: " << recon.deltael/M_PI*180 << std::endl
	      << "deltaew: " << recon.deltaew/M_PI*180 << std::endl
	      << "deltaRl: " << recon.deltaRl << std::endl
	      << "deltaRw: " << recon.deltaRw << std::endl << std::endl;
	  else
	    std::cout
	      << ' ' << recon.e
	      << ' ' << acos(-recon.e*epri)/M_PI*180
	      //<< ' ' << recon.R
	      << ' ' << R0
	      << ' ' << DR
	      << ' ' << recon.chi2e
	      << ' ' << recon.chi2R
	      << ' ' << recon.deltael/M_PI*180
	      << ' ' << recon.deltaew/M_PI*180
	      << ' ' << recon.deltaRl
	      << ' ' << recon.deltaRw;
#endif

	  if(recon.e.z() < 0)
	    for(unsigned iscope=0;iscope<param.scopes.size();iscope++)
	      {
		const VSAReconstruction::ScopeParameters& 
		  sp(param.scopes[iscope]);
		
		if(sp.Ni == 0)continue;
		
		Eigen3D delta2i_eigen;
		sp.delta2i.eigen(delta2i_eigen);
		
		if(single_event_mode)
		  std::cerr
		    << "i:         " << iscope << std::endl
		    << "Ni:        " << sp.Ni << std::endl
		    << "d1i:       " << sp.d1i << std::endl
		    << "d2i:       " << sp.d2i << std::endl
		    << "l1i:       " << sp.l1i << std::endl
		    << "l2i:       " << sp.l2i << std::endl
		    << "theta1i:   " << sp.theta1i/M_PI*180 << std::endl
		    << "theta2i:   " << sp.theta2i/M_PI*180 << std::endl
		    << "delta1i:   " << sp.delta1i << std::endl
		    //<< "delta2i:   " << sp.delta2i << std::endl
		    << "delta2i_1: " << sqrt(delta2i_eigen.val[0]) << " " 
		    << delta2i_eigen.vec[0] << std::endl
		    << "delta2i_2: " << sqrt(delta2i_eigen.val[1]) << " " 
		    << delta2i_eigen.vec[1] << std::endl
		    << "delta2i_3: " << sqrt(delta2i_eigen.val[2]) << " " 
		    << delta2i_eigen.vec[2] << std::endl
		    
		    << "D1i:       " << sp.D1i << std::endl
		    << "h1i:       " << sp.h1i << std::endl
		    << "h2i:       " << sp.h2i << std::endl
		    << "G1i:       " << sp.G1i << std::endl
		    << "G2i:       " << sp.G2i << std::endl
		    << "thetac1:   " << sp.thetac1/M_PI*180 << std::endl
		    
		    << "t01i:      " << sp.t01i << std::endl
		    << "t02i:      " << sp.t02i << std::endl
		    
		    << "lambdadi:  " << sp.lambdadi << std::endl
		    << "lambdaci:  " << sp.lambdaci << std::endl
		    << std::endl;
	      }
	  ntrig_ev++;
	}
      else 
	{
#ifndef NO_OUTPUT
	  if(!single_event_mode)
	    std::cout << ' ' << 0;
#endif
	}
      
#ifndef NO_OUTPUT
      if(!single_event_mode)
	std::cout << std::endl;
#endif      
    }
  
  std::cerr << progname << ": " << events.size() << " events read in "
	    << timer_query << " seconds" << std::endl;

  std::cerr << progname << ": " << ntrig_ev 
	    << " triggering events reconstructed in " 
	    << timer_reconstruction << " seconds" << std::endl;

}
