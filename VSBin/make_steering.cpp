//-*-mode:c++; mode:font-lock;-*-

/*! \file make_steering.cpp
  Make the steering file for CORSIKA

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    1.0
  \date       03/02/2005
*/

#ifndef NO_CEFFIC
#ifndef CEFFIC
#define CEFFIC
#endif
#endif

#ifndef NO_QGSJET
#ifndef QGSJET
#define QGSJET
#endif
#endif

#ifndef NO_FLUKA
#ifndef FLUKA
#define FLUKA
#endif
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <sys/time.h>
#include <unistd.h>

#include <VSOptions.hpp>
#include <VSDataConverter.hpp>
#include <xytohex.h>
#include <xytosquare.h>

#include <VSDBFactory.hpp>
#include <VSDBParameterTable.hpp>
#include <VSSimDB.hpp>
#include <VSOTelescopeArray.hpp>
#include <VSTargeting.hpp>
#include <VSSimDBCORSIKADatasets.hpp>

#include "make_steering_help.h"

using namespace VERITAS;

struct DataFile
{
  const char* filename;
  const char* data;
  unsigned data_len;
};

void writeCORSIKADataFile(bool test_only, const DataFile& df)
{
  if(!test_only)
    {
      FILE* fp = fopen(df.filename,"w");
      fwrite(df.data, 1, df.data_len, fp);
      fclose(fp);
    }
}

void writeDataOrDie(bool test_only,
		    const std::map<std::string, const DataFile*>& datafiles,
		    const std::string& choice,
		    const std::string& blurb)
{
  std::map<std::string, const DataFile*>::const_iterator idata = 
    datafiles.find(choice);

  if(idata == datafiles.end())
    {
      std::cerr << "Unknown " << blurb << " choice: \"" << choice << '"'
		<< std::endl
		<< "Valid choices are:";
      for(idata=datafiles.begin(); idata!=datafiles.end(); idata++)
	std::cerr << " \"" << idata->first << '"';
      std::cerr << std::endl;
      exit(EXIT_FAILURE);
    }
  else if(idata->second)
    {
      writeCORSIKADataFile(test_only,*idata->second);
    }
}

void writeAllCORSIKADataFiles(bool test_only,
			      const std::string& choice_modtran,
			      const std::string& choice_qeff,
			      const std::string& choice_atmo,
			      const std::string& choice_mirror)
{
#include "bern_atmprof1_dat.h"  
#include "bern_atmprof2_dat.h"  
#include "bern_atmprof3_dat.h"  
#include "bern_atmprof4_dat.h"  
#include "bern_atmprof5_dat.h"  
#include "bern_atmprof6_dat.h"  
#include "bern_atmprof9_dat.h"  
  std::map<std::string, const DataFile*> Modtran;
  Modtran["none"] = 0;
  Modtran["corsika1"] = &BERN_ATMPROF1_DAT;
  Modtran["corsika2"] = &BERN_ATMPROF2_DAT;
  Modtran["corsika3"] = &BERN_ATMPROF3_DAT;
  Modtran["corsika4"] = &BERN_ATMPROF4_DAT;
  Modtran["corsika5"] = &BERN_ATMPROF5_DAT;
  Modtran["corsika6"] = &BERN_ATMPROF6_DAT;
  Modtran["corsika9"] = &BERN_ATMPROF9_DAT;   
  writeDataOrDie(test_only,Modtran,choice_modtran,"paramterized atmosphere");

#ifdef CEFFIC
#include "cors_quanteff_dat.h"
  std::map<std::string, const DataFile*> Qeff;
  Qeff["none"] = 0;
  Qeff["corsika"] = &CORS_QUANTEFF_DAT;
  writeDataOrDie(test_only,Qeff,choice_qeff,"quantum efficiency");

#include "cors_atmabs_dat.h"
  std::map<std::string, const DataFile*> Atmo;
  Atmo["none"] = 0;
  Atmo["corsika"] = &CORS_ATMABS_DAT;
  writeDataOrDie(test_only,Atmo,choice_atmo,"atmospheric absorption");

#include "cors_mirreff_corsika_dat.h"
#include "cors_mirreff_falcone_dat.h"
#include "cors_mirreff_falcone_squared_dat.h"
  std::map<std::string, const DataFile*> Mirror;
  Mirror["none"] = 0;
  Mirror["corsika"] = &CORS_MIRREFF_CORSIKA_DAT;
  Mirror["double"] = &CORS_MIRREFF_FALCONE_SQUARED_DAT;
  Mirror["falcone"] = &CORS_MIRREFF_FALCONE_DAT;
  writeDataOrDie(test_only,Mirror,choice_mirror,"mirror reflectivity");
#else
  choice_qeff.empty(); // so that we don't get an unused argument warning.
  choice_atmo.empty(); // so that we don't get an unused argument warning.
  choice_mirror.empty(); // so that we don't get an unused argument warning.
#endif 

#include "cors_nucnuccs.h"
#include "cors_egsdat5_005.h"
#include "cors_egsdat5_015.h"
#include "cors_egsdat5_025.h"
#include "cors_egsdat5_040.h"
#include "cors_egsdat5_100.h"
#include "cors_egsdat5_300.h"
  writeCORSIKADataFile(test_only,CORS_NUCNUCCS);
  writeCORSIKADataFile(test_only,CORS_EGSDAT5_005);
  writeCORSIKADataFile(test_only,CORS_EGSDAT5_015);
  writeCORSIKADataFile(test_only,CORS_EGSDAT5_025);
  writeCORSIKADataFile(test_only,CORS_EGSDAT5_040);
  writeCORSIKADataFile(test_only,CORS_EGSDAT5_100);
  writeCORSIKADataFile(test_only,CORS_EGSDAT5_300);

#ifdef QGSJET
#include "cors_qgsdat.h"
#include "cors_sectnu.h"
  writeCORSIKADataFile(test_only,CORS_QGSDAT);
  writeCORSIKADataFile(test_only,CORS_SECTNU);
#endif

#ifdef DPMJET
#include "cors_glaubtar_dat.h"
#include "cors_nuclear_bin.h"
  writeCORSIKADataFile(test_only,CORS_GLAUBTAR_DAT);
  writeCORSIKADataFile(test_only,CORS_NUCLEAR_BIN);
#endif

#ifdef VENUS
#include "cors_venusdat.h"
  writeCORSIKADataFile(test_only,CORS_VENUSDAT);
#endif
  
#ifdef FLUKA
#include "fluka_elasct_bin.h"
#include "fluka_nuclear_bin.h"
#include "fluka_sigmapi_bin.h"
  writeCORSIKADataFile(test_only,FLUKA_ELASCT_BIN);
  writeCORSIKADataFile(test_only,FLUKA_NUCLEAR_BIN);
  writeCORSIKADataFile(test_only,FLUKA_SIGMAPI_BIN);
#endif
}

void write_card_option(std::ostream& stream, const std::string& opt)
{
  stream << std::setw(10) << std::left  << opt << std::endl;
}

template<typename T> void write_card_option(std::ostream& stream,
					    const std::string& opt,
					    const T& value)
{
  stream << std::setw(10) << std::left  << opt << ' ' 
	 << VSDataConverter::toString(value) << std::endl;
}

template<> void write_card_option<>(std::ostream& stream,
				    const std::string& opt,
				    const bool& value)
{
  stream << std::setw(10) << std::left  << opt << ' ' 
	 << (value?'T':'F') << std::endl;
}

template<> void write_card_option<>(std::ostream& stream,
				    const std::string& opt,
				    const std::string& value)
{
  stream << std::setw(10) << std::left  << opt << ' ' << '\''
	 << VSDataConverter::toString(value) << '\'' << std::endl;
}

struct CSCAT_Struct
{
  unsigned event_use;
  float sampling_radius;
};

template<> void write_card_option<>(std::ostream& stream,
				    const std::string& opt,
				    const CSCAT_Struct& value)
{
  stream << std::setw(10) << std::left  << opt << ' '
	 << VSDataConverter::toString(unsigned(value.event_use)) << ' ' 
	 << VSDataConverter::toString(value.sampling_radius) << ' '
	 << VSDataConverter::toString((float)0.0) << std::endl;
}

#ifdef QGSJET
struct QGSJET_Struct
{
  unsigned debug;
};

template<> void write_card_option<>(std::ostream& stream,
				    const std::string& opt,
				    const QGSJET_Struct& value)
{
  stream << std::setw(10) << std::left  << opt << ' '
	 << 'T' << ' ' << value.debug << std::endl;
}
#endif

struct Modtran_Struct
{
  unsigned atmosphere;
  bool refraction;
};

template<> void write_card_option<>(std::ostream& stream,
				    const std::string& opt,
				    const Modtran_Struct& value)
{
  stream << std::setw(10) << std::left  << opt << ' '
	 << VSDataConverter::toString(value.atmosphere) << ' ' 
	 << (value.refraction?'T':'F') << std::endl;
}

template<typename T> void write_card_option(std::ostream& stream,
					    const std::string& opt,
					    const std::vector<T>& value)
{
  stream << std::setw(10) << std::left  << opt;
  for(typename std::vector<T>::const_iterator i=value.begin(); 
      i!=value.end(); i++)
    stream << ' ' << VSDataConverter::toString(*i);
  stream << std::endl;
}

template<> void write_card_option<>(std::ostream& stream,
				    const std::string& opt,
				    const std::vector<bool>& value)
{
  stream << std::setw(10) << std::left  << opt;
  for(std::vector<bool>::const_iterator i=value.begin(); i!=value.end(); i++)
    stream << ' ' << ((*i)?'T':'F');
  stream << std::endl;
}

template<> void write_card_option<>(std::ostream& stream,
				    const std::string& opt,
				    const std::vector<std::string>& value)
{
  stream << std::setw(10) << std::left  << opt;
  for(std::vector<std::string>::const_iterator i=value.begin(); i
	!=value.end(); i++)
    stream << ' ' << '\'' << VSDataConverter::toString(*i) << '\'';
  stream << std::endl;
}

#define NUMOF(x) (sizeof(x)/sizeof(*(x)))

int main(int argc, char** argv)
{
  std::string progname(*argv);

  VSOptions options(argc, argv, true);

  bool usage_error = false;
  if(options.find("h")!=VSOptions::FS_NOT_FOUND)usage_error=true;
  if(options.find("help")!=VSOptions::FS_NOT_FOUND)usage_error=true;

  // --------------------------------------------------------------------------
  // DATABASE FACTORY CONFIGURATION
  // --------------------------------------------------------------------------

  VSDBFactory::configure(&options);

  // --------------------------------------------------------------------------
  // SOME INITIALIZATION OPTIONS
  // --------------------------------------------------------------------------

  bool ctl_init = false;
  bool ctl_load = true;
  bool ctl_save = false;

  if(options.find("INIT",HELP_INIT)!=VSOptions::FS_NOT_FOUND)ctl_init=true;
  if(options.find("no_load",HELP_NO_LOAD)!=VSOptions::FS_NOT_FOUND)
    ctl_load=false;
  if(options.find("SAVE",HELP_SAVE)!=VSOptions::FS_NOT_FOUND)ctl_save=true;

  // --------------------------------------------------------------------------
  // CHECK FOR OPTION TO USE DATABASE
  // --------------------------------------------------------------------------
  
  bool        ctl_db_use            = false;
  std::string ctl_db_name           = "";
  bool        ctl_table_use         = false;
  std::string ctl_table_name        = "";

  switch(options.findWithValue("db",ctl_db_name,HELP_DB))
    {
    case VSOptions::FS_FOUND:
      ctl_db_use = true;
      break;
    case VSOptions::FS_NOT_FOUND:
      ctl_db_use = false;
      break;
    case VSOptions::FS_FOUND_BUT_WITHOUT_VALUE:
      std::cerr << progname 
		<< ": database option used but no database specified"
		<< std::endl;
      usage_error = true;
      break;
    case VSOptions::FS_FOUND_BUT_COULD_NOT_CONVERT_VALUE:
    case VSOptions::FS_FOUND_WITH_UNDESIRED_VALUE:
    default:
      // Should not happen... we can always convert to string and we
      // explicitly asked for a value
      assert(0);
    };

  switch(options.findWithValue("table",ctl_table_name,HELP_TABLE))
    {
    case VSOptions::FS_FOUND:
      ctl_table_use = true;
      break;
    case VSOptions::FS_NOT_FOUND:
      ctl_table_use = false;
      break;
    case VSOptions::FS_FOUND_BUT_WITHOUT_VALUE:
      std::cerr << progname << ": table option used but no table specified"
		<< std::endl;
      usage_error = true;
      break;
    case VSOptions::FS_FOUND_BUT_COULD_NOT_CONVERT_VALUE:
    case VSOptions::FS_FOUND_WITH_UNDESIRED_VALUE:
    default:
      // Should not happen... we can always convert to string and we
      // explicitly asked for a value
      assert(0);
    };

  // --------------------------------------------------------------------------
  // CONNECT TO DATABASE
  // --------------------------------------------------------------------------

  VSDatabase* db(0);
  VSDBParameterTable* db_param(0);
  if(ctl_db_use)
    {
      db = VSDBFactory::getInstance()->createVSDB();
      if(db->useDatabase(ctl_db_name) < 0)
	{
	  std::cerr << progname << ": could not connect to database " 
		    << ctl_db_name << std::endl;
	  return EXIT_FAILURE;
	}
      db_param = new VSDBParameterTable(db);
    }

  VSSimDB* db_sim(0);
  if(ctl_table_use)
    {
      db_sim = new VSSimDB(db);
    }

  // --------------------------------------------------------------------------
  // LOAD COMMAND LINE OPTIONS FROM DATABASE
  // --------------------------------------------------------------------------

  if((ctl_db_use)&&(ctl_load))
    {
      VSDBParameterSet parameters;
      db_param->retrieveParameterSet("MakeSteering",parameters);
      for(VSDBParameterSet::const_iterator iopt = parameters.begin();
	  iopt != parameters.end(); iopt++)
	{
	  if(iopt->second.empty())options.addOption(iopt->first);
	  else options.addOptionWithValue(iopt->first,iopt->second);
	}
    }

  // --------------------------------------------------------------------------
  // SET UP SOME OF THE DEFAULTS FOR THE STEERING FILE
  // --------------------------------------------------------------------------

  char init_host_buffer[256];
  gethostname(init_host_buffer, NUMOF(init_host_buffer));
  for(unsigned i=0; i<NUMOF(init_host_buffer); i++)
    if((init_host_buffer[i]=='.')||(i==255))init_host_buffer[i]='\0';

  struct timeval init_tv;
  gettimeofday(&init_tv,0);
  srandom(init_tv.tv_usec^init_tv.tv_sec);
  for(unsigned i=0;i<100000;i++)random();

  unsigned    init_seed1 = random();
  unsigned    init_seed2 = random();
  unsigned    init_seed3 = random();
  unsigned    init_seed4 = random();
  
  while((init_seed1>900000000)||(init_seed1<1))init_seed1 = random();
  while((init_seed2>900000000)||(init_seed2<1))init_seed2 = random();
  while((init_seed3>900000000)||(init_seed3<1))init_seed3 = random();
  while((init_seed4>900000000)||(init_seed4<1))init_seed4 = random();

  char init_user_buffer[256];
  if(getlogin() == NULL)
    sprintf(init_user_buffer,"UID_%d",getuid());
  else
    strcpy(init_user_buffer,getlogin());

  // --------------------------------------------------------------------------
  // DEFINE ALL VARIABLES TO HOLD STEERING FILE OPTIONS
  // --------------------------------------------------------------------------

  unsigned    cor_debug             = 0;
  bool        cor_no_qgs            = false;
  bool        cor_no_ceffic         = false;
  bool        cor_no_fluka          = false;
  bool        cor_no_modtran        = false;

  std::string run_hostname          = init_host_buffer;
  std::string run_username          = init_user_buffer;
  unsigned    run_seed1             = init_seed1;
  unsigned    run_seed2             = init_seed2;
  unsigned    run_seed3             = init_seed3;
  unsigned    run_seed4             = init_seed4;
  unsigned    run_runno             = 1;
  unsigned    run_eventno           = 1;
  unsigned    run_num_shower        = 100;
	      
  uint32_t    tel_optics_id         = 0;
  unsigned    tel_rings             = 0;
  float       tel_spacing           = 80.0;
  float       tel_elevation         = 3500.0;
  float       tel_aperture          = 7.0;
  float       tel_max_dist          = std::numeric_limits<float>::infinity();
  std::string tel_layout            = "hex";

  unsigned    pri_species           = 1;
  float       pri_zenith_lo         = 20.0;
  float       pri_zenith_hi         = 20.0;
  float       pri_azimuth_lo        = 0.0;
  float       pri_azimuth_hi        = 360.0;
  float       pri_viewcone_lo       = 0.0;
  float       pri_viewcone_hi       = 0.0;
  float       pri_sampling_radius   = 0.0;
  float       pri_energy            = 100.0;
  float       pri_start_depth       = 0.0;

  unsigned    sim_atmosphere        = 6;
  std::string sim_modtran_data      = "corsika6";
  float       sim_transition_energy = 80.0;
  float       sim_bx                = 25.2;
  float       sim_bz                = 40.83;
  float       sim_elec_mul_scat_fac = 1.0;
#ifdef CORSIKA_THINNING
  float       sim_thinning_energy   = 0.01;
#endif
  float       sim_cutoff_hadron     = 0.10;
  float       sim_cutoff_muon       = 0.30;
  float       sim_cutoff_electron   = 0.020;
  float       sim_cutoff_gamma      = 0.020;

  unsigned    cer_event_use         = 1;
  float       cer_lambda_lo         = 180.0;
  float       cer_lambda_hi         = 700.0;
  unsigned    cer_bunch             = 1;
  bool        cer_qeff              = true;
  bool        cer_atmo              = true;
  bool        cer_mirror            = true;
  std::string cer_choice_qeff       = "corsika";
  std::string cer_choice_atmo       = "corsika";
  std::string cer_choice_mirror     = "falcone";

  unsigned    out_max_event_print   = 0;
  float       out_min_print_energy  = 0.0;
  std::string out_particle_file     = "corsika_particle_NNNNN.dat";
  bool        out_particle_enable   = false;
  bool        out_table_enable      = false;
  std::string out_cerenkov_file     = "corsika_cerenkov_NNNNN.dat";

  bool        iact_no_impact_corr   = false;

#if 0  
  int         iact_max_print_tel    = 10;
  int         iact_max_print_evt    = 100;
  int         iact_skip_print       = 1;
  int         iact_skip_print2      = 100;
  int         iact_skip_offset2     = 1;
  int         iact_max_int_bunches  = 1000000;
  int         iact_io_buffer        = 200000000;
#endif

  // --------------------------------------------------------------------------
  // PARSE COMMAND LINE OPTIONS AND SET VARIABLES
  // --------------------------------------------------------------------------

  // CORSIKA OPTIONS

  options.findWithValue("debug", cor_debug, HELP_DEBUG);
  if(options.find("no_qgs", HELP_NO_QGS)!=VSOptions::FS_NOT_FOUND)
    cor_no_qgs = true;
  if(options.find("no_fluka", HELP_NO_FLUKA)!=VSOptions::FS_NOT_FOUND)
    cor_no_fluka = true;
  if(options.find("no_ceffic", HELP_NO_CEFFIC)!=VSOptions::FS_NOT_FOUND)
    cor_no_ceffic = true;
  if(options.find("no_modtran", HELP_NO_MODTRAN)!=VSOptions::FS_NOT_FOUND)
    cor_no_modtran = true, sim_atmosphere=1, sim_modtran_data="none";

  // RUN OPTIONS

  options.findWithValue("host",run_hostname,HELP_HOST);
  options.findWithValue("user",run_username,HELP_USER);
  options.findWithValue("seed1",run_seed1,HELP_SEED1);
  options.findWithValue("seed2",run_seed2,HELP_SEED2);
  options.findWithValue("seed3",run_seed3,HELP_SEED3);
  options.findWithValue("seed4",run_seed4,HELP_SEED4);
  options.findWithValue("runno",run_runno,HELP_RUNNO);
  options.findWithValue("eventno",run_eventno,HELP_EVENTNO);
  options.findWithValue("n",run_num_shower,HELP_N);

  // TELESCOPE (OPTICS) OPTIONS

  options.findWithValue("tel_optics_id",tel_optics_id,HELP_TEL_OPTICS_ID);
  options.findWithValue("tel_rings",tel_rings,HELP_TEL_RINGS);
  options.findWithValue("tel_spacing",tel_spacing,HELP_TEL_SPACING);
  options.findWithValue("tel_elevation",tel_elevation,HELP_TEL_ELEVATION);
  options.findWithValue("tel_aperture",tel_aperture,HELP_TEL_APERTURE);
  options.findWithValue("tel_max_dist",tel_max_dist,HELP_TEL_MAX_DIST);
  options.findWithValue("tel_layout",tel_layout,HELP_TEL_LAYOUT);

  pri_sampling_radius = tel_spacing/sqrt(3.0)*1.01;
  if(tel_layout == "square")
    pri_sampling_radius = tel_spacing/sqrt(2.0)*1.01;

  // PRIMARY OPTIONS

#define NUCLEUS(A,Z) ((A)*100+(Z))
  if(options.find("gamma",HELP_GAMMA)!=VSOptions::FS_NOT_FOUND)
    pri_species=1;
  if(options.find("proton",HELP_PROTON)!=VSOptions::FS_NOT_FOUND)
    pri_species=14;
  if(options.find("electron",HELP_ELECTRON)!=VSOptions::FS_NOT_FOUND)
    pri_species=3;
  if(options.find("muon",HELP_MUON)!=VSOptions::FS_NOT_FOUND)
    pri_species=6;
  if(options.find("helium",HELP_HELIUM)!=VSOptions::FS_NOT_FOUND)
    pri_species=NUCLEUS(4,2);
  if(options.find("iron",HELP_IRON)!=VSOptions::FS_NOT_FOUND)
    pri_species=NUCLEUS(56,26);
  options.findWithValue("primary",pri_species,HELP_PRIMARY);
  if(options.findWithValue("zenith",pri_zenith_lo,HELP_ZENITH)
     !=VSOptions::FS_NOT_FOUND)pri_zenith_hi=pri_zenith_lo;
  options.findWithValue("zenith_lo",pri_zenith_lo,HELP_ZENITH_LO);
  options.findWithValue("zenith_hi",pri_zenith_hi,HELP_ZENITH_HI);
  if(options.findWithValue("azimuth",pri_azimuth_lo,HELP_AZIMUTH)
     !=VSOptions::FS_NOT_FOUND)
    pri_azimuth_hi=pri_azimuth_lo;
  options.findWithValue("azimuth_lo",pri_azimuth_lo,HELP_AZIMUTH_LO);
  options.findWithValue("azimuth_hi",pri_azimuth_hi,HELP_AZIMUTH_HI);
  if(options.findWithValue("viewcone",pri_viewcone_hi,HELP_VIEWCONE)
     !=VSOptions::FS_NOT_FOUND)pri_viewcone_lo=0;
  options.findWithValue("viewcone_lo",pri_viewcone_lo,HELP_VIEWCONE_LO);
  options.findWithValue("viewcone_hi",pri_viewcone_hi,HELP_VIEWCONE_HI);
  options.findWithValue("sampling_radius",pri_sampling_radius,
			HELP_SAMPLING_RADIUS);
  options.findWithValue("e",pri_energy,HELP_ENERGY);
  options.findWithValue("starting_depth",pri_start_depth,HELP_STARTING_DEPTH);

  // SIMULATION OPTIONS

  options.findWithValue("atmosphere",sim_atmosphere,HELP_ATMOSPHERE);
  sim_modtran_data = 
    std::string("corsika") + VSDataConverter::toString(sim_atmosphere);
  options.findWithValue("modtran_data",sim_modtran_data,HELP_MODTRAN_DATA);
  options.findWithValue("transition_energy",sim_transition_energy,
			HELP_TRANSITION_ENERGY);
  if(options.find("magnetic_arizona",HELP_MAGNETIC_ARIZONA)
     !=VSOptions::FS_NOT_FOUND)sim_bx = 25.2, sim_bz = 40.83;
  if(options.find("magnetic_argentina",HELP_MAGNETIC_ARGENTINA)
     !=VSOptions::FS_NOT_FOUND)sim_bx = 31.7, sim_bz = -14.8;
  if(options.find("magnetic_equator",HELP_MAGNETIC_EQUATOR)
     !=VSOptions::FS_NOT_FOUND)sim_bx = 40.0, sim_bz = 0.01;
  options.findWithValue("magnetic_bx",sim_bx,HELP_MAGNETIC_BX);
  options.findWithValue("magnetic_bz",sim_bz,HELP_MAGNETIC_BZ);
  options.findWithValue("elec_mult_scat_factor",sim_elec_mul_scat_fac,
			HELP_ELEC_MULT_SCAT_FACT);
#ifdef CORSIKA_THINNING
  options.findWithValue("thinning_energy",sim_thinning_energy,
			HELP_THINNING_ENERGY);
#endif
  options.findWithValue("cutoff_hadron",sim_cutoff_hadron,
			HELP_CUTOFF_HADRON);
  options.findWithValue("cutoff_muon",sim_cutoff_muon,
			HELP_CUTOFF_MUON);
  options.findWithValue("cutoff_electron",sim_cutoff_electron,
			HELP_CUTOFF_ELECTRON);
  options.findWithValue("cutoff_gamma",sim_cutoff_gamma,
			HELP_CUTOFF_GAMMA);

  // CERENKOV OPTIONS

  options.findWithValue("cerenkov_event_use", cer_event_use,
			HELP_CERENKOV_EVENT_USE);
  options.findWithValue("cerenkov_lambda_lo", cer_lambda_lo,
			HELP_CERENKOV_LAMBDA_LO);
  options.findWithValue("cerenkov_lambda_hi", cer_lambda_hi,
			HELP_CERENKOV_LAMBDA_HI);
  options.findWithValue("cerenkov_bunch", cer_bunch,
			HELP_CERENKOV_BUNCH);
#ifdef CEFFIC
  options.findWithValue("cerenkov_qeff",cer_qeff,
			HELP_CERENKOV_QEFF);
  options.findWithValue("cerenkov_atmo",cer_atmo,
			HELP_CERENKOV_ATMO);
  options.findWithValue("cerenkov_mirror",cer_mirror,
			HELP_CERENKOV_MIRROR);
  options.findWithValue("qeff_data",cer_choice_qeff,
			HELP_QEFF_DATA);
  options.findWithValue("atmo_data",cer_choice_atmo,
			HELP_ATMO_DATA);
  options.findWithValue("mirror_data",cer_choice_mirror,
			HELP_MIRROR_DATA);
#endif

  // OUTPUT OPTIONS

  std::ostringstream of_par_stream;
  of_par_stream << "corsika_particle_" << run_runno << ".dat";
  out_particle_file = of_par_stream.str();

  std::ostringstream of_cer_stream;
  of_cer_stream << "corsika_cerenkov_" << run_runno << ".dat";
  out_cerenkov_file = of_cer_stream.str();

  options.findWithValue("max_event_print",out_max_event_print,
			HELP_MAX_EVENT_PRINT);
  options.findWithValue("min_print_energy",out_min_print_energy,
			HELP_MIN_PRINT_ENERGY);
  options.findWithValue("particle_output_enable",out_particle_enable,
			HELP_PARTICLE_OUTPUT_ENABLE);
  options.findWithValue("table_output_enable",out_table_enable,
			HELP_TABLE_OUTPUT_ENABLE);
  options.findWithValue("o",out_cerenkov_file,
			HELP_O);

  if(out_particle_enable == false)out_particle_file="/dev/null";
  options.findWithValue("particle_output_file",out_particle_file,
			HELP_PARTICLE_OUTPUT_FILE);

  // IACT OPTIONS

  if(options.find("no_impact_correction",HELP_NO_IMPACT_CORRECTION)!=
     VSOptions::FS_NOT_FOUND)
    iact_no_impact_corr = true;

  // --------------------------------------------------------------------------
  // ALL KNOWN OPTIONS HAVE BEEN PROCESSED
  // --------------------------------------------------------------------------

  if(argc>1)
    {
      std::cerr << progname << ": Unknown command line options:";
      argc--,argv++;
      while(argc)
	{
	  std::cerr << ' ' << *argv;
	  argc--,argv++;
	}
      std::cerr << std::endl;
      usage_error=true;
    }

  if(usage_error)
    {
      std::cerr << "Usage: " << progname << " [options]" << std::endl
		<< "Options:" << std::endl;
      options.printUsage(std::cerr);
      return EXIT_FAILURE;
    }

  // --------------------------------------------------------------------------
  // OVERRIDE COMMAND LINE OPTIONS OR DEFAULTS WITH VALUES FROM DATABASE
  // --------------------------------------------------------------------------

  if((ctl_db_use)&&(ctl_table_use)&&(!ctl_save))
    {
      VSSimDBTableParam* param = db_sim->getDataTableByName(ctl_table_name);
      if(param==0)
	{
	  std::cerr << progname 
		    << ": could not find table " << ctl_table_name
		    << std::endl;
	  return EXIT_FAILURE;
	}
      
      pri_species         = param->fPrimaryID;
      pri_energy          = param->fEnergyGeV;
      pri_zenith_lo       = param->fZenithMinRad/M_PI*180.0;
      pri_zenith_hi       = param->fZenithMaxRad/M_PI*180.0;
      pri_azimuth_lo      = param->fAzimuthMinRad/M_PI*180.0;
      pri_azimuth_hi      = param->fAzimuthMaxRad/M_PI*180.0;
      pri_sampling_radius = param->fSamplingRadiusM;
      tel_optics_id       = param->fOpticsID;
      run_num_shower      = param->fWorkunitEventCount;

      if(cer_event_use > 1) run_num_shower /= cer_event_use;

      delete param;
    }

  if((ctl_db_use)&&(!ctl_save))
    {
      VSPrimaryArrivalDistribution* pad =
	VSPADFactory::getInstance()->getPrimaryArrivalDistribution(db_param);

      if(pad)pad->setCORSIKAParameters(pri_zenith_lo, pri_zenith_hi,
				       pri_azimuth_lo, pri_azimuth_hi,
				       pri_viewcone_lo, pri_viewcone_hi);

      delete pad;
    }

  // --------------------------------------------------------------------------
  // VALIDATE STEERING FILE PARAMETERS
  // --------------------------------------------------------------------------

  if(pri_zenith_hi < pri_zenith_lo)
    {
      std::cerr << progname << ": Upper bound on zenith range (" 
		<< pri_zenith_hi 
		<< ") must not be smaller than lower bound ("
		<< pri_zenith_lo 
		<< ")" << std::endl;
      return EXIT_FAILURE;
    }

  if(pri_azimuth_hi < pri_azimuth_lo)
    {
      std::cerr << progname << ": Upper bound on azimuth range (" 
		<< pri_azimuth_hi 
		<< ") must not be smaller than lower bound ("
		<< pri_azimuth_lo 
		<< ")" << std::endl;
      return EXIT_FAILURE;
    }

  // --------------------------------------------------------------------------
  // WRITE OPTIONS TO DATABASE IF REQUESTED
  // --------------------------------------------------------------------------

  if((ctl_save)&&(ctl_db_use))
    {
      VSDBParameterSet param;
#define SETPARAM(x,y) param[(x)]=VSDataConverter::toString(y)

      SETPARAM("debug",cor_debug);
      if(cor_no_qgs)param["no_qgs"]="";
      if(cor_no_fluka)param["no_fluka"]="";
      if(cor_no_ceffic)param["no_ceffic"]="";
      if(cor_no_modtran)param["no_modtran"]="";

      SETPARAM("tel_optics_id",tel_optics_id);
      SETPARAM("tel_rings",tel_rings);
      SETPARAM("tel_spacing",tel_spacing);
      SETPARAM("tel_elevation",tel_elevation);
      SETPARAM("tel_aperture",tel_aperture);
      SETPARAM("tel_max_dist",tel_max_dist);
      SETPARAM("tel_layout",tel_layout);
      
      SETPARAM("primary",pri_species);
      SETPARAM("zenith_lo",pri_zenith_lo);
      SETPARAM("zenith_hi",pri_zenith_hi);
      SETPARAM("azimuth_lo",pri_azimuth_lo);
      SETPARAM("azimuth_hi",pri_azimuth_hi);
      SETPARAM("viewcone_lo",pri_viewcone_lo);
      SETPARAM("viewcone_hi",pri_viewcone_hi);
      SETPARAM("sampling_radius",pri_sampling_radius);
      SETPARAM("e",pri_energy);
      SETPARAM("starting_depth",pri_start_depth);

      SETPARAM("atmosphere",sim_atmosphere);
      SETPARAM("modtran_data",sim_modtran_data);
      SETPARAM("transition_energy",sim_transition_energy);
      SETPARAM("magnetic_bx",sim_bx);
      SETPARAM("magnetic_bz",sim_bz);
      SETPARAM("elec_mult_scat_factor",sim_elec_mul_scat_fac);
#ifdef CORSIKA_THINNING
      SETPARAM("thinning_energy",sim_thinning_energy);
#endif
      SETPARAM("cutoff_hadron",sim_cutoff_hadron);
      SETPARAM("cutoff_muon",sim_cutoff_muon);
      SETPARAM("cutoff_electron",sim_cutoff_electron);
      SETPARAM("cutoff_gamma",sim_cutoff_gamma);

      SETPARAM("cerenkov_event_use", cer_event_use);
      SETPARAM("cerenkov_lambda_lo", cer_lambda_lo);
      SETPARAM("cerenkov_lambda_hi", cer_lambda_hi);
      SETPARAM("cerenkov_bunch", cer_bunch);
      SETPARAM("cerenkov_qeff",cer_qeff);
      SETPARAM("cerenkov_atmo",cer_atmo);
      SETPARAM("cerenkov_mirror",cer_mirror);
      SETPARAM("qeff_data",cer_choice_qeff);
      SETPARAM("atmo_data",cer_choice_atmo);
      SETPARAM("mirror_data",cer_choice_mirror);

      SETPARAM("max_event_print",out_max_event_print);
      SETPARAM("min_print_energy",out_min_print_energy);
      SETPARAM("particle_output_enable",out_particle_enable);
      SETPARAM("table_output_enable",out_table_enable);
      SETPARAM("particle_output_file",out_particle_file);

      if(iact_no_impact_corr)param["no_impact_correction"]="";

      db_param->deleteParameterSet("MakeSteering");
      db_param->storeParameterSet("MakeSteering",param);
    }

  if(ctl_save)return EXIT_SUCCESS;

  // --------------------------------------------------------------------------
  // SAMPLE UNIQUE TARGET DIRECTION IF NECESSARY
  // --------------------------------------------------------------------------
  
  if(fabs(pri_viewcone_hi - pri_viewcone_lo)>0)
    {
      float cos_zn_lo = cos(pri_zenith_lo*M_PI/180.0);
      float cos_zn_hi = cos(pri_zenith_hi*M_PI/180.0);
      float x = float(random())/float(RAND_MAX);
      float cos_zn = cos_zn_lo*(1.0-x) + cos_zn_hi*x;
      pri_zenith_hi=pri_zenith_lo=acos(cos_zn)/M_PI*180.0;

      float y = float(random())/float(RAND_MAX);
      float az = pri_azimuth_lo*(1.0-y) + pri_azimuth_hi*y;
      pri_azimuth_hi=pri_azimuth_lo=az;
    }
  
  // --------------------------------------------------------------------------
  // CONVERT UNITS TO THOSE USED BY CORSIKA AND IMPOSE LIMITS ON STRINGS
  // --------------------------------------------------------------------------

  pri_azimuth_hi = 180 - pri_azimuth_hi;
  pri_azimuth_lo = 180 - pri_azimuth_lo;
  std::swap(pri_azimuth_hi, pri_azimuth_lo);
  if(pri_azimuth_hi-pri_azimuth_lo >= 360.0)
    pri_azimuth_hi=180.0, pri_azimuth_lo=-180.0;
  if(pri_azimuth_lo < -360.0)
    pri_azimuth_hi += floor(fabs(pri_azimuth_lo)/360.0)*360.0,
      pri_azimuth_lo += floor(fabs(pri_azimuth_lo)/360.0)*360.0;
  if(pri_azimuth_hi > 360.0)
    pri_azimuth_lo -= floor(fabs(pri_azimuth_hi)/360.0)*360.0,
      pri_azimuth_hi -= floor(fabs(pri_azimuth_hi)/360.0)*360.0;  

  pri_sampling_radius *= 100.0;
  tel_spacing       *= 100.0;
  tel_elevation     *= 100.0;
  tel_aperture      *= 100.0;
  tel_max_dist      *= 100.0;

#define LIMIT(x,y) if(x.size()>y)x=x.substr(0,y)
  LIMIT(run_hostname,20);
  LIMIT(run_username,20);
  LIMIT(out_particle_file,100);
  LIMIT(out_cerenkov_file,100);

  // --------------------------------------------------------------------------
  // SET UP VARIOUS ARRAYS AND STRUCTURES USED IN WRITING CORSIKA OPTIONS
  // --------------------------------------------------------------------------

  std::vector<unsigned> _seed1(3,0); _seed1[0]=run_seed1; _seed1[1]=100000;
  std::vector<unsigned> _seed2(3,0); _seed2[0]=run_seed2; _seed2[1]=100000;
  std::vector<unsigned> _seed3(3,0); _seed3[0]=run_seed3; _seed3[1]=100000;
  std::vector<unsigned> _seed4(3,0); _seed4[0]=run_seed4; _seed4[1]=100000;
  std::vector<float> thetap(2);
  thetap[0]=pri_zenith_lo;
  thetap[1]=pri_zenith_hi;
  std::vector<float> phip(2);
  phip[0]=pri_azimuth_lo;
  phip[1]=pri_azimuth_hi;
  std::vector<float> vuecon(2);
  vuecon[0]=pri_viewcone_lo;
  vuecon[1]=pri_viewcone_hi;
  std::vector<float> erange(2,pri_energy);
  std::vector<float> magnet; 
  magnet.push_back(sim_bx); 
  magnet.push_back(sim_bz);
  std::vector<float> cutoff;
  cutoff.push_back(sim_cutoff_hadron);
  cutoff.push_back(sim_cutoff_muon);
  cutoff.push_back(sim_cutoff_electron);
  cutoff.push_back(sim_cutoff_gamma);
  std::vector<bool> parout; 
  parout.push_back(out_particle_enable);
  parout.push_back(out_table_enable);
  std::vector<float> cwavlg; 
  cwavlg.push_back(cer_lambda_lo); 
  cwavlg.push_back(cer_lambda_hi);
  std::vector<bool> cerqef;
  cerqef.push_back(cer_qeff);
  cerqef.push_back(cer_atmo);
  cerqef.push_back(cer_mirror);
  CSCAT_Struct cscat;
  cscat.event_use = cer_event_use;
  cscat.sampling_radius = pri_sampling_radius;
#ifdef QGSJET
  QGSJET_Struct qgs; 
  qgs.debug=(cor_debug>0)?(cor_debug-1):0; 
  if(qgs.debug>4)qgs.debug=4;
#endif
  Modtran_Struct modtran;
  modtran.atmosphere=sim_atmosphere;
  modtran.refraction=false;

  // --------------------------------------------------------------------------
  // SET UP TELESCOPE ARRAY (AUTOMATICALLY GENERATED OR FROM DATABASE)
  // --------------------------------------------------------------------------

  std::vector<std::vector<float> > telescopes;
  if(ctl_db_use)
    {
      VSOTelescopeArray array;
      array.readFromDatabase(db,tel_optics_id);

      assert(array.numTelescopes());
      tel_elevation = 
	array.telescope(0)->altitude()-array.telescope(0)->reflectorIP()/2.0;
      std::vector<float> telescope;
      for(unsigned iscope = 0; iscope<array.numTelescopes(); iscope++)
	{
	  const VSOTelescope* scope = array.telescope(iscope);
	  Physics::Vec3D pos = scope->position();
	  std::vector<float> telescope;	  
	  // CORSIKA coordinate system is rotated WRT optical simulation
	  telescope.push_back(pos.y);
	  telescope.push_back(-pos.x);
	  telescope.push_back(pos.z);
	  telescope.push_back(scope->reflectorIP()/2.0);
	  telescopes.push_back(telescope);
	  if((pos.z-scope->reflectorIP()/2.0) < tel_elevation)
	    tel_elevation = pos.z-scope->reflectorIP()/2.0;
	}

      for(std::vector<std::vector<float> >::iterator itel = 
	    telescopes.begin(); itel != telescopes.end(); itel++)
	(*itel)[2]-=tel_elevation;
    }
  else
    {
      unsigned ntel = 3*tel_rings*(tel_rings+1)+1;
      if(tel_layout=="square")ntel = (2*tel_rings+1)*(2*tel_rings+1);

      for(unsigned itel=0; itel<ntel; itel++)
	{
	  double x;
	  double y;
	  double z=0;

	  int n=itel+1;	
	  if(tel_layout=="square")ns_to_xy(n, &x, &y);
	  else nh_to_xy(&n, &x, &y);
	  
	  x *= tel_spacing;
	  y *= tel_spacing;
	  double r2 = x*x+y*y;

	  if(std::isfinite(tel_max_dist) && r2>tel_max_dist*tel_max_dist)
	    continue;
	  
	  std::vector<float> telescope;
	  telescope.push_back(x);
	  telescope.push_back(y);
	  telescope.push_back(z);
	  telescope.push_back(tel_aperture/2.0);
	  telescopes.push_back(telescope);
	}
    }

  // --------------------------------------------------------------------------
  // WRITE CORSIKA FILES
  // --------------------------------------------------------------------------

  if((ctl_init)&&(db))
    {
      VSSimDBCORSIKADatasets db_datasets(db);

      VSSimDBWavelengthDataset*          ds_quantum_efficiency(0);
      VSSimDBWavelengthDataset*          ds_mirror_reflectivity(0);
      VSSimDBWavelengthAltitudeDataset*  ds_atmospheric_absorption(0);
      VSSimDBModtranProfileDataset*      ds_modtran_profile(0);

      ds_quantum_efficiency =
	db_datasets.retrieveWavelengthDataset(VSIMDB_TABLE_NAME_QUANEFF);

      ds_mirror_reflectivity = 
	db_datasets.retrieveWavelengthDataset(VSIMDB_TABLE_NAME_MIRRREF);

      ds_atmospheric_absorption =
	db_datasets.
	retrieveWavelengthAltitudeDataset(VSIMDB_TABLE_NAME_ATMOABS);

      ds_modtran_profile =
	db_datasets.retrieveModtranProfileDataset(VSIMDB_TABLE_NAME_MODTRAN);

      if(ds_quantum_efficiency)
	{
	  ds_quantum_efficiency->writeToCORSIKA("quanteff.dat");
	  cer_choice_qeff="none";
	  delete ds_quantum_efficiency;
	}

      if(ds_mirror_reflectivity)
	{
	  ds_mirror_reflectivity->writeToCORSIKA("mirreff.dat");
	  cer_choice_mirror="none";
	  delete ds_mirror_reflectivity;
	}

      if(ds_atmospheric_absorption)
	{
	  ds_atmospheric_absorption->writeToCORSIKA("atmabs.dat");
	  cer_choice_atmo="none";
	  delete ds_atmospheric_absorption;
	}

      if(ds_modtran_profile)
	{
	  ds_modtran_profile->writeToCORSIKA("atmprof99.dat");
	  sim_modtran_data="none";
	  cor_no_modtran=false;
	  modtran.atmosphere=99;
	  delete ds_modtran_profile;
	}
    }

  writeAllCORSIKADataFiles(!ctl_init,
			   sim_modtran_data,
			   cer_choice_qeff, 
			   cer_choice_atmo,
			   cer_choice_mirror);
  
  // --------------------------------------------------------------------------
  // WRITE CORSIKA STEERING FILE
  // --------------------------------------------------------------------------

  std::ostream* stream = &std::cout;

  write_card_option(*stream,   "RUNNR",        run_runno);              // 4.1
  write_card_option(*stream,   "EVTNR",        run_eventno) ;           // 4.2
  write_card_option(*stream,   "SEED",         _seed1);                 // 4.3
  write_card_option(*stream,   "SEED",         _seed2);                 // 4.3
  write_card_option(*stream,   "SEED",         _seed3);                 // 4.3
  write_card_option(*stream,   "SEED",         _seed4);                 // 4.3
  write_card_option(*stream,   "NSHOW",        run_num_shower);         // 4.4
  write_card_option(*stream,   "PRMPAR",       pri_species);            // 4.5
  write_card_option(*stream,   "ERANGE",       erange);                 // 4.6
  write_card_option(*stream,   "THETAP",       thetap);                 // 4.8
  write_card_option(*stream,   "PHIP",         phip);                   // 4.9
  write_card_option(*stream,   "VIEWCONE",     vuecon);                 // 4.10
  write_card_option(*stream,   "FIXCHI",       pri_start_depth);        // 4.11
  if(cor_no_modtran)
    write_card_option(*stream, "ATMOD",        sim_atmosphere);         // 4.15
  else
    write_card_option(*stream, "ATMOSPHERE",   modtran);                // 4.20
  write_card_option(*stream,   "MAGNET",       magnet);                 // 4.21

#ifdef QGSJET
  if(!cor_no_qgs)
    write_card_option(*stream, "QGSJET",       qgs);                    // 4.28
#endif

  write_card_option(*stream,   "HILOW",        sim_transition_energy);  // 4.36
  write_card_option(*stream,   "ELMFLG",       "F T");                  // 4.37
  write_card_option(*stream,   "STEPFC",       sim_elec_mul_scat_fac);  // 4.38
#ifdef CORSIKA_THINNING
  write_card_option(*stream,   "THIN", sim_thinning_energy/erange[1]);  // 4.40
#endif
  write_card_option(*stream,   "ECUTS",        cutoff);                 // 4.43
  write_card_option(*stream,   "OBSLEV",       tel_elevation);          // 4.47
  write_card_option(*stream,   "MAXPRT",       out_max_event_print);    // 4.50
  write_card_option(*stream,   "ECTMAP",       out_min_print_energy);   // 4.51
  write_card_option(*stream,   "DIRECT",       out_particle_file);      // 4.52
  write_card_option(*stream,   "PAROUT",       parout);                 // 4.53
  write_card_option(*stream,   "CWAVLG",       cwavlg);                 // 4.57
  write_card_option(*stream,   "CERSIZ",       (float)cer_bunch);       // 4.58
  write_card_option(*stream,   "CERFIL",       true);                   // 4.59

#ifdef CEFFIC
  if(!cor_no_ceffic)
    write_card_option(*stream, "CERQEF",       cerqef);                 // 4.60
#endif

  write_card_option(*stream,   "CSCAT",        cscat);                  // 4.61

  for(std::vector<std::vector<float> >::const_iterator itel = 
	telescopes.begin(); itel != telescopes.end(); itel++)
    write_card_option(*stream, "TELESCOPE",    *itel);                  // 4.62

  write_card_option(*stream,   "TELFIL",       out_cerenkov_file);      // 4.63
  write_card_option(*stream,   "USER",         run_username);           // 4.65
  write_card_option(*stream,   "HOST",         run_hostname);           // 4.66

  if(cor_debug>0)
    write_card_option(*stream, "DEBUG",        "T 6 F 0");              // 4.67
  if(cor_debug>1)
    write_card_option(*stream, "EGSDEB",       unsigned(0));            // 4.68

#ifdef FLUKA
  if((!cor_no_fluka)&&(cor_debug>1))
    write_card_option(*stream, "FLUDBG",       true);                   // 4.69
#endif  

  if(iact_no_impact_corr)
    write_card_option(*stream, "IACT","impact_correction F");           // IACT

  write_card_option(*stream,   "EXIT");                                 // 4.79

  // --------------------------------------------------------------------------
  // CLEAN UP AND WE ARE FINISHED
  // --------------------------------------------------------------------------

  delete db_sim;
  delete db_param;
  delete db;

  return(EXIT_SUCCESS);
}

