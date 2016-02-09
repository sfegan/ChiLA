//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBFactory.cpp

  Factory class to magically create database connection extracting the
  DB location user name and password from the config file, command
  line or environment variables.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       21/02/2005
*/

#include "VSOptions.hpp"
#include "VSDBFactory.hpp"

using namespace VERITAS;

std::auto_ptr<VSDBFactory> VSDBFactory::sInstance;

#ifndef VSCONFIG_NO_MYSQL
std::string     VSDBFactory::sDefaultDBI("MYSQL");
#endif

bool            VSDBFactory::sDefaultLoud(false);
std::string     VSDBFactory::sDefaultHost("");
unsigned short  VSDBFactory::sDefaultPort(0);
std::string     VSDBFactory::sDefaultUser("");
std::string     VSDBFactory::sDefaultPass("");
std::string     VSDBFactory::sDefaultUnix("");

VSDBFactory::~VSDBFactory()
{
  // nothing to see here
}

VSDatabase* VSDBFactory::createVSDB()
{
  VSDatabase* db(0);
  switch(fDBI)
    {
    case VSDBFactory::DBI_MYSQL:
#ifndef VSCONFIG_NO_MYSQL
#if(MYSQL_VERSION_ID>=40100)
      db = new VSDBMySQL41(fLoud,false,fHost,fUser,fPass,fPort,fUnix,0,true);
#else
      db = new VSDBMySQL3X(fLoud,false,fHost,fUser,fPass,fPort,fUnix,0,true);
#endif
#endif
      break;

    case VSDBFactory::DBI_MYSQL41:
#ifndef VSCONFIG_NO_MYSQL
#if(MYSQL_VERSION_ID>=40100)
      db = new VSDBMySQL41(fLoud,false,fHost,fUser,fPass,fPort,fUnix,0,true);
#endif
#endif
      break;
      
    case VSDBFactory::DBI_MYSQL3X:
#ifndef VSCONFIG_NO_MYSQL
      db = new VSDBMySQL3X(fLoud,false,fHost,fUser,fPass,fPort,fUnix,0,true);
#endif
      break;
    }
  return db;
}

VSDBFactory* VSDBFactory::getInstance()
{
  if(sInstance.get() == 0)sInstance.reset(new VSDBFactory());
  return sInstance.get();
}

template<typename T> static bool stringToData(T& x, const std::string& s)
{
  std::istringstream stream(s);
  return(stream >> x);
}

template<> bool stringToData<>(std::string& x, const std::string& s)
{
  x=s;
  return true;
}

template<typename T> static bool 
getOptEnvDef(T& x, VSOptions* options, const char* opt, const char* env,
	     const T& def, const char* help, const std::string catagory)
{
  if((options)&&
     (options->findWithValue(opt, x, help, catagory)==VSOptions::FS_FOUND))
    {
      return true;
    }
  else if(char* env_val_ptr = getenv(env))
    {
      std::string env_val(env_val_ptr);
      if(stringToData(x,env_val.substr(env_val.find('=')+1)))return true;
      else { x = def; return false; }
    }
  else
    {
      x=def;
      return false;
    }
}

template<typename T>
static inline void
getEnvDef(T& x, const std::string& env)
{
  if(!env.empty())
    {
      char* env_val_ptr = getenv(env.c_str());
      if(env_val_ptr)
	{
	  std::string env_val(env_val_ptr);
	  stringToData(x,env_val.substr(env_val.find('=')+1));
	}
    }
}

#define HELP_VSDBRDBMS "Set relational database manager. Valid choices are: " \
  "\"MYSQL\", use MySQL 3.X or 4.1 database driver (default); " \
  "\"MYSQL3X\", use MySQL 3.X interface; " \
  "\"MYSQL41\", use MySQL 4.1 interface which takes advantage of server " \
  "side prepared statements. This option overrides the VSDB_RDBMS " \
  "environment variable"
  
#define HELP_VSDBLOUD "Write all database transactions to the console. " \
  "This option overrides the VSDB_LOUD environment variable"

#define HELP_VSDBHOST "Set the hostname to connect to. " \
  "This option overrides the VSDB_HOST environment variable"

#define HELP_VSDBUSER "Set the username to use when connecting to the " \
  "database. This option overrides the VSDB_USER environment variable"

#define HELP_VSDBPASS "Set the password to use when connecting to the " \
  "database. This option overrides the VSDB_PASS environment variable"

#define HELP_VSDBPORT "Set to port number to connect to. " \
  "This option overrides the VSDB_PORT environment variable"

#define HELP_VSDBUNIX "Set to UNIX socket to connect to. This is probably " \
  "only meaningful if the RDBMS is on the local machine. " \
  "This option overrides the VSDB_UNIX environment variable"

void VSDBFactory::configure(VSOptions* options, const std::string catagory)
{
  getOptEnvDef(sDefaultDBI, options, "VSDBRdbms", "VSDB_RDBMS", 
	       std::string("MYSQL"),HELP_VSDBRDBMS,catagory);

  if(((options)&&(options->find("VSDBLoud",HELP_VSDBLOUD,catagory)
		  !=VSOptions::FS_NOT_FOUND))||(getenv("VSDB_LOUD")))
    sDefaultLoud=true;

  getOptEnvDef(sDefaultHost, options, "VSDBHost", "VSDB_HOST", 
	       std::string(""),HELP_VSDBHOST,catagory);

  getOptEnvDef(sDefaultUser, options, "VSDBUser", "VSDB_USER", 
	       std::string(""),HELP_VSDBUSER,catagory);

  getOptEnvDef(sDefaultPass, options, "VSDBPass", "VSDB_PASS", 
	       std::string(""),HELP_VSDBPASS,catagory);

  getOptEnvDef(sDefaultPort, options, "VSDBPort", "VSDB_PORT", 
	       (unsigned short)0,HELP_VSDBPORT,catagory);

  getOptEnvDef(sDefaultUnix, options, "VSDBUnix", "VSDB_UNIX", 
	       std::string(""),HELP_VSDBUNIX,catagory);
}

std::string VSDBFactory::getDefaultDBI()
{
  return sDefaultDBI;
}

bool VSDBFactory::getDefaultLoud()
{
  return sDefaultLoud;
}

std::string VSDBFactory::getDefaultHost()
{
  return sDefaultHost;
}

unsigned short VSDBFactory::getDefaultPort()
{ 
  return sDefaultPort;
}

std::string VSDBFactory::getDefaultUser()
{
  return sDefaultUser;
}

std::string VSDBFactory::getDefaultPass()
{
  return sDefaultPass;
}

std::string VSDBFactory::getDefaultUnix()
{
  return sDefaultUnix;
}

VSDBFactory::VSDBFactory():
  fDBI(DBI_MYSQL), fLoud(sDefaultLoud), 
  fHost(sDefaultHost), fPort(sDefaultPort), 
  fUser(sDefaultUser), fPass(sDefaultPass), fUnix(sDefaultUnix)
{
  if(sDefaultDBI=="MYSQL")fDBI=DBI_MYSQL;
  else if(sDefaultDBI=="MYSQL3X")fDBI=DBI_MYSQL3X;
  else if(sDefaultDBI=="MYSQL41")fDBI=DBI_MYSQL41;
}

