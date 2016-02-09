//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDBFactory.hpp

  Factory class to magically create database connection extracting the
  DB location user name and password from the config file, command
  line or environment variables.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       21/02/2005
*/

#ifndef VSDBFACTORY_HPP
#define VSDBFACTORY_HPP

#include <string>
#include <memory>

#include "VSDatabase.hpp"

#include "VSDBMySQL3X.hpp"
#include "VSDBMySQL41.hpp"

namespace VERITAS
{

  //! Magically create database connection
  class VSDBFactory
  {
  public:
    virtual ~VSDBFactory();
    VSDatabase* createVSDB();
    static VSDBFactory* getInstance();
    static void configure(VSOptions* options = 0, 
			  const std::string catagory="");
    static void setDefaultDBI(const std::string& dbi) { sDefaultDBI=dbi; }
    static void setDefaultLoud(bool loud) { sDefaultLoud=loud; }
    static void setDefaultHost(const std::string& host) { sDefaultHost=host; }
    static void setDefaultPort(unsigned short port) { sDefaultPort=port; }
    static void setDefaultUser(const std::string& user) { sDefaultUser=user; }
    static void setDefaultPass(const std::string& pass) { sDefaultPass=pass; }
    static void setDefaultUnix(const std::string& _unx) { sDefaultUnix=_unx; }
    static std::string getDefaultDBI();
    static bool getDefaultLoud();
    static std::string getDefaultHost();
    static unsigned short getDefaultPort();
    static std::string getDefaultUser();
    static std::string getDefaultPass();
    static std::string getDefaultUnix();
  private:
    VSDBFactory();
    VSDBFactory(const VSDatabase&);
    VSDBFactory& operator =(const VSDatabase&);

    enum DBInterface { DBI_MYSQL, DBI_MYSQL3X, DBI_MYSQL41 };

    DBInterface                             fDBI;
    bool                                    fLoud;
    std::string                             fHost;
    unsigned short                          fPort;
    std::string                             fUser;
    std::string                             fPass;
    std::string                             fUnix;
    
    static std::auto_ptr<VSDBFactory>       sInstance;

    static std::string                      sDefaultDBI;
    static bool                             sDefaultLoud;
    static std::string                      sDefaultHost;
    static unsigned short                   sDefaultPort;
    static std::string                      sDefaultUser;
    static std::string                      sDefaultPass;
    static std::string                      sDefaultUnix;
  };

};

#endif // VSDBFACTORY_HPP
