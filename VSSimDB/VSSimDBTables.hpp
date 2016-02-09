//-*-mode:c++; mode:font-lock;-*-

/*! \file VSSimDBTables.hpp

  Definitions of TABLE names in the simulations database

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       09/01/2005
  \note
*/

#ifndef VSIMDBTABLES_HPP
#define VSIMDBTABLES_HPP

#define VSIMDB_TABLE_NAME_PARAMETERS  "VS_Parameters"

#define VSIMDB_TABLE_PREFIX_OPTICS    "VSO_"
#define VSIMDB_TABLE_NAME_ARRAY       VSIMDB_TABLE_PREFIX_OPTICS "Array"
#define VSIMDB_TABLE_NAME_TELESCOPE   VSIMDB_TABLE_PREFIX_OPTICS "Telescopes"
#define VSIMDB_TABLE_NAME_MIRROR      VSIMDB_TABLE_PREFIX_OPTICS "Mirrors"
#define VSIMDB_TABLE_NAME_PIXEL       VSIMDB_TABLE_PREFIX_OPTICS "Pixels"

#define VSIMDB_TABLE_PREFIX_DATA      "VSD_"
#define VSIMDB_TABLE_NAME_DIRECTORY   VSIMDB_TABLE_PREFIX_DATA "TableDirectory"
#define VSIMDB_TABLE_NAME_WORKUNIT    VSIMDB_TABLE_PREFIX_DATA "WorkunitRun"
#define VSIMDB_TABLE_POSTFIX_EVENTS   "_Events"
#define VSIMDB_TABLE_POSTFIX_SCOPES   "_Scopes"
#define VSIMDB_TABLE_POSTFIX_PES      "_PEs"

#define VSIMDB_TABLE_PREFIX_SHOWER    "VSS_"
#define VSIMDB_TABLE_NAME_QUANEFF     VSIMDB_TABLE_PREFIX_SHOWER "QuanEff"
#define VSIMDB_TABLE_NAME_ATMOABS     VSIMDB_TABLE_PREFIX_SHOWER "AtmoAbs"
#define VSIMDB_TABLE_NAME_MIRRREF     VSIMDB_TABLE_PREFIX_SHOWER "MirrRef"
#define VSIMDB_TABLE_NAME_MODTRAN     VSIMDB_TABLE_PREFIX_SHOWER "Modtran"

#endif // VSIMDBTABLES_HPP
