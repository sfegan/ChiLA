//-*-mode:c++; mode:font-lock;-*-

#include <sstream>
#include <cassert>

#include <VSDataConverter.hpp>
#include <VSSimDBTables.hpp>

#include "VSOPixel.hpp"
#include "VSOTelescope.hpp"

using namespace Physics;
using namespace VERITAS;

VSOPixel::VSOPixel():
  fTelescope(0), fID(0), fHexID(0), fRemoved(false), fPos()
{
  // nothing to see here
}

VSOPixel::VSOPixel(const VSOTelescope* T, 
		   unsigned PID, unsigned PHID, bool REM,
		   const Physics::Vec3D& P):
  fTelescope(T), fID(PID), fHexID(PHID), fRemoved(REM), fPos(P)
{
  // nothing to see here
}

VSOPixel::~VSOPixel()
{
  // nothing to see here
}

Vec3D VSOPixel::incomingSkyVectorAtZenith(double plate_scale) const
{
  double x = fPos.x*plate_scale;
  double y = -fPos.z*plate_scale;
  double z = -(fTelescope->focalPlanePosition().y+fPos.y);

  Vec3D p(x,y,z);
  p/=p.Norm();

  return p;
}

VSDBStatement* VSOPixel::createInsertQuery(VSDatabase* db)
{
  return db->createInsertQuery(VSIMDB_TABLE_NAME_PIXEL,8);
}

VSDBStatement* VSOPixel::createSelectQuery(VSDatabase* db)
{
  return db->createSelectQuery(VSIMDB_TABLE_NAME_PIXEL,
			       "OpticsID=? AND TelescopeID=?");
}

void VSOPixel::writeToDatabase(VSDBStatement *stmt, uint32_t optics_id) const
{
  stmt->bindToParam(optics_id);

  stmt->bindToParam(fTelescope->id());
  stmt->bindToParam(fID);
  stmt->bindToParam(fHexID);
  stmt->bindToParam(fRemoved);
  stmt->bindToParam(fPos.x);

  stmt->bindToParam(fPos.y);
  stmt->bindToParam(fPos.z);
  
  assert(stmt->execute() == 1);
}

VSOPixel* VSOPixel::createFromDatabaseRow(VSDBStatement* stmt, 
					  const VSOTelescope* T)
{
  VSOPixel* pixel = new VSOPixel;
  unsigned OID;
  unsigned TID;
  
  stmt->bindToResult(OID);

  stmt->bindToResult(TID);
  stmt->bindToResult(pixel->fID);
  stmt->bindToResult(pixel->fHexID);
  stmt->bindToResult(pixel->fRemoved);
  stmt->bindToResult(pixel->fPos.x);

  stmt->bindToResult(pixel->fPos.y);
  stmt->bindToResult(pixel->fPos.z);

  if(stmt->retrieveNextRow() == 1)
    {
      pixel->fTelescope=T;      
      assert(TID == pixel->fTelescope->id());
    }
  else
    {
      delete pixel;
      pixel = 0;
    }

  return pixel;
}

void VSOPixel::createPixelTable(VSDatabase* db)
{
  // stupid to create this --- C++ should have a "typeof" operator
  // (g++ does but it is not portable .. and doesn't work in a static
  // function anyway)
  VSOPixel temp; 
  unsigned id;
  
  db->createTable(VSIMDB_TABLE_NAME_PIXEL,

		  db->sqlSpecOf("OpticsID",id,true,"NOT NULL")+

		  db->sqlSpecOf("TelescopeID",id,false,"NOT NULL")+
		  db->sqlSpecOf("PixelID",temp.fID,false,"NOT NULL")+
		  db->sqlSpecOf("PixelHexID",temp.fHexID,false,"NOT NULL")+
		  db->sqlSpecOf("Removed",temp.fRemoved,false,"NOT NULL")+
		  db->sqlSpecOf("PosX",temp.fPos.x,false,"NOT NULL")+

		  db->sqlSpecOf("PosY",temp.fPos.y,false,"NOT NULL")+
		  db->sqlSpecOf("PosZ",temp.fPos.z,false,"NOT NULL")+

		  ", PRIMARY KEY ( OpticsID, TelescopeID, PixelID )",

		  VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
}

void VSOPixel::dumpShort(std::ostream& stream) const
{
  stream
    << "PIXEL "

    << VSDataConverter::toString(fTelescope->id()) << ' '
    << VSDataConverter::toString(fID) << ' '
    << VSDataConverter::toString(fHexID) << ' '
    << VSDataConverter::toString(fRemoved) << ' '
    << VSDataConverter::toString(fPos.x) << ' '

    << VSDataConverter::toString(fPos.y) << ' '
    << VSDataConverter::toString(fPos.z) << std::endl;
}

void VSOPixel::dump(std::ostream& stream, unsigned l) const
{
  stream
    << FDAV("Telescope ID", fTelescope->id(), "", 30, l) << std::endl
    << FDAV("Pixel ID", fID, "", 30, l) << std::endl
    << FDAV("Pixel Hex ID", fHexID, "", 30, l) << std::endl
    << FDAV("Removed", fRemoved, "", 30, l) << std::endl
    << FDAV("Position X", fPos.x, "", 30, l) << std::endl

    << FDAV("Position Y", fPos.y, "", 30, l) << std::endl
    << FDAV("Position Z", fPos.z, "", 30, l) << std::endl;
}

VSOPixel* VSOPixel::createFromShortDump(std::istream& stream,
					  const VSOTelescope* T)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return 0;

  std::istringstream linestream(line);

  VSOPixel* pixel = new VSOPixel;
  pixel->fTelescope=T;

  std::string keyword;
  linestream >> keyword;
  assert(keyword==std::string("PIXEL"));

  unsigned TID;

  linestream
    >> TID
    >> pixel->fID
    >> pixel->fHexID
    >> pixel->fRemoved
    >> pixel->fPos.x

    >> pixel->fPos.y
    >> pixel->fPos.z;

  if(!linestream)
    {
      delete pixel;
      return 0;
    }

  return pixel;
}
