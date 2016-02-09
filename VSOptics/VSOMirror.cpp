//-*-mode:c++; mode:font-lock;-*-

#include <sstream>
#include <cassert>

#include <VSDataConverter.hpp>
#include <VSOMirror.hpp>
#include <VSOTelescope.hpp>

using namespace Physics;
using namespace VERITAS;

VSOMirror::VSOMirror():
  fTelescope(0), fID(0), fHexID(0), fRemoved(false), fPos(), fAlign(), 
  fFocalLength(0), fSpotSize(0), fDegradingFactor(0), fRotationVector()
{
  // nothing to see here
}

VSOMirror::VSOMirror(const VSOTelescope* T, 
		     unsigned ID, unsigned MHID, bool REM,
		     const Physics::Vec3D& P, const Physics::Vec3D& A, 
		     double FL, double SS, double DF):
  fTelescope(T), fID(ID), fHexID(MHID), fRemoved(REM), fPos(P), fAlign(A), 
  fFocalLength(FL), fSpotSize(SS), fDegradingFactor(DF), fRotationVector()
{
  calculateRotationVector();
}

VSOMirror::~VSOMirror()
{
  // nothing to see here
}

void VSOMirror::calculateRotationVector()
{
  Vec3D rrot(0,-fTelescope->reflectorRotation(),0);
  Vec3D align(fAlign);
  align.Rotate(rrot);
  Vec3D normalize = align^Vec3D(0,1,0);
  double sintheta = normalize.Norm();
  if(sintheta == 0)
    {
      fRotationVector = rrot;
    }
  else 
    {
      double costheta = align*Vec3D(0,1,0);
      normalize *= atan(sintheta/costheta)/sintheta;
      fRotationVector = rrot & normalize;
    }
}

void VSOMirror::reflectorToMirror(Physics::Particle& p) const
{
  //  std::cerr << fRotationVector << std::endl;
  // First: Translate from reflector to mirror mount point
  p.TranslateOrigin(Vec4D(0,fPos));
  // Second: Rotate coordinate system so mirror normal is y-axis
  p.Rotate(fRotationVector);
  // Third: Fix parity
  // if(fTelescope->mirrorParity())mom3.x=-mom3.x;
  //  std::cout << "A: " << p.Momenta().r << ' ' << fRotationVector << ' ';
  //  std::cout << p.Momenta().r << std::endl;
#warning clean me up
}

void VSOMirror::mirrorToReflector(Physics::Particle& p) const
{
#warning Fix parity
  // First: Fix parity
  // if(fTelescope->mirrorParity())v.x=-v.x;
  // Second: Rotate coordinate system back to the reflector
  p.Rotate(-fRotationVector);
  // Third: Translate from mirror mount point to the reflector
  p.TranslateOrigin(Vec4D(0,-fPos));
}

VSDBStatement* VSOMirror::
createInsertQuery(VSDatabase* db)
{
  return db->createInsertQuery(VSIMDB_TABLE_NAME_MIRROR,14);
}

VSDBStatement* VSOMirror::
createSelectQuery(VSDatabase* db)
{
  return db->createSelectQuery(VSIMDB_TABLE_NAME_MIRROR,
			       "OpticsID=? AND TelescopeID=?");
}

void VSOMirror::writeToDatabase(VSDBStatement* stmt, uint32_t optics_id) const
{
  stmt->bindToParam(optics_id);

  stmt->bindToParam(fTelescope->id());
  stmt->bindToParam(fID);
  stmt->bindToParam(fHexID);
  stmt->bindToParam(fRemoved);
  stmt->bindToParam(fPos.x);

  stmt->bindToParam(fPos.y);
  stmt->bindToParam(fPos.z);
  stmt->bindToParam(fAlign.x);
  stmt->bindToParam(fAlign.y);
  stmt->bindToParam(fAlign.z);

  stmt->bindToParam(fFocalLength);
  stmt->bindToParam(fSpotSize);
  stmt->bindToParam(fDegradingFactor);

  assert(stmt->execute() == 1);
}

VSOMirror* VSOMirror::createFromDatabaseRow(VSDBStatement* stmt, 
					    const VSOTelescope* T)
{
  VSOMirror* mirror = new VSOMirror;
  unsigned OID;
  unsigned TID;

  stmt->bindToResult(OID);

  stmt->bindToResult(TID);
  stmt->bindToResult(mirror->fID);
  stmt->bindToResult(mirror->fHexID);
  stmt->bindToResult(mirror->fRemoved);
  stmt->bindToResult(mirror->fPos.x);

  stmt->bindToResult(mirror->fPos.y);
  stmt->bindToResult(mirror->fPos.z);
  stmt->bindToResult(mirror->fAlign.x);
  stmt->bindToResult(mirror->fAlign.y);
  stmt->bindToResult(mirror->fAlign.z);

  stmt->bindToResult(mirror->fFocalLength);
  stmt->bindToResult(mirror->fSpotSize);
  stmt->bindToResult(mirror->fDegradingFactor);
  
  if(stmt->retrieveNextRow() == 1)
    {
      mirror->fTelescope=T;
      assert(TID == mirror->fTelescope->id());
      mirror->calculateRotationVector();
    }
  else
    {
      delete mirror;
      mirror = 0;
    }

  return mirror;
}

void VSOMirror::createMirrorTable(VSDatabase* db)
{
  unsigned id;    // stupid to create this --- C++ should have a "typeof" 
  VSOMirror temp; // operator (g++ does but it is not portable .. and doesn't 
                  // work in a static function anyway

  db->createTable(VSIMDB_TABLE_NAME_MIRROR,

		  db->sqlSpecOf("OpticsID",id,true,"NOT NULL")+

		  db->sqlSpecOf("TelescopeID",id,false,"NOT NULL")+
		  db->sqlSpecOf("MirrorID",temp.fID,false,"NOT NULL")+
		  db->sqlSpecOf("MirrorHexID",temp.fHexID,false,"NOT NULL")+
		  db->sqlSpecOf("Removed",temp.fRemoved,false,"NOT NULL")+
		  db->sqlSpecOf("PosX",temp.fPos.x,false,"NOT NULL")+

		  db->sqlSpecOf("PosY",temp.fPos.y,false,"NOT NULL")+
		  db->sqlSpecOf("PosZ",temp.fPos.z,false,"NOT NULL")+
		  db->sqlSpecOf("AlignmentX",temp.fAlign.x,false,"NOT NULL")+
		  db->sqlSpecOf("AlignmentY",temp.fAlign.y,false,"NOT NULL")+
		  db->sqlSpecOf("AlignmentZ",temp.fAlign.z,false,"NOT NULL")+

		  db->sqlSpecOf("FocalLength",temp.fFocalLength,false,"NOT NULL")+
		  db->sqlSpecOf("SpotSize",temp.fSpotSize,false,"NOT NULL")+
		  db->sqlSpecOf("DegradingFactor",temp.fDegradingFactor,false,
				"NOT NULL")+
		  
		  ", PRIMARY KEY ( OpticsID, TelescopeID, MirrorID )",

		  VSDatabase::FLAG_NO_ERROR_ON_EXIST_OR_NOT_EXIST);
}

void VSOMirror::dumpShort(std::ostream& stream) const
{
  stream
    << "MIRROR "
    << VSDataConverter::toString(fTelescope->id()) << ' '
    << VSDataConverter::toString(fID) << ' '
    << VSDataConverter::toString(fHexID) << ' '
    << VSDataConverter::toString(fRemoved) << ' '
    << VSDataConverter::toString(fPos.x) << ' '

    << VSDataConverter::toString(fPos.y) << ' '
    << VSDataConverter::toString(fPos.z) << ' '
    << VSDataConverter::toString(fAlign.x) << ' '
    << VSDataConverter::toString(fAlign.y) << ' '
    << VSDataConverter::toString(fAlign.z) << ' '

    << VSDataConverter::toString(fFocalLength) << ' '
    << VSDataConverter::toString(fSpotSize) << ' '
    << VSDataConverter::toString(fDegradingFactor) << std::endl;
}

void VSOMirror::dump(std::ostream& stream, unsigned l) const
{
  stream
    << FDAV("Telescope ID", fTelescope->id(), "", 30, l) << std::endl
    << FDAV("Mirror ID", fID, "", 30, l) << std::endl
    << FDAV("Mirror Hex ID", fHexID, "", 30, l) << std::endl
    << FDAV("Removed", fRemoved, "", 30, l) << std::endl
    << FDAV("Position X", fPos.x, "", 30, l) << std::endl

    << FDAV("Position Y", fPos.y, "", 30, l) << std::endl
    << FDAV("Position Z", fPos.z, "", 30, l) << std::endl
    << FDAV("AlignmentX", fAlign.x, "", 30, l) << std::endl
    << FDAV("AlignmentY", fAlign.y, "", 30, l) << std::endl
    << FDAV("AlignmentZ", fAlign.z, "", 30, l) << std::endl

    << FDAV("FocalLength", fFocalLength, "", 30, l) << std::endl
    << FDAV("SpotSize", fSpotSize, "", 30, l) << std::endl
    << FDAV("DegradingFactor", fDegradingFactor, "", 30, l) << std::endl;
}

VSOMirror* VSOMirror::createFromShortDump(std::istream& stream,
					  const VSOTelescope* T)
{
  std::string line;
  std::getline(stream,line);
  if(line.empty())return 0;

  std::istringstream linestream(line);

  VSOMirror* mirror = new VSOMirror;
  mirror->fTelescope=T;

  std::string keyword;
  linestream >> keyword;
  assert(keyword==std::string("MIRROR"));

  unsigned TID;

  linestream
    >> TID
    >> mirror->fID
    >> mirror->fHexID
    >> mirror->fRemoved
    >> mirror->fPos.x

    >> mirror->fPos.y
    >> mirror->fPos.z
    >> mirror->fAlign.x
    >> mirror->fAlign.y
    >> mirror->fAlign.z

    >> mirror->fFocalLength
    >> mirror->fSpotSize
    >> mirror->fDegradingFactor;

  if(!linestream)
    {
      delete mirror;
      return 0;
    }

  mirror->calculateRotationVector();

  return mirror;
}
