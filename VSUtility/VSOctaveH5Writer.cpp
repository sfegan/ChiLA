//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOctaveH5Writer.cpp

  Class to write to GNU/Octave file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       05/20/2006
*/

#include<iomanip>
#include<cctype>
#include<vsassert>
#include<stdint.h>
#include<string>
#include<sstream>

#include"VSOctaveH5Writer.hpp"

using namespace VERITAS;

// ============================================================================
// WRITER BASE
// ============================================================================

bool VSOctaveH5WriterBase::s_default_compress = false;

VSOctaveH5WriterBase::~VSOctaveH5WriterBase()
{
  if(m_parent)m_parent->notifyOfChildDestruction(this);
  for(std::set<VSOctaveH5WriterBase*>::iterator iparty = 
	m_interested_parties.begin(); iparty != m_interested_parties.end();
      iparty++)(*iparty)->notifyOfChildDestruction(this);
}

void VSOctaveH5WriterBase::makeEmptyMatrix(hid_t gid, 
					   unsigned rows, unsigned cols)
{
  makeAttribute<uint8_t>(gid, OCTAVE_H5_ATTR_EMPTY_MATRIX, 1);

  int double_zero[2] = { cols, rows };
  hsize_t dim = 2;
  hid_t sid = H5Screate_simple(1, &dim, NULL);
  hid_t did = H5Dcreate(gid, OCTAVE_H5_EL_VALUE, H5T_STD_I32LE, 
			sid, H5P_DEFAULT);
  if(did<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not create dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
    
  if(H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	      double_zero)<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not write to dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
  
  H5Dclose(did);
  H5Sclose(sid);
}

void VSOctaveH5WriterBase::registerInterestedParty(VSOctaveH5WriterBase* x)
{
  m_interested_parties.insert(x);
}

void VSOctaveH5WriterBase::
notifyOfChildDestruction(VSOctaveH5WriterBase* child)
{
  m_children.erase(child);
}

void VSOctaveH5WriterBase::flush()
{
  // nothing to see here
}

// ============================================================================
// WRITER COMPOSITE VECTOR ELEMENT
// ============================================================================

VSOH5WCVElementWriter::~VSOH5WCVElementWriter()
{
  // nothing to see here
}

VSOH5WCVEWString::~VSOH5WCVEWString()
{
  // nothing to see here
}

bool VSOH5WCVEWString::
write(const void* data, const VSOctaveH5CompositeDefinition::Member& field)
{
  return 
    m_ecv->appendString(*((const std::string*)((const char*)(data) + field.offset))); 
}

// ============================================================================
// WRITER CELL ARRAY
// ============================================================================

VSOctaveH5WriterCellArray::
VSOctaveH5WriterCellArray(hid_t gid, unsigned rows, unsigned cols,
			  VSOctaveH5WriterBase* parent):
  VSOctaveH5WriterBase(parent), m_struct(), m_rows(rows), m_cols(cols),
  m_name_precision(0), m_written(rows*cols,false)
{
  if((cols!=0)&&(rows!=0))
    {
      m_struct = new VSOctaveH5WriterStruct(gid, 0);
      m_name_precision = VSOctaveH5Common::cellPrecision(rows*cols);
    }
}

VSOctaveH5WriterCellArray::~VSOctaveH5WriterCellArray()
{
  for(unsigned irow=0;irow<m_rows;irow++)
    for(unsigned icol=0;icol<m_cols;icol++)
      if(m_written[m_cols*irow+icol]==false)
	writeEmptyMatrix(irow, icol);
  if(m_struct)delete m_struct;
}

bool VSOctaveH5WriterCellArray::writeEmptyMatrix(unsigned row, unsigned col)
{
  return m_struct->writeEmptyMatrix(name(row,col));
}

bool VSOctaveH5WriterCellArray::
writeString(unsigned row, unsigned col, const std::string& s)
{
  return m_struct->writeString(name(row,col), s);
}

VSOctaveH5WriterStruct* VSOctaveH5WriterCellArray::
writeStruct(unsigned row, unsigned col)
{
  return m_struct->writeStruct(name(row,col));
}

VSOctaveH5WriterCellArray* VSOctaveH5WriterCellArray::
writeCellArray(unsigned row, unsigned col, unsigned rows, unsigned cols)
{
  return m_struct->writeCellArray(name(row,col), rows, cols);
}

VSOctaveH5WriterCellVector* VSOctaveH5WriterCellArray::
writeCellVector(unsigned row, unsigned col, unsigned nel,
		VSOctaveH5VecOrientation orient)
{
  return m_struct->writeCellVector(name(row,col), nel, orient);
}

VSOctaveH5WriterExpandableCellVector* VSOctaveH5WriterCellArray::
writeExpandableCellVector(unsigned row, unsigned col,
			  VSOctaveH5VecOrientation orient)
{
  return m_struct->writeExpandableCellVector(name(row,col), orient);
}

std::string VSOctaveH5WriterCellArray::name(unsigned row, unsigned col)
{
  if((row>=m_rows)||(col>=m_cols))throw std::out_of_range(__PRETTY_FUNCTION__);
  unsigned index = m_cols*row+col;
  m_written[index]=true;
  std::ostringstream stream;
  stream << '_' << std::setfill('0') << std::setw(m_name_precision) << index;
  return stream.str();
}

void VSOctaveH5WriterCellArray::flush()
{
  if(m_struct)m_struct->flush();
}

// ============================================================================
// WRITER CELL VECTOR
// ============================================================================

VSOctaveH5WriterCellVector::
VSOctaveH5WriterCellVector(hid_t gid, unsigned nel,
			   VSOctaveH5WriterBase* parent):
  VSOctaveH5WriterBase(parent), m_struct(), m_nel(nel), 
  m_name_precision(0), m_written(nel,false)
{
  if(nel!=0)
    {
      m_struct = new VSOctaveH5WriterStruct(gid, 0);
      m_name_precision = VSOctaveH5Common::cellPrecision(nel);
    }
}

VSOctaveH5WriterCellVector::~VSOctaveH5WriterCellVector()
{
  for(unsigned iel=0;iel<m_nel;iel++)
    if(m_written[iel]==false)
      writeEmptyMatrix(iel);
  if(m_struct)delete m_struct;
}

bool VSOctaveH5WriterCellVector::writeEmptyMatrix(unsigned iel)
{
  return m_struct->writeEmptyMatrix(name(iel));
}

bool VSOctaveH5WriterCellVector::
writeString(unsigned iel, const std::string& s)
{
  return m_struct->writeString(name(iel), s);
}

VSOctaveH5WriterStruct* VSOctaveH5WriterCellVector::
writeStruct(unsigned iel)
{
  return m_struct->writeStruct(name(iel));
}

VSOctaveH5WriterCellArray* VSOctaveH5WriterCellVector::
writeCellArray(unsigned iel, unsigned rows, unsigned cols)
{
  return m_struct->writeCellArray(name(iel), rows, cols);
}

VSOctaveH5WriterCellVector* VSOctaveH5WriterCellVector::
writeCellVector(unsigned iel, unsigned nel,
		VSOctaveH5VecOrientation orient)
{
  return m_struct->writeCellVector(name(iel), nel, orient);
}

VSOctaveH5WriterExpandableCellVector* VSOctaveH5WriterCellVector::
writeExpandableCellVector(unsigned iel, VSOctaveH5VecOrientation orient)
{
  return m_struct->writeExpandableCellVector(name(iel), orient);
}

std::string VSOctaveH5WriterCellVector::name(unsigned iel)
{
  if(iel>=m_nel)throw std::out_of_range(__PRETTY_FUNCTION__);
  m_written[iel]=true;
  std::ostringstream stream;
  stream << '_' << std::setfill('0') << std::setw(m_name_precision) << iel;
  return stream.str();
}

void VSOctaveH5WriterCellVector::flush()
{
  if(m_struct)m_struct->flush();
}

// ============================================================================
// WRITER EXPANGABLE CELL VECTOR
// ============================================================================

VSOctaveH5WriterExpandableCellVector::
VSOctaveH5WriterExpandableCellVector(hid_t gid, VSOctaveH5WriterBase* parent,
				     VSOctaveH5VecOrientation orient):
  VSOctaveH5WriterBase(parent), 
  m_nel(0), m_orient(orient), m_gid(gid), m_value_gid(-1), m_struct()
{
  // nothing to see here
}

VSOctaveH5WriterExpandableCellVector::~VSOctaveH5WriterExpandableCellVector()
{
  if(m_struct)
    {
      // Write the dimensions dataset
      int rows = (m_orient == VO_COLUMN)?m_nel:1;
      int cols = (m_orient == VO_COLUMN)?1:m_nel;
      int double_zero[2] = { rows, cols }; // inconsistant col,row order!!
      hsize_t dim = 2;
      hid_t sid = H5Screate_simple(1, &dim, NULL);
      hid_t pid = H5P_DEFAULT;
#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
      pid = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_layout(pid, H5D_COMPACT);
#endif
      hid_t did = H5Dcreate(m_value_gid, OCTAVE_H5_EL_DIMS, 
			    H5T_STD_I32LE, sid, pid);
      if(did<0)
	{
	  VSOctaveH5Exception 
	    e(std::string("Could not create dataset: " OCTAVE_H5_EL_VALUE));
	  throw(e);
	}
    
      if(H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		  double_zero)<0)
	{
	  VSOctaveH5Exception 
	    e(std::string("Could not write to dataset: " OCTAVE_H5_EL_VALUE));
	  throw(e);
	}
  
#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
      H5Pclose(pid);
#endif
      H5Dclose(did);
      H5Sclose(sid);

      // Clean up struct
      delete m_struct;

      // Rename variables
      unsigned name_precision = VSOctaveH5Common::cellPrecision(m_nel);

      m_value_gid = H5Gopen(m_gid, OCTAVE_H5_EL_VALUE);
      if(m_value_gid<0)
	{
	  VSOctaveH5Exception 
	    e(std::string("Could not re-open group: " OCTAVE_H5_EL_VALUE));
	  throw(e);
	}

      for(unsigned iel=0;iel<m_nel;iel++)
	{
	  std::ostringstream fr_stream;
	  fr_stream << '_' << iel;

	  std::ostringstream to_stream;
	  to_stream << '_' << std::setfill('0') 
		    << std::setw(name_precision) << iel;

	  std::string fr_name = fr_stream.str();
	  std::string to_name = to_stream.str();
	  if(fr_name != to_name)
	    H5Gmove(m_value_gid, 
		    fr_stream.str().c_str(), to_stream.str().c_str());
	}

      H5Gclose(m_value_gid);
    }
  else
    {
      if(m_orient == VO_COLUMN)makeEmptyMatrix(m_gid,0,1);
      else makeEmptyMatrix(m_gid,1,0);
    }
}

bool VSOctaveH5WriterExpandableCellVector::appendEmptyMatrix()
{
  return m_struct->writeEmptyMatrix(name());
}

bool VSOctaveH5WriterExpandableCellVector::appendString(const std::string& s)
{
  return m_struct->writeString(name(), s);
}

VSOctaveH5WriterStruct* VSOctaveH5WriterExpandableCellVector::appendStruct()
{
  return m_struct->writeStruct(name());
}

VSOctaveH5WriterCellArray* VSOctaveH5WriterExpandableCellVector::
appendCellArray(unsigned rows, unsigned cols)
{
  return m_struct->writeCellArray(name(), rows, cols);
}

VSOctaveH5WriterCellVector* VSOctaveH5WriterExpandableCellVector::
appendCellVector(unsigned nel, VSOctaveH5VecOrientation orient)
{
  return m_struct->writeCellVector(name(), nel, orient);
}

VSOctaveH5WriterExpandableCellVector* VSOctaveH5WriterExpandableCellVector::
appendExpandableCellVector(VSOctaveH5VecOrientation orient)
{
  return m_struct->writeExpandableCellVector(name(), orient);
}

std::string VSOctaveH5WriterExpandableCellVector::name()
{
  if(m_nel==0)
    {
      // Create the value group
      m_value_gid = H5Gcreate(m_gid, OCTAVE_H5_EL_VALUE, 0);
      if(m_value_gid<0)
	{
	  VSOctaveH5Exception 
	    e(std::string("Could not create group: " OCTAVE_H5_EL_VALUE));
	  throw(e);
	}

      m_struct = new VSOctaveH5WriterStruct(m_value_gid, 0);
    }

  std::ostringstream stream;
  stream << '_' << m_nel++;
  return stream.str();
}

void VSOctaveH5WriterExpandableCellVector::flush()
{
  if(m_struct)m_struct->flush();
}

// ============================================================================
// WRITER STRUCT
// ============================================================================

VSOctaveH5WriterStruct::
VSOctaveH5WriterStruct(const std::string& filename, bool overwrite,
		       const std::string& comment):
  VSOctaveH5WriterBase(0),
  m_my_file(true), m_gid(0), m_vars(), m_gid_map()
{
  m_gid = H5Fcreate(filename.c_str(), overwrite?H5F_ACC_TRUNC:H5F_ACC_EXCL,
		    H5P_DEFAULT, H5P_DEFAULT);
  if(m_gid < 0)
    {
      VSOctaveH5Exception e(std::string("Could not open file ")+filename);
      throw e;
    }

  if(!comment.empty())H5Gset_comment(m_gid, "/", comment.c_str());
}

VSOctaveH5WriterStruct::
VSOctaveH5WriterStruct(hid_t gid, VSOctaveH5WriterBase* parent):
  VSOctaveH5WriterBase(parent),
  m_my_file(false), m_gid(gid), m_vars(), m_gid_map()
{
  // nothing to see here
}

VSOctaveH5WriterStruct::~VSOctaveH5WriterStruct()
{
  unsigned nchildren = m_children.size();
  while(!m_children.empty())
    {
      delete *m_children.begin();
      unsigned mchildren = m_children.size();
      vsassert(mchildren<nchildren);
      nchildren = mchildren;
    }
//  std::set<VSOctaveH5WriterBase*> delete_us = m_children;
//   for(std::set<VSOctaveH5WriterBase*>::iterator idel = delete_us.begin();
//       idel != delete_us.end(); idel++)delete *idel;
//   vsassert(m_children.empty());
  vsassert(m_gid_map.empty());
  if(m_my_file)
    {
      H5Fflush(m_gid,H5F_SCOPE_GLOBAL);
      H5Fclose(m_gid);
    }
  else H5Gclose(m_gid);
}

void VSOctaveH5WriterStruct::
registerChild(VSOctaveH5WriterBase* child, hid_t gid)
{
  m_gid_map[child]=gid;
  m_children.insert(child);
}

void VSOctaveH5WriterStruct::
notifyOfChildDestruction(VSOctaveH5WriterBase* child)
{
  std::map<VSOctaveH5WriterBase*, hid_t>::iterator
    ifind = m_gid_map.find(child);
  vsassert(ifind != m_gid_map.end());
  if(ifind->second >= 0)H5Gclose(ifind->second);
  m_gid_map.erase(ifind);
  m_children.erase(child);
}

hid_t VSOctaveH5WriterStruct::
makeVariable(const std::string& name, 
	     const std::string& type, const std::string& eltype)
{
  std::string _name = name;
  //if(isdigit(name[0]))
  
  // Check that variable doesn't exist already
  
  std::set<std::string>::iterator ivar = m_vars.find(_name);
  if(ivar != m_vars.end())
    {
      VSOctaveH5Exception e(std::string("Variable \"")+_name
			    +std::string("\" already exists"));
      throw(e);
    }

  // Open the group

  hid_t gid = H5Gcreate(m_gid, name.c_str(), 0);
  if(gid<0)
    {
      VSOctaveH5Exception e(std::string("Could not create group: ")+_name);
      throw(e);
    }

  hid_t tid;
  hid_t sid;
  hid_t did;
  hid_t pid = H5P_DEFAULT;

  // Write the OCTAVE_NEW_FORMAT attribute

  makeAttribute<uint8_t>(gid, OCTAVE_H5_ATTR_NEW_FORMAT, 1);

  // Write the TYPE dataset

  std::string full_type;
  if(!eltype.empty())full_type = eltype + std::string(" ");
  full_type += type;

  if(full_type == "bool scalar")full_type="bool"; // Bool exception from octave

  tid = H5Tcopy (H5T_C_S1);
  H5Tset_size(tid, full_type.size()+1);
  sid = H5Screate(H5S_SCALAR);

#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
  pid = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_layout(pid, H5D_COMPACT);
#endif

  did = H5Dcreate(gid, OCTAVE_H5_EL_TYPE, tid, sid, pid);

  if(did<0)
    {
      VSOctaveH5Exception 
	e(name+std::string(": could not create dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
  

  if(H5Dwrite(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, full_type.c_str()) < 0)
    {
      VSOctaveH5Exception 
	e(name
	  + std::string(": could not write to dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }

#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
  H5Pclose(pid);
#endif
  H5Dclose(did);
  H5Sclose(sid);
  H5Tclose(tid);

  m_vars.insert(name);
  return gid;
}

bool VSOctaveH5WriterStruct::writeEmptyMatrix(const std::string& name)
{
  return writeMatrix(name, 0, 0, (double*)0);
}

bool VSOctaveH5WriterStruct::
writeString(const std::string& name, const std::string& s)
{
  hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_STRING);
  
  if(s.size() == 0)
    {
      makeEmptyMatrix(gid,1,0);
      return true;
    }

  hsize_t dims[2] = { s.size(), 1 };
  hid_t sid = H5Screate_simple(2, dims, NULL);
  hid_t did = H5Dcreate(gid, OCTAVE_H5_EL_VALUE, 
			H5T_STD_I8LE, sid, H5P_DEFAULT);
  if(did<0)
    {
      VSOctaveH5Exception 
	e(name+std::string(": could not create dataset: type"));
      throw(e);
    }
  
  if(H5Dwrite(did, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	      s.c_str()) < 0)
    {
      VSOctaveH5Exception 
	e(name+std::string(": could not write to dataset: type"));
      throw(e);
    }

  H5Dclose(did);
  H5Sclose(sid);
  H5Gclose(gid);

  return true;
}

VSOctaveH5WriterStruct* 
VSOctaveH5WriterStruct::writeStruct(const std::string& name)
{
  hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_STRUCT);
  
  hid_t value_gid = H5Gcreate(gid, OCTAVE_H5_EL_VALUE, 0);
  if(value_gid<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not create group: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
  
  VSOctaveH5WriterStruct* s = new VSOctaveH5WriterStruct(value_gid,this);
  registerChild(s,gid);

  return s;
}

VSOctaveH5WriterCellArray* 
VSOctaveH5WriterStruct::writeCellArray(const std::string& name, 
				       unsigned rows, unsigned cols)
{
  hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_CELL);

  // If its empty then 

  if((rows==0)||(cols==0))
    {
      makeEmptyMatrix(gid, rows, cols);
      VSOctaveH5WriterCellArray* c = 
	new VSOctaveH5WriterCellArray(-1,rows,cols,this);
      registerChild(c,gid);
      return c;
    }

  // Create the value group

  hid_t value_gid = H5Gcreate(gid, OCTAVE_H5_EL_VALUE, 0);
  if(value_gid<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not create group: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }

  // Write the dimensions dataset

  int double_zero[2] = { rows, cols }; // inconsistant col,row order!!
  hsize_t dim = 2;
  hid_t sid = H5Screate_simple(1, &dim, NULL);
      hid_t pid = H5P_DEFAULT;
#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
      pid = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_layout(pid, H5D_COMPACT);
#endif
  hid_t did = H5Dcreate(value_gid, OCTAVE_H5_EL_DIMS, 
			H5T_STD_I32LE, sid, pid);
  if(did<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not create dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
    
  if(H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	      double_zero)<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not write to dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
  
#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
      H5Pclose(pid);
#endif
  H5Dclose(did);
  H5Sclose(sid);

  // Create the cell object
  
  VSOctaveH5WriterCellArray* c = 
    new VSOctaveH5WriterCellArray(value_gid,rows,cols,this);
  registerChild(c,gid);

  return c;
}

VSOctaveH5WriterCellVector* VSOctaveH5WriterStruct::
writeCellVector(const std::string& name, unsigned nel,
		VSOctaveH5VecOrientation orient)
{
  hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_CELL);

  // If its empty then 

  if(nel==0)
    {
      if(orient == VO_COLUMN)makeEmptyMatrix(gid,nel,1);
      else makeEmptyMatrix(gid,1,nel);
      VSOctaveH5WriterCellVector* c = 
	new VSOctaveH5WriterCellVector(-1,nel,this);
      registerChild(c,gid);
      return c;
    }

  // Create the value group

  hid_t value_gid = H5Gcreate(gid, OCTAVE_H5_EL_VALUE, 0);
  if(value_gid<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not create group: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }

  // Write the dimensions dataset

  int rows = (orient == VO_COLUMN)?nel:1;
  int cols = (orient == VO_COLUMN)?1:nel;
  int double_zero[2] = { rows, cols }; // inconsistant col,row order!!
  hsize_t dim = 2;
  hid_t sid = H5Screate_simple(1, &dim, NULL);
      hid_t pid = H5P_DEFAULT;
#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
      pid = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_layout(pid, H5D_COMPACT);
#endif
  hid_t did = H5Dcreate(value_gid, OCTAVE_H5_EL_DIMS, 
			H5T_STD_I32LE, sid, pid);
  if(did<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not create dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
    
  if(H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
	      double_zero)<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not write to dataset: " OCTAVE_H5_EL_VALUE));
      throw(e);
    }
  
#if (H5_VERS_MAJOR>1)||(H5_VERS_MAJOR==1 && H5_VERS_MINOR>=6)
      H5Pclose(pid);
#endif
  H5Dclose(did);
  H5Sclose(sid);

  // Create the cell object
  
  VSOctaveH5WriterCellVector* c = 
    new VSOctaveH5WriterCellVector(value_gid,nel,this);
  registerChild(c,gid);

  return c;
}

VSOctaveH5WriterExpandableCellVector* VSOctaveH5WriterStruct::
writeExpandableCellVector(const std::string& name,
			  VSOctaveH5VecOrientation orient)
{
  hid_t gid = makeVariable(name, OCTAVE_H5_TYPE_CELL);

  // Create the cell object
  
  VSOctaveH5WriterExpandableCellVector* c = 
    new VSOctaveH5WriterExpandableCellVector(gid,this,orient);
  registerChild(c,gid);

  return c;
}

void VSOctaveH5WriterStruct::flush()
{
  H5Fflush(m_gid,H5F_SCOPE_LOCAL);
  for(std::set<VSOctaveH5WriterBase*>::iterator ichild = m_children.begin();
      ichild != m_children.end(); ichild++)(*ichild)->flush();
  H5Fflush(m_gid,H5F_SCOPE_LOCAL);
}

// ****************************************************************************
// TEST MAIN FUNCTION
// ****************************************************************************

#ifdef TEST_MAIN_H5WRITER

// g++ -DTEST_MAIN_H5WRITER -I. -I../VSCommon -L../VSCommon -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -g -Wall -DTEST_MAIN -o test VATime.o VSOctaveH5Writer.cpp VSOctaveIO.cpp -lVSCommon -lhdf5

#include<iostream>
#include<VSTime.hpp>

struct test_class
{
  unsigned i;
  float f;
  double d;

  static void _compose(VSOctaveH5CompositeDefinition& c)
  {
    H5_ADDMEMBER(c,test_class,i);
    H5_ADDMEMBER(c,test_class,f);
    H5_ADDMEMBER(c,test_class,d);
  }
};

struct test_class_with_string
{
  test_class tc;
  std::string s;

  static void _compose(VSOctaveH5CompositeDefinition& c)
  {
    H5_ADDSIMPLECOMPOSITE(c,test_class_with_string,tc);
    H5_ADDMEMBER(c,test_class_with_string,s);
  }
};

int main(int argc, char** argv)
{
  try
    {
      VSOctaveH5Writer writer("test.h5", true, "# produced by test code");
      writer.writeString("a","Hello there!!");
      double x = 10.0;
      double y = 11.0;
      writer.writeScalar("x",x);
      writer.writeScalar("y",y);
      std::vector<double> v(10);
      for(unsigned i=0;i<v.size();i++)v[i]=double(i)*double(i);
      writer.writeVector("r",v,VO_ROW);
      writer.writeVector("c",v,VO_COLUMN);

      uint32_t ix = 1010;
      uint32_t iy = 1011;
      writer.writeScalar("ix",ix);
      writer.writeScalar("iy",iy);
      std::vector<uint32_t> iv(10);
      for(unsigned i=0;i<iv.size();i++)iv[i]=i*i;
      writer.writeVector("ir",iv,VO_ROW);
      writer.writeVector("ic",iv,VO_COLUMN);
      double m[3][4] = { { 11,12,13,14 } , { 21,22,23,24} ,{ 31,32,33,34 } };
      writer.writeMatrix("m",3,4,(double*)m);

      VSOctaveH5WriterVector<double>* vv = 
	writer.writeExpandableVector<double>("vv");
      vv->append(100.0);

      VSOctaveH5WriterStruct *s = writer.writeStruct("s");
      s->writeVector("c",v,VO_COLUMN);
      s->writeVector("r",v,VO_ROW);

      VSOctaveH5WriterStruct *t = writer.writeStruct("t");
      t->writeScalar("x",x);
      t->writeScalar("y",y);
      
      vv->append(200.0);

      VSOctaveH5WriterTable<double>* mm = 
	writer.writeExpandableTable<double>("table",3);
      std::vector<double> row(3);
      row[0] = 1; row[1] = 2; row[2] = 3; 
      mm->append(row);
      row[0] = 10; row[1] = 11; row[2] = 12; 
      mm->append(row);

      writer.writeScalar("u",10U);

      writer.writeString("empty_s","");
      writer.writeExpandableVector<double>("empty_ev");
      writer.writeExpandableTable<double>("empty_t",3);
      writer.writeVector("empty_rv",std::vector<double>(),VO_ROW);
      writer.writeVector("empty_cv",std::vector<double>(),VO_COLUMN);
      writer.writeMatrix("empty_m",10,0,(double*)NULL);

      VSOctaveH5WriterCellArray* array = writer.writeCellArray("ca",3,2);
      array->writeScalar(0,0,(double)11.0);
      array->writeScalar(0,1,(double)12.0);
      array->writeScalar(1,0,(double)21.0);
      array->writeScalar(1,1,(double)22.0);
      array->writeScalar(2,0,(double)31.0);
      VSOctaveH5WriterStruct *u = array->writeStruct(2,1);

      VSOctaveH5WriterCellVector* cvector = writer.writeCellVector("cv",7);
      cvector->writeScalar(0,(double)70.0);
      cvector->writeScalar(1,(double)71.0);
      cvector->writeScalar(2,(double)72.0);
      cvector->writeScalar(3,(double)73.0);
      cvector->writeScalar(4,(double)74.0);
      cvector->writeScalar(5,(double)75.0);
      cvector->writeScalar(6,(double)76.0);
      
      u->writeScalar("x",x);
      u->writeScalar("y",y);
      u->writeVector("c",v,VO_COLUMN);
      u->writeVector("r",v,VO_ROW);      

      bool b = true;
      writer.writeScalar("b",b);
      bool mb[12];
      for(unsigned i=0;i<sizeof(mb)/sizeof(*mb);i++)mb[i]=false;
      mb[6]=true;
      writer.writeMatrix("bm",3,4,mb);
      std::vector<bool> bv;
      bv.push_back(false);
      bv.push_back(false);
      bv.push_back(false);
      bv.push_back(true);
      bv.push_back(true);
      bv.push_back(true);
      bv.push_back(false);
      bv.push_back(false);
      bv.push_back(false);
      writer.writeVector("bv",bv);

      test_class tc1 = { 1, 2.0, 0.5 };
      writer.writeComposite("tc1",tc1);

      std::vector<test_class> tc2;
      for(unsigned i=0;i<1000;i++)
	{
	  tc2.push_back(tc1);
	  tc1.i += 2;
	  tc1.f *= 1.001;
	  tc1.d = 3.7*tc1.d*(1.0-tc1.d); // CHAOS!!
	}
      writer.writeCompositeVector("tc2",tc2);

      VSOctaveH5WriterCompositeVector<test_class>* tc3 = 
	writer.writeCompositeExpandableVector<test_class>("tc3");
      
      for(unsigned i=0;i<2000;i++)
	{
	  tc3->append(tc1);
	  tc1.i += 2;
	  tc1.f *= 1.001;
	  tc1.d = 3.7*tc1.d*(1.0-tc1.d);
	}

      test_class_with_string tc4;
      tc4.tc.i = 1000;
      tc4.tc.f = 2000;
      tc4.tc.d = 3000;
      tc4.s = "4000";
      writer.writeComposite("tc4",tc4);

      std::vector<VSTime> vstime;
      for(unsigned i=0;i<2000;i++)vstime.push_back(VSTime::now());
      writer.writeCompositeVector("vstime",vstime);

      std::pair<double, int> pc1(1.234,567);
      writer.writeComposite("pc1",pc1);

      VSOctaveH5WriterExpandableCellVector* ecv = 
	writer.writeExpandableCellVector("ecv");
      for(unsigned i=0;i<100;i++)
	ecv->appendScalar(i);

      writer.writeCellArray("ca-r",1,10);
      writer.writeCellArray("ca-c",10,1);
      writer.writeCellVector("cv-r",10,VO_ROW);
      writer.writeCellVector("cv-c",10,VO_COLUMN);
      ecv = writer.writeExpandableCellVector("ecv-r",VO_ROW);
      for(unsigned i=0;i<10;i++)ecv->appendEmptyMatrix();
      ecv = writer.writeExpandableCellVector("ecv-c",VO_COLUMN);
      for(unsigned i=0;i<10;i++)ecv->appendEmptyMatrix();
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << e.message() << std::endl;
    }

  H5close();
}
#endif
