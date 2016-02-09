//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOctaveH5Reader.cpp

  Class to read to GNU/Octave file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/12/2006
*/

#include<iomanip>
#include<iostream>
#include<cctype>
#include<vsassert>
#include<stdint.h>
#include<string>
#include<sstream>

#include"VSFileUtility.hpp"
#include"VSOctaveH5Reader.hpp"

using namespace VERITAS;

// ============================================================================
//
// READER BASE
//
// ============================================================================

VSOctaveH5ReaderBase::~VSOctaveH5ReaderBase()
{
  std::set<VSOctaveH5ReaderBase*> delete_us = m_children;
  for(std::set<VSOctaveH5ReaderBase*>::iterator idel = delete_us.begin();
      idel != delete_us.end(); idel++)delete *idel;
  vsassert(m_children.empty());

  if(m_parent)m_parent->notifyOfChildDestruction(this);
}

bool VSOctaveH5ReaderBase::
getEmpty(Variable* var, unsigned& rows, unsigned& cols)
{
  uint8_t is_empty = 0;
  var->getAttribute(OCTAVE_H5_ATTR_EMPTY_MATRIX, is_empty);
  if(!is_empty)return false;

  hid_t did = H5Dopen(var->gid, OCTAVE_H5_EL_VALUE);
  if(did < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not open empty dataset: " OCTAVE_H5_EL_VALUE));
      throw e;
    }

  hid_t sid = H5Dget_space(did);
  if(sid < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not get space for empty dataset"));
      throw e;
    }

  if(H5Sget_simple_extent_type(sid) != H5S_SIMPLE)
    {
      VSOctaveH5Exception 
	e(std::string("Space for empty dataset is not SIMPLE"));
      throw e;	  
    }

  if(H5Sget_simple_extent_ndims(sid) != 1)
    {
      VSOctaveH5Exception 
	e(std::string("Empty dataset does not have 1 dimension"));
      throw e;
    }

  if(H5Sget_simple_extent_npoints(sid) != 2)
    {
      VSOctaveH5Exception 
	e(std::string("Empty dataset does not have 2 elements"));
      throw e;
    }

  unsigned dims[2];
  if(H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,dims) < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not read dimensions of empty dataset"));
      throw e;
    }
      
  cols = dims[0];
  rows = dims[1];

  H5Sclose(sid);
  H5Dclose(did);

  return true;
}

bool VSOctaveH5ReaderBase::
getCString(hid_t gid, const std::string& name, std::string& s)
{
  hid_t did = H5Dopen(gid,name.c_str());
  if(did < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not open dataset: ")+std::string(name));
      throw e;
    }

  hid_t tid = H5Dget_type(did);
  if(tid < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not get dataset type: ")+std::string(name));
      throw e;
    }
    
  if(H5Tget_class(tid) != H5T_STRING)
    {
      VSOctaveH5Exception 
	e(std::string("Dataset type does not have string class: ")
	  + std::string(name));
      throw e;
    }

  hid_t sid = H5Dget_space(did);
  if(sid < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not get dataspace: ")+std::string(name));
      throw e;
    }

  if(H5Sget_simple_extent_type(sid) != H5S_SCALAR)
    {
      VSOctaveH5Exception 
	e(std::string("String dataspace is not scalar: ")
	  + std::string(name));
      throw e;
    }
  
  unsigned nchar = H5Tget_size(tid);
  char* buffer = new char[nchar];

  if(H5Dread(did, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not read string dataspace: ")
	  + std::string(name));
      throw e;
    }
  
  s = buffer;

  delete[] buffer;
  H5Sclose(sid);
  H5Tclose(tid);
  H5Dclose(did);

  return true;
}

void VSOctaveH5ReaderBase::
registerChild(VSOctaveH5ReaderBase *reader)
{
  m_children.insert(reader);
}


void VSOctaveH5ReaderBase::
notifyOfChildDestruction(VSOctaveH5ReaderBase* child)
{
  m_children.erase(child);
}

// ============================================================================
// READER BASE / VARIABLE
// ============================================================================

herr_t VSOctaveH5ReaderBase::Variable::
iterateOverAtr(hid_t gid, const char*  name, void *object)
{
  static_cast<VSOctaveH5ReaderStruct::Variable*>(object)->
    registerAttribute(name);
  return 0;
}

void VSOctaveH5ReaderBase::Variable::registerAttribute(const char* name)
{
  attr.insert(name);
}

void VSOctaveH5ReaderBase::Variable::readAttributes()
{
  if(H5Aiterate(gid, NULL, &iterateOverAtr, this) != 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not iterate over members of group"));
      throw(e);
    }
}

// ==========================================================================
//
// READER COMPOSITE VECTOR
//
// ==========================================================================

VSOH5RCVElementReader::~VSOH5RCVElementReader()
{
  // nothing to see here
}

bool VSOH5RCVElementReader::isNullReader() const
{
  return false;
}

VSOH5RCVERString::~VSOH5RCVERString()
{
  delete m_cv;
}

bool VSOH5RCVERString::
read(void* data, unsigned index,
     const VSOctaveH5CompositeDefinition::Member& field)
{
  return
    m_cv->readString(index, *((std::string*)((char*)(data)+field.offset)));
}

void VSOH5RCVERString::releaseReader() 
{ 
  m_cv = 0;
}

// ============================================================================
//
// READER MATRIX
//
// ============================================================================

VSOctaveH5ReaderMatrix::VSOctaveH5ReaderMatrix(Variable* var):
  VSOctaveH5ReaderBase(var->parent), m_did(), m_cols(), m_rows()
{
  if(getEmpty(var, m_rows, m_cols))return;

  m_did = H5Dopen(var->gid, OCTAVE_H5_EL_VALUE);
  if(m_did < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not open dataset: " OCTAVE_H5_EL_VALUE));
      throw e;
    }

  hid_t sid = H5Dget_space(m_did);
  if(sid < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not get space for dataset"));
      throw e;
    }

  H5S_class_t type = H5Sget_simple_extent_type(sid);
  if(type == H5S_SCALAR)
    {
      m_rows = 1;
      m_cols = 1;
    }
  else if(type == H5S_SIMPLE)
    {
      hsize_t ndim = H5Sget_simple_extent_ndims(sid);
      if(ndim <= 0)
	{
	  VSOctaveH5Exception 
	    e(std::string("Could not get dimensions of dataset"));
	  throw e;
	}
      else if(ndim > 2)
	{
	  VSOctaveH5Exception 
	    e(std::string("Dataset has more than 2 dimensions"));
	  throw e;
	}

      hsize_t dims[2];
      if(H5Sget_simple_extent_dims(sid, dims, NULL) < 0)
	{
	  VSOctaveH5Exception 
	    e(std::string("Could not get extent of dataset"));
	  throw e;
	}
      
      if(ndim == 1)m_rows=1, m_cols=dims[0];
      else m_rows=dims[0], m_cols=dims[1];
    }
  else
    {
      VSOctaveH5Exception 
	e(std::string("Dataset is neither SCALER or SIMPLE"));
      throw e;
    }

  H5Sclose(sid);
}

VSOctaveH5ReaderMatrix::~VSOctaveH5ReaderMatrix()
{
  std::set<VSOctaveH5ReaderBase*> delete_us = m_children;
  for(std::set<VSOctaveH5ReaderBase*>::iterator idel = delete_us.begin();
      idel != delete_us.end(); idel++)delete *idel;
  vsassert(m_children.empty());
  if(m_did)H5Dclose(m_did);
}

bool VSOctaveH5ReaderMatrix::readString(std::string& s)
{
  if(!isVector())return false;
  if(isEmpty()){ s.clear(); return true; }
  if(m_rows!=1)s.resize(m_rows);
  else s.resize(m_cols);
  hid_t core_tid = VSOctaveH5Type<uint8_t>::h5_type();
  if(H5Dread(m_did,core_tid,H5S_ALL,H5S_ALL,H5P_DEFAULT,&(*s.begin())) < 0)
    {
      VSOctaveH5Exception e(std::string("Could not read from dataset"));
      throw(e);
    }
  H5Tclose(core_tid);
  return true;
}

// ============================================================================
//
// READER CELL ARRAY
//
// ============================================================================

VSOctaveH5ReaderCellArray::VSOctaveH5ReaderCellArray(Variable* var):
  VSOctaveH5ReaderBase(var->parent),
  m_struct(), m_rows(), m_cols(), m_name_precision()
{
  if(getEmpty(var, m_rows, m_cols))return;

  hid_t gid = H5Gopen(var->gid, OCTAVE_H5_EL_VALUE);
  if(gid < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not open group " OCTAVE_H5_EL_VALUE));
      throw e;
    }

  hid_t did = H5Dopen(gid, OCTAVE_H5_EL_DIMS);
  if(did < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not open cell dataset: " OCTAVE_H5_EL_VALUE));
      throw e;
    }

  hid_t sid = H5Dget_space(did);
  if(sid < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not get space for cell dataset"));
      throw e;
    }

  if(H5Sget_simple_extent_type(sid) != H5S_SIMPLE)
    {
      VSOctaveH5Exception 
	e(std::string("Space for cell dataset is not SIMPLE"));
      throw e;	  
    }

  if(H5Sget_simple_extent_ndims(sid) != 1)
    {
      VSOctaveH5Exception 
	e(std::string("Cell dataset does not have 1 dimension"));
      throw e;
    }

  if(H5Sget_simple_extent_npoints(sid) != 2)
    {
      VSOctaveH5Exception 
	e(std::string("Cell dataset does not have 2 elements"));
      throw e;
    }

  unsigned dims[2];
  if(H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,dims) < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not read dimensions of cell dataset"));
      throw e;
    }
      
  H5Sclose(sid);
  H5Dclose(did);  
  H5Gclose(gid);  

  m_rows = dims[0];
  m_cols = dims[1];

  if((m_rows!=0)&&(m_cols!=0))
    m_name_precision = VSOctaveH5Common::cellPrecision(m_rows*m_cols);
  else
    m_name_precision = 0;

  Variable bogus_var = *var;
  bogus_var.parent = 0; // struct does not inform anyone of its deletion

  m_struct = new VSOctaveH5ReaderStruct(&bogus_var);  
}

VSOctaveH5ReaderCellArray::~VSOctaveH5ReaderCellArray()
{
  if(m_struct)delete m_struct;
}

std::string VSOctaveH5ReaderCellArray::variableType(unsigned row, unsigned col)
{
  return m_struct->variableType(name(row,col));
}

std::string VSOctaveH5ReaderCellArray::elementType(unsigned row, unsigned col)
{
  return m_struct->elementType(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isValid(unsigned row, unsigned col)
{
  return (row<m_rows)&&(col<m_cols);
}

bool VSOctaveH5ReaderCellArray::isEmpty(unsigned row, unsigned col)
{
  return m_struct->isEmpty(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isScalar(unsigned row, unsigned col)
{
  return m_struct->isScalar(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isVector(unsigned row, unsigned col)
{
  return m_struct->isVector(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isMatrix(unsigned row, unsigned col)
{
  return m_struct->isMatrix(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isString(unsigned row, unsigned col)
{
  return m_struct->isString(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isCellArray(unsigned row, unsigned col)
{
  return m_struct->isCellArray(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isCellVector(unsigned row, unsigned col)
{
  return m_struct->isCellVector(name(row,col));
}

bool VSOctaveH5ReaderCellArray::isStruct(unsigned row, unsigned col)
{
  return m_struct->isStruct(name(row,col));
}

bool VSOctaveH5ReaderCellArray::dimensions(unsigned row, unsigned col,
					   unsigned& rows, unsigned& cols)
{
  return m_struct->dimensions(name(row,col), rows, cols);  
}

bool VSOctaveH5ReaderCellArray::
readString(unsigned row, unsigned col, std::string& s)
{
  return m_struct->readString(name(row,col), s);
}

VSOctaveH5ReaderStruct* VSOctaveH5ReaderCellArray::
readStruct(unsigned row, unsigned col)
{
  return m_struct->readStruct(name(row,col));
}

VSOctaveH5ReaderCellArray* VSOctaveH5ReaderCellArray::
readCellArray(unsigned row, unsigned col)
{
  return m_struct->readCellArray(name(row,col));
}

VSOctaveH5ReaderCellVector* VSOctaveH5ReaderCellArray::
readCellVector(unsigned row, unsigned col)
{
  return m_struct->readCellVector(name(row,col));
}

bool VSOctaveH5ReaderCellArray::
copyTo(unsigned srow, unsigned scol,
       VSOctaveH5WriterStruct* s, std::string dname)
{
  return m_struct->copyTo(name(srow,scol), s, dname);
}

bool VSOctaveH5ReaderCellArray::
copyTo(unsigned srow, unsigned scol,
       VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol)
{
  return m_struct->copyTo(name(srow,scol), c, dirow, dicol);
}

bool VSOctaveH5ReaderCellArray::
copyTo(unsigned srow, unsigned scol,
       VSOctaveH5WriterCellVector* c, unsigned diel)
{
  return m_struct->copyTo(name(srow,scol), c, diel);
}

bool VSOctaveH5ReaderCellArray::
copyTo(unsigned srow, unsigned scol,
       VSOctaveH5WriterExpandableCellVector* c)
{
  return m_struct->copyTo(name(srow,scol), c);
}

std::string VSOctaveH5ReaderCellArray::name(unsigned row, unsigned col)
{
  if(!isValid(row,col))throw std::out_of_range(__PRETTY_FUNCTION__);
  unsigned index = m_cols*row+col;
  std::ostringstream stream;
  stream << '_' << std::setfill('0') << std::setw(m_name_precision) << index;
  return stream.str();
}

// ============================================================================
//
// READER CELL VECTOR
//
// ============================================================================

VSOctaveH5ReaderCellVector::
VSOctaveH5ReaderCellVector(VSOctaveH5ReaderCellArray* cell):
  VSOctaveH5ReaderBase(cell), m_cell(cell), m_nel(), m_column(false)
{
  unsigned row;
  unsigned col;
  m_cell->dimensions(row,col);
  vsassert(row==1 || col==1);
  if((col==1)&&(row!=1))m_column=true,m_nel=row;
  else m_nel = col;
}

VSOctaveH5ReaderCellVector::~VSOctaveH5ReaderCellVector()
{
  // nothing to see here
}

std::string VSOctaveH5ReaderCellVector::variableType(unsigned iel)
{
  return m_cell->variableType(row(iel),col(iel));
}

std::string VSOctaveH5ReaderCellVector::elementType(unsigned iel)
{
  return m_cell->elementType(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isValid(unsigned iel)
{
  return iel<m_nel;
}

bool VSOctaveH5ReaderCellVector::isEmpty(unsigned iel)
{
  return m_cell->isEmpty(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isScalar(unsigned iel)
{
  return m_cell->isScalar(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isVector(unsigned iel)
{
  return m_cell->isVector(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isMatrix(unsigned iel)
{
  return m_cell->isMatrix(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isString(unsigned iel)
{
  return m_cell->isString(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isCellArray(unsigned iel)
{
  return m_cell->isCellArray(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isCellVector(unsigned iel)
{
  return m_cell->isCellVector(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::isStruct(unsigned iel)
{
  return m_cell->isStruct(row(iel),col(iel));
}

bool VSOctaveH5ReaderCellVector::dimensions(unsigned iel,
					    unsigned& rows, unsigned& cols)
{
  return m_cell->dimensions(row(iel),col(iel), rows, cols);  
}

bool VSOctaveH5ReaderCellVector::
readString(unsigned iel, std::string& s)
{
  return m_cell->readString(row(iel),col(iel), s);
}

VSOctaveH5ReaderStruct* VSOctaveH5ReaderCellVector::
readStruct(unsigned iel)
{
  return m_cell->readStruct(row(iel),col(iel));
}

VSOctaveH5ReaderCellArray* VSOctaveH5ReaderCellVector::
readCellArray(unsigned iel)
{
  return m_cell->readCellArray(row(iel),col(iel));
}

VSOctaveH5ReaderCellVector* VSOctaveH5ReaderCellVector::
readCellVector(unsigned iel)
{
  return m_cell->readCellVector(row(iel),col(iel));
}

// ============================================================================
//
// READER STRUCT
//
// ============================================================================

VSOctaveH5ReaderStruct::
VSOctaveH5ReaderStruct(const std::string& filename):
  VSOctaveH5ReaderBase(0), 
  m_my_file(true), m_gid(0), m_vars(), m_index(), m_assistants()
{
  if(!VSFileUtility::isFile(filename))
    {
      VSOctaveH5Exception e(std::string("Not a file: ")+filename);
      throw e;
    }

  if(!isHDF5(filename))
    {
      VSOctaveH5Exception e(std::string("File is not HDF5: ")+filename);
      throw e;
    }

  m_gid = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(m_gid < 0)
    {
      VSOctaveH5Exception e(std::string("Could not open file: ")+filename);
      throw e;
    }

  readGroup();
}

VSOctaveH5ReaderStruct::
VSOctaveH5ReaderStruct(Variable* var):
  VSOctaveH5ReaderBase(var->parent),
  m_my_file(false), m_gid(), m_vars(), m_index(), m_assistants()
{
  m_gid = H5Gopen(var->gid, OCTAVE_H5_EL_VALUE);
  if(m_gid < 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not open group: " OCTAVE_H5_EL_VALUE));
      throw e;
    }  
  readGroup();
}

VSOctaveH5ReaderStruct::~VSOctaveH5ReaderStruct()
{
  std::set<VSOctaveH5ReaderBase*> delete_us;

  delete_us = m_assistants;
  for(std::set<VSOctaveH5ReaderBase*>::iterator idel = delete_us.begin();
      idel != delete_us.end(); idel++)delete *idel;
  vsassert(m_assistants.empty());

  delete_us.clear();
  for(std::map<std::string, Variable*>::iterator ivar = m_vars.begin();
      ivar != m_vars.end(); ivar++)delete_us.insert(ivar->second->reader);
  for(std::set<VSOctaveH5ReaderBase*>::iterator idel = delete_us.begin();
      idel != delete_us.end(); idel++)delete *idel;
  for(std::map<std::string, Variable*>::iterator ivar = m_vars.begin();
      ivar != m_vars.end(); ivar++)
    {
      vsassert(ivar->second->reader == 0);
      H5Gclose(ivar->second->gid);
      delete ivar->second;
      ivar->second = 0;
    }
  vsassert(m_children.empty());

  if(m_my_file)H5Fclose(m_gid);
  else H5Gclose(m_gid);
}

void VSOctaveH5ReaderStruct::readGroup()
{
  m_index = 0;
  if(H5Giterate(m_gid, ".", NULL, &iterateOverGrp, this) != 0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not iterate over members of group"));
      throw(e);
    }
}

herr_t VSOctaveH5ReaderStruct::
iterateOverGrp(hid_t gid, const char*  name, void *object)
{
  static_cast<VSOctaveH5ReaderStruct*>(object)->registerVariable(name);
  return 0;
}

void VSOctaveH5ReaderStruct::registerVariable(const char*  name)
{
  // Skip datasets and what-not (like the cell array "dims"!)
  if(H5Gget_objtype_by_idx(m_gid,m_index++) != H5G_GROUP)
    return;

  hid_t gid = H5Gopen(m_gid,name);
  if(gid<0)
    {
      VSOctaveH5Exception 
	e(std::string("Could not open group: ")+std::string(name));
      throw(e);
    }

  Variable* var = new Variable;
  var->name = name;
  var->parent = this;
  var->gid = gid;
  var->readAttributes();

  uint8_t is_octave = false;
  if((!var->getAttribute(OCTAVE_H5_ATTR_NEW_FORMAT, is_octave))||(!is_octave))
    {
      H5Gclose(gid);
      return;
    }

  std::string super_type;
  getCString(gid, OCTAVE_H5_EL_TYPE, super_type);
  if(super_type == "bool")super_type = "bool scalar";

  std::string::iterator ifind = super_type.begin();
  while((ifind != super_type.end())&&(!isspace(*ifind)))ifind++;
  
  if(ifind == super_type.end())
    {
      var->type = super_type;
      var->eltype = "";
    }
  else
    {
      var->eltype = std::string(super_type.begin(), ifind);
      ifind++;
      var->type = std::string(ifind, super_type.end());
    }

  m_vars[name]=var;
}

void VSOctaveH5ReaderStruct::
registerAssistant(VSOctaveH5ReaderBase *reader)
{
  m_assistants.insert(reader);
}

void VSOctaveH5ReaderStruct::
notifyOfChildDestruction(VSOctaveH5ReaderBase* child)
{
  if(m_assistants.find(child) != m_assistants.end())
    {
      m_assistants.erase(child);
      return;
    }

  std::map<std::string, Variable*>::iterator ivar = m_vars.begin();
  while(ivar != m_vars.end())
    {
      if(ivar->second->reader == child)break;
      ivar++;
    }
  vsassert(ivar != m_vars.end());
  ivar->second->reader=0;
  m_children.erase(child);
}

VSOctaveH5ReaderBase::Variable* VSOctaveH5ReaderStruct::
resolveVariable(const std::string& name)
{
  std::string::const_iterator ifind = name.begin();
  while((ifind != name.end())&&(*ifind != '.')&&(*ifind != '{'))ifind++;
  
  std::string mine = std::string(name.begin(), ifind);
  
  std::map<std::string, Variable*>::iterator ivar = m_vars.find(mine);
  if(ivar == m_vars.end())return 0;
  else if(ifind == name.end())return ivar->second;
  else if(*ifind == '.')
    {
      ifind++;
      std::string rest = std::string(ifind,name.end());
      if(rest.empty())return 0; // parse error
      VSOctaveH5ReaderBase* r = loadVariable(ivar->second);
      VSOctaveH5ReaderStruct* rs = 
	dynamic_cast<VSOctaveH5ReaderStruct*>(r);
      if(rs == 0)return 0;
      return rs->resolveVariable(rest);
    }
  else /* (*ifind == '{') */
    {
      ifind++;
      std::string::const_iterator ifind2 = ifind;
      while((ifind2 != name.end())&&(*ifind2 != '}'))ifind2++;
      if(ifind2 == name.end())return 0; // parse error
      VSOctaveH5ReaderBase* r = loadVariable(ivar->second);
      VSOctaveH5ReaderCellArray* rc = 
	dynamic_cast<VSOctaveH5ReaderCellArray*>(r);

      std::string address(ifind,ifind2);
      unsigned row = 0;
      unsigned col = 0;
      unsigned n = sscanf(address.c_str(),"%u,%u",&row,&col);
      if(n == 0)return 0; // parse error
      else if(n == 1)
	{
	  unsigned nrow = 0;
	  unsigned ncol = 0;
	  rc->dimensions(nrow,ncol);
	  if(nrow==1){ col=row; row=0; }
	}
      if(!rc->isValid(row,col))return 0;
      std::string rest = rc->name(row,col);
      rest.append(++ifind2,name.end());
      return rc->m_struct->resolveVariable(rest);
    }
  return 0;
}

VSOctaveH5ReaderBase* VSOctaveH5ReaderStruct::loadVariable(Variable* var)
{
  if(var->reader)return var->reader;
  if(var->type == OCTAVE_H5_TYPE_STRUCT)
    {
      var->reader = new VSOctaveH5ReaderStruct(var);
      var->parent->registerChild(var->reader);
    }
  else if(var->type == OCTAVE_H5_TYPE_CELL)
    {
      var->reader = new VSOctaveH5ReaderCellArray(var);
      var->parent->registerChild(var->reader);
    }
  else if(var->type == OCTAVE_H5_TYPE_STRING)
    {
      var->reader = new VSOctaveH5ReaderMatrix(var);
      var->parent->registerChild(var->reader);      
    }
  else if(var->type == OCTAVE_H5_TYPE_MATRIX)
    {
      var->reader = new VSOctaveH5ReaderMatrix(var);
      var->parent->registerChild(var->reader);      
    }
  else if(var->type == OCTAVE_H5_TYPE_SCALAR)
    {
      var->reader = new VSOctaveH5ReaderMatrix(var);
      var->parent->registerChild(var->reader);      
    }
  else
    {
      VSOctaveH5Exception 
	e(std::string("Unknown variable type: ")+var->type);
      throw(e);
    }
  return var->reader;
}

std::vector<std::string> VSOctaveH5ReaderStruct::variables() const
{
  std::vector<std::string> names;
  names.reserve(m_vars.size());
  for(std::map<std::string, Variable*>::const_iterator ivar = m_vars.begin();
      ivar != m_vars.end(); ivar++)
    names.push_back(ivar->first);
  return names;
}

bool VSOctaveH5ReaderStruct::hasVariable(const std::string& name)
{
  return resolveVariable(name) != 0;
}

std::string VSOctaveH5ReaderStruct::variableType(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if(var == 0)return std::string();
  else return var->type;
}

std::string VSOctaveH5ReaderStruct::elementType(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if(var == 0)return std::string();
  else return var->eltype;
}

bool VSOctaveH5ReaderStruct::isValid(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if(var == 0)return false;
  else return true;
}

bool VSOctaveH5ReaderStruct::isEmpty(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if(var == 0)return false;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  VSOctaveH5ReaderMatrix* mtx = dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
  if(mtx != 0)return mtx->isEmpty();
  VSOctaveH5ReaderStruct* str = dynamic_cast<VSOctaveH5ReaderStruct*>(reader);
  if(str != 0)return str->isEmpty();
  VSOctaveH5ReaderCellArray* cel = 
    dynamic_cast<VSOctaveH5ReaderCellArray*>(reader);
  if(cel != 0)return cel->isEmpty();
  vsassert(0);
  return false;
}

bool VSOctaveH5ReaderStruct::isScalar(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)
     ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	(var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  VSOctaveH5ReaderMatrix* mtx = dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
  vsassert(mtx != 0); // since type is TYPE_MATRIX or TYPE_SCALAR
  return mtx->isScalar();
}

bool VSOctaveH5ReaderStruct::isVector(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)
     ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	(var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  VSOctaveH5ReaderMatrix* mtx = dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
  vsassert(mtx != 0); // since type is TYPE_MATRIX or TYPE_SCALAR
  return mtx->isVector();
}

bool VSOctaveH5ReaderStruct::isMatrix(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)
     ||((var->type != OCTAVE_H5_TYPE_MATRIX)&&
	(var->type != OCTAVE_H5_TYPE_SCALAR)))return false;
  else return true;
}

bool VSOctaveH5ReaderStruct::isString(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_STRING))return false;
  else return true;
}

bool VSOctaveH5ReaderStruct::isCellArray(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_CELL))return false;
  else return true; 
}

bool VSOctaveH5ReaderStruct::isCellVector(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_CELL))return false;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  if(reader==0)return 0;
  VSOctaveH5ReaderCellArray* rca = 
    dynamic_cast<VSOctaveH5ReaderCellArray*>(reader);
  if(rca==0)return 0;
  unsigned rows;
  unsigned cols;
  rca->dimensions(rows,cols);
  if(rows!=1 && cols!=1)return false;
  else return true; 
}

bool VSOctaveH5ReaderStruct::isStruct(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_STRUCT))return false;
  else return true;
}

bool VSOctaveH5ReaderStruct::
dimensions(const std::string& name, unsigned& rows, unsigned& cols)
{
  Variable* var = resolveVariable(name);
  if(var == 0)return 0;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  if(reader==0)return 0;

  VSOctaveH5ReaderStruct* rs = dynamic_cast<VSOctaveH5ReaderStruct*>(reader);
  if(rs)return rs->dimensions(rows,cols);

  VSOctaveH5ReaderMatrix* rm = dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
  if(rm)return rm->dimensions(rows,cols);

  VSOctaveH5ReaderCellArray* rc = 
    dynamic_cast<VSOctaveH5ReaderCellArray*>(reader);
  if(rc)return rc->dimensions(rows,cols);

  return false;
}

bool VSOctaveH5ReaderStruct::
readString(const std::string& name, std::string& s)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_STRING))return false;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  VSOctaveH5ReaderMatrix* mtx = 
    dynamic_cast<VSOctaveH5ReaderMatrix*>(reader);
  vsassert(mtx != 0); // since type is TYPE_STRING
  return mtx->readString(s);
}

VSOctaveH5ReaderStruct* VSOctaveH5ReaderStruct::
readStruct(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_STRUCT))return 0;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  if(reader==0)return 0;
  return dynamic_cast<VSOctaveH5ReaderStruct*>(reader);
}

VSOctaveH5ReaderCellArray* VSOctaveH5ReaderStruct::
readCellArray(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_CELL))return 0;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  if(reader==0)return 0;
  return dynamic_cast<VSOctaveH5ReaderCellArray*>(reader);
}

VSOctaveH5ReaderCellVector* VSOctaveH5ReaderStruct::
readCellVector(const std::string& name)
{
  Variable* var = resolveVariable(name);
  if((var == 0)||(var->type != OCTAVE_H5_TYPE_CELL))return 0;
  VSOctaveH5ReaderBase* reader = loadVariable(var);
  if(reader==0)return 0;
  VSOctaveH5ReaderCellArray* rca = 
    dynamic_cast<VSOctaveH5ReaderCellArray*>(reader);
  if(rca==0)return 0;
  unsigned rows;
  unsigned cols;
  rca->dimensions(rows,cols);
  if(rows!=1 && cols!=1)return 0;
  VSOctaveH5ReaderCellVector* rcv =  new VSOctaveH5ReaderCellVector(rca);
  rca->registerChild(rcv);
  return rcv;
}

bool VSOctaveH5ReaderStruct::isHDF5(const std::string& filename)
{
  return H5Fis_hdf5(filename.c_str());
}

bool VSOctaveH5ReaderStruct::
copyTo(const std::string& sname,
       VSOctaveH5WriterStruct* s, std::string dname)
{
  if(dname.empty())dname=sname;
  if(!isValid(sname))return false;

  if(isStruct(sname))
    {
      VSOctaveH5ReaderStruct* rs = readStruct(sname); 
      if(rs==0)return false;
      VSOctaveH5WriterStruct* ws = s->writeStruct(dname);
      if(ws==0)return false;
      if(!rs->copyAllTo(ws))return false;
      delete ws;
      return true;
    }
  else if(isCellArray(sname))
    {
      VSOctaveH5ReaderCellArray* rc = readCellArray(sname); 
      if(rc==0)return false;
      unsigned nrow;
      unsigned ncol;
      rc->dimensions(nrow, ncol);
      VSOctaveH5WriterCellArray* wc = s->writeCellArray(dname,nrow,ncol);
      if(wc==0)return false;
      bool _ret = true;
      for(unsigned irow=0;irow<nrow;irow++)
	for(unsigned icol=0;icol<ncol;icol++)
	  if(!rc->copyTo(irow,icol,wc,irow,icol))_ret = false;
      delete wc;
      return _ret;
    }
  else if(isString(sname))
    {
      std::string str;
      if(!readString(sname, str))return false;
      return(s->writeString(dname, str));
    }
  else if(isScalar(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealScalar<bool>(sname, s, dname));
      else if(eltype == "uint8")
	return(copyRealScalar<uint8_t>(sname, s, dname));
      else if(eltype == "uint16")
	return(copyRealScalar<uint16_t>(sname, s, dname));
      else if(eltype == "uint32")
	return(copyRealScalar<uint32_t>(sname, s, dname));
      else if(eltype == "uint64")
	return(copyRealScalar<uint64_t>(sname, s, dname));
      else if(eltype == "int8")
	return(copyRealScalar<int8_t>(sname, s, dname));
      else if(eltype == "int16")
	return(copyRealScalar<int16_t>(sname, s, dname));
      else if(eltype == "int32")
	return(copyRealScalar<int32_t>(sname, s, dname));
      else if(eltype == "int64")
	return(copyRealScalar<int64_t>(sname, s, dname));
      else 
	return(copyRealScalar<double>(sname, s, dname));
    }
  else if(isMatrix(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealMatrix<bool>(sname, s, dname));
      else if(eltype == "uint8")
	return(copyRealMatrix<uint8_t>(sname, s, dname));
      else if(eltype == "uint16")
	return(copyRealMatrix<uint16_t>(sname, s, dname));
      else if(eltype == "uint32")
	return(copyRealMatrix<uint32_t>(sname, s, dname));
      else if(eltype == "uint64")
	return(copyRealMatrix<uint64_t>(sname, s, dname));
      else if(eltype == "int8")
	return(copyRealMatrix<int8_t>(sname, s, dname));
      else if(eltype == "int16")
	return(copyRealMatrix<int16_t>(sname, s, dname));
      else if(eltype == "int32")
	return(copyRealMatrix<int32_t>(sname, s, dname));
      else if(eltype == "int64")
	return(copyRealMatrix<int64_t>(sname, s, dname));
      else 
	return(copyRealMatrix<double>(sname, s, dname));
    }
  else
    {
      return false;
    }
  
  assert(0);
}

bool VSOctaveH5ReaderStruct::
copyTo(const std::string& sname,
       VSOctaveH5WriterCellArray* c, unsigned dirow, unsigned dicol)
{
  if(!isValid(sname))return false;

  if(isStruct(sname))
    {
      VSOctaveH5ReaderStruct* rs = readStruct(sname); 
      if(rs==0)return false;
      VSOctaveH5WriterStruct* ws = c->writeStruct(dirow,dicol);
      if(ws==0)return false;
      if(!rs->copyAllTo(ws))return false;
      delete ws;
      return true;
    }
  else if(isCellArray(sname))
    {
      VSOctaveH5ReaderCellArray* rc = readCellArray(sname); 
      if(rc==0)return false;
      unsigned nrow;
      unsigned ncol;
      rc->dimensions(nrow, ncol);
      VSOctaveH5WriterCellArray* wc = c->writeCellArray(dirow,dicol,nrow,ncol);
      if(wc==0)return false;
      bool _ret = true;
      for(unsigned irow=0;irow<nrow;irow++)
	for(unsigned icol=0;icol<ncol;icol++)
	  _ret &= rc->copyTo(irow,icol,wc,irow,icol);
      delete wc;
      return _ret;
    }
  else if(isString(sname))
    {
      std::string str;
      if(!readString(sname, str))return false;
      return(c->writeString(dirow, dicol, str));
    }
  else if(isScalar(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealScalar<bool>(sname, c, dirow, dicol));
      else if(eltype == "uint8")
	return(copyRealScalar<uint8_t>(sname, c, dirow, dicol));
      else if(eltype == "uint16")
	return(copyRealScalar<uint16_t>(sname, c, dirow, dicol));
      else if(eltype == "uint32")
	return(copyRealScalar<uint32_t>(sname, c, dirow, dicol));
      else if(eltype == "uint64")
	return(copyRealScalar<uint64_t>(sname, c, dirow, dicol));
      else if(eltype == "int8")
	return(copyRealScalar<int8_t>(sname, c, dirow, dicol));
      else if(eltype == "int16")
	return(copyRealScalar<int16_t>(sname, c, dirow, dicol));
      else if(eltype == "int32")
	return(copyRealScalar<int32_t>(sname, c, dirow, dicol));
      else if(eltype == "int64")
	return(copyRealScalar<int64_t>(sname, c, dirow, dicol));
      else 
	return(copyRealScalar<double>(sname, c, dirow, dicol));
    }
  else if(isMatrix(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealMatrix<bool>(sname, c, dirow, dicol));
      else if(eltype == "uint8")
	return(copyRealMatrix<uint8_t>(sname, c, dirow, dicol));
      else if(eltype == "uint16")
	return(copyRealMatrix<uint16_t>(sname, c, dirow, dicol));
      else if(eltype == "uint32")
	return(copyRealMatrix<uint32_t>(sname, c, dirow, dicol));
      else if(eltype == "uint64")
	return(copyRealMatrix<uint64_t>(sname, c, dirow, dicol));
      else if(eltype == "int8")
	return(copyRealMatrix<int8_t>(sname, c, dirow, dicol));
      else if(eltype == "int16")
	return(copyRealMatrix<int16_t>(sname, c, dirow, dicol));
      else if(eltype == "int32")
	return(copyRealMatrix<int32_t>(sname, c, dirow, dicol));
      else if(eltype == "int64")
	return(copyRealMatrix<int64_t>(sname, c, dirow, dicol));
      else 
	return(copyRealMatrix<double>(sname, c, dirow, dicol));
    }
  else
    {
      return false;
    }
  
  assert(0);
}

bool VSOctaveH5ReaderStruct::
copyTo(const std::string& sname,
       VSOctaveH5WriterCellVector* c, unsigned diel)
{
  if(!isValid(sname))return false;

  if(isStruct(sname))
    {
      VSOctaveH5ReaderStruct* rs = readStruct(sname); 
      if(rs==0)return false;
      VSOctaveH5WriterStruct* ws = c->writeStruct(diel);
      if(ws==0)return false;
      if(!rs->copyAllTo(ws))return false;
      delete ws;
      return true;
    }
  else if(isCellArray(sname))
    {
      VSOctaveH5ReaderCellArray* rc = readCellArray(sname); 
      if(rc==0)return false;
      unsigned nrow;
      unsigned ncol;
      rc->dimensions(nrow, ncol);
      VSOctaveH5WriterCellArray* wc = c->writeCellArray(diel,nrow,ncol);
      if(wc==0)return false;
      bool _ret = true;
      for(unsigned irow=0;irow<nrow;irow++)
	for(unsigned icol=0;icol<ncol;icol++)
	  _ret &= rc->copyTo(irow,icol,wc,irow,icol);
      delete wc;
      return _ret;
    }
  else if(isString(sname))
    {
      std::string str;
      if(!readString(sname, str))return false;
      return(c->writeString(diel, str));
    }
  else if(isScalar(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealScalar<bool>(sname, c, diel));
      else if(eltype == "uint8")
	return(copyRealScalar<uint8_t>(sname, c, diel));
      else if(eltype == "uint16")
	return(copyRealScalar<uint16_t>(sname, c, diel));
      else if(eltype == "uint32")
	return(copyRealScalar<uint32_t>(sname, c, diel));
      else if(eltype == "uint64")
	return(copyRealScalar<uint64_t>(sname, c, diel));
      else if(eltype == "int8")
	return(copyRealScalar<int8_t>(sname, c, diel));
      else if(eltype == "int16")
	return(copyRealScalar<int16_t>(sname, c, diel));
      else if(eltype == "int32")
	return(copyRealScalar<int32_t>(sname, c, diel));
      else if(eltype == "int64")
	return(copyRealScalar<int64_t>(sname, c, diel));
      else 
	return(copyRealScalar<double>(sname, c, diel));
    }
  else if(isMatrix(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealMatrix<bool>(sname, c, diel));
      else if(eltype == "uint8")
	return(copyRealMatrix<uint8_t>(sname, c, diel));
      else if(eltype == "uint16")
	return(copyRealMatrix<uint16_t>(sname, c, diel));
      else if(eltype == "uint32")
	return(copyRealMatrix<uint32_t>(sname, c, diel));
      else if(eltype == "uint64")
	return(copyRealMatrix<uint64_t>(sname, c, diel));
      else if(eltype == "int8")
	return(copyRealMatrix<int8_t>(sname, c, diel));
      else if(eltype == "int16")
	return(copyRealMatrix<int16_t>(sname, c, diel));
      else if(eltype == "int32")
	return(copyRealMatrix<int32_t>(sname, c, diel));
      else if(eltype == "int64")
	return(copyRealMatrix<int64_t>(sname, c, diel));
      else 
	return(copyRealMatrix<double>(sname, c, diel));
    }
  else
    {
      return false;
    }
  
  assert(0);
}

bool VSOctaveH5ReaderStruct::
copyTo(const std::string& sname,VSOctaveH5WriterExpandableCellVector* c)
{
  if(!isValid(sname))return false;

  if(isStruct(sname))
    {
      VSOctaveH5ReaderStruct* rs = readStruct(sname); 
      if(rs==0)return false;
      VSOctaveH5WriterStruct* ws = c->appendStruct();
      if(ws==0)return false;
      if(!rs->copyAllTo(ws))return false;
      delete ws;
      return true;
    }
  else if(isCellArray(sname))
    {
      VSOctaveH5ReaderCellArray* rc = readCellArray(sname); 
      if(rc==0)return false;
      unsigned nrow;
      unsigned ncol;
      rc->dimensions(nrow, ncol);
      VSOctaveH5WriterCellArray* wc = c->appendCellArray(nrow,ncol);
      if(wc==0)return false;
      bool _ret = true;
      for(unsigned irow=0;irow<nrow;irow++)
	for(unsigned icol=0;icol<ncol;icol++)
	  _ret &= rc->copyTo(irow,icol,wc,irow,icol);
      delete wc;
      return _ret;
    }
  else if(isString(sname))
    {
      std::string str;
      if(!readString(sname, str))return false;
      return(c->appendString(str));
    }
  else if(isScalar(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealScalar<bool>(sname, c));
      else if(eltype == "uint8")
	return(copyRealScalar<uint8_t>(sname, c));
      else if(eltype == "uint16")
	return(copyRealScalar<uint16_t>(sname, c));
      else if(eltype == "uint32")
	return(copyRealScalar<uint32_t>(sname, c));
      else if(eltype == "uint64")
	return(copyRealScalar<uint64_t>(sname, c));
      else if(eltype == "int8")
	return(copyRealScalar<int8_t>(sname, c));
      else if(eltype == "int16")
	return(copyRealScalar<int16_t>(sname, c));
      else if(eltype == "int32")
	return(copyRealScalar<int32_t>(sname, c));
      else if(eltype == "int64")
	return(copyRealScalar<int64_t>(sname, c));
      else 
	return(copyRealScalar<double>(sname, c));
    }
  else if(isMatrix(sname))
    {
      std::string eltype = elementType(sname);
      if(eltype == "bool")
	return(copyRealMatrix<bool>(sname, c));
      else if(eltype == "uint8")
	return(copyRealMatrix<uint8_t>(sname, c));
      else if(eltype == "uint16")
	return(copyRealMatrix<uint16_t>(sname, c));
      else if(eltype == "uint32")
	return(copyRealMatrix<uint32_t>(sname, c));
      else if(eltype == "uint64")
	return(copyRealMatrix<uint64_t>(sname, c));
      else if(eltype == "int8")
	return(copyRealMatrix<int8_t>(sname, c));
      else if(eltype == "int16")
	return(copyRealMatrix<int16_t>(sname, c));
      else if(eltype == "int32")
	return(copyRealMatrix<int32_t>(sname, c));
      else if(eltype == "int64")
	return(copyRealMatrix<int64_t>(sname, c));
      else 
	return(copyRealMatrix<double>(sname, c));
    }
  else
    {
      return false;
    }
  
  assert(0);
}

bool VSOctaveH5ReaderStruct::
copyAllTo(VSOctaveH5WriterStruct* s)
{
  std::vector<std::string> fieldnames = variables();
  unsigned nfield = fieldnames.size();
  // Loop over fields
  bool _ret = true;
  for(unsigned ifield=0;ifield<nfield;ifield++)
    _ret &= copyTo(fieldnames[ifield],s);
  return _ret;
}


// ****************************************************************************
//
// TEST MAIN FUNCTION
//
// ****************************************************************************

#ifdef TEST_MAIN_H5READER

// g++ -DTEST_MAIN_H5READER -I. -I../VSCommon -L../VSCommon -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -g -Wall  -o test VATime.o VSOctaveH5Reader.cpp VSOctaveH5Writer.cpp VSFileUtility.cpp VSOctaveIO.cpp -lVSCommon -lhdf5

#include<iostream>
#include<iterator>
#include<algorithm>

#include<VSTime.hpp>

int main(int argc, char** argv)
{
  VSOctaveH5Reader* reader;
  try
    {
      reader = new VSOctaveH5Reader("test.h5");

      std::vector<std::string> vars = reader->variables();
      for(std::vector<std::string>::const_iterator ivar = vars.begin();
	  ivar != vars.end(); ivar++)
	{
	  std::cerr << std::setw(30) << std::left << *ivar << ' '
		    << std::setw(10) << reader->variableType(*ivar) << ' '
		    << std::setw(10) << reader->elementType(*ivar)
		    << std::endl;
	}

      VSOctaveH5ReaderStruct* reader2 = reader->readStruct("s");
      vars = reader2->variables();
#if 1
      for(std::vector<std::string>::const_iterator ivar = vars.begin();
	  ivar != vars.end(); ivar++)
	{
	  std::cerr << "s." << std::setw(28) << std::left << *ivar << ' '
		    << std::setw(10) << reader->variableType(*ivar) << ' '
		    << std::setw(10) << reader->elementType(*ivar)
		    << std::endl;
	}
#endif
      std::cerr << std::endl;

      unsigned rows;
      unsigned cols;
      bool good;
      good = reader->dimensions("c",rows,cols);
      std::cerr << good << ' ' << rows << ' ' << cols << std::endl;
      good = reader->dimensions("r",rows,cols);
      std::cerr << good << ' ' << rows << ' ' << cols << std::endl;
      good = reader->dimensions("m",rows,cols);
      std::cerr << good << ' ' << rows << ' ' << cols << std::endl;
      good = reader->dimensions("empty_cv",rows,cols);
      std::cerr << good << ' ' << rows << ' ' << cols << std::endl;
      good = reader->dimensions("u",rows,cols);
      std::cerr << good << ' ' << rows << ' ' << cols << std::endl;
      std::cerr << std::endl;
      
      double x;
      double y;
      unsigned u;
      double ud;

      unsigned ix;
      unsigned iy;
      double fix;
      double fiy;

      good = reader->readScalar("x",x);
      std:: cerr << good << ' ' << x << std::endl;
      good = reader->readScalar("y",y);
      std:: cerr << good << ' ' << y << std::endl;

      good = reader->readScalar("ix",ix);
      std:: cerr << good << ' ' << ix << std::endl;
      good = reader->readScalar("iy",iy);
      std:: cerr << good << ' ' << iy << std::endl;
      good = reader->readScalar("ix",fix);
      std:: cerr << good << ' ' << fix << std::endl;
      good = reader->readScalar("iy",fiy);
      std:: cerr << good << ' ' << fiy << std::endl;

      good = reader->readScalar("u",u);
      std:: cerr << good << ' ' << u << std::endl;
      good = reader->readScalar("u",ud);
      std:: cerr << good << ' ' << ud << std::endl;
      std::cerr << std::endl;

      std::vector<double> c;
      good = reader->readVector("c",c);
      std::copy(c.begin(), c.end(),
		std::ostream_iterator<double>(std::cerr," "));
      std::cerr << std::endl;

      std::vector<float> cu;
      good = reader->readVector("c",cu);
      std::copy(cu.begin(), cu.end(),
		std::ostream_iterator<double>(std::cerr," "));
      std::cerr << std::endl;
      std::cerr << std::endl;
      
      std::vector<unsigned> ic;
      good = reader->readVector("ic",ic);
      std::copy(ic.begin(), ic.end(),
		std::ostream_iterator<double>(std::cerr," "));
      std::cerr << std::endl;

      std::vector<double> cic;
      good = reader->readVector("ic",cic);
      std::copy(cic.begin(), cic.end(),
		std::ostream_iterator<double>(std::cerr," "));
      std::cerr << std::endl;

      std::vector<double> empty_cv;
      good = reader->readVector("empty_cv",empty_cv);
      std::copy(empty_cv.begin(), empty_cv.end(),
		std::ostream_iterator<double>(std::cerr," "));
      std::cerr << std::endl;

      unsigned mr;
      unsigned mc;
      double* m;
      good = reader->readMatrix("m",mr,mc,m);
      for(unsigned ir=0;ir<mr;ir++)
	{
	  std::copy(m+ir*mc,m+(ir+1)*mc,
		    std::ostream_iterator<double>(std::cerr," "));
	  std::cerr << std::endl;
	}
      delete[] m;
      std::cerr << std::endl;

      double m2[3][4];
      good = reader->readMatrix("m",(double*)m2);
      for(unsigned ir=0;ir<3;ir++)
	{
	  for(unsigned ic=0;ic<4;ic++)
	    std::cerr << m2[ir][ic] << ' ';
	  std::cerr << std::endl;
	}
      std::cerr << std::endl;

      std::string a;
      good = reader->readString("a",a);
      std::cerr << a << std::endl;
      std::cerr << std::endl;

      good = reader->readVector("a",cu);
      std::cerr << good << std::endl;
      good = reader->readScalar("a",x);
      std::cerr << good << std::endl;
      good = reader->readScalar("c",x);
      std::cerr << good << std::endl;
      std::cerr << std::endl;

      VSOctaveH5ReaderCellArray* rc = reader->readCellArray("ca");
      rc->dimensions(rows,cols);
      std::cerr << rows << ' ' << cols << std::endl;
      std::cerr << std::endl;

      for(unsigned ir=0;ir<rows;ir++)
	for(unsigned ic=0;ic<cols;ic++)
	  {
	    unsigned r=0;
	    unsigned c=0;
	    rc->dimensions(ir,ic,r,c);
	    std::cerr << rc->variableType(ir,ic) << ' ' << r << ' ' << c
		      << std::endl;
	  }
      std::cerr << std::endl;

      VSOctaveH5ReaderCellVector* rcv = reader->readCellVector("cv");
      std::cerr << "cell vec: " << rcv->dimensions() << std::endl;
      std::cerr << std::endl;

      for(unsigned iel=0;iel<rcv->dimensions();iel++)
	  {
	    unsigned r=0;
	    unsigned c=0;
	    rcv->dimensions(iel,r,c);
	    std::cerr << rcv->variableType(iel) << ' ' << r << ' ' << c
		      << std::endl;
	  }
      std::cerr << std::endl;

      rcv = reader->readCellVector("cv-r");
      std::cerr << "cell vec (row): " << rcv->dimensions() << std::endl;
      std::cerr << std::endl;

      for(unsigned iel=0;iel<rcv->dimensions();iel++)
	  {
	    unsigned r=0;
	    unsigned c=0;
	    rcv->dimensions(iel,r,c);
	    std::cerr << rcv->variableType(iel) << ' ' << r << ' ' << c
		      << std::endl;
	  }
      std::cerr << std::endl;

      rcv = reader->readCellVector("cv-c");
      std::cerr << "cell vec (col): " << rcv->dimensions() << std::endl;
      std::cerr << std::endl;

      for(unsigned iel=0;iel<rcv->dimensions();iel++)
	  {
	    unsigned r=0;
	    unsigned c=0;
	    rcv->dimensions(iel,r,c);
	    std::cerr << rcv->variableType(iel) << ' ' << r << ' ' << c
		      << std::endl;
	  }
      std::cerr << std::endl;

      VSOctaveH5ReaderVector<double>* vec =
	reader->readVectorElements<double>("c",3);
      for(unsigned irow=0;irow<vec->rows();irow++)
	std::cerr << (*vec)[irow] << std::endl;
      std::cerr << std::endl;

      std::vector<bool> bv;
      reader->readVector("bv",bv);
      std::copy(bv.begin(), bv.end(),
		std::ostream_iterator<bool>(std::cerr," "));
      std::cerr << std::endl;

      std::vector<double> bvd;
      reader->readVector("bv",bvd);
      std::copy(bvd.begin(), bvd.end(),
		std::ostream_iterator<double>(std::cerr," "));
      std::cerr << std::endl;

      std::vector<unsigned> bvu;
      reader->readVector("bv",bvu);
      std::copy(bvu.begin(), bvu.end(),
		std::ostream_iterator<unsigned>(std::cerr," "));
      std::cerr << std::endl;
      std::cerr << std::endl;

      //std::vector<VSTime> vstime;
      //reader->readCompositeVector("vstime",vstime);
      VSOctaveH5Reader* temps = reader->readStruct("vstime");
      vsassert(temps);
      //temps->readCompositeVectorHere(vstime);
      VSOctaveH5ReaderCompositeVector<VSTime>* v =
	temps->readCompositeVectorElementsHere<VSTime>();
      for(unsigned i=0;i<v->rows();i+=200)
	std::cout << v->at(i).toString() << std::endl;
      delete v;
      delete temps;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << e.message() << std::endl;
    }

  try
    {
      delete reader;
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << e.message() << std::endl;
    }
}
#endif
