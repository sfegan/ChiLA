//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOctaveH5Library.cpp

  Functions to implement library of objects indexed by a tree of "paths" in
  the octave H5 format 

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/21/2007

  $Id: VSH5Library.cpp,v 1.12 2009/04/29 01:16:49 matthew Exp $

*/

#include<VSH5Library.hpp>
#include<VSDataConverter.hpp>

using namespace VERITAS;

// ============================================================================
//
// Class: VSH5LibraryPath
//
// ============================================================================

VSH5LibraryPath::~VSH5LibraryPath()
{
  // nothing to see here
}

const VSH5LibraryDefaultPath* VSH5LibraryPath::castToDefaultPath() const
{
  return 0;
}

const VSH5LibraryNamedPath* VSH5LibraryPath::castToNamedPath() const
{
  return 0;
}

const VSH5LibraryIndexedPath* VSH5LibraryPath::castToIndexedPath() const
{
  return 0;
}

const VSH5LibraryRangedPath* VSH5LibraryPath::castToRangedPath() const
{
  return 0;
}

bool VSH5LibraryPath::isValidForDefinition() const
{
  return false;
}

bool VSH5LibraryPath::isValidForReading(const VSH5LibraryPath* val) const
{
  return false;
}

bool VSH5LibraryPath::isValidForWriting(const VSH5LibraryPath* val) const
{
  return false;
}

VSH5LibraryPath* VSH5LibraryPath::
readPath(VSOctaveH5ReaderStruct* s)
{
  std::string _type;
  std::string _name;
  s->readString("type",_type);
  s->readString("name",_name);

  if(_type == VSH5LibraryDefaultPath::typeName())
    {
      vsassert(false);
    }
  else if(_type == VSH5LibraryNamedPath::typeName())
    {
      return VSH5LibraryNamedPath::readPath(_name,s);
    }
  else if(_type == VSH5LibraryIndexedPath::typeName())
    {
      return VSH5LibraryIndexedPath::readPath(_name,s);
    }
  else if(_type == VSH5LibraryRangedPath::typeName())
    {
      return VSH5LibraryRangedPath::readPath(_name,s);
    }
  else
    {
      vsassert(false);
    }
}

void VSH5LibraryPath::writePath(VSOctaveH5WriterStruct* s) const
{
  s->writeString("type",this->type());
  s->writeString("name",this->name());
  this->writePathDetails(s);
}
    
void VSH5LibraryPath::writePathDetails(VSOctaveH5WriterStruct* s) const
{
  // nothing to see here
}

VSH5LibraryPath::Writer::~Writer()
{
  // nothing to see here
}

// ============================================================================
//
// Class: VSH5LibraryDefaultPath
//
// ============================================================================

VSH5LibraryDefaultPath::~VSH5LibraryDefaultPath()
{
  // nothing to see here
}

bool VSH5LibraryDefaultPath::allowsDefault() const
{
  return false;
}

std::string VSH5LibraryDefaultPath::type() const
{
  return typeName();
}

const VSH5LibraryDefaultPath* VSH5LibraryDefaultPath::castToDefaultPath() const
{
  return this;
}

bool VSH5LibraryDefaultPath::isValidForDefinition() const
{
  return false;
}

bool VSH5LibraryDefaultPath::isValidForReading(const VSH5LibraryPath* val) const
{
  // should never get here
  vsassert(0);
}

bool VSH5LibraryDefaultPath::isValidForWriting(const VSH5LibraryPath* val) const
{
  // should never get here
  vsassert(0);
}

VSH5LibraryPath* VSH5LibraryDefaultPath::copy() const
{
  return new VSH5LibraryDefaultPath(*this);
}

VSH5LibraryPath::Writer* VSH5LibraryDefaultPath::
getWriter(VSOctaveH5WriterStruct* s, const VSH5LibraryPathSet* remainder) const
{
  vsassert(false);
}

VSOctaveH5ReaderStruct* VSH5LibraryDefaultPath::
read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const
{
  vsassert(false);
}

VSH5LibraryPathSet* VSH5LibraryDefaultPath::
list(VSOctaveH5ReaderStruct* s) const
{
  vsassert(false);
}

// ============================================================================
//
// Class: VSH5LibraryNamedPath
//
// ============================================================================

VSH5LibraryNamedPath::~VSH5LibraryNamedPath()
{
  // nothing to see here
}

bool VSH5LibraryNamedPath::allowsDefault() const
{
  return true;
}

std::string VSH5LibraryNamedPath::type() const
{
  return typeName();
}

const VSH5LibraryNamedPath* VSH5LibraryNamedPath::castToNamedPath() const
{
  return this;
}

bool VSH5LibraryNamedPath::isValidForDefinition() const
{
  return fNamedValue.empty();
}

bool VSH5LibraryNamedPath::isValidForReading(const VSH5LibraryPath* val) const
{
  const VSH5LibraryNamedPath* cpath = val->castToNamedPath();
  if((val->name() == name())
     &&((val->castToDefaultPath())
	||((cpath)&&(!cpath->namedValue().empty()))))
    return true;
  else return false;
}

bool VSH5LibraryNamedPath::isValidForWriting(const VSH5LibraryPath* val) const
{
  return isValidForReading(val);
}

VSH5LibraryNamedPath* VSH5LibraryNamedPath::
readPath(const std::string& _name, VSOctaveH5ReaderStruct* s)
{
  return new VSH5LibraryNamedPath(_name);
}

VSH5LibraryPath* VSH5LibraryNamedPath::copy() const
{
  return new VSH5LibraryNamedPath(*this);
}

VSH5LibraryNamedPath::Writer::
Writer(VSOctaveH5WriterStruct* s, const VSH5LibraryNamedPath* path,
       const VSH5LibraryPathSet *remainder):
  VSH5LibraryPath::Writer(),
  fPath(path), fRemainder(remainder), 
  fBaseStruct(s), fDataStruct(s->writeStruct("data")), fData(), fDefault()
{
  // nothing to see here
}

VSH5LibraryNamedPath::Writer::~Writer()
{
  delete fDefault;
  for(std::map<std::string,VSH5LibraryPath::Writer*>::iterator 
	iwriter=fData.begin(); iwriter!=fData.end(); iwriter++)
    delete iwriter->second;
  delete fDataStruct;
  delete fBaseStruct;
  delete fRemainder;
  delete fPath;
}

VSOctaveH5WriterStruct* 
VSH5LibraryNamedPath::Writer::write(const VSH5LibraryPathSet* pathset)
{
  VSOctaveH5WriterStruct* return_struct = 0;

  vsassert(pathset);
  const VSH5LibraryPath* path = (*pathset)[0];  
  vsassert(fPath->isValidForWriting(path));
  VSH5LibraryPathSet* sub_pathset = pathset->shift();

  if(fRemainder)
    {
      vsassert(sub_pathset);
      VSH5LibraryPath::Writer* sub_writer = 0;
      if(path->castToDefaultPath())
	{
	  if(fDefault == 0)
	    {
	      VSH5LibraryPathSet* sub_remainder = fRemainder->shift();
	      VSOctaveH5WriterStruct* s = fBaseStruct->writeStruct("default");
	      fDefault = (*fRemainder)[0]->getWriter(s, sub_remainder);
	    }

	  sub_writer = fDefault;
	}
      else
	{
	  const VSH5LibraryNamedPath* cpath = path->castToNamedPath();
	  vsassert(cpath);
	  if(fData.find(cpath->namedValue()) == fData.end())
	    {
	      VSH5LibraryPathSet* sub_remainder = fRemainder->shift();
	      VSOctaveH5WriterStruct* s = 
		fDataStruct->writeStruct(cpath->namedValue());
	      fData[cpath->namedValue()] = 
		(*fRemainder)[0]->getWriter(s, sub_remainder);
	    }

	  sub_writer = fData[cpath->namedValue()];
	}

      return_struct = sub_writer->write(sub_pathset);
    }
  else
    { 
      vsassert(sub_pathset==0);
      if(path->castToDefaultPath())
	{
	  return_struct = fBaseStruct->writeStruct("default");
	}
      else
	{
	  const VSH5LibraryNamedPath* cpath = path->castToNamedPath();
	  vsassert(cpath);
	  return_struct = fDataStruct->writeStruct(cpath->namedValue());
	}
    }

  delete pathset;
  return return_struct;
}

VSH5LibraryPath::Writer* 
VSH5LibraryNamedPath::getWriter(VSOctaveH5WriterStruct* s,
				const VSH5LibraryPathSet* remainder) const
{
  return new Writer(s, new VSH5LibraryNamedPath(*this), remainder);
}

VSOctaveH5ReaderStruct* VSH5LibraryNamedPath::
read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const
{
  vsassert(path);
  vsassert(this->isValidForReading(path));
  
  const VSH5LibraryNamedPath* cpath = path->castToNamedPath();
  if(cpath)
    {
      VSOctaveH5ReaderStruct* sdata = s->readStruct("data");
      vsassert(sdata);
      if(sdata->isStruct(cpath->namedValue()))
	return sdata->readStruct(cpath->namedValue());
    }

  // If we get here then we need to check for the default since we
  // were either asked for it explicitly (cpath==0) or the value we
  // were asked for is not present

  if(s->isStruct("default"))return s->readStruct("default");
  
  return 0;
}

VSH5LibraryPathSet* VSH5LibraryNamedPath::
list(VSOctaveH5ReaderStruct* s) const
{
  VSOctaveH5ReaderStruct* sdata = s->readStruct("data");
  vsassert(sdata);
  std::vector<std::string> fn = sdata->fieldNames();
  bool has_default = s->isStruct("default");
  VSH5LibraryPathSet* lps = 
    new VSH5LibraryPathSet(fn.size()+(has_default?1:0));
  for(unsigned ifn = 0;ifn<fn.size();ifn++)
    lps->set(ifn, new VSH5LibraryNamedPath(name(),fn[ifn]));
  if(has_default)
    lps->set(fn.size(),new VSH5LibraryDefaultPath(name()));
  return lps;
}

// ============================================================================
//
// Class: VSH5LibraryIndexedPath
//
// ============================================================================

VSH5LibraryIndexedPath::~VSH5LibraryIndexedPath()
{
  // nothing to see here
}

bool VSH5LibraryIndexedPath::allowsDefault() const
{
  return true;
}

std::string VSH5LibraryIndexedPath::type() const
{
  return typeName();
}

const VSH5LibraryIndexedPath* VSH5LibraryIndexedPath::castToIndexedPath() const
{
  return this;
}

bool VSH5LibraryIndexedPath::isValidForDefinition() const
{
  return fIndexedValue>0;
}

bool VSH5LibraryIndexedPath::isValidForReading(const VSH5LibraryPath* val) const
{
  if(val->name() != name())return false;
  if(val->castToDefaultPath())return true;
  const VSH5LibraryIndexedPath* cval = val->castToIndexedPath();
  if(cval)return true;
  else return false;
}

bool VSH5LibraryIndexedPath::isValidForWriting(const VSH5LibraryPath* val) const
{
  return isValidForReading(val);
}

void VSH5LibraryIndexedPath::writePathDetails(VSOctaveH5WriterStruct* s) const
{
  s->writeScalar("count",fIndexedValue);
}

VSH5LibraryIndexedPath* VSH5LibraryIndexedPath::
readPath(const std::string& _name, VSOctaveH5ReaderStruct* s)
{
  unsigned count;
  s->readScalar("count",count);
  return new VSH5LibraryIndexedPath(_name,count);
}

VSH5LibraryPath* VSH5LibraryIndexedPath::copy() const
{
  return new VSH5LibraryIndexedPath(*this);
}

VSH5LibraryIndexedPath::Writer::
Writer(VSOctaveH5WriterStruct* s, const VSH5LibraryIndexedPath* path,
       const VSH5LibraryPathSet *remainder):
  VSH5LibraryPath::Writer(),
  fPath(path), fRemainder(remainder), fBaseStruct(s),
  fDataCellVector(), fData(path->indexedValue()), fDefault()
{
  fDataCellVector = fBaseStruct->writeCellVector("data",fData.size());
}

VSH5LibraryIndexedPath::Writer::~Writer()
{
  delete fDefault;
  for(std::vector<VSH5LibraryPath::Writer*>::iterator 
	iwriter=fData.begin(); iwriter!=fData.end(); iwriter++)
    delete (*iwriter);
  delete fDataCellVector;
  delete fBaseStruct;
  delete fRemainder;
  delete fPath;
}

VSOctaveH5WriterStruct* 
VSH5LibraryIndexedPath::Writer::write(const VSH5LibraryPathSet* pathset)
{
  VSOctaveH5WriterStruct* return_struct = 0;

  vsassert(pathset);
  const VSH5LibraryPath* path = (*pathset)[0];  
  vsassert(fPath->isValidForWriting(path));
  VSH5LibraryPathSet* sub_pathset = pathset->shift();

  if(fRemainder)
    {
      vsassert(sub_pathset);
      VSH5LibraryPath::Writer* sub_writer = 0;
      if(path->castToDefaultPath())
	{
	  if(fDefault == 0)
	    {
	      VSH5LibraryPathSet* sub_remainder = fRemainder->shift();
	      VSOctaveH5WriterStruct* s = fBaseStruct->writeStruct("default");
	      fDefault = (*fRemainder)[0]->getWriter(s, sub_remainder);
	    }

	  sub_writer = fDefault;
	}
      else
	{
	  const VSH5LibraryIndexedPath* cpath = path->castToIndexedPath();
	  vsassert(cpath);
	  vsassert(cpath->indexedValue() < fData.size());
	  if(!fData[cpath->indexedValue()])
	    {
	      VSH5LibraryPathSet* sub_remainder = fRemainder->shift();
	      VSOctaveH5WriterStruct* s = 
		fDataCellVector->writeStruct(cpath->indexedValue());
	      fData[cpath->indexedValue()] = 
		(*fRemainder)[0]->getWriter(s, sub_remainder);
	    }

	  sub_writer = fData[cpath->indexedValue()];
	}

      return_struct = sub_writer->write(sub_pathset);
    }
  else
    { 
      vsassert(sub_pathset==0);
      if(path->castToDefaultPath())
	{
	  return_struct = fBaseStruct->writeStruct("default");
	}
      else
	{
	  const VSH5LibraryIndexedPath* cpath = path->castToIndexedPath();
	  vsassert(cpath);
	  vsassert(cpath->indexedValue() < fData.size());
	  return_struct = fDataCellVector->writeStruct(cpath->indexedValue());
	}
    }

  delete pathset;
  return return_struct;
}

VSH5LibraryPath::Writer* 
VSH5LibraryIndexedPath::getWriter(VSOctaveH5WriterStruct* s,
				  const VSH5LibraryPathSet* remainder) const
{
  return new Writer(s, new VSH5LibraryIndexedPath(*this), remainder);
}

VSOctaveH5ReaderStruct* VSH5LibraryIndexedPath::
read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const
{
  vsassert(path);
  vsassert(this->isValidForReading(path));
  
  const VSH5LibraryIndexedPath* cpath = path->castToIndexedPath();
  if(cpath)
    {
      VSOctaveH5ReaderCellVector* cdata = s->readCellVector("data");

      vsassert(cdata);
      unsigned nkey = cdata->dimensions();
      vsassert(nkey==indexedValue());
      if(cpath->indexedValue() >= nkey && s->isStruct("default"))
	return s->readStruct("default");
      else if(cdata->isStruct(cpath->indexedValue()))
	return cdata->readStruct(cpath->indexedValue());
    }

  // If we get here then we need to check for the default since we
  // were either asked for it explicitly (cpath==0) or the value we
  // were asked for is not present

  if(s->isStruct("default"))return s->readStruct("default");
  
  return 0;
}

VSH5LibraryPathSet* VSH5LibraryIndexedPath::
list(VSOctaveH5ReaderStruct* s) const
{
  VSOctaveH5ReaderCellVector* cdata = s->readCellVector("data");
  vsassert(cdata);
  std::vector<unsigned> el;
  for(unsigned iel=0;iel<cdata->dimensions();iel++)
    if(!cdata->isEmpty(iel))el.push_back(iel);
  bool has_default = s->isStruct("default");
  VSH5LibraryPathSet* lps = 
    new VSH5LibraryPathSet(el.size()+(has_default?1:0));
  for(unsigned iel = 0;iel<el.size();iel++)
    lps->set(iel, new VSH5LibraryIndexedPath(name(),el[iel]));
  if(has_default)
    lps->set(el.size(),new VSH5LibraryDefaultPath(name()));
  return lps;
}

// ============================================================================
//
// Class: VSH5LibraryRangedPath
//
// ============================================================================

VSH5LibraryRangedPath::~VSH5LibraryRangedPath()
{
  // nothing to see here
}

bool VSH5LibraryRangedPath::allowsDefault() const
{
  return true;
}

std::string VSH5LibraryRangedPath::type() const
{
  return typeName();
}

const VSH5LibraryRangedPath* VSH5LibraryRangedPath::castToRangedPath() const
{
  return this;
}

bool VSH5LibraryRangedPath::isValidForDefinition() const
{
  return (fRangedValueLo == 0)&&(fRangedValueHi == 0);
}

bool VSH5LibraryRangedPath::isValidForReading(const VSH5LibraryPath* val) const
{
  const VSH5LibraryRangedPath* cval = val->castToRangedPath();
  if((val->name() == name())
     &&((val->castToDefaultPath())||(cval)))return true;
  else return false;
}

bool VSH5LibraryRangedPath::isValidForWriting(const VSH5LibraryPath* val) const
{
#if 0
  const VSH5LibraryRangedPath* cval = val->castToRangedPath();
  if((val->name() == name())
     &&((val->castToDefaultPath())
	||((cval)&&(cval->rangedValueLo() < cval->rangedValueHi()))))
     return true;
  else return false;
#endif
  return isValidForReading(val);
}

VSH5LibraryRangedPath* VSH5LibraryRangedPath::
readPath(const std::string& _name, VSOctaveH5ReaderStruct* s)
{
  return new VSH5LibraryRangedPath(_name);
}

VSH5LibraryPath* VSH5LibraryRangedPath::copy() const
{
  return new VSH5LibraryRangedPath(*this);
}

VSH5LibraryRangedPath::Writer::
Writer(VSOctaveH5WriterStruct* s, const VSH5LibraryRangedPath* path,
       const VSH5LibraryPathSet *remainder):
  VSH5LibraryPath::Writer(),
  fPath(path), fRemainder(remainder), fBaseStruct(s),
  fDataECellVector(), fKeyList(), fData(), fDefault()
{
  fDataECellVector = fBaseStruct->writeExpandableCellVector("data");
  fKeyList = fBaseStruct->writeCompositeExpandableVector<key_type>("keys");
}

VSH5LibraryRangedPath::Writer::~Writer()
{
  delete fDefault;
  for(std::map<key_type,VSH5LibraryPath::Writer*>::iterator 
	iwriter=fData.begin(); iwriter!=fData.end(); iwriter++)
    delete iwriter->second;
  delete fDataECellVector;
  delete fBaseStruct;
  delete fRemainder;
  delete fPath;
}

VSOctaveH5WriterStruct* 
VSH5LibraryRangedPath::Writer::write(const VSH5LibraryPathSet* pathset)
{
  VSOctaveH5WriterStruct* return_struct = 0;

  vsassert(pathset);
  const VSH5LibraryPath* path = (*pathset)[0];  
  vsassert(fPath->isValidForWriting(path));
  VSH5LibraryPathSet* sub_pathset = pathset->shift();

  if(fRemainder)
    {
      vsassert(sub_pathset);
      VSH5LibraryPath::Writer* sub_writer = 0;
      if(path->castToDefaultPath())
	{
	  if(fDefault == 0)
	    {
	      VSH5LibraryPathSet* sub_remainder = fRemainder->shift();
	      VSOctaveH5WriterStruct* s = fBaseStruct->writeStruct("default");
	      fDefault = (*fRemainder)[0]->getWriter(s, sub_remainder);
	    }

	  sub_writer = fDefault;
	}
      else
	{
	  const VSH5LibraryRangedPath* cpath = path->castToRangedPath();
	  vsassert(cpath);
	  key_type key(cpath->rangedValueLo(),cpath->rangedValueHi());
	  if(fData.find(key) == fData.end())
	    {
	      VSH5LibraryPathSet* sub_remainder = fRemainder->shift();
	      fKeyList->append(key);
	      VSOctaveH5WriterStruct* s = fDataECellVector->appendStruct();
	      fData[key] = (*fRemainder)[0]->getWriter(s, sub_remainder);
	    }

	  sub_writer = fData[key];
	}

      return_struct = sub_writer->write(sub_pathset);
    }
  else
    { 
      vsassert(sub_pathset==0);
      if(path->castToDefaultPath())
	{
	  return_struct = fBaseStruct->writeStruct("default");
	}
      else
	{
	  const VSH5LibraryRangedPath* cpath = path->castToRangedPath();
	  vsassert(cpath);
	  key_type key(cpath->rangedValueLo(),cpath->rangedValueHi());
	  fKeyList->append(key);
	  return_struct = fDataECellVector->appendStruct();
	}
    }

  delete pathset;
  return return_struct;
}

VSH5LibraryPath::Writer* 
VSH5LibraryRangedPath::getWriter(VSOctaveH5WriterStruct* s,
				const VSH5LibraryPathSet* remainder) const
{
  return new Writer(s, new VSH5LibraryRangedPath(*this), remainder);
}

VSOctaveH5ReaderStruct* VSH5LibraryRangedPath::
read(VSOctaveH5ReaderStruct* s, const VSH5LibraryPath* path) const
{
  vsassert(path);
  vsassert(this->isValidForReading(path));
  
  const VSH5LibraryRangedPath* cpath = path->castToRangedPath();
  if(cpath)
    {
      std::vector<Writer::key_type> keys;
      vsassert(s->readCompositeVector("keys",keys));
      for(unsigned ikey=0;ikey<keys.size();ikey++)
	{
	  if((cpath->rangedValueLo() >= keys[ikey].first)
	     &&((cpath->rangedValueHi() < keys[ikey].second)
		||((keys[ikey].first == keys[ikey].second)
		   &&(cpath->rangedValueHi() == keys[ikey].second))))
	    {
	      VSOctaveH5ReaderCellVector* cdata = s->readCellVector("data");
	      vsassert(cdata);
	      if(cdata->isStruct(ikey))return cdata->readStruct(ikey);
	    }
	}
    }

  // If we get here then we need to check for the default since we
  // were either asked for it explicitly (cpath==0) or the value we
  // were asked for is not present

  if(s->isStruct("default"))return s->readStruct("default");
  
  return 0;
}

VSH5LibraryPathSet* VSH5LibraryRangedPath::
list(VSOctaveH5ReaderStruct* s) const 
{
  std::vector<Writer::key_type> keys;
  vsassert(s->readCompositeVector("keys",keys));
  bool has_default = s->isStruct("default");
  VSH5LibraryPathSet* lps = 
    new VSH5LibraryPathSet(keys.size()+(has_default?1:0));
  for(unsigned ikey = 0;ikey<keys.size();ikey++)
    lps->set(ikey, new VSH5LibraryRangedPath(name(),keys[ikey].first,
					     keys[ikey].second));
  if(has_default)
    lps->set(keys.size(),new VSH5LibraryDefaultPath(name()));
  return lps;
}

// ============================================================================
//
// VSH5LibraryPathSet
//
// ============================================================================

VSH5LibraryPathSet::VSH5LibraryPathSet(const VSH5LibraryPathSet& o):
  fPaths(o.fPaths.size()) 
{
  for(unsigned ipath=0;ipath<fPaths.size();ipath++)
    fPaths[ipath] = o.fPaths[ipath]->copy();
}

VSH5LibraryPathSet& VSH5LibraryPathSet::operator=(const VSH5LibraryPathSet& o)
{
  for(unsigned ipath=0;ipath<fPaths.size();ipath++)delete fPaths[ipath];
  fPaths.resize(o.fPaths.size());
  for(unsigned ipath=0;ipath<fPaths.size();ipath++)
    fPaths[ipath] = o.fPaths[ipath]->copy();
  return *this;
}

VSH5LibraryPathSet::~VSH5LibraryPathSet()
{
  for(unsigned ipath=0;ipath<fPaths.size();ipath++)delete fPaths[ipath];
}

VSH5LibraryPathSet* VSH5LibraryPathSet::copy() const
{
  VSH5LibraryPathSet* path_set(new VSH5LibraryPathSet(fPaths.size()));
  for(unsigned ipath=0;ipath<fPaths.size();ipath++)
    path_set->fPaths[ipath] = fPaths[ipath]->copy();
  return path_set;
}

VSH5LibraryPathSet* VSH5LibraryPathSet::shift() const
{  
  if(size()==1)return 0;
  VSH5LibraryPathSet* path_set(new VSH5LibraryPathSet(fPaths.size()-1));
  for(unsigned ipath=1;ipath<fPaths.size();ipath++)
    path_set->fPaths[ipath-1] = fPaths[ipath]->copy();
  return path_set;
}

bool VSH5LibraryPathSet::validate() const
{
  for(unsigned ipath=0;ipath<fPaths.size();ipath++)
    if(fPaths[ipath] == 0)return false;
  return true;
}

void VSH5LibraryPathSet::resize(unsigned path_count)
{
  vsassert(path_count);
  if(path_count<fPaths.size())
    for(unsigned ipath=path_count;ipath<fPaths.size();ipath++)
      delete fPaths[ipath];
  fPaths.resize(path_count);
}

bool VSH5LibraryPathSet::isValidForDefinition() const
{
  for(unsigned ipath=0;ipath<fPaths.size();ipath++)
    if((fPaths[ipath] == 0)||(fPaths[ipath]->isValidForDefinition() == false))
      return false;
  return true;
}

void VSH5LibraryPathSet::
writePathSet(VSOctaveH5WriterStruct* s, const std::string& name) const
{
  VSOctaveH5WriterCellVector* c = s->writeCellVector(name,size());
  for(unsigned ipath=0;ipath<size();ipath++)
    {
      VSOctaveH5WriterStruct* cs = c->writeStruct(ipath);
      fPaths[ipath]->writePath(cs);
      delete cs;
    }
  delete c;
}

VSH5LibraryPathSet* VSH5LibraryPathSet::
readPathSet(VSOctaveH5ReaderStruct* s, const std::string& name)
{
  VSOctaveH5ReaderCellVector* c = s->readCellVector(name);
  vsassert(c);
  unsigned npath = c->dimensions();
  VSH5LibraryPathSet* pathset = new VSH5LibraryPathSet(npath);
  for(unsigned ipath=0;ipath<npath;ipath++)
    {
      VSOctaveH5ReaderStruct* cs = c->readStruct(ipath);
      pathset->set(ipath, VSH5LibraryPath::readPath(cs));
      delete cs;
    }
  delete c;
  return pathset;
}

// ============================================================================
//
// VSH5LibraryWriter
//
// ============================================================================

VSH5LibraryWriter::
VSH5LibraryWriter(VSOctaveH5WriterStruct* base,
		  const VSH5LibraryPathSet& paths):
  fBase(base), fPaths(), fWriter()
{
  vsassert(paths.isValidForDefinition());
  fPaths = paths.copy();
  fPaths->writePathSet(fBase,"path_spec");

  VSH5LibraryPathSet* remainder = fPaths->shift();
  fWriter = (*fPaths)[0]->getWriter(fBase->writeStruct("library"), remainder);
}

VSH5LibraryWriter::~VSH5LibraryWriter()
{
  delete fWriter;
  delete fPaths;
}

VSOctaveH5WriterStruct* 
VSH5LibraryWriter::write(const VSH5LibraryPathSet& pathset)
{
  return fWriter->write(pathset.copy());
}

// ============================================================================
//
// VSH5LibraryReader
//
// ============================================================================

VSH5LibraryReader::
VSH5LibraryReader(VSOctaveH5ReaderStruct* base):
  fBase(base), fPaths()
{
  fPaths = VSH5LibraryPathSet::readPathSet(fBase,"path_spec");
  vsassert(!fPaths->empty());
}

VSH5LibraryReader::~VSH5LibraryReader()
{
  delete fPaths;
}

VSOctaveH5ReaderStruct* 
VSH5LibraryReader::read(const VSH5LibraryPathSet& pathset) const
{
  vsassert(pathset.size() == fPaths->size());
  VSOctaveH5ReaderStruct* s = fBase->readStruct("library");
  for(unsigned ipath=0;ipath<fPaths->size();ipath++)
    {
      vsassert(pathset[ipath]);
      s = (*fPaths)[ipath]->read(s, pathset[ipath]);
      if(s==0)return 0;
    }
  return s;
}

VSH5LibraryPathSet* 
VSH5LibraryReader::list() const
{
  VSOctaveH5ReaderStruct* s = fBase->readStruct("library");
  return (*fPaths)[0]->list(s);
}

VSH5LibraryPathSet* 
VSH5LibraryReader::list(const VSH5LibraryPathSet& pathset) const
{
  vsassert(pathset.size() < fPaths->size());
  VSOctaveH5ReaderStruct* s = fBase->readStruct("library");
  for(unsigned ipath=0;ipath<pathset.size();ipath++)
    {
      vsassert(pathset[ipath]);
      s = (*fPaths)[ipath]->read(s, pathset[ipath]);
      if(s==0)return 0;
    }
  return (*fPaths)[pathset.size()]->list(s);
}

// ****************************************************************************
// TEST MAIN FUNCTION
// ****************************************************************************

#ifdef TEST_MAIN_LIBRARY_WRITER

// g++ -I. -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -g -Wall -DTEST_MAIN_LIBRARY_WRITER -o test VSH5Library.cpp VSOctaveH5Writer.o VSOctaveH5Reader.o VSOctaveIO.o VSFileUtility.o -lhdf5

#include<cmath>
#include<iostream>
#include<VSTime.hpp>

int main(int argc, char** argv)
{
  try
    {
      VSOctaveH5Writer* writer = 
	new VSOctaveH5Writer("test.h5", true, "# produced by test code");

      VSH5LibraryPathSet pathset(3);
      pathset.set(0,new VSH5LibraryNamedPath("dataset"));
      pathset.set(1,new VSH5LibraryRangedPath("zenith"));
      pathset.set(2,new VSH5LibraryIndexedPath("scope",4));
      
      VSH5LibraryWriter library(writer, pathset);
      unsigned n=0;
      const char* sets[] = { "mscw", "mscw_sigma", "mscl", "mscl_sigma" };
      for(unsigned iset=0; iset<sizeof(sets)/sizeof(*sets); iset++)
	{
	  for(double costheta=0;costheta<0.5;costheta+=0.1)
	    {
	      double lo = acos(1.0-costheta)/M_PI*180;
	      double hi = acos(1.0-costheta-0.1)/M_PI*180;
	      for(unsigned iscope=0;iscope<4;iscope++)
		{
		  VSH5LibraryPathSet pathset(3);
		  pathset.set(0,new VSH5LibraryNamedPath("dataset",sets[iset]));
		  pathset.set(1,new VSH5LibraryRangedPath("zenith",lo,hi));
		  pathset.set(2,new VSH5LibraryIndexedPath("scope",iscope));
		  if(n%7 != 0)library.write(pathset)->writeScalar("n",n);
		  n++;
		}
	    }
	  
	  VSH5LibraryPathSet pathset(3);
	  pathset.set(0,new VSH5LibraryNamedPath("dataset",sets[iset]));
	  pathset.set(1,new VSH5LibraryDefaultPath("zenith"));
	  pathset.set(2,new VSH5LibraryDefaultPath("scope"));

	  library.write(pathset)->writeScalar("iset",iset);
	}
    }
  catch(const VSOctaveH5Exception& e)
    {
      std::cerr << e.message() << std::endl;
    }

  H5close();
  return 0;
}

#endif
