/*! \file VSNSpaceOctaveH5IO.hpp

  NSpace IO to Octave H5 file

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       06/16/2006

*/

#include<VSOctaveIO.hpp>
#include<VSNSpaceOctaveH5IO.hpp>

using namespace VERITAS;

VSNSpaceOctaveH5IO::~VSNSpaceOctaveH5IO()
{
  // nothing to see here
}

bool VSNSpaceOctaveH5IO::
readSpace(const std::string& filename, VSNSpace::Space& space)
{
  VSOctaveH5Reader file(filename);
  return readSpace(&file, space);
}

bool VSNSpaceOctaveH5IO::
writeSpace(const std::string& filename, const VSNSpace::Space& space)
{
  VSOctaveH5Writer file(filename,true,"# written by VSNSpaveOctaveH5IO");
  return writeSpace(&file, space);
}
bool VSNSpaceOctaveH5IO::
readVolume(const std::string& filename, VSNSpace::Volume& vol)
{
  VSOctaveH5Reader file(filename);
  return readVolume(&file, vol);
}

bool VSNSpaceOctaveH5IO::
writeVolume(const std::string& filename, const VSNSpace::Volume& vol,
	    bool sparse)
{
  VSOctaveH5Writer file(filename,true,"# written by VSNSpaveOctaveH5IO");
  return writeVolume(&file,vol);
}

bool VSNSpaceOctaveH5IO::
readOrdering(const std::string& filename, VSNSpace::Ordering& ordering)
{
  VSOctaveH5Reader file(filename);
  VSNSpace::Space space;
  if(!readSpace(&file, space))return false;

  ordering.m_space = space;
  if(!file.readVector("index", ordering.m_index))return false;
  if(!file.readVector("counts_on", ordering.m_counts_on))return false;
  if(!file.readVector("counts_off", ordering.m_counts_off))return false;
  if(!file.readVector("excess", ordering.m_excess))return false;
  if(!file.readVector("Q", ordering.m_Q))return false;
  if(!file.readVector("sigma", ordering.m_sigma))return false;

  return true;
}

bool VSNSpaceOctaveH5IO::
writeOrdering(const std::string& filename, const VSNSpace::Ordering& ordering)
{
  VSOctaveH5Writer file(filename,true,"# written by VSNSpaveOctaveH5IO");
  if(!writeSpace(&file, ordering.m_space))return false;

  if(!file.writeVector("index", ordering.m_index))return false;
  if(!file.writeVector("counts_on",ordering.m_counts_on))return false;
  if(!file.writeVector("counts_off",ordering.m_counts_off))return false;
  if(!file.writeVector("excess",ordering.m_excess))return false;
  if(!file.writeVector("Q",ordering.m_Q))return false;
  if(!file.writeVector("sigma",ordering.m_sigma))return false;

  return true;
}

bool VSNSpaceOctaveH5IO::
readHistogram(const std::string& filename, VSNSpace& hist)
{
  VSOctaveH5Reader file(filename);
  return readHistogram(&file, hist);
}

bool VSNSpaceOctaveH5IO::
writeHistogram(const std::string& filename, const VSNSpace& hist)
{
  VSOctaveH5Writer file(filename,true,"# written by VSNSpaveOctaveH5IO");
  return writeHistogram(&file, hist);
}

bool VSNSpaceOctaveH5IO::
readSpace(VSOctaveH5ReaderStruct* s, VSNSpace::Space& space)
{
  if(!s) return false;    
  if(!isSpace(s)) return false;

  VSOctaveH5ReaderCellVector* c = s->readCellVector("dims");
  if(c==0)return false;
  unsigned nel = c->dimensions();
  space.resize(nel);
  assert(nel == space.ndim);
  for(unsigned idim=0;idim<space.ndim;idim++)
    {
      VSOctaveH5ReaderStruct* s2 = c->readStruct(idim);
      if(s2==0)return false;
      VSNSpace::Coord lo;
      VSNSpace::Coord hi;
      VSNSpace::Index nbin;
      std::string name;
      if((!s2->readScalar("lo_bound",lo))&&(!s2->readScalar("lo_bound{0}",lo)))
	return false;
      if((!s2->readScalar("hi_bound",hi))&&(!s2->readScalar("hi_bound{0}",hi)))
	return false;
      if((!s2->readScalar("nbin",nbin))&&(!s2->readScalar("nbin{0}",nbin)))
	return false;
      if((!s2->readString("name",name))&&(!s2->readString("name{0}",name)))
	return false;
      space.axes[idim] = VSNSpace::Axis(lo, hi, nbin, name);
    }
  return true;
}

bool VSNSpaceOctaveH5IO::isSpace(VSOctaveH5ReaderStruct* s)
{
  if(!s) return false;    
  std::string nspace_type;
  if(!s->readString("nspace_type",nspace_type)) return false;
  if(nspace_type != "space")
    return false;
  else
    return true;
}

bool VSNSpaceOctaveH5IO::
writeSpace(VSOctaveH5WriterStruct* s, const VSNSpace::Space& space)
{
  s->writeString("nspace_type","space");
  VSOctaveH5WriterCellVector* c = s->writeCellVector("dims",space.ndim);
  if(c==0)return false;
  for(unsigned idim=0;idim<space.ndim;idim++)
    {
      VSOctaveH5WriterStruct* s2 = c->writeStruct(idim);
      if(s2==0)return false;
      if(!s2->writeScalar("lo_bound",space.axes[idim].lo_bound))return false;
      if(!s2->writeScalar("hi_bound",space.axes[idim].hi_bound))return false;
      if(!s2->writeScalar("nbin",space.axes[idim].nbin))return false;
      if(!s2->writeString("name",space.axes[idim].name))return false;
    }
  return true;
}

bool VSNSpaceOctaveH5IO::
readVolume(VSOctaveH5ReaderStruct* s, VSNSpace::Volume& vol)
{
  if(!s) return false;    
  if(!isVolume(s)) return false;

  VSNSpace::Space space;
  readSpace(s->readStruct("space"), space);
  bool sparse;
  s->readScalar("sparse", sparse);
  if(sparse)
    {
      std::vector<VSNSpace::Index> index;
      if(!s->readVector("volume", index))return false;
      return(vol.loadSparse(space,index));
    }
  else
    {
      std::vector<bool> volume;
      if(!s->readVector("volume", volume))return false;
      return(vol.load(space,volume));
    }
}

bool VSNSpaceOctaveH5IO::isVolume(VSOctaveH5ReaderStruct* s)
{
  if(!s) return false;    
  std::string nspace_type;
  if(!s->readString("nspace_type",nspace_type)) return false;

  if(nspace_type == "volume")
    return true;
  else
    return false;
}

bool VSNSpaceOctaveH5IO::
writeVolume(VSOctaveH5WriterStruct* s, const VSNSpace::Volume& vol,
	    bool sparse)
{
  s->writeString("nspace_type","volume");
  if(!writeSpace(s->writeStruct("space"), vol.space()))return false;
  if(!s->writeScalar("sparse", sparse))return false;
  if(sparse)
    {
      std::vector<VSNSpace::Index> index;
      vol.sparse(index);
      if(!s->writeVector("volume", index))return false;
    }
  else
    {
      if(!s->writeVector("volume", vol.volume()))return false;
    }
  return true;
}


bool VSNSpaceOctaveH5IO::
readHistogram(VSOctaveH5ReaderStruct* s, VSNSpace& hist)
{
  if(!s) return false;    
  if(!isHistogram(s)) return false;

  std::string comment;
  s->readString("comment", comment);
  VSNSpace::Space space;
  readSpace(s->readStruct("space"), space);
  std::vector<VSNSpace::Weight> weight;
  s->readVector("weight", weight);
  return(hist.load(space, weight, comment));
}

bool VSNSpaceOctaveH5IO::isHistogram(VSOctaveH5ReaderStruct* s)
{
  if(!s) return false;    
  std::string nspace_type;
  if(!s->readString("nspace_type",nspace_type)) return false;
  if(nspace_type != "hist")
    return false;
  else
    return true;
}

bool VSNSpaceOctaveH5IO::
writeHistogram(VSOctaveH5WriterStruct* s, const VSNSpace& hist)
{
  s->writeString("nspace_type","hist");
  s->writeString("comment", hist.comment());
  if(!writeSpace(s->writeStruct("space"), hist.space()))return false;
  if(!s->writeVector("weight", hist.weight()))return false;
  return true;
}
