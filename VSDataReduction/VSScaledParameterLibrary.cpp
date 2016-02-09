//-*-mode:c++; mode:font-lock;-*-

/*! \file VSScaledParameterLibrary.cpp

  Library of scaled parameters (MSCW, MSCL etc..)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/23/2007
*/

#include<vsassert>

#include<VSScaledParameterLibrary.hpp>
#include<VSNSpaceOctaveH5IO.hpp>

using namespace VERITAS;

// ============================================================================
//
// WRITER
//
// ============================================================================

unsigned VSScaledParameterLibraryWriter::s_default_ndim = 2;

VSScaledParameterLibraryWriter::
VSScaledParameterLibraryWriter(VSOctaveH5WriterStruct* base):
  m_writer(base,defaultPathSet())
{
  // nothing to see here
}

VSScaledParameterLibraryWriter::~VSScaledParameterLibraryWriter()
{
  // nothing to see here
}

VSOctaveH5WriterStruct* 
VSScaledParameterLibraryWriter::
writePath(ScaledParameterSet sps,
	  bool default_zn, double zn_lo, double zn_hi,
	  bool default_az, double az_lo, double az_hi,
	  bool default_scope, unsigned iscope,
	  bool default_ped_rms, double ped_rms_lo, double ped_rms_hi)
{
  VSH5LibraryPathSet def_path_set = defaultPathSet();
  VSH5LibraryPathSet path_set(def_path_set.size());

  VSH5LibraryPath* p0 =
    new VSH5LibraryNamedPath(def_path_set[0]->name(), spsName(sps));

  VSH5LibraryPath* p1 = 0;
  if(default_zn)
    p1 = new VSH5LibraryDefaultPath(def_path_set[1]->name());
  else
    p1 = new VSH5LibraryRangedPath(def_path_set[1]->name(), zn_lo, zn_hi);

  VSH5LibraryPath* p2 = 0;
  if(default_az)
    p2 = new VSH5LibraryDefaultPath(def_path_set[2]->name());
  else
    p2 = new VSH5LibraryRangedPath(def_path_set[2]->name(), az_lo, az_hi);

  VSH5LibraryPath* p3 = 0;
  if(default_scope)
    p3 = new VSH5LibraryDefaultPath(def_path_set[3]->name());
  else
    p3 = new VSH5LibraryIndexedPath(def_path_set[3]->name(), iscope);

  VSH5LibraryPath* p4 = 0;
  if(default_ped_rms)
    p4 = new VSH5LibraryDefaultPath(def_path_set[4]->name());
  else
    p4 = new VSH5LibraryRangedPath(def_path_set[4]->name(), 
				   ped_rms_lo, ped_rms_hi);

  path_set.set(0,p0);
  path_set.set(1,p1);
  path_set.set(2,p2);
  path_set.set(3,p3);
  path_set.set(4,p4);

  VSOctaveH5WriterStruct* s = m_writer.write(path_set);
  vsassert(s);

  return s;
}

bool VSScaledParameterLibraryWriter::
write(const VSNSpace& data,
      const VSNSpace& mask,
      ScaledParameterSet sps,
      bool default_zn, double zn_lo, double zn_hi,
      bool default_az, double az_lo, double az_hi,
      bool default_scope, unsigned iscope,
      bool default_ped_rms, double ped_rms_lo, double ped_rms_hi) 
{
  VSOctaveH5WriterStruct* s = 
    writePath(sps,default_zn,zn_lo,zn_hi,
	      default_az,az_lo,az_hi,
	      default_scope,iscope,default_ped_rms,ped_rms_lo,ped_rms_hi);

  data.save(s->writeStruct("data"));
  mask.save(s->writeStruct("mask"));

  return true;
}

bool VSScaledParameterLibraryWriter::
write(const VSNSpace& data,
      ScaledParameterSet sps,
      bool default_zn, double zn_lo, double zn_hi,
      bool default_az, double az_lo, double az_hi,
      bool default_scope, unsigned iscope,
      bool default_ped_rms, double ped_rms_lo, double ped_rms_hi) 
{
  VSOctaveH5WriterStruct* s = 
    writePath(sps,default_zn,zn_lo,zn_hi,
	      default_az,az_lo,az_hi,
	      default_scope,iscope,default_ped_rms,ped_rms_lo,ped_rms_hi);
  

  VSNSpace mask = data;
  mask.clear(1.0);

  data.save(s->writeStruct("data"));
  mask.save(s->writeStruct("mask"));
  return true;
  // VSNSpaceOctaveH5IO io;
  // return io.writeHistogram(s, data);
}

VSNSpace::Space VSScaledParameterLibraryWriter::defaultSpace(unsigned ndim)
{
  VSNSpace::Space space(ndim);
  if(ndim==2)
    {
space.axes[0] = VSNSpace::Axis(-10.0, 800.0, 20.0, 0, "Impact distance [m]");
space.axes[1] = VSNSpace::Axis(  1.0,   6.0,  0.1, 0, "ADC signal [log10 DC]");
    }
  else if(ndim==3)
    {
space.axes[0] = VSNSpace::Axis(-20.0, 800.0, 40.0, 0, "Impact distance [m]");
space.axes[1] = VSNSpace::Axis(  1.0,   6.0,  0.2, 0, "ADC signal [log10 DC]");
space.axes[2] = VSNSpace::Axis( -0.1,   3.6,  0.2, 0, "FP displacement [deg]");
    }
  else
    {
      vsassert(0);
    }
  return space;
}

VSH5LibraryPathSet VSScaledParameterLibraryWriter::defaultPathSet()
{
  VSH5LibraryPathSet path_set(5);
  path_set.set(0,new VSH5LibraryNamedPath("Parameter set"));
  path_set.set(1,new VSH5LibraryRangedPath("Zenith angle [deg]"));
  path_set.set(2,new VSH5LibraryRangedPath("Azimuth angle (E from N) [deg]"));
  path_set.set(3,new VSH5LibraryIndexedPath("Telescope number",16));
  path_set.set(4,new VSH5LibraryRangedPath("Ped RMS"));
  return path_set;
}

std::string VSScaledParameterLibraryWriter::spsName(ScaledParameterSet sps)
{
  switch(sps)
    {
    case SPS_WIDTH_EXPECTED:  return "width_expected";
    case SPS_WIDTH_RMS:       return "width_rms";
    case SPS_LENGTH_EXPECTED: return "length_expected";
    case SPS_LENGTH_RMS:      return "length_rms";
    case SPS_DISP_EXPECTED:   return "disp_expected";
    case SPS_DISP_RMS:        return "disp_rms";
    case SPS_ENERGY_EXPECTED: return "energy_expected";
    case SPS_ENERGY_RMS:      return "energy_rms";
    case SPS_KERNEL_SIGMA1:   return "energy_kernel_sigma1";
    case SPS_KERNEL_BIAS1:    return "energy_kernel_bias1";
    case SPS_KERNEL_SIGMA2:   return "energy_kernel_sigma2";
    case SPS_KERNEL_BIAS2:    return "energy_kernel_bias2";
    case SPS_KERNEL_ALPHA:    return "energy_kernel_alpha";
    case SPS_EFFECTIVE_AREA:  return "effective_area";
    case SPS_PSF_SIGMA1:      return "psf_sigma1";
    case SPS_PSF_SIGMA2:      return "psf_sigma2";
    case SPS_PSF_ALPHA:       return "psf_alpha";
    }
  vsassert(0);
}

#define MATCHENUM(x) if(s==spsName(x))return x

VSScaledParameterLibraryWriter::ScaledParameterSet 
VSScaledParameterLibraryWriter::spsEnum(const std::string& s)
{
  MATCHENUM(SPS_WIDTH_EXPECTED);
  MATCHENUM(SPS_WIDTH_RMS);
  MATCHENUM(SPS_LENGTH_EXPECTED);
  MATCHENUM(SPS_LENGTH_RMS);
  MATCHENUM(SPS_DISP_EXPECTED);
  MATCHENUM(SPS_DISP_RMS);
  MATCHENUM(SPS_ENERGY_EXPECTED);
  MATCHENUM(SPS_ENERGY_RMS);
  MATCHENUM(SPS_KERNEL_SIGMA1);
  MATCHENUM(SPS_KERNEL_BIAS1);
  MATCHENUM(SPS_KERNEL_SIGMA2);
  MATCHENUM(SPS_KERNEL_BIAS2);
  MATCHENUM(SPS_KERNEL_ALPHA);
  MATCHENUM(SPS_EFFECTIVE_AREA);
  vsassert(0);
}

// ============================================================================
//
// READER
//
// ============================================================================

VSScaledParameterLibraryReader::
VSScaledParameterLibraryReader(VSOctaveH5ReaderStruct* base):
  m_reader(base)
{
  // nothing to see here
}

VSScaledParameterLibraryReader::~VSScaledParameterLibraryReader()
{
  // nothing to see here
}

VSOctaveH5ReaderStruct* VSScaledParameterLibraryReader::
readPath(ScaledParameterSet sps,
	 bool force_default_zn, double zn,
	 bool force_default_az, double az,
	 bool force_default_scope, unsigned iscope,
	 bool force_default_ped_rms, double ped_rms)
{
  VSH5LibraryPathSet def_path_set = defaultPathSet();
  VSH5LibraryPathSet path_set(def_path_set.size());

  VSH5LibraryPath* p0 =
    new VSH5LibraryNamedPath(def_path_set[0]->name(), spsName(sps));

  VSH5LibraryPath* p1 = 0;
  if(force_default_zn)
    p1 = new VSH5LibraryDefaultPath(def_path_set[1]->name());
  else
    p1 = new VSH5LibraryRangedPath(def_path_set[1]->name(), zn);

  VSH5LibraryPath* p2 = 0;
  if(force_default_az)
    p2 = new VSH5LibraryDefaultPath(def_path_set[2]->name());
  else
    p2 = new VSH5LibraryRangedPath(def_path_set[2]->name(), az);

  VSH5LibraryPath* p3 = 0;
  if(force_default_scope)
    p3 = new VSH5LibraryDefaultPath(def_path_set[3]->name());
  else
    p3 = new VSH5LibraryIndexedPath(def_path_set[3]->name(), iscope);

  VSH5LibraryPath* p4 = 0;
  if(force_default_ped_rms)
    p4 = new VSH5LibraryDefaultPath(def_path_set[4]->name());
  else
    p4 = new VSH5LibraryRangedPath(def_path_set[4]->name(), ped_rms);

  path_set.set(0,p0);
  path_set.set(1,p1);
  path_set.set(2,p2);
  path_set.set(3,p3);
  path_set.set(4,p4);

  VSOctaveH5ReaderStruct* s = m_reader.read(path_set);
  if(s==0) return 0; else return s;
}

VSNSpace* VSScaledParameterLibraryReader::
read(ScaledParameterSet sps,
     bool force_default_zn, double zn,
     bool force_default_az, double az,
     bool force_default_scope, unsigned iscope,
     bool force_default_ped_rms, double ped_rms)
{
  VSOctaveH5ReaderStruct* s = 
    readPath(sps, force_default_zn, zn, force_default_az, az,
	     force_default_scope, iscope, force_default_ped_rms, ped_rms);

  return loadData(s);
}

VSNSpace* VSScaledParameterLibraryReader::loadData(VSOctaveH5ReaderStruct* s)
{
  if(s == 0) return 0;
  VSNSpace* data = new VSNSpace;

  if(VSNSpaceOctaveH5IO::isHistogram(s))
    {
      if(data->load(s)) return data;
    }
  else if(s->isStruct("data"))
    {
      VSOctaveH5ReaderStruct* s_data = s->readStruct("data");
      if(data->load(s_data)) return data;
    }

  delete data;
  return 0;
}

VSNSpace* VSScaledParameterLibraryReader::
readMask(ScaledParameterSet sps,
     bool force_default_zn, double zn,
     bool force_default_az, double az,
     bool force_default_scope, unsigned iscope,
     bool force_default_ped_rms, double ped_rms)
{
  VSOctaveH5ReaderStruct* s = 
    readPath(sps, force_default_zn, zn, force_default_az, az,
	     force_default_scope, iscope, force_default_ped_rms, ped_rms);

  return loadMask(s);
}

VSNSpace* VSScaledParameterLibraryReader::loadMask(VSOctaveH5ReaderStruct* s)
{
  if(s == 0) return 0;
  VSNSpace* mask = new VSNSpace;

  if(s->isStruct("mask"))
    {
      VSOctaveH5ReaderStruct* s_mask = s->readStruct("mask");
      if(mask->load(s_mask)) return mask;
    }
  else if(VSNSpaceOctaveH5IO::isHistogram(s))
    {
      if(mask->load(s)) 
	{
	  mask->clear(1.0);
	  mask->setComment("");
	  return mask;
	}
    }

  delete mask;
  return 0; 
}

void VSScaledParameterLibraryReader::
getRangedPaths(const VSH5LibraryPathSet& base_path_set,
	       std::vector<VSH5LibraryPathSet>& paths, double par)
{
  VSH5LibraryPathSet def_path_set = defaultPathSet();

  VSH5LibraryPathSet* path_set = m_reader.list(base_path_set);
  vsassert(path_set);
  std::vector<double> par_list;

  for(unsigned iset=0;iset<path_set->size();iset++)
    {
      if((*path_set)[iset]->castToRangedPath())
	par_list.push_back((*path_set)[iset]->castToRangedPath()->
			   rangedValueLo());
    }
  delete path_set;
  
  
  VSH5LibraryPathSet p1(base_path_set);
  VSH5LibraryPathSet p2(base_path_set);

  const unsigned ipar = base_path_set.size();
  const std::string par_name = def_path_set[ipar]->name();

  p1.resize(base_path_set.size()+1);
  p2.resize(base_path_set.size()+1);

  if(par_list.empty())
    {
      p1.set(ipar,new VSH5LibraryDefaultPath(par_name));
      p2.set(ipar,new VSH5LibraryDefaultPath(par_name));
    }
  else if(par_list.size()==1)
    {
      double par1 = par_list.front();
      p1.set(ipar,new VSH5LibraryRangedPath(par_name, par1));
      p2.set(ipar,new VSH5LibraryRangedPath(par_name, par1));
    }
  else
    {
      std::sort(par_list.begin(),par_list.end());
      std::vector<double>::const_iterator itr = par_list.begin();
      while(itr != par_list.end()) 
	{
	  if(*itr > par) break;
	  else itr++;
	}      

      double par1, par2;

      if((itr == par_list.begin() || itr == par_list.end()) && ipar == 2)
	{
	  par1 = par_list.back();
	  par2 = par_list.front();
	}
      else if(itr == par_list.begin())
	{
	  par1 = *itr;
	  par2 = *(itr+1);
	}
      else if(itr == par_list.end())
	{
	  par1 = *(itr-2);
	  par2 = *(itr-1);
	}
      else
	{
	  par1 = *(itr-1);
	  par2 = *itr;
	}

      p1.set(ipar,new VSH5LibraryRangedPath(par_name, par1));
      p2.set(ipar,new VSH5LibraryRangedPath(par_name, par2));
    }

  paths.push_back(p1);
  paths.push_back(p2);
}

void VSScaledParameterLibraryReader::
interpolate(std::vector< VSNSpace* >& data,
	    std::vector<VSH5LibraryPathSet>::iterator path_itr, double par)
{
  for(std::vector< VSNSpace* >::iterator data_itr = data.begin();
      data_itr != data.end(); )
    {
      const VSH5LibraryPath* path1 = *((*path_itr).end()-1);
      const VSH5LibraryPath* path2 = *((*(path_itr+1)).end()-1);
     
      if(path1->castToRangedPath() && path2->castToRangedPath())
	{
	  double p1 = path1->castToRangedPath()->rangedValueLo();
	  double p2 = path2->castToRangedPath()->rangedValueLo();

	  VSNSpace& d1 = *(*data_itr);
	  VSNSpace& d2 = *(*(data_itr+1));

	  double dx = 1;
	  if(p1 > p2) dx = fmod(p2-par+360,360)/fmod(p2-p1+360,360);
	  else if(p1 != p2) dx = (p2-par)/(p2-p1);

	  d1 *= dx;
	  d2 *= (1-dx);
	  d1 += d2;
	}

      delete *(data_itr+1);
      data_itr = data.erase(data_itr+1,data_itr+2);
      path_itr+=2;
    }
}

void VSScaledParameterLibraryReader::
readAndInterpolate(VSNSpace& data_nspace, VSNSpace& mask_nspace,
		   ScaledParameterSet sps, 
		   double zn, double az, unsigned iscope, double ped_rms)
{
  VSH5LibraryPathSet def_path_set = defaultPathSet();

  VSH5LibraryPath* p0 =
    new VSH5LibraryNamedPath(def_path_set[0]->name(), spsName(sps));

  VSH5LibraryPath* p3 = 
    new VSH5LibraryIndexedPath(def_path_set[3]->name(), iscope);

  VSH5LibraryPathSet path_set(1);
  path_set.set(0,p0);

  std::vector<VSH5LibraryPathSet> zn_paths;
  std::vector<VSH5LibraryPathSet> az_paths;
  std::vector<VSH5LibraryPathSet> ped_paths;

  // Get list of zenith angles in file ----------------------------------------
  getRangedPaths(path_set,zn_paths,zn);

  // Get list of azimuth angles in file ---------------------------------------
  for(std::vector<VSH5LibraryPathSet>::iterator zn_itr = 
	zn_paths.begin(); zn_itr != zn_paths.end(); ++zn_itr)
    getRangedPaths(*zn_itr,az_paths,az);

  // Get list of ped rms in file ----------------------------------------------
  for(std::vector<VSH5LibraryPathSet>::iterator az_itr = 
	az_paths.begin(); az_itr != az_paths.end(); ++az_itr)
    {
      VSH5LibraryPathSet path(*az_itr);

      path.resize(path.size()+1);
      path.set(3,p3->copy());
      getRangedPaths(path,ped_paths,ped_rms);
    }

  std::vector< VSNSpace* > data;
  VSNSpace* mask = 0;

  for(std::vector<VSH5LibraryPathSet>::iterator path_itr = 
	ped_paths.begin(); path_itr != ped_paths.end(); ++path_itr)
    {     
      VSOctaveH5ReaderStruct* s = m_reader.read(*path_itr);

      if(s==0)return;

      VSNSpace* d = loadData(s);
      VSNSpace* m = loadMask(s);
      vsassert(d != 0 && m != 0);

      data.push_back(d);

      if(mask == 0) mask = m;
      else
	{
	  *mask *= *m;
	  delete m;
	}
    }

  interpolate(data,ped_paths.begin(),ped_rms);  
  interpolate(data,az_paths.begin(),az);  
  interpolate(data,zn_paths.begin(),zn);  

  data_nspace = *data[0];
  mask_nspace = *mask;

  delete data[0];
  delete mask;
}

VSNSpace* VSScaledParameterLibraryReader::
readAndInterpolate(ScaledParameterSet sps, 
		   double zn, double az, unsigned iscope, double ped_rms)
{
  VSNSpace* data = new VSNSpace;
  VSNSpace* mask = new VSNSpace;

  readAndInterpolate(*data,*mask,sps,zn,az,iscope,ped_rms);
  delete mask;

  return data;
}

#define NUMEL(x) (sizeof(x)/sizeof(*x))

void VSScaledParameterLibraryReader::
copyTables(VSScaledParameterLibraryWriter& writer, double zn, double az,
	   double ped_rms)
{
  const VSH5LibraryPathSet def_path_set = defaultPathSet();

  VSH5LibraryPathSet* p0set = m_reader.list();
  vsassert(p0set);
  
  unsigned nset = p0set->size();
  for(unsigned iset=0;iset<nset;iset++)
    {
      if(!(*p0set)[iset]->castToNamedPath())continue;

      ScaledParameterSet sps = 
	spsEnum((*p0set)[iset]->castToNamedPath()->namedValue());

      VSH5LibraryPathSet path_set(3);

      VSH5LibraryPath* p1 = 
	new VSH5LibraryRangedPath(def_path_set[1]->name(), zn);

      VSH5LibraryPath* p2 = 0;
      if(az < 0)
	p2 = new VSH5LibraryDefaultPath(def_path_set[2]->name());
      else
	p2 = new VSH5LibraryRangedPath(def_path_set[2]->name(), az);

      path_set.set(0,(*p0set)[iset]->copy());
      path_set.set(1,p1);
      path_set.set(2,p2);

      VSH5LibraryPathSet* p3set = m_reader.list(path_set);
      vsassert(p3set);

      path_set.resize(def_path_set.size());
      
      unsigned nscope = p3set->size();
      for(unsigned iscope=0;iscope<nscope;iscope++)
	{
	  path_set.set(3,(*p3set)[iscope]->copy());

	  VSH5LibraryPath* p4 = 0;
	  if(ped_rms < 0)
	    p4 = new VSH5LibraryDefaultPath(def_path_set[4]->name());
	  else
	    p4 = new VSH5LibraryRangedPath(def_path_set[4]->name(), ped_rms);
	  
	  path_set.set(4,p4);

	  VSOctaveH5ReaderStruct* s = m_reader.read(path_set);
	  if(s==0)continue;

	  VSNSpace* data = new VSNSpace;
	  VSNSpaceOctaveH5IO io;
	  if(io.readHistogram(s, *data) == false)
	    {
	      delete data;
	      continue;
	    }

	  bool def = (*p3set)[iscope]->castToDefaultPath();
	  unsigned iiscope = 0;
	  if(!def)
	    {
	      const VSH5LibraryIndexedPath* cpath = 
		(*p3set)[iscope]->castToIndexedPath();
	      vsassert(cpath);
	      iiscope = cpath->indexedValue();
	    }

	  writer.write(*data, sps, false, zn, zn, az<0, az, az, def, iiscope,
		       ped_rms<0, ped_rms,ped_rms);
	  
	  delete data;
	}

      delete p3set;
    }

  delete p0set;
}
