//-*-mode:c++; mode:font-lock;-*-

/*! \file VSScaledParameterLibrary.hpp

  Library of scaled parameters (MSCW, MSCL etc..)

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       08/23/2007
*/

#ifndef VSSCALEDPARAMETERLIBRARY_HPP
#define VSSCALEDPARAMETERLIBRARY_HPP

#include<cmath>

#include<VSH5Library.hpp>
#include<VSNSpace.hpp>

namespace VERITAS
{

  class VSScaledParameterLibraryWriter
  {
  public:
    VSScaledParameterLibraryWriter(VSOctaveH5WriterStruct* base);
    ~VSScaledParameterLibraryWriter();

    enum ScaledParameterSet { SPS_WIDTH_EXPECTED,  SPS_WIDTH_RMS,
			      SPS_LENGTH_EXPECTED, SPS_LENGTH_RMS,
			      SPS_DISP_EXPECTED,   SPS_DISP_RMS,
			      SPS_ENERGY_EXPECTED, SPS_ENERGY_RMS,
			      SPS_KERNEL_SIGMA1,   SPS_KERNEL_BIAS1,   
			      SPS_KERNEL_SIGMA2,   SPS_KERNEL_BIAS2,   
			      SPS_KERNEL_ALPHA,    SPS_EFFECTIVE_AREA,
			      SPS_PSF_SIGMA1,      SPS_PSF_SIGMA2,
			      SPS_PSF_ALPHA };

    VSOctaveH5WriterStruct* 
    writePath(ScaledParameterSet sps,
	      bool default_zn, double zn_lo, double zn_hi,
	      bool default_az, double az_lo, double az_hi,
	      bool default_scope, unsigned iscope,
	      bool default_ped_rms, double ped_rms_lo, double ped_rms_hi);
    
    bool write(const VSNSpace& data,
	       const VSNSpace& mask,
	       ScaledParameterSet sps,
	       bool default_zn = true, double zn_lo = 0, double zn_hi = 0,
	       bool default_az = true, double az_lo = 0, double az_hi = 0,
	       bool default_scope = true, unsigned iscope = 0,
	       bool default_ped_rms = true, 
	       double ped_rms_lo = 0, double ped_rms_hi = 0);

    bool write(const VSNSpace& data,
	       ScaledParameterSet sps,
	       bool default_zn = true, double zn_lo = 0, double zn_hi = 0,
	       bool default_az = true, double az_lo = 0, double az_hi = 0,
	       bool default_scope = true, unsigned iscope = 0,
	       bool default_ped_rms = true, 
	       double ped_rms_lo = 0, double ped_rms_hi = 0);

    static VSNSpace::Space defaultSpace(unsigned ndim = defaultNDim());
    static VSH5LibraryPathSet defaultPathSet();
    static std::string spsName(ScaledParameterSet sps);
    static ScaledParameterSet spsEnum(const std::string& s);

    static void point(VSNSpace::Point& p, double Ri, double Ni, double disp,
		      unsigned ndim = defaultNDim())
    {
      p.resize(ndim);
      p.x[0] = Ri;
      p.x[1] = log10(Ni);
      if(ndim>2)p.x[2] = disp;
    }

    static unsigned defaultNDim() { return s_default_ndim; }
    static unsigned setDefaultNDim(unsigned d) 
    { unsigned od=s_default_ndim; s_default_ndim=d; return od; }

  private:
    VSScaledParameterLibraryWriter(const VSScaledParameterLibraryWriter&);
    VSScaledParameterLibraryWriter& 
    operator=(const VSScaledParameterLibraryWriter&);

    VSH5LibraryWriter m_writer;

    static unsigned   s_default_ndim;
  };

  class VSScaledParameterLibraryReader
  {
  public:
    VSScaledParameterLibraryReader(VSOctaveH5ReaderStruct* base);
    ~VSScaledParameterLibraryReader();
    
    typedef VSScaledParameterLibraryWriter::ScaledParameterSet 
    ScaledParameterSet;

    VSOctaveH5ReaderStruct* 
    readPath(ScaledParameterSet sps,
	     bool force_default_zn, double zn,
	     bool force_default_az, double az,
	     bool force_default_scope, unsigned iscope,
	     bool force_default_ped_rms, double ped_rms);

    VSNSpace* loadData(VSOctaveH5ReaderStruct* s);
    VSNSpace* loadMask(VSOctaveH5ReaderStruct* s);

    VSNSpace* read(ScaledParameterSet sps,
		   bool force_default_zn = true, double zn = 0,
		   bool force_default_az = true, double az = 0,
		   bool force_default_scope = true, unsigned iscope = 0,
		   bool force_default_ped_rms = true, double ped_rms = 0);

    VSNSpace* readMask(ScaledParameterSet sps,
		       bool force_default_zn = true, double zn = 0,
		       bool force_default_az = true, double az = 0,
		       bool force_default_scope = true, unsigned iscope = 0,
		       bool force_default_ped_rms = true, double ped_rms = 0);

    VSNSpace* readAndInterpolate(ScaledParameterSet sps,
				 double zn, double az, 
				 unsigned iscope, double ped_rms);

    void readAndInterpolate(VSNSpace& data, VSNSpace& mask,
			    ScaledParameterSet sps, 
			    double zn, double az, unsigned iscope, 
			    double ped_rms);

    VSNSpace* read(ScaledParameterSet sps,
		   double zn, double az, unsigned iscope, double ped_rms)
    {
      return read(sps,false,zn,false,az,false,iscope,false,ped_rms);
    }

    void copyTables(VSScaledParameterLibraryWriter& writer,
		    double zn, double az, double ped_rms);

    static VSNSpace::Space defaultSpace() 
    { return VSScaledParameterLibraryWriter::defaultSpace(); }
    static VSH5LibraryPathSet defaultPathSet()
    { return VSScaledParameterLibraryWriter::defaultPathSet(); }
    static std::string spsName(ScaledParameterSet sps)
    { return VSScaledParameterLibraryWriter::spsName(sps); }
    static ScaledParameterSet spsEnum(const std::string& s)
    { return VSScaledParameterLibraryWriter::spsEnum(s); }

  private:
    void getRangedPaths(const VSH5LibraryPathSet& base_path_set,
			std::vector<VSH5LibraryPathSet>& paths, 
			double par);

    void interpolate(std::vector< VSNSpace* >& data,
		     std::vector<VSH5LibraryPathSet>::iterator path_itr,
		     double par);
    

    VSScaledParameterLibraryReader(const VSScaledParameterLibraryReader&);
    VSScaledParameterLibraryReader& 
    operator=(const VSScaledParameterLibraryReader&);

    VSH5LibraryReader m_reader;      
  };

}

#endif // defined VSSCALEDPARAMETERLIBRARY_HPP
