//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstruction.hpp
  Various reconstruction methods

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/12/2005
*/

#ifndef VSARECONSTRUCTION_HPP
#define VSARECONSTRUCTION_HPP

#define VSA_ASSERT

#include<vector>
#include<memory>
#include<limits>

#include<VSAAlgebra.hpp>
#include<VSAReconstructionCommon.hpp>
#include<VSAReconstructionMethod1.hpp>
#include<VSAReconstructionMethod2.hpp>
#include<VSAReconstructionMethod3.hpp>
#include<VSAAtmosphere.hpp>

namespace VERITAS 
{

#define VSAR_M2_GRID 500.0   /* Radius for Method 2 search               [m] */
#define VSAR_M3_GRID 3.0     /* Radius for Method 3 search         [degrees] */

  class VSAReconstruction
  {
  public:

    typedef VSAReconstructionInternals::ArrayInfo ArrayInfo;
    typedef VSAReconstructionInternals::ScopeInfo ScopeInfo;
    typedef VSAReconstructionInternals::PixelImage PixelImage;
    typedef VSAReconstructionInternals::ScopeImage ScopeImage;
    typedef VSAReconstructionInternals::ArrayMoments ArrayMoments;
    typedef VSAReconstructionInternals::ScopeMoments ScopeMoments;

    enum Method { M_1,
		  M_2,
		  M_3,
		  M_2_GRIDSEARCH,
		  M_3_GRIDSEARCH };

    enum ScopeWeighting { SW_ONE,
			  SW_SIZE,
			  SW_ELLIPTICITY_MEMO,
			  SW_SIZE_ELLIPTICITY_MEMO,
			  SW_ELLIPTICITY_TRAD,
			  SW_WIDTH };

    struct Reconstruction
    {
      Reconstruction(): 
	moments(), 
	e(), R(), chi2e(), chi2R(), deltael(), deltaew(), deltaRl(), deltaRw()
      { /* nothing to see here */ }

      // Multitude of intermediate parameters
      ArrayMoments      moments;  

      // Reconstruction parameters
      VSAAlgebra::Vec3D e;        // Arrival direction
      VSAAlgebra::Vec3D R;        // Impact point
      double            chi2e;    // Chi^2 of e-functional at minimum point
      double            chi2R;    // Chi^2 of R-functional at minimum point
      double            deltael;  // Length of e-error contour;
      double            deltaew;  // Width of e-error contour;  
      double            deltaRl;  // Length of R-error contour;
      double            deltaRw;  // Width of R-error contour;
    };

    struct ScopeParameters
    {
      ScopeParameters():
	Ni(), Ri(), d1i(), d2i(), l1i(), l2i(), theta1i(), theta2i(), 
	delta1i(), delta2i(), D1i(), h1i(), h2i(), G1i(), G2i(), thetac1(),
	t01i(), t02i(), 
	lambdadi(), lambdaci() { /* nothing to see here */ }
      
      // Geometric properties
      double            Ni;       // Total signal size
      double            Ri;       // Impact parameter WRT telescope
      double            d1i;      // Mean distance from impact to emission
      double            d2i;      // RMS distance from impact to emission
      double            l1i;      // Mean distance from scope to emission
      double            l2i;      // RMS distance from scope to emission
      double            theta1i;  // Mean emission angle
      double            theta2i;  // RMS emission angle
      VSAAlgebra::Vec3D delta1i;  // Mean vector from emission to photons
      VSAAlgebra::Symmetric3D delta2i;  // RMS vector from emission to photons
      VSAAlgebra::Vec3D D1i;      // Vector from impact to mean emission point
      double            h1i;      // Mean emission height ASL
      double            h2i;      // Mean emission height ASL
      double            G1i;      // Mean depth in atmosphere at emission point
      double            G2i;      // RMS depth in atmosphere at emission point
      double            thetac1;  // Cherenkov angle at mean depth emission pt

      // Timing parameters
      double            t01i;     // Mean arrival time of photons at impact
      double            t02i;     // RMS arrival time of photons at impact

      // Luminosity parameters
      double            lambdadi; // Density of emmitters in diffusive regime
      double            lambdaci; // Density of emmitters in coherent regime
    };

    struct ArrayParameters
    {
      ArrayParameters(unsigned nscopes=0)
	: N(0), scopes(nscopes)
      { /* nothing to see here */ }

      double                       N;
      std::vector<ScopeParameters> scopes;      
    };

    struct ScopeQualityCuts
    {
      ScopeQualityCuts(unsigned _nimage_min = 0,
		       double _N_min = 0, double _width_min = -1,
		       double _dist2_max = 
		       std::numeric_limits<double>::infinity()):
	nimage_min(_nimage_min), 
	N_min(_N_min), width_min(_width_min), dist2_max(_dist2_max)
      { /* nothing to see here */ }

      unsigned                     nimage_min;
      double                       N_min;
      double                       width_min;
      double                       dist2_max;
    };

    struct ArrayQualityCuts
    {
      ArrayQualityCuts(unsigned nscope = 0, 
		       unsigned _nuse_in_reconstruction_min = 0):
	scopes(nscope),
	nuse_in_reconstruction_min(_nuse_in_reconstruction_min)
      { /* nothing to see here */ }

      std::vector<ScopeQualityCuts> scopes;
      unsigned                      nuse_in_reconstruction_min;
    };

    VSAReconstruction(const ArrayInfo& array_info, 
		      const VSAAtmosphere& atmosphere,
		      const ArrayQualityCuts& cuts = ArrayQualityCuts())
      : fArray(array_info), fAtmo(atmosphere),
	fQualityCuts(cuts),
	//fObsThick(atmosphere.thickness(array_info.elevation/1000.0)),
	fObsDelay(atmosphere.delay(array_info.elevation/1000.0))
    { /* nothing to see here */ }
    ~VSAReconstruction();

    bool reconstruct(Reconstruction& recon,
		     Method method, ScopeWeighting weighting,
		     const std::vector<ScopeImage>& images) const;
    
    void calculateArrayMoments(ArrayMoments& eventmoments,
			       ScopeWeighting weighting,
			       const std::vector<ScopeImage>& image,
			       const ArrayQualityCuts* qc = 0) const;
    
    void calculateArrayParameters(ArrayParameters& param,
				  const std::vector<ScopeImage>& images,
				  const Reconstruction& reconstruction) const;

    static VSAReconstructionInternals::Method1* method1();
    static VSAReconstructionInternals::Method2* method2();
    static VSAReconstructionInternals::Method3* method3();

    unsigned nscopes() const { return fArray.scopes.size(); }
    const ArrayInfo& array() const { return fArray; }

  private:
    ArrayInfo                        fArray;
    VSAAtmosphere                    fAtmo;
    ArrayQualityCuts                 fQualityCuts;
    //double                           fObsThick;
    double                           fObsDelay;
    static std::auto_ptr<VSAReconstructionInternals::Method1> sMethod1;
    static std::auto_ptr<VSAReconstructionInternals::Method2> sMethod2;
    static std::auto_ptr<VSAReconstructionInternals::Method3> sMethod3;
  };

}

#endif // VSARECONSTRUCTION_HPP
