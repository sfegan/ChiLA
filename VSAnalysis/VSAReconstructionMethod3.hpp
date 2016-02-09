//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstructionMethod3.hpp

  Implementation of reconstruction 3, 3D simultaneously in real
  space

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       12/09/2005
*/

#ifndef VSARECONSTRUCTIONMETHOD3_HPP
#define VSARECONSTRUCTIONMETHOD3_HPP

#define VSA_ASSERT

#include<vector>

#include<VSAAlgebra.hpp>
#include<VSAReconstructionCommon.hpp>

namespace VERITAS 
{
  namespace VSAReconstructionInternals
  {

    class Method3
    {
    public:
      Method3() { /* nothing to see here */ }
      virtual ~Method3();
      
      double 
      calculateChi2ForHypothesis(const VSAAlgebra::Vec3D& e,
				 VSAAlgebra::Vec3D& R,
				 VSAAlgebra::Symmetric3D& Y,
				 const ArrayMoments& event,
				 const std::vector<ScopeImage>& images,
				 const std::vector<ScopeInfo>& scopes,
				 bool minimize_R = true) const;
      
      struct MinimizationResults
      {
	MinimizationResults()
	  : chi2(0), e(), R(), at_edge(false), 
	    chi2_e_curv_l(), chi2_e_curv_w(), chi2_e_curv_axis(),
	    chi2_R_curv_l(), chi2_R_curv_w(), chi2_R_curv_axis() { }
	double              chi2;
	VSAAlgebra::Vec3D   e;
	VSAAlgebra::Vec3D   R;
	bool                at_edge;
	double              chi2_e_curv_l;
	double              chi2_e_curv_w;
	VSAAlgebra::Vec3D   chi2_e_curv_axis;
	double              chi2_R_curv_l;
	double              chi2_R_curv_w;
	VSAAlgebra::Vec3D   chi2_R_curv_axis;
      };
      
      template<class OutputIterator> OutputIterator 
      gridSearch(MinimizationResults& results, const ArrayMoments& event,
		 const std::vector<ScopeImage>& images,
		 const std::vector<ScopeInfo>& scopes,
		 const std::vector<GridSearchSpecs>& grid,
		 OutputIterator map_output);
      
      void gridSearch(MinimizationResults& results,
		      const ArrayMoments& event,
		      const std::vector<ScopeImage>& images,
		      const std::vector<ScopeInfo>& scopes,
		      const std::vector<GridSearchSpecs>& grid);
      
    private:
      
      class GridCalculator
      {
      public:
	GridCalculator(const Method3& reconstruction, 
		       const ArrayMoments& event,
		       const std::vector<ScopeImage>& images,
		       const std::vector<ScopeInfo>& scopes)
	  : fReconstruction(reconstruction),
	    fEvent(event), fImages(images), fScopes(scopes) { }
	
	double operator () (double x, double y) const;
	double calculate(double x, double y,
			 VSAAlgebra::Vec3D& e, VSAAlgebra::Vec3D& R,
			 VSAAlgebra::Symmetric3D& Y,
			 bool minimize = true) const;
	
      private:
	const Method3& fReconstruction;
	const ArrayMoments& fEvent;
	const std::vector<ScopeImage>& fImages;
	const std::vector<ScopeInfo>& fScopes;
      };
      
    };
    
    template<class OutputIterator> OutputIterator Method3::
    gridSearch(MinimizationResults& results, const ArrayMoments& event,
	       const std::vector<ScopeImage>& images,
	       const std::vector<ScopeInfo>& scopes,
	       const std::vector<GridSearchSpecs>& grid,
	       OutputIterator map_output)
    {
      GridCalculator calculator(*this,event,images,scopes);
      GridSearchResults search_results;
      runGridSearch(search_results,calculator,grid,map_output);
      GridCurvatureResults curv;
      calculateGridCurvature(curv,search_results,calculator);

      results.chi2             = search_results.min_val;
      results.e                = curv.min_e;
      results.R                = curv.min_R;
      results.at_edge          = search_results.min_at_edge;
      results.chi2_R_curv_l    = curv.curv_Y_l;
      results.chi2_R_curv_w    = curv.curv_Y_w;
      results.chi2_R_curv_axis = curv.curv_Y_axis;
      results.chi2_e_curv_l    = curv.curv_grid_l;
      results.chi2_e_curv_w    = curv.curv_grid_w;
      results.chi2_e_curv_axis.set(cos(curv.curv_grid_theta),
				   sin(curv.curv_grid_theta), 0);
      
#warning Recheck this calculation
      double x = 
	search_results.min_x+search_results.min_step*cos(curv.curv_grid_theta);
      double y = 
	search_results.min_y+search_results.min_step*sin(curv.curv_grid_theta);
      
      results.chi2_e_curv_axis.set(0,0,-1);
      results.chi2_e_curv_axis.rotate(VSAAlgebra::Vec3D(0,-sqrt(x*x+y*y),0));
      results.chi2_e_curv_axis.rotate(VSAAlgebra::Vec3D(0,0,atan2(y,x)));
      
      x = search_results.min_x;
      y = search_results.min_y;
      results.chi2_e_curv_axis.rotate(VSAAlgebra::Vec3D(0,0,-atan2(y,x)));
      results.chi2_e_curv_axis.rotate(VSAAlgebra::Vec3D(0,sqrt(x*x+y*y),0));

      results.chi2_e_curv_axis = 
	VSAAlgebra::Vec3D(0,0,-1) - results.chi2_e_curv_axis;
      results.chi2_e_curv_axis.normalize();
      
      results.chi2_e_curv_axis.rotate(VSAAlgebra::Vec3D(0,-sqrt(x*x+y*y),0));
      results.chi2_e_curv_axis.rotate(VSAAlgebra::Vec3D(0,0,atan2(y,x)));
      
      return map_output;
    }

  } // namespace VSAReconstructionInternals


} // namespace VERITAS

#endif // VSARECONSTRUCTIONMETHOD3_HPP
